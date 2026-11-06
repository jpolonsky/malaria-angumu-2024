# Builds the de-identified analysis files in data/ (df_ind.rds, df_hh.rds,
# entomology_*.rds, pnlp_health_area_monthly.rds, dhis2_monthly.rds,
# climate_weekly.rds and the map layers in data/maps/) from the study team's
# internal data store: the Kobo survey snapshots, the PNLP and DHIS2
# surveillance exports and the entomology workbooks. It cannot be run outside
# that store; it is included to document how the files used by the analysis
# scripts are derived. It runs from the root of the study team's project and
# sources that project's per-round survey preparation scripts, which clean each
# round and assign the as-treated zone. Names and device, submission and
# household identifiers are dropped, and household and individual keys are
# replaced by integers that keep the original sort order. Household locations
# (lon, lat) stay in df_ind.rds, which is not released.

pacman::p_load(janitor, sf, tidyverse)

out_dir <- "../malaria-angumu-2024/data"
dir.create(file.path(out_dir, "maps"), recursive = TRUE, showWarnings = FALSE)

# Survey rounds ----
# Each round script leaves df_hh and df_ind in the global environment, so the
# objects are captured before the next script overwrites them.
round_scripts <- c(
  "code/prep/02_data_prep_2022-03.R", "code/prep/03_data_prep_2022-11.R",
  "code/prep/04_data_prep_2023-03.R", "code/prep/05_data_prep_2023-10.R",
  "code/prep/06_data_prep_2024-04.R"
)
rounds <- map(round_scripts, \(path) {
  source(path)
  list(hh = df_hh, ind = df_ind)
})

# Coding of questionnaire items ----
ascii <- \(x) stringi::stri_trans_general(as.character(x), "Latin-ASCII")
item_chr <- \(x) str_to_lower(str_squish(ascii(x)))

# Net item to any net / intact / holed / none; "don't know" and missing are NA.
# Survey 1 codes the item -1 / 1 / 2 / 3, surveys 3 and 5 carry text labels.
net_any_code <- \(x) {
  x <- item_chr(x)
  case_when(x %in% c("1", "2") | str_detect(x, "^oui") ~ 1L, x == "3" | x == "non" ~ 0L, .default = NA_integer_)
}
net_cat_code <- \(x) {
  x <- item_chr(x)
  case_when(
    x == "1" | str_detect(x, "bonne") ~ "intact", x == "2" | str_detect(x, "mauvais") ~ "holed",
    x == "3" | x == "non" ~ "none", .default = NA_character_
  )
}

# MDA receipt item to card / recall / no / dk; missing stays NA.
receipt_code <- \(x) {
  x <- item_chr(x)
  case_when(
    str_detect(x, "^oui selon la carte") ~ "card", str_detect(x, "^oui selon le rappel") ~ "recall",
    x == "non" ~ "no", str_detect(x, "^ne sait pas") ~ "dk", .default = NA_character_
  )
}
yes_no_code <- \(x) {
  x <- item_chr(x)
  case_when(x %in% c("oui", "yes") ~ "yes", x %in% c("non", "no") ~ "no", .default = NA_character_)
}

# Household net ownership, surveys 1, 3 and 5, from the raw household sheets.
hh_net <- function(path, sheet, survey) {
  h <- read_rds(path)[[sheet]] |> clean_names()
  tibble(
    survey = survey,
    index_hh = as.character(h$index),
    hh_own_net = case_when(
      as.character(haven::as_factor(h$moustiquaire_own)) %in% c("1", "Oui") ~ 1L,
      as.character(haven::as_factor(h$moustiquaire_own)) %in% c("0", "Non") ~ 0L,
      .default = NA_integer_
    ),
    hh_net_count = suppressWarnings(as.numeric(as.character(haven::as_factor(h$moustiquaire_simple))))
  )
}
hh_net_lookup <- bind_rows(
  hh_net("data/srvy/srvy_raw_2022-03.rds", "HH", 1L),
  hh_net("data/srvy/srvy_raw_2023-03.rds", "main", 3L),
  hh_net("data/srvy/srvy_raw_2024-04.rds", "main", 5L)
)
stopifnot(!anyDuplicated(hh_net_lookup[c("survey", "index_hh")]))

# Survey 2: the prepared rows carry no receipt items or the file's own weight
# (pw), so they are matched to the raw sheet on household key and roster position.
raw2 <- read_rds("data/srvy/srvy_raw_2022-11.rds")
raw2_key <- paste(paste0(raw2$`_parent_table_name`, "|", raw2$v001, "-", raw2$index_v2), raw2$RangPresent)
match2 <- match(paste(rounds[[2]]$ind$index_hh, rounds[[2]]$ind$hh_num), raw2_key)
stopifnot(!anyDuplicated(raw2_key), !anyNA(match2))

# Individuals ----
# Survey 1 codes sex 0 = Femme, 1 = Homme; every other round has the labels.
rounds[[1]]$ind <- rounds[[1]]$ind |>
  mutate(sex = as.character(sex) |> case_match("0" ~ "Femme", "1" ~ "Homme", .default = NA))

rounds[[2]]$ind <- rounds[[2]]$ind |>
  mutate(
    mda_r1 = receipt_code(raw2$C5EstcequePrenomarecu[match2]),
    mda_r2 = receipt_code(raw2$C8EstcequePrenomarecu[match2]),
    irs_hh = yes_no_code(raw2$B5Uninsecticideatilétép[match2]),
    pw = raw2$pw[match2]
  )

rounds[[3]]$ind <- rounds[[3]]$ind |>
  mutate(mda_r3 = receipt_code(mda_round_3), mda_r4 = receipt_code(mda_round_4), irs_hh = yes_no_code(irs))

# Survey 4 has no usable net item, so its net columns stay missing.
net_rounds <- c(1, 3, 5)
for (i in net_rounds) {
  rounds[[i]]$ind <- rounds[[i]]$ind |>
    mutate(net_any = net_any_code(bednet_sleep), net_cat = net_cat_code(bednet_sleep))
}

df_ind <-
  map2(rounds, 1:5, \(r, survey) {
    xy <- st_coordinates(r$ind$geometry)
    r$ind |>
      as_tibble() |>
      mutate(
        survey = survey,
        index_hh = as.character(index_hh),
        index_ind = as.character(index_ind),
        zone = factor(as.character(zone), levels = c("Routine", "IRS", "MDA")),
        across(any_of(c("sex", "status", "rdt_result", "died_cause")), as.character),
        lon = xy[, 1],
        lat = xy[, 2]
      ) |>
      left_join(hh_net_lookup, by = c("survey", "index_hh")) |>
      select(any_of(c(
        "survey", "index_hh", "index_ind", "cluster", "area", "zone", "age", "age_grp_few", "sex", "status",
        "hh_size", "weight", "rdt_result", "net_any", "net_cat", "hh_own_net", "hh_net_count",
        "mda_r1", "mda_r2", "mda_r3", "mda_r4", "irs_hh", "pw",
        "died", "died_cause", "date_died", "pt_start", "pt_end", "pt", "lon", "lat"
      )))
  }) |>
  list_rbind()

# Households ----
df_hh <-
  map2(rounds, 1:5, \(r, survey) {
    r$hh |>
      as_tibble() |>
      transmute(
        survey = survey,
        index_hh = as.character(index_hh),
        cluster, area,
        zone = factor(as.character(zone), levels = c("Routine", "IRS", "MDA")),
        consent = as.character(consent),
        date_int = as_date(date_int)
      )
  }) |>
  list_rbind()

# De-identification ----
# Integer codes follow the C-locale sort order of the original keys, and all have
# the same number of digits, so they sort the same way as numbers and as text
# (mixed-effects fits order their grouping levels as text). Household keys repeat
# across surveys (the same number is a different household), so they are coded
# on one scale: the code is as unique as the original key.
code <- function(x, universe = x) {
  keys <- sort(unique(universe), method = "radix")
  match(x, keys) + as.integer(10^nchar(length(keys)))
}

hh_keys <- c(df_ind$index_hh, df_hh$index_hh)
df_ind <- df_ind |>
  mutate(index_hh = code(index_hh, hh_keys), index_ind = code(index_ind))
df_hh <- df_hh |>
  mutate(index_hh = code(index_hh, hh_keys))

write_rds(df_ind, file.path(out_dir, "df_ind.rds"))
write_rds(df_hh, file.path(out_dir, "df_hh.rds"))

# Entomology ----
# The preparation script cleans the trap-night and mosquito sheets and audits
# the cleaning; its own file writes are not wanted here.
ento_env <- new.env()
ento_env$write_csv <- ento_env$write_rds <- \(...) invisible(NULL)
source("code/report/manuscript/entomology_prep.R", local = ento_env)

# Trap-nights and mosquitoes are linked on the sheet index, which is recoded.
tn_keys <- ento_env$ento_trapnights$tn_index
ento_env$ento_trapnights |>
  mutate(tn_index = code(tn_index, tn_keys)) |>
  write_rds(file.path(out_dir, "entomology_trapnights.rds"))
ento_env$ento_mosq |>
  mutate(tn_index = code(tn_index, tn_keys)) |>
  select(-code_beem) |>
  write_rds(file.path(out_dir, "entomology_mosquitoes.rds"))
write_rds(ento_env$ento_audit, file.path(out_dir, "entomology_audit.rds"))

# Routine surveillance ----
source("code/misc/surveillance_zones.R")

pnlp <- read_rds("data/misc/pnlp_secondary.rds")

pnlp$indicators_aire |>
  filter(indicateur_name %in% c("Positivite Tdr", "Cas Suspects Teste", "Cas Confirme Traité")) |>
  transmute(
    area = clean_surveillance_area(aire_de_sante),
    period, indicateur_name, n_fosa, taux_moyen, total_absolu
  ) |>
  write_rds(file.path(out_dir, "pnlp_health_area_monthly.rds"))

pnlp$climate |>
  select(year, epi_week, precipitation, tmax_mean, tmin_mean) |>
  write_rds(file.path(out_dir, "climate_weekly.rds"))

read_rds("data/misc/surveillance.rds")$monthly_dhis2 |>
  filter(metric %in% c("cas_suspect", "tdr_positif")) |>
  transmute(area = clean_surveillance_area(health_area), date, metric, value) |>
  write_rds(file.path(out_dir, "dhis2_monthly.rds"))

# Map layers ----
boundaries <- read_rds("data/misc/boundaries.rds")

boundaries$shp_health_areas |>
  select(pcode, name_fr, name_en) |>
  write_rds(file.path(out_dir, "maps/health_areas.rds"))
boundaries$lake_albert |>
  select(osm_id, name) |>
  write_rds(file.path(out_dir, "maps/lake_albert.rds"))

# Province outlines and health areas (GRID3, CC BY), the national boundary (geoBoundaries,
# CC BY 4.0) and lakes (Natural Earth, public domain) are copied unchanged.
file.copy(
  file.path(
    "data/misc/maps",
    c("drc_provinces_grid3_v9.gpkg", "grid3_cod_health_areas_v9_angumu.geojson", "geoBoundaries-COD-ADM0_simplified.geojson", "ne_10m_lakes.rds")
  ),
  file.path(out_dir, "maps"),
  overwrite = TRUE
)
