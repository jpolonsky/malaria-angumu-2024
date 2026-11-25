# Paths below are relative to the study team's project root, where prep.R runs.
if (!identical(basename(here::here()), "angumu_malaria")) stop("here() resolved to '", here::here(), "' instead of the expected project root 'angumu_malaria' — check the working directory before running this script.")
## Entomology data preparation and audit: the cleaning steps behind
## data/entomology_trapnights.rds and data/entomology_mosquitoes.rds.
## Sourced by prep.R from the study team's project root (the workbooks are
## not public); identical to the script in the study team's pipeline.
##
## Input: data/misc/entomology/Result Angumu/Angumu_Entomo_Results.xlsx
## (Arnaud Desjardins's cleaned workbook of 2026-10-15, sheet 'Trap_night' =
## one row per trap-night, sheet 'Mosquitos' = one row per collected
## Anopheles; laboratory results re-entered, so a mosquito with a result in
## 'Sporozoite positif ?' was tested and an empty cell was not). It is a
## superset of the 2026-10-07 workbook (Analyse Entomo.xlsx, same _index),
## which supplies the three fields the cleaned file dropped (submission time,
## indoor/outdoor, trap number) and the trap dates that were overwritten
## during the 2025 re-entry. Written so a refreshed workbook re-runs through
## the same code: nothing below depends on the 2026-10-07 row counts, which
## are only compared against in the audit table. Analysis window: to
## 2024-12-31 (Arnaud's recommendation; later rows are read and dropped).
##
## Outputs (objects in the calling session):
##   ento_trapnights  one row per analysed trap-night (duplicates, outdoor
##                    traps and bad dates handled per the audit rules)
##   ento_mosq        one row per Anopheles linked to an analysed trap-night
##   ento_audit       the audit table (also written to CSV)
##   ento_periods     half-year period lookup
## and, on disk, output/data/tables/entomology_audit.csv and
## output/data/entomology_analysis.rds (no GPS, house codes only: house
## coordinates are identifying).
##
## Zone is derived from health area and date, never from the workbook's own
## zone column (11 Awasi trap-nights are coded Nord there). Underlying zone
## values follow the project convention (Routine, IRS, MDA); reader-facing
## relabelling to "MDA+IRS" happens in the table/figure script.

pacman::p_load(here, readxl, sf, tidyverse)

ento_xlsx <- here("data/misc/entomology/Result Angumu/Angumu_Entomo_Results.xlsx")
ento_xlsx_ref <- here("data/misc/entomology/Analyse Entomo.xlsx")
ento_window_end <- as.Date("2024-12-31")

## ── Fixed definitions (from vector_control_recap.xlsx and the plan) ─────────

## Original assignment of each health area to a study zone
ento_area_zone <- tribble(
  ~area,          ~zone_assigned,
  "Abia",         "Routine",
  "Ara",          "Routine",
  "Awasi",        "Routine",
  "Cawa",         "Routine",
  "Linga",        "Routine",
  "Anzika",       "IRS",
  "Pakala",       "IRS",
  "Panyadong",    "IRS",
  "Udongo Abira", "IRS",
  "Apala",        "MDA",
  "Gengere",      "MDA",
  "Lanyi",        "MDA",
  "Ugudo-Zii",    "MDA",
  "Uyandu",       "MDA"
)

## Date from which a Routine area leaves the zone comparison (as-treated
## rule: dropped from the date of its first intervention, never merged into
## the zone it crossed into). Abia and Cawa: IRS Aug-Sep 2023; Ara, Linga
## (and Cawa): MDA Jan-Feb 2024; Awasi: IRS phase III from Dec 2024.
## Month starts are used as the conservative cut-off.
ento_crossover_date <- tribble(
  ~area,   ~out_from,
  "Abia",  as.Date("2023-08-01"),
  "Cawa",  as.Date("2023-08-01"),
  "Ara",   as.Date("2024-01-01"),
  "Linga", as.Date("2024-01-01"),
  "Awasi", as.Date("2024-12-01")
)

## Half-year periods, each containing one survey (plan section 3)
ento_periods <- tribble(
  ~period, ~p_start,              ~p_end,                ~survey_in_period,
  "P0",    as.Date("2022-05-01"), as.Date("2022-06-30"), NA_real_,
  "P1",    as.Date("2022-07-01"), as.Date("2022-12-31"), 2,
  "P2",    as.Date("2023-01-01"), as.Date("2023-06-30"), 3,
  "P3",    as.Date("2023-07-01"), as.Date("2023-12-31"), 4,
  "P4",    as.Date("2024-01-01"), as.Date("2024-06-30"), 5,
  "P5",    as.Date("2024-07-01"), as.Date("2024-11-30"), NA_real_
)
irs_round2_end <- as.Date("2023-05-31")

assign_period <- function(d) {
  out <- rep(NA_character_, length(d))
  for (i in seq_len(nrow(ento_periods))) out[!is.na(d) & d >= ento_periods$p_start[i] & d <= ento_periods$p_end[i]] <- ento_periods$period[i]
  out
}

## Zone as treated on a given date; NA when the area has left the comparison
zone_as_treated <- function(area, date) {
  base <- ento_area_zone$zone_assigned[match(area, ento_area_zone$area)]
  out_from <- ento_crossover_date$out_from[match(area, ento_crossover_date$area)]
  if_else(!is.na(out_from) & date >= out_from, NA_character_, base)
}

## ── Read ─────────────────────────────────────────────────────────────────────

stopifnot("entomology workbook not found" = file.exists(ento_xlsx), "reference workbook (Analyse Entomo.xlsx) not found" = file.exists(ento_xlsx_ref))
raw_tn <- suppressWarnings(read_excel(ento_xlsx, sheet = "Trap_night", .name_repair = "minimal"))
raw_mq <- suppressWarnings(read_excel(ento_xlsx, sheet = "Mosquitos", .name_repair = "minimal"))
stopifnot("Index and _index differ" = all(raw_tn$Index == raw_tn$`_index`))

## Fields the cleaned workbook dropped, from the 2026-10-07 workbook by _index
ref_tn <-
  suppressWarnings(read_excel(ento_xlsx_ref, sheet = "Enquete Angumu", .name_repair = "minimal")) |>
  transmute(
    tn_index = `_index`,
    date_ref_recorded = as.Date(`Enter a date`),
    date_ref_submitted = as.Date(`_submission_time`),
    indoor_ref = `Interieur ou Exterieur ?` == "Interieur",
    trap_no = `Numero du Piege`
  )

area_crosswalk <- c(
  "ABIA" = "Abia", "ANZIKA" = "Anzika", "APALA" = "Apala", "ARA" = "Ara", "AWASI" = "Awasi",
  "CAWA ANJU" = "Cawa", "GENGERE" = "Gengere", "LANYI" = "Lanyi", "LINGA" = "Linga",
  "PAKALA" = "Pakala", "PANYADONG" = "Panyadong", "UDONGO ABIRA" = "Udongo Abira",
  "UGUDO ZII" = "Ugudo-Zii", "UYANDU" = "Uyandu"
)

tn <-
  raw_tn |>
  transmute(
    tn_index = `_index`,
    date_recorded = as.Date(`Enter a date`),
    area_raw = `Aire De Santé`,
    zone_workbook = `Zone d'intervention`,
    site = `Site de capture`,
    village_code = `Code Village`,
    random_house = `Est une maison selectionnée Aleatoirement (symbole Cercle) ou Volontairement ( symbole Carré)` == "Cercle",
    night = `Nuit de capture`,
    working = `Piege a fonctionné ?` == "oui",
    persons = `Nombre de personnes qui habitent dans la maison`,
    n_culicid = `Nombre de Culicidés`,
    n_anoph = `Nombre d'Anopheles`,
    walls = `Nature des murs`,
    roof = `Nature du toit`,
    nets_owned = `Nombre de moustiquaires possedées par le menage`,
    nets_used = `Nombre de personnes qui ont dormi sous moustiquaire la nuit derniere`,
    gps = `Record your current location`
  ) |>
  left_join(ref_tn, by = "tn_index") |>
  mutate(
    lat = suppressWarnings(as.numeric(str_split_fixed(gps, " ", 4)[, 1])),
    lon = suppressWarnings(as.numeric(str_split_fixed(gps, " ", 4)[, 2])),
    indoor = coalesce(indoor_ref, TRUE),
    area = unname(area_crosswalk[area_raw]),
    site_id = str_c(area, " / ", site)
  )

## ── Audit ────────────────────────────────────────────────────────────────────
## Each check records what the rule did and the value seen in the 2026-10-07
## workbook (`expected`), so a refreshed workbook shows where it differs.

audit <- tibble(id = character(), check = character(), n = numeric(), expected = numeric(), action = character())
add_audit <- function(id, check, n, expected = NA_real_, action = "") {
  audit <<- bind_rows(audit, tibble(id = id, check = check, n = as.numeric(n), expected = as.numeric(expected), action = action))
}

## A1 date range and missing dates
stopifnot("unmapped health-area names (extend area_crosswalk)" = !anyNA(tn$area))
add_audit("A1a", "trap-nights with a missing recorded date", sum(is.na(tn$date_recorded)), 0, "stop if > 0")
stopifnot("missing recorded dates" = !anyNA(tn$date_recorded))
add_audit("A1b", "trap-nights with recorded date before 2022-05-01 or after 2025-12-31", sum(tn$date_recorded < as.Date("2022-05-01") | tn$date_recorded > as.Date("2025-12-31")), 0, "stop if > 0")

## A2 dates. Two defects, both repaired from the 2026-10-07 workbook by _index:
## (a) recorded date later than submission date: the date is set to the
## submission date; (b) in the cleaned workbook 47 trap-nights carry a date in
## May-June 2025 (the re-entry), 5 to 27 months after the trap-night: the
## 2026-10-07 date is restored. Trap-nights not in the 2026-10-07 workbook have
## no second date to check against: those dated inside that workbook's coverage
## are held out (A2d), later ones are outside the analysis window anyway.
early <- !is.na(tn$date_ref_submitted) & tn$date_ref_submitted < tn$date_ref_recorded
overwritten <- !is.na(tn$date_ref_recorded) & tn$date_recorded != tn$date_ref_recorded
add_audit("A2a", "trap-nights submitted before their recorded date (2026-10-07 workbook)", sum(early), 135, "date set to the submission date; the night offset (up to 3 days) cannot be recovered and moves at most 0.9% of any estimate")
add_audit("A2b", "trap-nights whose date in the cleaned workbook differs from the 2026-10-07 workbook (all moved to 2025)", sum(overwritten), 47, "2026-10-07 date restored; ask Arnaud to confirm")
add_audit("A2c", "trap-nights not in the 2026-10-07 workbook", sum(is.na(tn$date_ref_recorded)), 6246, "date and indoor status taken as recorded; trap number unknown")
tn <-
  tn |>
  mutate(date = if_else(!is.na(date_ref_recorded), pmin(date_ref_recorded, date_ref_submitted, na.rm = TRUE), date_recorded), date_fixed = date != date_recorded)
## A2d trap-nights absent from the 2026-10-07 workbook but dated inside its coverage
## (52, all Awasi: 38 dated 2023-08-27 whose activity-month counter, night
## sequence and neighbouring submissions point to December 2024, and 14 dated
## 2024-10-08 that look genuine). Held out until Arnaud confirms their dates;
## the 38 would fall outside Routine if they are from December 2024.
unverified <- is.na(tn$date_ref_recorded) & tn$date <= max(tn$date_ref_recorded, na.rm = TRUE)
add_audit("A2d", "trap-nights absent from the 2026-10-07 workbook but dated within its coverage", sum(unverified), 52, "held out until Arnaud confirms the dates (_index 22313-22357 look like December 2024)")
tn <- tn |> filter(!unverified)
add_audit("A1c", "trap-nights after the analysis window (2024-12-31)", sum(tn$date > ento_window_end), NA, "dropped (Arnaud: stop at December 2024)")
tn <- tn |> filter(date <= ento_window_end)
stopifnot("dates outside the study window after correction" = !any(tn$date < as.Date("2022-05-01") | tn$date > ento_window_end))
stopifnot("trap-nights without GPS" = !anyNA(tn$lat) && !anyNA(tn$lon))

## GPS in UTM 36N (metres), used by the duplicate rules and the house grouping
tn_sf <- st_as_sf(tn, coords = c("lon", "lat"), crs = 4326, remove = FALSE) |> st_transform(32636)
tn$x <- st_coordinates(tn_sf)[, 1]
tn$y <- st_coordinates(tn_sf)[, 2]

## A3 duplicates: identical on every recorded field and within 50 m of an
## earlier record (two records identical in every field but 350 m apart are
## separate houses, not a double upload)
dup_key <- c("date", "area", "site", "night", "trap_no", "n_anoph", "n_culicid", "persons", "random_house", "working", "walls", "roof", "nets_owned", "nets_used")
tn <- tn |> mutate(dupkey = do.call(paste, c(across(all_of(dup_key)), sep = "|")))
dup3 <- rep(FALSE, nrow(tn))
for (g in unique(tn$dupkey[duplicated(tn$dupkey)])) {
  idx <- which(tn$dupkey == g)
  for (j in idx[-1]) {
    earlier <- idx[idx < j & !dup3[idx]]
    if (any(sqrt((tn$x[earlier] - tn$x[j])^2 + (tn$y[earlier] - tn$y[j])^2) <= 50)) dup3[j] <- TRUE
  }
}
add_audit("A3", "duplicate trap-nights (identical fields, within 50 m of an earlier record)", sum(dup3), 15, "dropped")
tn <- tn[!dup3, ]

## A4 same date/area/site/trap/night recorded more than once
tn <- tn |> mutate(grp4 = paste(date, area, site, night, trap_no, sep = "|"))
grp4_n <- tn |> count(grp4)
add_audit("A4a", "date/area/site/trap/night groups recorded more than once", sum(grp4_n$n > 1), 148, "see A4b")
drop4 <- rep(FALSE, nrow(tn))
for (g in grp4_n$grp4[grp4_n$n > 1]) {
  idx <- which(tn$grp4 == g)
  for (j in idx[-1]) {
    earlier <- idx[idx < j & !drop4[idx]]
    close_same <- earlier[sqrt((tn$x[earlier] - tn$x[j])^2 + (tn$y[earlier] - tn$y[j])^2) <= 50 & tn$n_anoph[earlier] %in% tn$n_anoph[j] & tn$n_culicid[earlier] %in% tn$n_culicid[j]]
    if (length(close_same) > 0) drop4[j] <- TRUE
  }
}
add_audit("A4b", "repeat records within 50 m with identical counts", sum(drop4), NA, "dropped; other repeats kept as separate houses and listed for Arnaud")
tn <- tn[!drop4, ]

## A5 trap did not work
add_audit("A5a", "trap-nights where the trap did not work", sum(!tn$working), 465, "excluded from the trap-night denominator")
add_audit("A5b", "non-working trap-nights that recorded Anopheles", sum(!tn$working & tn$n_anoph > 0, na.rm = TRUE), 66, "mosquitoes kept for the sporozoite denominator once tested")

## A6 implausible counts
add_audit("A6a", "negative Anopheles counts", sum(tn$n_anoph < 0, na.rm = TRUE), 1, "set to missing")
add_audit("A6b", "largest single-trap-night Anopheles count", max(tn$n_anoph, na.rm = TRUE), 433, "kept; winsorised sensitivity")
tn <- tn |> mutate(n_anoph = if_else(n_anoph < 0, NA_real_, n_anoph))

## A7 household size
add_audit("A7", "household size 0 or above 20", sum(tn$persons == 0 | tn$persons > 20, na.rm = TRUE), 4, "set to missing (only used in the per-person sensitivity)")
tn <- tn |> mutate(persons = if_else(persons == 0 | persons > 20, NA_real_, persons))

## A8 indoor/outdoor
add_audit("A8", "traps set outdoors", sum(!tn$indoor), 5, "dropped")
tn <- tn |> filter(indoor)

## A9 workbook zone vs zone implied by the health area
zone_label_wb <- c("Zone Centre" = "Routine", "Zone Nord" = "IRS", "Zone Sud" = "MDA")
zone_implied <- ento_area_zone$zone_assigned[match(tn$area, ento_area_zone$area)]
add_audit("A9", "trap-nights whose workbook zone differs from the zone implied by the health area", sum(unname(zone_label_wb[tn$zone_workbook]) != zone_implied, na.rm = TRUE), 11, "zone taken from the health area")

## ── Zone as treated, period, house ──────────────────────────────────────────

tn <-
  tn |>
  mutate(
    zone_assigned = ento_area_zone$zone_assigned[match(area, ento_area_zone$area)],
    zone = zone_as_treated(area, date),
    period = assign_period(date),
    month = floor_date(date, "month")
  )

## A15 houses: no identifier in the workbook; group GPS points within 25 m
## inside a health area (complete linkage, so a house cannot chain along a
## row of neighbouring compounds)
tn$house_id <- NA_character_
for (a in unique(tn$area)) {
  idx <- which(tn$area == a)
  if (length(idx) == 1) {
    tn$house_id[idx] <- paste0(a, "_h1")
    next
  }
  cl <- cutree(hclust(dist(cbind(tn$x[idx], tn$y[idx])), method = "complete"), h = 25)
  tn$house_id[idx] <- paste0(a, "_h", cl)
}
tn <- tn |> mutate(house_id = str_replace_all(house_id, " ", "_"))
add_audit("A15a", "houses after grouping GPS points within 25 m (all areas)", n_distinct(tn$house_id), NA, "complete linkage splits some real houses; resampling uses sites, houses are a sensitivity")
add_audit("A15c", "sites (health area / site name) after cleaning", n_distinct(tn$site_id), NA, "primary resampling unit")
add_audit("A15b", "distinct area/site/trap-number combinations (alternative house key)", n_distinct(tn$area, tn$site, tn$trap_no), NA, "compare with A15a")

## A16 trapping effort (reported in the CSV, not acted on)
effort <- tn |> filter(!is.na(period)) |> group_by(period, zone_assigned) |> summarise(trap_nights = n(), houses = n_distinct(house_id), .groups = "drop")
for (i in seq_len(nrow(effort))) add_audit(paste0("A16_", effort$period[i], "_", effort$zone_assigned[i]), paste0("trap-nights by originally assigned zone ", effort$period[i], " ", effort$zone_assigned[i]), effort$trap_nights[i], NA, "effort")

## ── Mosquito records ────────────────────────────────────────────────────────

mq <-
  raw_mq |>
  transmute(
    tn_index = `_parent_index`,
    code_beem = `CODE BEEM`,
    species_field = `Especes Anopheles`,
    abdomen = `Etat de l'abdomen`,
    spz_workbook = `Sporozoite positif ?`
  )

## A10 mosquito rows against the trap-night's own Anopheles count (all raw trap-nights)
mq_per_tn <- mq |> count(tn_index, name = "n_rows")
cmp <- tn |> select(tn_index, n_anoph) |> left_join(mq_per_tn, by = "tn_index") |> mutate(n_rows = replace_na(n_rows, 0))
add_audit("A10", "trap-nights (within the window) whose mosquito rows differ from the recorded Anopheles count", sum(cmp$n_rows != cmp$n_anoph, na.rm = TRUE), 8, "counts taken from the trap-night sheet")

## A11 mosquito code uniqueness
add_audit("A11", "mosquito codes (CODE BEEM) used more than once", sum(table(mq$code_beem) > 1), 777, "not a unique ID: mosquitoes are linked to trap-nights by _parent_index")

mq <-
  mq |>
  inner_join(tn |> select(tn_index, date, area, site_id, zone, zone_assigned, period, month, house_id, working, random_house), by = "tn_index") |>
  mutate(
    spz_workbook = factor(spz_workbook, levels = c("Negative", "Positive")),
    species_group = case_when(species_field == "gambiae s.l" ~ "gambiae", species_field == "funestus s.l" ~ "funestus", .default = "other")
  )

## A13 species field
add_audit("A13", "mosquito rows with a species other than gambiae/funestus s.l.", sum(mq$species_group == "other"), 199, "counted as Anopheles; excluded in a sensitivity")

## ── Tested flag ──────────────────────────────────────────────────────────────
## In the cleaned workbook (2026-10-15) only mosquitoes with a laboratory
## result carry a value in 'Sporozoite positif ?'; an empty cell means not
## tested. (The 2026-10-07 workbook coded untested mosquitoes 'Negative'.)
## Arnaud: the 2022 results from the first laboratory (South Africa) were
## almost all negative (1 positive in 61) and he does not trust them (technique or sample storage), so
## sporozoite rate and EIR are reported from January 2023 only.
ento_tested_source <- "Angumu_Entomo_Results.xlsx, column 'Sporozoite positif ?': a result means tested, empty means not tested; 2022 results (61 mosquitoes) not used (laboratory reliability doubted by Arnaud)"
ento_spz_start <- as.Date("2023-01-01")
mq <-
  mq |>
  mutate(
    tested = !is.na(spz_workbook) & date >= ento_spz_start,
    spz_result = if_else(tested, as.character(spz_workbook), NA_character_)
  )

## A12 sporozoite coding by half-year
spz_by_period <- mq |> filter(!is.na(period)) |> group_by(period) |> summarise(rows = n(), with_result = sum(!is.na(spz_workbook)), positive = sum(spz_workbook == "Positive", na.rm = TRUE), .groups = "drop")
for (i in seq_len(nrow(spz_by_period))) add_audit(paste0("A12_", spz_by_period$period[i]), paste0("mosquito rows ", spz_by_period$period[i]), spz_by_period$rows[i], NA, paste0("with a laboratory result: ", spz_by_period$with_result[i], "; positives: ", spz_by_period$positive[i]))
add_audit("B1", "mosquito rows used as tested (result recorded, from 2023-01-01)", sum(mq$tested), NA, "sporozoite rate and EIR denominators")
add_audit("B1b", "sporozoite positives among them", sum(mq$spz_result %in% "Positive"), NA, "")

## ── Analysis frames ─────────────────────────────────────────────────────────

ento_trapnights <- tn |> select(tn_index, date, month, period, area, site_id, zone_assigned, zone, house_id, random_house, night, trap_no, working, persons, n_culicid, n_anoph, date_fixed)
ento_mosq <- mq |> select(tn_index, code_beem, date, month, period, area, site_id, zone_assigned, zone, house_id, random_house, working, species_group, abdomen, spz_workbook, tested, spz_result)
ento_audit <- audit

write_csv(ento_audit, here("output/data/tables/entomology_audit.csv"))
dir.create(here("output/data"), showWarnings = FALSE, recursive = TRUE)
write_rds(list(trapnights = ento_trapnights, mosq = ento_mosq, periods = ento_periods, tested_source = ento_tested_source), here("output/data/entomology_analysis.rds"))
