# Study area, design and participants: Figure 1 (study-area map), Figure 2 (timeline and zone
# grid), S1 Table (baseline characteristics by zone), S4 Table (clusters and households),
# S7 Table (participant characteristics by survey) and S8 Table (MDA and IRS coverage).

pacman::p_load(sf, ggspatial, cowplot, rnaturalearth, ggforce, patchwork, showtext, sysfonts, gtsummary, srvyr, tidyverse)

source("code/functions.R")

sf_use_s2(FALSE)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

df_ind <- read_rds("data/df_ind.rds")
df_hh <- read_rds("data/df_hh.rds")

# Helpers ----

# Share of a 0/1 outcome (`y`) in percent with a logit confidence interval, single-stage
# cluster design stratified by zone. `d` needs y, cluster, zone and weight.
prop_logit <- function(d) {
  d <- d |> filter(!is.na(cluster), !is.na(zone)) |> mutate(y = as.numeric(y))
  des <- survey::svydesign(ids = ~cluster, strata = ~zone, weights = ~weight, data = d)
  ci <- survey::svyciprop(~y, des, method = "logit", level = 0.95)
  tibble(n = nrow(d), prop = 100 * as.numeric(ci), prop_low = 100 * attr(ci, "ci")[[1]], prop_upp = 100 * attr(ci, "ci")[[2]])
}

# Figure 1 ----

# Zone assignment at survey 1 (before any crossover), keyed on the health-area names in the map.
sampled_zones <- study_areas |> select(name_en = map_name, zone = baseline_zone)

provinces <- st_read("data/maps/drc_provinces_grid3_v9.gpkg", layer = "provinces", quiet = TRUE) |> st_make_valid()
ituri <- provinces |> filter(province == "Ituri")

# Health areas: GRID3 COD Health Areas v9.0 (CC BY), the 22 areas of Angumu.
health_areas <-
  st_read("data/maps/grid3_cod_health_areas_v9_angumu.geojson", quiet = TRUE) |>
  rename(name_en = airesante) |>
  st_transform(4326) |>
  st_make_valid() |>
  left_join(sampled_zones, by = "name_en") |>
  mutate(zone = fct_relevel(factor(zone), "Routine", "IRS", "MDA"))
stopifnot(sum(!is.na(health_areas$zone)) == nrow(sampled_zones))

lake_albert <- read_rds("data/maps/lake_albert.rds") |> st_transform(4326) |> st_make_valid()

# Merged outline of Angumu; the buffer out and back in closes slivers between neighbouring areas.
angumu <-
  health_areas |>
  st_union() |>
  st_buffer(0.0008) |>
  st_buffer(-0.0008) |>
  st_make_valid() |>
  st_sf(geometry = _)

lakes <- read_rds("data/maps/ne_10m_lakes.rds") |> st_make_valid()

countries <-
  ne_countries(scale = "medium", returnclass = "sf") |>
  filter(name %in% c("Uganda", "South Sudan", "Rwanda", "Burundi", "Central African Rep.")) |>
  select(name)

drc <- st_read("data/maps/geoBoundaries-COD-ADM0_simplified.geojson", quiet = TRUE)

ituri_col <- "#fbf3b0"
angumu_col <- "#f0c800"
lake_col <- "#a9d6f0"
land_col <- "grey94"
border_col <- "grey55"
ink <- "grey15"

anchor <- function(x) {
  xy <- x |> st_geometry() |> st_point_on_surface() |> st_coordinates()
  x |> st_drop_geometry() |> mutate(x = xy[, 1], y = xy[, 2])
}

# Left panel: Ituri in its regional context.
bb <- st_bbox(angumu)
pad <- 0.06
zoom_box <-
  st_bbox(
    c(xmin = bb[["xmin"]] - pad, xmax = bb[["xmax"]] + pad, ymin = bb[["ymin"]] - pad, ymax = bb[["ymax"]] + pad),
    crs = 4326
  ) |>
  st_as_sfc()

province_labels <-
  provinces |>
  filter(province %in% c("Haut-Uele", "Nord-Kivu", "Tshopo", "Bas-Uele")) |>
  anchor() |>
  mutate(x = if_else(province == "Haut-Uele", x - 0.5, x))

country_labels <- tribble(
  ~label, ~x, ~y,
  "Uganda", 31.75, 1.35,
  "South Sudan", 29.9, 4.2
)

ituri_centre <- ituri |> st_geometry() |> st_point_on_surface() |> st_coordinates()
angumu_centre <- angumu |> st_geometry() |> st_point_on_surface() |> st_coordinates()

map_context <-
  ggplot() +
  geom_sf(data = countries, fill = "grey90", colour = border_col, linewidth = 0.2) +
  geom_sf(data = provinces, fill = land_col, colour = border_col, linewidth = 0.2) +
  geom_sf(data = lakes, fill = lake_col, colour = NA) +
  geom_sf(data = ituri, fill = ituri_col, colour = ink, linewidth = 0.6) +
  geom_sf(data = health_areas, aes(fill = zone), colour = NA) +
  geom_sf(data = angumu, fill = NA, colour = ink, linewidth = 0.4) +
  scale_fill_manual(values = zone_cols, na.value = "grey80", guide = "none") +
  geom_sf(data = zoom_box, fill = NA, colour = ink, linewidth = 0.5, linetype = "dashed") +
  geom_text(data = province_labels, aes(x, y, label = province), colour = "grey45", fontface = "italic", size = 3) +
  geom_text(data = country_labels, aes(x, y, label = label), colour = "grey45", fontface = "italic", size = 3) +
  annotate(
    "segment",
    x = ituri_centre[1] - 0.1, y = ituri_centre[2] + 0.35, xend = 28.9, yend = 3.55, colour = ink, linewidth = 0.3
  ) +
  annotate("text", x = 28.9, y = 3.65, label = "Ituri", fontface = "bold", size = 4.5, colour = ink) +
  annotate(
    "segment",
    x = angumu_centre[1], y = angumu_centre[2], xend = 31.75, yend = 2.75, colour = ink, linewidth = 0.3
  ) +
  annotate("text", x = 31.75, y = 2.88, label = "Angumu", fontface = "bold", size = 4.5, colour = ink) +
  annotation_scale(location = "br", width_hint = 0.25, text_cex = 0.8, line_width = 0.5) +
  coord_sf(xlim = c(26.6, 32.2), ylim = c(0.1, 4.4), expand = FALSE) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", colour = NA))

# Inset: Ituri within the DRC.
map_drc <-
  ggplot() +
  geom_sf(data = drc, fill = land_col, colour = border_col, linewidth = 0.2) +
  geom_sf(data = ituri, fill = angumu_col, colour = ink, linewidth = 0.2) +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = ink, linewidth = 0.4),
    plot.margin = margin(4, 4, 4, 4)
  )

# Right panel: the 14 surveyed health areas by zone.
area_labels <-
  health_areas |>
  filter(!is.na(zone)) |>
  anchor() |>
  select(name_en, zone, x, y)

# Leader lines start near the shore inside each zone and end in the empty lake.
zone_labels <-
  tribble(
    ~zone, ~x, ~y, ~lx, ~ly,
    "IRS", 31.160, 2.052, 31.20, 2.00,
    "Routine", 31.092, 2.018, 31.17, 1.955,
    "MDA", 31.050, 1.990, 31.10, 1.915
  ) |>
  mutate(label = unname(zone_display_labels[zone]))

map_zoom <-
  ggplot() +
  geom_sf(data = ituri, fill = ituri_col, colour = NA) +
  geom_sf(data = lake_albert, fill = lake_col, colour = NA) +
  geom_sf(data = health_areas, fill = "grey88", colour = "white", linewidth = 0.2) +
  geom_sf(data = health_areas |> filter(!is.na(zone)), aes(fill = zone), colour = "white", linewidth = 0.25) +
  geom_sf(data = angumu, fill = NA, colour = ink, linewidth = 0.6) +
  scale_fill_manual(values = zone_cols, name = "2022 intervention zone") +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = ink), guide = "none") +
  geom_text(
    data = area_labels,
    aes(x, y, label = if_else(name_en == "Awasi Luga", "Awasi", name_en), colour = zone != "Routine"),
    size = 2.2, lineheight = 0.85
  ) +
  geom_segment(data = zone_labels, aes(x = x, y = y, xend = lx, yend = ly), colour = ink, linewidth = 0.3) +
  geom_text(data = zone_labels, aes(lx, ly, label = label), fontface = "bold", size = 4.5, colour = ink) +
  annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.8, line_width = 0.5) +
  coord_sf(
    xlim = c(bb[["xmin"]] - 0.02, bb[["xmax"]] + 0.05),
    ylim = c(bb[["ymin"]] - 0.02, bb[["ymax"]] + 0.04),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = ink, fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white", colour = NA)
  )

figure1 <-
  ggdraw() +
  draw_plot(map_context, x = 0, y = 0.02, width = 0.5, height = 0.96) +
  draw_plot(map_drc, x = 0.01, y = 0.66, width = 0.12, height = 0.3) +
  draw_plot(map_zoom, x = 0.52, y = 0.03, width = 0.47, height = 0.94)

ggsave("output/figures/figure1_study_area_map.png", figure1, width = 12, height = 4.8, dpi = 300, bg = "white")

# Figure 2 ----

# Inter is fetched from Google Fonts, one family per weight (a variable font cannot be
# reached with fontface under showtext).
font_add_google("Inter", "inter", regular.wt = 400)
font_add_google("Inter", "inter_md", regular.wt = 500)
font_add_google("Inter", "inter_sb", regular.wt = 600)
font_add_google("Inter", "inter_bold", regular.wt = 800)
showtext_auto()
showtext_opts(dpi = 300)

# Panel A: timeline. Campaign and survey dates come from the programme records and the
# protocol, not from the survey data.
surveys <-
  tribble(
    ~survey, ~label, ~start, ~end,
    1, "Survey 1 (baseline)", "2022-03-15", "2022-04-15",
    2, "Survey 2", "2022-11-01", "2022-12-15",
    3, "Survey 3", "2023-04-01", "2023-05-15",
    4, "Survey 4", "2023-10-01", "2023-10-31",
    5, "Survey 5", "2024-04-01", "2024-05-15"
  ) |>
  mutate(start = as.Date(start), end = as.Date(end), mid = start + (end - start) / 2)

interventions <-
  tribble(
    ~category, ~label, ~start, ~end, ~zone,
    "MDA", "Precursor (3 rounds)", "2020-09-01", "2021-01-31", "4 precursor HAs (2020-21)",
    "IRS", "Round 1", "2022-06-27", "2022-11-10", "IRS + MDA+IRS zones",
    "MDA", "Rounds 1-2 (Pyramax)", "2022-07-01", "2022-08-31", "MDA+IRS zone",
    "MDA", "Rounds 3-4 (DHA-piperaquine)", "2023-01-01", "2023-02-28", "MDA+IRS zone",
    "IRS", "Round 2", "2023-02-02", "2023-05-31", "IRS + MDA+IRS zones",
    "Crossover", "IRS (Abia, Cawa)", "2023-08-01", "2023-09-30", "Ex-Routine → IRS",
    "ITN", "National distribution", "2023-12-01", "2023-12-31", "All zones",
    "Crossover", "MDA (Cawa, Ara, Linga)", "2024-01-01", "2024-02-29", "Ex-Routine → MDA+IRS"
  ) |>
  mutate(start = as.Date(start), end = as.Date(end))

category_levels <- c("ITN", "Crossover", "IRS", "MDA", "Mortality recall", "Survey")
category_cols <- c(
  "MDA" = "#D55E00", "IRS" = "#0072B2", "Crossover" = "#8a8a8a", "ITN" = "#4a3aa7",
  "Survey" = "#3d3b38", "Mortality recall" = "#c9bfae"
)
ink_primary <- "#0b0b0b"
ink_secondary <- "#52514e"
ink_muted <- "#898781"
line_baseline <- "#c3c2b7"
band_fill <- "#f3f2ee"

row_h <- 1.25
lane_h <- 0.6
category_y <- set_names(seq_along(category_levels) * row_h, category_levels)
guide_ymin <- min(category_y) - lane_h / 2 - 0.05

# Mortality recall windows run from recall start to the end of each mortality survey's
# fieldwork; surveys 3 and 5 overlap, so survey 3 sits on its own sub-row.
recall_h <- 0.42
recall <-
  tribble(
    ~label, ~start, ~end, ~y_off,
    "Survey 1", "2021-01-01", "2022-04-15", recall_h / 2 + 0.03,
    "Survey 3", "2022-04-18", "2023-05-15", -(recall_h / 2 + 0.03),
    "Survey 5", "2023-04-09", "2024-05-15", recall_h / 2 + 0.03
  ) |>
  mutate(start = as.Date(start), end = as.Date(end))

interventions <-
  bind_rows(
    interventions,
    surveys |> transmute(category = "Survey", label, start, end),
    recall |> transmute(category = "Mortality recall", label, start, end, y_off)
  ) |>
  mutate(
    category = factor(category, levels = category_levels),
    y_off = replace_na(y_off, 0),
    y = category_y[as.character(category)] + y_off,
    mid = start + (end - start) / 2,
    id = row_number()
  ) |>
  arrange(category, start) |>
  mutate(.by = category, label_below = row_number() %% 2 == 1) |>
  mutate(bar_h = if_else(category == "Mortality recall", recall_h, lane_h))

rounded_bars <-
  interventions |>
  reframe(
    .by = c(id, category, label),
    x = c(start, end, end, start),
    y = c(y - bar_h / 2, y - bar_h / 2, y + bar_h / 2, y + bar_h / 2)
  )

survey_guides <- interventions |> filter(category == "Survey")

lane_labels <- tibble(
  category = factor(category_levels, levels = category_levels),
  y = category_y[category_levels],
  label = toupper(category_levels)
)

x_limits <- c(as.Date("2020-07-01"), as.Date("2024-08-01"))
y_limits <- c(guide_ymin - 0.55, max(category_y) + row_h / 2 + 0.55)

fig_timeline <-
  ggplot() +
  geom_rect(
    data = lane_labels |> filter((length(category_levels) - as.integer(category)) %% 2 == 0),
    aes(ymin = y - row_h / 2, ymax = y + row_h / 2),
    xmin = x_limits[1], xmax = x_limits[2], fill = band_fill
  ) +
  geom_segment(
    data = survey_guides, aes(x = mid, xend = mid, y = y - lane_h / 2, yend = guide_ymin),
    colour = line_baseline, linewidth = 0.35, linetype = "22"
  ) +
  geom_shape(
    data = rounded_bars, aes(x = x, y = y, group = id, fill = category),
    radius = unit(1.6, "mm"), colour = "white", linewidth = 0.6
  ) +
  geom_text(
    data = interventions |> filter(!label_below, category != "Mortality recall"),
    aes(x = mid, y = y + lane_h / 2 + 0.09, label = label),
    family = "inter_md", size = 3.05, colour = ink_secondary, vjust = 0
  ) +
  geom_text(
    data = interventions |> filter(label_below, category != "Mortality recall"),
    aes(x = mid, y = y - lane_h / 2 - 0.09, label = label),
    family = "inter_md", size = 3.05, colour = ink_secondary, vjust = 1
  ) +
  geom_text(
    data = interventions |> filter(category == "Mortality recall"),
    aes(x = start + 8, y = y, label = label),
    family = "inter_md", size = 3.05, colour = ink_primary, hjust = 0
  ) +
  geom_tile(data = lane_labels, aes(x = x_limits[1] - 8, y = y, fill = category), width = 5, height = lane_h * 0.7) +
  geom_text(
    data = lane_labels, aes(x = x_limits[1] - 22, y = y, label = label, colour = category),
    family = "inter_bold", size = 3.5, hjust = 1
  ) +
  scale_fill_manual(values = category_cols, guide = "none") +
  # The recall bar fill is too pale for text, so its lane label is a darker brown.
  scale_colour_manual(values = replace(category_cols, "Mortality recall", "#7a6c55"), guide = "none") +
  scale_x_date(
    date_breaks = "6 months", date_labels = "%b %Y",
    limits = c(x_limits[1] - 90, x_limits[2]), expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 12, base_family = "inter") +
  theme(
    axis.text.x = element_text(family = "inter_md", colour = ink_muted, size = 9.5, margin = margin(t = 5)),
    axis.line.x = element_line(colour = line_baseline, linewidth = 0.4),
    axis.ticks.x = element_line(colour = line_baseline, linewidth = 0.4),
    axis.ticks.length.x = unit(1.5, "mm"),
    plot.margin = margin(t = 10, r = 16, b = 6, l = 130)
  )

# Panel B: zone of each health area at each survey. Areas that crossed over have no zone
# from the survey they crossed over onward.
zone_grid <-
  df_ind |>
  filter(!is.na(area)) |>
  distinct(survey, area, zone) |>
  mutate(
    area = if_else(area == "Ugudo-Zii", "Ugudo Zii", area),
    zone = factor(
      if_else(is.na(zone), "Excluded (crossover)", as.character(zone)),
      levels = c("Routine", "MDA", "IRS", "Excluded (crossover)")
    ),
    round_label = factor(
      survey,
      levels = 1:5,
      labels = c(
        "Survey 1\n(Mar 2022)", "Survey 2\n(Nov 2022)", "Survey 3\n(Apr 2023)",
        "Survey 4\n(Oct 2023)", "Survey 5\n(Apr 2024)"
      )
    )
  )
stopifnot(!anyDuplicated(zone_grid[c("survey", "area")]), nrow(zone_grid) == 14 * 5)

# Areas are grouped by their survey-1 zone.
area_order <- zone_grid |> filter(survey == 1) |> arrange(zone, area) |> pull(area)
zone_grid <- zone_grid |> mutate(area = factor(area, levels = rev(area_order)))

grid_cols <- c(zone_cols[c("Routine", "MDA", "IRS")], "Excluded (crossover)" = "grey60")

fig_zone_grid <-
  zone_grid |>
  ggplot(aes(x = round_label, y = area, fill = zone)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  scale_fill_manual(
    "Zone (as treated)", values = grid_cols, drop = FALSE,
    breaks = names(grid_cols), labels = c("Routine", "MDA+IRS", "IRS", "Excluded (crossover)")
  ) +
  labs(x = NULL, y = "Health area") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold"),
    legend.position = "right"
  )

# wrap_elements() makes each half one leaf, so tag_levels gives one tag per panel.
figure2 <-
  wrap_elements(full = fig_timeline) / wrap_elements(full = fig_zone_grid) +
  plot_layout(heights = c(6.4, 5.9)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

# dpi must equal showtext_opts(dpi = 300) or the timeline text misaligns.
ggsave("output/figures/figure2_timeline_zone_grid.png", figure2, width = 9, height = 12.3, dpi = 300)

# S1 Table ----

# Survey 1, health areas in the zone comparison.
df_baseline <-
  df_ind |>
  filter(survey == 1, !is.na(zone), !is.na(cluster)) |>
  mutate(
    zone = zone_display(zone),
    sex = case_when(sex == "Femme" ~ "Female", sex == "Homme" ~ "Male"),
    rdt_positive = case_when(rdt_result == "Positif" ~ TRUE, rdt_result == "Negatif" ~ FALSE, .default = NA),
    died_yn = died,
    died_malaria = died_cause %in% "Malaria et/ou anemie"
  ) |>
  select(zone, area, sex, age, age_grp_few, status, hh_size, rdt_positive, died_yn, died_malaria, cluster, weight)

svy_baseline <-
  df_baseline |>
  select(-area) |>
  as_survey(ids = cluster, strata = zone, weights = weight)

tbl_s1 <-
  svy_baseline |>
  tbl_svysummary(
    by = zone,
    include = -c(cluster, weight),
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n_unweighted} ({p}%)"
    ),
    label = list(
      sex ~ "Sex",
      age ~ "Age (years)",
      age_grp_few ~ "Age group",
      status ~ "Displacement status",
      hh_size ~ "Household size",
      rdt_positive ~ "Malaria prevalence (RDT positive)",
      died_yn ~ "Deaths in recall period",
      died_malaria ~ "Malaria/anaemia deaths"
    ),
    missing = "ifany",
    missing_text = "Missing",
    missing_stat = "{N_miss_unweighted}"
  ) |>
  add_overall()

write_csv(tbl_s1$table_body, "output/tables/s1_table_baseline_characteristics.csv")

s1_areas_by_zone <-
  df_baseline |>
  distinct(zone, area) |>
  arrange(zone, area) |>
  summarise(areas = str_c(area, collapse = ", "), n = n(), .by = zone)
write_csv(s1_areas_by_zone, "output/tables/s1_table_areas_by_zone.csv")
print(s1_areas_by_zone)

# S4 Table ----

survey_dates <- tribble(
  ~survey, ~survey_label, ~dates,
  1L, "1 (baseline)", "Mar-Apr 2022",
  2L, "2", "Nov-Dec 2022",
  3L, "3", "Apr-May 2023",
  4L, "4", "Oct 2023",
  5L, "5", "Apr-May 2024"
)

# Households that declined consent are dropped at survey 1 only (26 households).
cluster_counts <-
  df_hh |>
  filter(survey != 1 | consent %in% "Oui", !is.na(zone), !is.na(cluster)) |>
  summarise(clusters = n_distinct(cluster), households = n(), .by = c(survey, zone)) |>
  left_join(survey_dates, by = "survey") |>
  left_join(survey_timing |> select(survey = survey_n, months_from_baseline) |> mutate(survey = as.integer(survey)), by = "survey") |>
  mutate(zone = zone_display(zone))

widen_counts <- function(counts, value) {
  counts |>
    select(Survey = survey_label, Dates = dates, `Months from baseline` = months_from_baseline, zone, all_of(value)) |>
    pivot_wider(names_from = zone, values_from = all_of(value), values_fill = 0) |>
    mutate(Total = Routine + `MDA+IRS` + IRS) |>
    relocate(Routine, `MDA+IRS`, IRS, .after = `Months from baseline`)
}

s4_clusters <- widen_counts(cluster_counts, "clusters")
s4_households <- widen_counts(cluster_counts, "households")
write_csv(s4_clusters, "output/tables/s4_table_clusters.csv")
write_csv(s4_households, "output/tables/s4_table_households.csv")
print(s4_clusters)
print(s4_households)

# S7 Table ----

participant_data <-
  survey_dates |>
  mutate(data_collected = c("RDT, mortality, CoD", "RDT only", "RDT, mortality, CoD", "RDT only", "RDT, mortality, CoD")) |>
  pmap(\(survey, survey_label, dates, data_collected) {
    s <- survey
    hh <- df_hh |> filter(survey == s)
    ind <- df_ind |> filter(survey == s)
    if (s == 1L) {
      keep_hh <- hh |> filter(consent %in% "Oui") |> pull(index_hh)
      hh <- hh |> filter(index_hh %in% keep_hh)
      ind <- ind |> filter(index_hh %in% keep_hh)
    }
    tibble(
      survey = survey_label,
      survey_n = s,
      dates = dates,
      households = nrow(hh),
      individuals = nrow(ind),
      female_pct = round(100 * mean(ind$sex %in% "Femme", na.rm = TRUE), 0),
      age_median = median(ind$age, na.rm = TRUE),
      age_q1 = quantile(ind$age, 0.25, na.rm = TRUE),
      age_q3 = quantile(ind$age, 0.75, na.rm = TRUE),
      u5_pct = round(100 * mean(ind$age_grp_few %in% "0-4", na.rm = TRUE), 1),
      displaced_pct = round(100 * mean(ind$status %in% "Deplace", na.rm = TRUE), 0),
      # `==` keeps untested individuals missing, so the denominator is those tested.
      rdt_positive_pct = round(100 * mean(ind$rdt_result == "Positif", na.rm = TRUE), 1),
      data_collected = data_collected
    )
  }) |>
  list_rbind() |>
  left_join(survey_timing |> select(survey_n, months_from_baseline) |> mutate(survey_n = as.integer(survey_n)), by = "survey_n") |>
  mutate(age_iqr = str_c(round(age_median), " (", round(age_q1), "-", round(age_q3), ")"))

s7_table <-
  participant_data |>
  select(
    Survey = survey, `Months from baseline` = months_from_baseline, Dates = dates, Households = households,
    Individuals = individuals, `Female (%)` = female_pct, `Median age, yrs (IQR)` = age_iqr,
    `Children <5 (%)` = u5_pct, `Displaced (%)` = displaced_pct, `Malaria prevalence (%)` = rdt_positive_pct,
    `Data collected` = data_collected
  )
write_csv(s7_table, "output/tables/s7_table_participant_characteristics.csv")
print(s7_table, width = 200)

# S8 Table ----

coverage_row <- function(label, round, prop_tbl) {
  prop_tbl |>
    transmute(intervention = label, round = round, zone, pct = round(prop, 1), lcl = round(prop_low, 1), ucl = round(prop_upp, 1))
}

received <- c("card", "recall")

# MDA rounds 1 and 2 and IRS round 1 were asked at survey 2.
s2_pop <- df_ind |> filter(survey == 2, !is.na(zone), !is.na(cluster)) |> mutate(zone = as.character(zone))

s2_mda <- function(var) {
  s2_pop |>
    filter(zone == "MDA") |>
    mutate(y = as.numeric(.data[[var]] %in% received)) |>
    prop_logit() |>
    mutate(zone = "MDA")
}
tbl_mda_1 <- s2_mda("mda_r1")
tbl_mda_2 <- s2_mda("mda_r2")

# IRS round 1 is a household report: one row per household with an answer, weighted by
# the person weight times the people listed, which removes the household-size factor.
s2_hh <-
  s2_pop |>
  filter(zone %in% c("IRS", "MDA")) |>
  summarise(irs_hh = first(irs_hh), cluster = first(cluster), zone = first(zone), weight = first(weight) * n(), .by = index_hh) |>
  filter(!is.na(irs_hh)) |>
  mutate(y = as.numeric(irs_hh == "yes"))
tbl_irs_1 <- c("IRS", "MDA") |>
  map_dfr(\(z) prop_logit(s2_hh |> filter(zone == z)) |> mutate(zone = z))

# MDA rounds 3 and 4 and IRS round 2 were asked at survey 3.
s3_mda <- df_ind |> filter(survey == 3, zone %in% "MDA")
tbl_mda_3 <- prop_logit(s3_mda |> mutate(y = mda_r3 %in% received)) |> mutate(zone = "MDA")
tbl_mda_4 <- prop_logit(s3_mda |> mutate(y = mda_r4 %in% received)) |> mutate(zone = "MDA")

# Apala and Uyandu are excluded from the MDA+IRS denominator of IRS round 2.
s3_irs <- df_ind |> filter(survey == 3, zone %in% c("IRS", "MDA"), !area %in% c("Apala", "Uyandu"))
tbl_irs_2 <- c("IRS", "MDA") |>
  map_dfr(\(z) prop_logit(s3_irs |> filter(zone == z) |> mutate(y = irs_hh %in% "yes")) |> mutate(zone = z))

zone_name <- \(z) if_else(z == "MDA", "MDA+IRS", z)

s8_table <-
  bind_rows(
    coverage_row("MDA (Pyramax)", "Round 1 (Jul 2022)", tbl_mda_1),
    coverage_row("MDA (Pyramax)", "Round 2 (Aug 2022)", tbl_mda_2),
    coverage_row("IRS (SumiShield 50WG)", "Round 1 (Jul-Aug 2022)", tbl_irs_1),
    coverage_row("MDA (DHA-piperaquine)", "Round 3 (2023)", tbl_mda_3),
    coverage_row("MDA (DHA-piperaquine)", "Round 4 (2023)", tbl_mda_4),
    coverage_row("IRS (SumiShield 50WG)", "Round 2 (Apr 2023)", tbl_irs_2)
  ) |>
  mutate(zone = zone_name(zone), ci = str_c(pct, " (", lcl, "-", ucl, ")")) |>
  select(intervention, round, zone, ci) |>
  mutate(intervention = factor(intervention, levels = c("MDA (Pyramax)", "MDA (DHA-piperaquine)", "IRS (SumiShield 50WG)"))) |>
  arrange(intervention, round, zone) |>
  mutate(intervention = as.character(intervention)) |>
  rename(Zone = zone, `% (95% CI)` = ci)
write_csv(s8_table, "output/tables/s8_table_coverage.csv")
print(s8_table, n = 20)

# Households without an answer to the IRS round 1 item, for the table note.
s8_notes <-
  s2_pop |>
  filter(zone %in% c("IRS", "MDA")) |>
  summarise(irs_hh = first(irs_hh), area = first(area), zone = first(zone), .by = index_hh) |>
  summarise(
    n_households = n(), n_irs_missing = sum(is.na(irs_hh)),
    n_households_gengere = sum(area %in% "Gengere"), n_irs_missing_gengere = sum(is.na(irs_hh) & area %in% "Gengere"),
    .by = zone
  ) |>
  mutate(zone = zone_name(zone)) |>
  arrange(desc(zone))
write_csv(s8_notes, "output/tables/s8_table_coverage_notes.csv")
print(s8_notes)
