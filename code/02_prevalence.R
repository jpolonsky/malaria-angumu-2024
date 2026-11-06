# Malaria prevalence by zone: the pooled baseline-adjusted model and its within-survey contrasts.
# Reproduces Table 1 (prevalence by zone), Figure 3 (prevalence trend), Figure 4 (prevalence maps),
# S4 Figure and S9 Table (prevalence by health area) and S6 Table (Gengere).

pacman::p_load(lme4, emmeans, srvyr, survey, sf, spatstat, mgcv, nngeo, stars, ggnewscale, ggrepel, patchwork, tidyverse)

source("code/functions.R")

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/models", recursive = TRUE, showWarnings = FALSE)

# Analysis data ----

# Valid RDT results with a zone (crossover-excluded areas carry no zone)
dat_rdt <-
  read_rds("data/df_ind.rds") |>
  transmute(survey, cluster, index_hh = as.character(index_hh), area, zone, age_years = age, sex, status, weight, rdt_result) |>
  filter(rdt_result %in% c("Positif", "Negatif"), !is.na(zone)) |>
  mutate(
    survey = factor(survey),
    outcome_rdt = as.integer(rdt_result == "Positif"),
    sex = factor(sex, levels = c("Homme", "Femme")),
    status = factor(status, levels = c("Autochtone", "Deplace")),
    age_years = if_else(age_years > 100, 100, age_years),
    cluster_uid = paste(survey, cluster, sep = "_")
  )

# Area-level baseline prevalence (survey 1), winsorised at the 95th percentile and standardised
baseline_prev_area <-
  dat_rdt |>
  filter(survey == "1") |>
  summarise(baseline_prev = mean(outcome_rdt), .by = area)

overall_s1_prev <- mean(baseline_prev_area$baseline_prev)
p95_prev <- quantile(baseline_prev_area$baseline_prev, 0.95)

dat_rdt <-
  dat_rdt |>
  left_join(baseline_prev_area, by = "area") |>
  mutate(
    baseline_prev = replace_na(baseline_prev, overall_s1_prev),
    baseline_prev_w = pmin(baseline_prev, p95_prev),
    baseline_prev_z = as.numeric(scale(baseline_prev_w))
  )

# Pooled model ----

# Baseline-adjusted model (primary) and the same model without the baseline covariate,
# which supplies the survey 1 row of Table 1
glmer_control <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

m4_durability <- glmer(
  outcome_rdt ~ zone * survey + age_years + sex + status + (1 | cluster_uid/index_hh),
  data = dat_rdt, family = binomial, control = glmer_control
)

m4_durability_adj <- glmer(
  outcome_rdt ~ zone * survey + age_years + sex + status + baseline_prev_z + (1 | cluster_uid/index_hh),
  data = dat_rdt, family = binomial, control = glmer_control
)

# Fits and analysis data for 03_prevalence_sensitivity.R
saveRDS(list(dat_rdt = dat_rdt, m4_durability = m4_durability, m4_durability_adj = m4_durability_adj), "output/models/prevalence_model.rds")

# Within-survey contrasts against the concurrent Routine zone
set.seed(1)

emm_prev <- emmeans(m4_durability_adj, ~ zone | survey, type = "response")
contrasts_prev <- tidy_contrasts(emm_prev) |> mutate(survey_n = as.integer(as.character(survey)))
contrasts_prev_unadj <- tidy_contrasts(emmeans(m4_durability, ~ zone | survey, type = "response")) |>
  mutate(survey_n = as.integer(as.character(survey)))

# Table 1 ----

# Survey-weighted crude prevalence and crude odds ratios, one design per survey
build_prev_design <- \(df) as_survey(df, ids = cluster_uid, strata = zone, weights = weight)

dat_by_survey <-
  dat_rdt |>
  mutate(survey_n = as.integer(as.character(survey))) |>
  group_split(survey_n)

crude_prev <- prevalence_by_zone(dat_rdt)

crude_or <-
  dat_by_survey |>
  map(\(df) {
    svyglm(outcome_rdt ~ zone, design = build_prev_design(df), family = quasibinomial()) |>
      broom::tidy(exponentiate = TRUE, conf.int = TRUE) |>
      filter(term != "(Intercept)") |>
      transmute(
        zone = str_remove(term, "^zone"),
        crude_or = estimate, crude_or_lcl = conf.low, crude_or_ucl = conf.high,
        survey_n = unique(df$survey_n)
      )
  }) |>
  list_rbind()

# Adjusted prevalence difference (percentage points) at the covariate means, random effects at zero
prev_diff_pp <-
  contrast(regrid(emm_prev), method = "trt.vs.ctrl", ref = "Routine", adjust = "mvt") |>
  summary(infer = TRUE) |>
  as_tibble() |>
  transmute(
    zone = str_remove(as.character(contrast), " - Routine"),
    survey_n = as.integer(as.character(survey)),
    pd_ci = sprintf("%.1f (%.1f to %.1f)", 100 * estimate, 100 * asymp.LCL, 100 * asymp.UCL)
  )

# Ratio of odds ratios versus survey 1 (the zone:survey interaction)
ror_table <-
  contrast(
    emmeans(m4_durability_adj, ~ zone * survey),
    interaction = c(zone = "trt.vs.ctrl", survey = "trt.vs.ctrl"), adjust = "none"
  ) |>
  summary(infer = TRUE) |>
  as_tibble() |>
  transmute(
    zone = str_remove(as.character(zone_trt.vs.ctrl), " - Routine"),
    survey_n = as.integer(str_extract(as.character(survey_trt.vs.ctrl), "^\\d")),
    ror_ci = fmt_ci(exp(estimate), exp(asymp.LCL), exp(asymp.UCL))
  )

table1 <-
  crude_prev |>
  mutate(crude_ci = sprintf("%.1f (%.1f-%.1f)", prev * 100, lcl * 100, ucl * 100)) |>
  select(zone, survey_n, n_total, crude_ci) |>
  left_join(
    bind_rows(
      crude_or |> mutate(or_ci = fmt_ci(crude_or, crude_or_lcl, crude_or_ucl)) |> select(zone, survey_n, or_ci),
      tibble(zone = "Routine", survey_n = 1:5, or_ci = "Reference")
    ),
    by = c("zone", "survey_n")
  ) |>
  # Survey 1 comes from the model without the baseline covariate, which is derived from survey 1 itself
  left_join(
    bind_rows(
      contrasts_prev_unadj |> filter(survey_n == 1),
      contrasts_prev |> filter(survey_n != 1)
    ) |>
      mutate(adj_or_ci = fmt_ci(or, lcl, ucl)) |>
      select(zone, survey_n, adj_or_ci, p.value) |>
      bind_rows(tibble(zone = "Routine", survey_n = 1:5, adj_or_ci = "Reference", p.value = NA_real_)),
    by = c("zone", "survey_n")
  ) |>
  left_join(bind_rows(prev_diff_pp, tibble(zone = "Routine", survey_n = 1:5, pd_ci = "Reference")), by = c("zone", "survey_n")) |>
  left_join(bind_rows(ror_table, tibble(zone = "Routine", survey_n = 2:5, ror_ci = "Reference")), by = c("zone", "survey_n")) |>
  left_join(season_lu, by = "survey_n") |>
  mutate(
    pd_ci = if_else(survey_n == 1 & zone != "Routine", "---", pd_ci),
    ror_ci = replace_na(ror_ci, "---"),
    zone = zone_display(zone),
    p_fmt = case_when(survey_n == 1 ~ "---", zone == "Routine" ~ "Reference", p.value < 0.001 ~ "<0.001", .default = sprintf("%.3f", p.value))
  ) |>
  arrange(survey_n, zone) |>
  select(
    Zone = zone, Survey = survey_n, `Months from baseline` = months_from_baseline, Season = season,
    `N (participants)` = n_total,
    `Crude malaria prevalence (95% CI)` = crude_ci, `Crude OR (95% CI)` = or_ci,
    `OR (95% CI)` = adj_or_ci, `p-value` = p_fmt,
    `Prevalence difference, pp (95% CI)` = pd_ci, `DiD ratio of ORs (95% CI)` = ror_ci
  )

print(table1, n = Inf, width = Inf)
write_csv(table1, "output/tables/table1_prevalence.csv")

# Figure 3 ----

zone_facet_labels <- c(Routine = "Routine", IRS = "IRS", MDA = "MDA+IRS")

figure3 <-
  crude_prev |>
  left_join(season_lu, by = "survey_n") |>
  mutate(zone = factor(zone, levels = zone_display_levels)) |>
  ggplot(aes(x = months_from_baseline, y = prev * 100, colour = zone, fill = zone)) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = lcl * 100, ymax = ucl * 100), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = season), size = 3) +
  scale_shape_manual(values = c("Low transmission" = 16, "High transmission" = 4)) +
  scale_colour_manual(values = zone_cols) +
  scale_fill_manual(values = zone_cols) +
  scale_x_continuous(breaks = survey_timing$months_from_baseline) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~zone, nrow = 1, labeller = labeller(zone = zone_facet_labels)) +
  labs(x = "Months from baseline", y = "Malaria prevalence (%)", shape = "Preceding season") +
  guides(colour = "none", fill = "none") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"), legend.position = "bottom")

ggsave("output/figures/figure3_prevalence_trend.png", figure3, width = 9, height = 6, dpi = 300)

# Figure 4 ----

areas <- read_rds("data/maps/health_areas.rds")

# Land and Lake Albert underneath each panel
ext_bbox <-
  st_bbox(areas) |>
  (\(b) st_bbox(c(xmin = b[["xmin"]] - 4000, xmax = b[["xmax"]] + 4000, ymin = b[["ymin"]] - 3000, ymax = b[["ymax"]] + 3000), crs = st_crs(areas)))()

lake_albert <-
  read_rds("data/maps/lake_albert.rds") |>
  st_transform(st_crs(areas)) |>
  st_crop(ext_bbox)

lake_label <-
  lake_albert |>
  st_cast("POLYGON", warn = FALSE) |>
  mutate(area = st_area(geometry)) |>
  slice_max(area, n = 1) |>
  st_centroid() |>
  st_coordinates()

# Areas that crossed over to another intervention are outlined in grey from the survey they crossed
crossed_over <- map(1:5, \(s) study_areas$map_name[!is.na(study_areas$crossover_survey) & study_areas$crossover_survey <= s]) |> set_names(1:5)
area_zone <- set_names(study_areas$baseline_zone, study_areas$map_name)

union_no_holes <- \(x) x |> st_union() |> nngeo::st_remove_holes()

# Outline of the accessible areas (the smoothing window) and zone outlines for one survey
zone_overlay <- function(survey) {
  shp <-
    areas |>
    mutate(
      zone = case_when(
        name_en %in% crossed_over[[survey]] ~ "Excluded (crossover)",
        name_en %in% names(area_zone) ~ unname(area_zone[name_en]),
        .default = "Inaccessible"
      ) |>
        factor(levels = c("Routine", "IRS", "MDA", "Excluded (crossover)", "Inaccessible"))
    )
  outlines <-
    levels(shp$zone) |>
    set_names() |>
    keep(\(z) any(shp$zone == z)) |>
    map(\(z) shp |> filter(zone == z) |> union_no_holes())
  list(
    window = shp |> filter(zone != "Inaccessible") |> union_no_holes(),
    layers = c(
      list(geom_sf(data = outlines[["Inaccessible"]], fill = "lightgrey", colour = "black")),
      outlines[names(outlines) != "Inaccessible"] |>
        imap(\(x, z) geom_sf(
          data = x, fill = NA,
          colour = if (z == "Excluded (crossover)") "grey60" else "black",
          linewidth = if (z == "Excluded (crossover)") 0.6 else 0.4
        ))
    )
  )
}

# Cluster-level prevalence, shrunk towards the survey mean with a binomial random-intercept model
cluster_prevalence <- function(dat) {
  by_cluster <-
    dat |>
    group_by(cluster) |>
    summarise(num_ppl = n(), prev = 100 * sum(case) / n())
  fit <- gam(case ~ s(cluster_fct, bs = "re"), data = dat |> mutate(case = as.numeric(case), cluster_fct = factor(cluster)), family = binomial())
  shrunk <-
    dat |>
    distinct(cluster) |>
    mutate(cluster_fct = factor(cluster), prev = 100 * predict(fit, newdata = tibble(cluster_fct), type = "response")) |>
    select(cluster, prev)
  by_cluster |>
    select(-prev) |>
    left_join(shrunk, by = "cluster") |>
    left_join(
      dat |> group_by(cluster) |> filter(!is.na(lon), !is.na(lat)) |> slice(1) |> select(cluster, lon, lat),
      by = "cluster"
    ) |>
    filter(!is.na(lon), !is.na(lat)) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
    st_transform(3857)
}

# Kernel-smoothed surface, weighted by cluster size, with the bandwidth chosen per survey
smooth_prevalence <- function(clusters, window) {
  win <- owin(poly = list(x = rev(st_coordinates(window)[, 1]), y = rev(st_coordinates(window)[, 2])))
  ppp(x = st_coordinates(clusters)[, 1], y = st_coordinates(clusters)[, 2], marks = clusters$prev, window = win) |>
    Smooth(weights = clusters$num_ppl) |>
    stars::st_as_stars() |>
    st_as_sf() |>
    st_set_crs(3857)
}

fill_limits <- c(0, 85)

df_ind <- read_rds("data/df_ind.rds")

map_panels <-
  1:5 |>
  map(\(s) {
    overlay <- zone_overlay(as.character(s))
    clusters <-
      df_ind |>
      filter(survey == s, rdt_result %in% c("Positif", "Negatif")) |>
      mutate(case = rdt_result == "Positif") |>
      cluster_prevalence()
    surface <- smooth_prevalence(clusters, overlay$window)
    print(c(survey = s, n_clusters = nrow(clusters), max_surface = max(surface$v)))
    ggplot() +
      geom_rect(aes(xmin = ext_bbox[["xmin"]], xmax = ext_bbox[["xmax"]], ymin = ext_bbox[["ymin"]], ymax = ext_bbox[["ymax"]], fill = "Land")) +
      geom_sf(data = lake_albert, aes(fill = "Sea (Lake Albert)"), colour = "#4a90c4", linewidth = 0.3) +
      annotate("text", x = lake_label[, "X"], y = lake_label[, "Y"], label = "Lake Albert", colour = "#2a5d8a", fontface = "italic", size = 3.5) +
      scale_fill_manual(name = "Area type", values = c("Sea (Lake Albert)" = "#a6cee3", "Land" = "#f2ede3"), limits = c("Sea (Lake Albert)", "Land")) +
      ggnewscale::new_scale_fill() +
      geom_sf(data = surface, aes(fill = v), colour = NA) +
      geom_sf(data = st_boundary(overlay$window), linewidth = 1) +
      scale_fill_viridis_c(option = "plasma", direction = -1, name = "Percentage (%)", limits = fill_limits) +
      overlay$layers +
      coord_sf(xlim = ext_bbox[c("xmin", "xmax")], ylim = ext_bbox[c("ymin", "ymax")], crs = st_crs(areas), expand = FALSE) +
      theme_void()
  })

figure4 <-
  wrap_plots(c(map_panels, list(guide_area())), ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right")

ggsave("output/figures/figure4_prevalence_maps.png", figure4, width = 15, height = 9, dpi = 300)

# S6 Table ----

# Gengere supplies most MDA+IRS mosquitoes: prevalence in Gengere and in the other four MDA+IRS areas
dat_prev <-
  dat_rdt |>
  mutate(
    survey_n = as.integer(as.character(survey)),
    group = case_when(zone == "MDA" & area == "Gengere" ~ "Gengere", zone == "MDA" ~ "Other four MDA+IRS areas", .default = NA_character_)
  )

# Zone rows as in Table 1 (one design per survey, grouped by zone); area groups each in their own design
zone_est <-
  dat_prev |>
  group_split(survey_n) |>
  map(\(d) {
    d |>
      as_survey(ids = cluster, strata = zone, weights = weight) |>
      group_by(zone) |>
      summarise(prop = survey_mean(outcome_rdt * 100, vartype = "ci", deff = TRUE)) |>
      mutate(survey = d$survey_n[1], group = as.character(zone_display(zone)), row_type = "zone (context)")
  }) |>
  list_rbind() |>
  select(survey, row_type, group, prop, prop_low, prop_upp, prop_deff)

area_est <-
  dat_prev |>
  filter(zone == "MDA") |>
  group_split(survey_n, group) |>
  map(\(d) {
    est <- calculate_prop(d |> mutate(zone = droplevels(zone)), key = rdt_result, value = "Positif")
    tibble(survey = d$survey_n[1], row_type = "MDA+IRS area group", group = d$group[1], prop = est$prop, prop_low = est$prop_low, prop_upp = est$prop_upp, prop_deff = est$prop_deff)
  }) |>
  list_rbind()

dat_groups <- bind_rows(
  dat_prev |> filter(zone == "MDA") |> mutate(row_type = "MDA+IRS area group"),
  dat_prev |> mutate(group = as.character(zone_display(zone)), row_type = "zone (context)")
)

s6_table <-
  dat_groups |>
  mutate(survey = survey_n) |>
  summarise(
    areas = str_c(sort(unique(area)), collapse = ", "), n_clusters = n_distinct(cluster), n_tested = n(), n_positive = sum(outcome_rdt),
    .by = c(survey, row_type, group)
  ) |>
  left_join(bind_rows(zone_est, area_est), by = c("survey", "row_type", "group")) |>
  mutate(unweighted_pct = 100 * n_positive / n_tested) |>
  rename(prevalence_pct = prop, lcl_pct = prop_low, ucl_pct = prop_upp, deff = prop_deff) |>
  left_join(survey_timing |> select(survey = survey_n, months_from_baseline), by = "survey") |>
  mutate(group = factor(group, levels = c("Gengere", "Other four MDA+IRS areas", "Routine", "MDA+IRS", "IRS"))) |>
  arrange(survey, match(row_type, c("MDA+IRS area group", "zone (context)")), group) |>
  mutate(group = as.character(group)) |>
  select(survey, months_from_baseline, row_type, group, areas, n_clusters, n_tested, n_positive, prevalence_pct, lcl_pct, ucl_pct, unweighted_pct, deff)

print(s6_table, n = Inf, width = Inf)
write_csv(s6_table, "output/tables/s6_table_gengere_prevalence.csv")

# S9 Table and S4 Figure ----

area_lookup <- study_areas |> select(area, baseline_zone, crossover_survey)

# Area estimates come from the same design as Table 1, grouped by area; a single-cluster area gets no interval
area_prev <-
  dat_rdt |>
  mutate(survey_n = as.integer(as.character(survey))) |>
  group_split(survey_n) |>
  map(\(df) {
    clusters_per_area <- df |> summarise(n_clusters = n_distinct(cluster_uid), .by = area)
    df |>
      as_survey(ids = cluster_uid, strata = zone, weights = weight) |>
      group_by(area) |>
      summarise(
        n_tested = unweighted(n()), n_positive = unweighted(sum(outcome_rdt)),
        prop = survey_mean(outcome_rdt, vartype = "ci")
      ) |>
      left_join(clusters_per_area, by = "area") |>
      mutate(
        survey_n = df$survey_n[1],
        prevalence_pct = 100 * prop,
        ci_ok = n_clusters > 1 & !is.na(prop_low),
        lcl_pct = if_else(ci_ok, 100 * prop_low, NA_real_),
        ucl_pct = if_else(ci_ok, 100 * prop_upp, NA_real_)
      ) |>
      select(area, survey_n, n_clusters, n_tested, n_positive, prevalence_pct, lcl_pct, ucl_pct)
  }) |>
  list_rbind()

# Every area at every survey, so that excluded areas keep a row
s9_long <-
  expand_grid(area = area_lookup$area, survey_n = 1:5) |>
  left_join(area_lookup, by = "area") |>
  left_join(area_prev, by = c("area", "survey_n")) |>
  left_join(survey_timing |> select(survey_n, months_from_baseline), by = "survey_n") |>
  mutate(
    status = case_when(
      !is.na(n_tested) ~ "in comparison",
      !is.na(crossover_survey) & survey_n >= crossover_survey ~ "excluded after crossover",
      .default = "not sampled"
    ),
    unweighted_pct = 100 * n_positive / n_tested,
    baseline_zone_display = zone_display(baseline_zone)
  ) |>
  arrange(baseline_zone_display, area, survey_n)

s9_table <-
  s9_long |>
  transmute(
    baseline_zone = as.character(baseline_zone_display), area, survey = survey_n, months_from_baseline, status,
    n_clusters, n_tested, n_positive, unweighted_pct, prevalence_pct, lcl_pct, ucl_pct
  )

s9_wide <-
  s9_long |>
  mutate(
    cell = case_when(
      status == "excluded after crossover" ~ "excluded",
      status == "not sampled" ~ "not sampled",
      is.na(lcl_pct) ~ str_c(sprintf("%.1f", prevalence_pct), " (CI not estimable)"),
      .default = str_c(sprintf("%.1f", prevalence_pct), " (", sprintf("%.1f", lcl_pct), " to ", sprintf("%.1f", ucl_pct), ")")
    ),
    survey_col = str_c("survey_", survey_n)
  ) |>
  select(baseline_zone = baseline_zone_display, area, survey_col, cell) |>
  pivot_wider(names_from = survey_col, values_from = cell) |>
  mutate(baseline_zone = as.character(baseline_zone))

print(s9_wide, n = Inf, width = Inf)
write_csv(s9_table, "output/tables/s9_table_prevalence_by_area.csv", na = "")
write_csv(s9_wide, "output/tables/s9_table_prevalence_by_area_wide.csv", na = "")

# A line stops at the survey before the area's crossover
fig_dat <- s9_long |> filter(status == "in comparison")
fig_labels <- fig_dat |> slice_max(survey_n, n = 1, by = area, with_ties = FALSE)

s4_figure <-
  fig_dat |>
  ggplot(aes(x = months_from_baseline, y = prevalence_pct, group = area, colour = baseline_zone)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1) +
  geom_text_repel(
    data = fig_labels, aes(label = area),
    direction = "y", hjust = 0, nudge_x = 1, size = 2.8, segment.size = 0.2, min.segment.length = 0.1,
    box.padding = 0.15, seed = 1, show.legend = FALSE
  ) +
  facet_wrap(vars(baseline_zone_display), nrow = 1) +
  scale_colour_manual(values = zone_cols, guide = "none") +
  scale_x_continuous(breaks = survey_timing$months_from_baseline, expand = expansion(mult = c(0.03, 0.4))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Months from baseline", y = "Malaria prevalence (%), survey-weighted") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(), strip.text = element_text(face = "bold"))

ggsave("output/figures/s4_figure_prevalence_by_area.png", s4_figure, width = 9, height = 5, dpi = 300)
