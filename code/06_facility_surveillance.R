# Routine facility surveillance: Figure 6 (facility TPR against survey malaria prevalence), S1 Figure
# (monthly facility TPR and survey prevalence by zone), S11 Table (monthly TPR by zone) and the
# concordance numbers quoted in the Results and S1 Text.

pacman::p_load(srvyr, survey, ggforce, patchwork, scales, tidyverse)

source("code/functions.R")

dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)

# Helpers ----

# Kendall's coefficient of concordance for k rankings (rows) of n items (columns)
kendall_w <- function(mat) {
  k <- nrow(mat)
  n <- ncol(mat)
  ranks <- t(apply(mat, 1, rank))
  r <- colSums(ranks)
  s <- sum((r - mean(r))^2)
  w <- 12 * s / (k^2 * (n^3 - n))
  chi2 <- k * (n - 1) * w
  list(w = w, chi2 = chi2, df = n - 1, p = pchisq(chi2, df = n - 1, lower.tail = FALSE))
}

# Lowest zone of one survey as its display label; a tie for the lowest value counts as no agreement
lowest_zone_label <- function(zone, value) {
  is_min <- value == min(value)
  if (sum(is_min) == 1) unname(zone_display_labels[as.character(zone[is_min])]) else NA_character_
}

# Per-survey agreement between survey prevalence and TPR on the lowest zone and the full ranking
agreement_by_survey <- function(df) {
  df |>
    filter(!is.na(prev_pct), !is.na(tpr_pct)) |>
    filter(n() == 3, .by = survey) |>
    summarise(
      lowest_zone_tpr = lowest_zone_label(zone, tpr_pct),
      lowest_zone_survey = lowest_zone_label(zone, prev_pct),
      agree_lowest = !is.na(lowest_zone_tpr) & !is.na(lowest_zone_survey) & lowest_zone_tpr == lowest_zone_survey,
      agree_full_ranking = all(rank(tpr_pct) == rank(prev_pct)),
      .by = survey
    ) |>
    arrange(survey)
}

# Pairs each zone-by-survey prevalence with the zone TPR in the survey's calendar month
pair_prev_tpr <- function(tpr_df) {
  survey_prev |>
    left_join(tpr_df |> select(zone, date, tpr), by = c("zone", "date")) |>
    mutate(prev_pct = prev * 100, tpr_pct = tpr * 100)
}

# Kendall's W over the survey-prevalence and TPR rankings of the three zones (ten rankings
# when all five surveys are used), plus the two agreement counts
concordance_stats <- function(df) {
  rank_mat <-
    bind_rows(
      df |> select(survey, zone, value = prev_pct),
      df |> select(survey, zone, value = tpr_pct)
    ) |>
    mutate(stream = rep(c("prev", "tpr"), each = nrow(df))) |>
    pivot_wider(id_cols = c(stream, survey), names_from = zone, values_from = value) |>
    select(all_of(zone_levels)) |>
    as.matrix()
  kw <- kendall_w(rank_mat[complete.cases(rank_mat), ])
  agree <- agreement_by_survey(df)
  tibble(
    n_obs = nrow(df |> filter(!is.na(prev_pct), !is.na(tpr_pct))),
    kendall_w = kw$w, kendall_chi2 = kw$chi2, kendall_df = kw$df, kendall_p = kw$p,
    agree_lowest_n = sum(agree$agree_lowest),
    agree_full_n = sum(agree$agree_full_ranking),
    n_surveys = nrow(agree)
  )
}

# Data ----

# Survey-weighted malaria prevalence by zone and survey
survey_prev <-
  read_rds("data/df_ind.rds") |>
  filter(rdt_result %in% c("Positif", "Negatif"), !is.na(zone)) |>
  mutate(outcome_rdt = as.integer(rdt_result == "Positif"), cluster_uid = paste(survey, cluster, sep = "_")) |>
  prevalence_by_zone() |>
  rename(survey = survey_n)

# Surveys are placed in their calendar month on the monthly series
survey_dates <- tibble(
  survey = 1:5,
  date = as_date(c("2022-04-01", "2022-11-01", "2023-04-01", "2023-10-01", "2024-04-01")),
  label = str_c("Survey ", 1:5)
)

survey_prev <-
  survey_prev |>
  left_join(survey_dates, by = "survey") |>
  mutate(zone = factor(zone, levels = zone_display_levels))

# PNLP monthly health-area indicators, with the Angumu zone-referral row dropped
df_pnlp <-
  read_rds("data/pnlp_health_area_monthly.rds") |>
  filter(area != "Angumu") |>
  mutate(date = ym(period), zone = zone_at_date(area, date))

# Zone TPR: mean of the health-area monthly rates
df_tdr <-
  df_pnlp |>
  filter(indicateur_name == "Positivite Tdr", !is.na(zone)) |>
  summarise(taux_moyen = mean(taux_moyen, na.rm = TRUE), .by = c(date, zone))

# Zone case counts: per-health-area mean, so that the shrinking Routine zone is not read as a fall
df_cascade <-
  df_pnlp |>
  filter(indicateur_name %in% c("Cas Suspects Teste", "Cas Confirme Traité"), !is.na(zone)) |>
  summarise(total_absolu = sum(total_absolu, na.rm = TRUE) / n_distinct(area), .by = c(date, zone, indicateur_name))

df_zone_monthly <-
  df_cascade |>
  pivot_wider(names_from = indicateur_name, values_from = total_absolu) |>
  rename(cas_suspect = `Cas Suspects Teste`, tdr_positif = `Cas Confirme Traité`) |>
  left_join(df_tdr |> rename(tpr = taux_moyen), by = c("zone", "date"))

# DHIS2 sensitivity source: zone TPR is positive RDTs over suspected cases, summed over areas
df_zone_monthly_dhis2 <-
  read_rds("data/dhis2_monthly.rds") |>
  mutate(zone = zone_at_date(area, date)) |>
  filter(!is.na(zone)) |>
  summarise(value = sum(value, na.rm = TRUE), .by = c(zone, date, metric)) |>
  pivot_wider(names_from = metric, values_from = value) |>
  mutate(tpr = if_else(cas_suspect > 0, tdr_positif / cas_suspect, NA_real_))

start_date <- as_date("2021-10-01")
end_date <- as_date("2024-10-01")

# S1 Figure ----

s1_fig_tpr <-
  df_zone_monthly |>
  filter(date >= start_date, date <= end_date) |>
  ggplot(aes(x = date, colour = zone)) +
  geom_vline(xintercept = survey_dates$date, linetype = "dotted", colour = "grey60", linewidth = 0.4) +
  annotate("text", x = survey_dates$date, y = Inf, label = survey_dates$label, vjust = 1.4, size = 3, colour = "grey40") +
  geom_line(aes(y = tpr * 100), linewidth = 0.5, alpha = 0.3) +
  geom_smooth(aes(y = tpr * 100), method = "loess", span = 0.25, se = FALSE, linewidth = 1.1) +
  geom_pointrange(
    data = survey_prev |> filter(date >= start_date, date <= end_date),
    aes(y = prev * 100, ymin = lcl * 100, ymax = ucl * 100),
    size = 0.9, linewidth = 1.0, shape = 21, fill = "white"
  ) +
  scale_colour_manual(values = zone_cols, breaks = zone_display_levels, labels = zone_display_labels) +
  scale_x_date(limits = c(start_date, end_date), date_breaks = "6 months", date_labels = "%b %Y") +
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0, NA)) +
  labs(x = NULL, y = "Survey malaria prevalence (points) and\ntest positivity rate (lines), %", colour = "Zone") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")

s1_fig_cases <-
  df_zone_monthly |>
  filter(date >= start_date, date <= end_date) |>
  select(zone, date, cas_suspect, tdr_positif) |>
  pivot_longer(cols = c(cas_suspect, tdr_positif), names_to = "case_type", values_to = "n") |>
  mutate(case_type = fct_recode(case_type, "Suspected" = "cas_suspect", "Confirmed (RDT+)" = "tdr_positif")) |>
  ggplot(aes(x = date, y = n, colour = zone, linetype = case_type)) +
  geom_vline(xintercept = survey_dates$date, linetype = "dotted", colour = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = zone_cols, breaks = zone_display_levels, labels = zone_display_labels) +
  scale_linetype_manual(values = c("Suspected" = "dashed", "Confirmed (RDT+)" = "solid")) +
  scale_x_date(limits = c(start_date, end_date), date_breaks = "6 months", date_labels = "%b %Y") +
  labs(x = NULL, y = "Cases (monthly, mean per health area)", colour = "Zone", linetype = "Case type") +
  guides(colour = "none") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")

s1_fig_monthly <- s1_fig_tpr / s1_fig_cases + plot_layout(guides = "collect", heights = c(2, 1))
ggsave("output/figures/s1_figure_monthly_facility_tpr.png", s1_fig_monthly, width = 11, height = 9, dpi = 300)

survey_prev |>
  transmute(survey, zone = zone_display(zone), prevalence_pct = round(prev * 100, 1)) |>
  arrange(survey, zone) |>
  print(n = Inf)

# S11 Table ----

n_areas_tdr <-
  df_pnlp |>
  filter(indicateur_name == "Positivite Tdr", !is.na(zone), !is.na(taux_moyen)) |>
  distinct(date, zone, area) |>
  count(date, zone, name = "n_areas")

tpr_monthly <-
  df_tdr |>
  filter(date >= start_date, date <= end_date) |>
  left_join(n_areas_tdr, by = c("date", "zone")) |>
  mutate(zone = zone_display(zone), tpr_pct = 100 * taux_moyen) |>
  arrange(date, zone) |>
  select(date, zone, tpr_pct, n_areas)

zone_names <- as.character(zone_display(zone_display_levels))

s11_table <-
  tpr_monthly |>
  mutate(month = format(date, "%Y-%m"), tpr_pct = round(tpr_pct, 1)) |>
  pivot_wider(
    id_cols = month, names_from = zone, values_from = c(tpr_pct, n_areas),
    names_glue = "{zone} {.value}"
  ) |>
  select(month, all_of(str_c(zone_names, " tpr_pct")), all_of(str_c(zone_names, " n_areas"))) |>
  rename_with(\(x) str_replace(x, " tpr_pct$", " TPR (%)")) |>
  rename_with(\(x) str_replace(x, " n_areas$", " health areas (n)"))

# TPR is written as text so that a value such as 78.0 keeps its decimal
s11_table |>
  mutate(across(ends_with("TPR (%)"), \(x) sprintf("%.1f", x))) |>
  write_csv("output/tables/s11_table_tpr_monthly_by_zone.csv")

# Figure 6 ----

concord <- pair_prev_tpr(df_zone_monthly)

concordance_results <-
  bind_rows(
    pair_prev_tpr(df_zone_monthly) |> concordance_stats() |> mutate(source = "pnlp_superset", dataset = "full_5_surveys"),
    pair_prev_tpr(df_zone_monthly) |> filter(survey != 4) |> concordance_stats() |> mutate(source = "pnlp_superset", dataset = "excl_survey_4"),
    pair_prev_tpr(df_zone_monthly_dhis2) |> concordance_stats() |> mutate(source = "dhis2_ids", dataset = "full_5_surveys"),
    pair_prev_tpr(df_zone_monthly_dhis2) |> filter(survey != 4) |> concordance_stats() |> mutate(source = "dhis2_ids", dataset = "excl_survey_4")
  ) |>
  select(source, dataset, everything())

concordance_results |> print(width = Inf)
write_csv(concordance_results, "output/tables/figure6_concordance_stats.csv")

concord_by_survey <- agreement_by_survey(concord)
print(concord_by_survey)
write_csv(concord_by_survey, "output/tables/figure6_concordance_by_survey.csv")

agree_lowest_n <- sum(concord_by_survey$agree_lowest)
agree_full_n <- sum(concord_by_survey$agree_full_ranking)
agree_n_surveys <- nrow(concord_by_survey)

fig6 <-
  ggplot(concord, aes(x = prev_pct, y = tpr_pct)) +
  geom_mark_hull(
    aes(colour = zone, fill = zone),
    concavity = 8, expand = unit(4, "mm"), radius = unit(3, "mm"),
    alpha = 0.08, linewidth = 0.4, linetype = "dashed", show.legend = FALSE
  ) +
  geom_errorbar(aes(xmin = lcl * 100, xmax = ucl * 100, colour = zone), width = 0.8, linewidth = 0.4, alpha = 0.6) +
  geom_point(aes(colour = zone, shape = factor(survey)), size = 3.5) +
  scale_colour_manual(values = zone_cols, breaks = zone_display_levels, labels = zone_display_labels) +
  scale_fill_manual(values = zone_cols) +
  scale_shape_manual(values = c("1" = 15, "2" = 16, "3" = 17, "4" = 18, "5" = 8), labels = str_c("Survey ", 1:5), name = NULL) +
  annotate(
    "text",
    x = min(concord$prev_pct, na.rm = TRUE) + 5, y = max(concord$tpr_pct, na.rm = TRUE) - 5,
    label = sprintf(
      "Lowest zone agrees in %d of %d surveys\nFull ranking agrees in %d of %d",
      agree_lowest_n, agree_n_surveys, agree_full_n, agree_n_surveys
    ),
    hjust = 0, size = 3
  ) +
  labs(x = "Survey malaria prevalence (%)", y = "Surveillance TPR (%)", colour = "Zone") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("output/figures/figure6_tpr_prevalence_concordance.png", fig6, width = 8, height = 6, dpi = 300)
