# Mortality: Table 2 (U5MR and CMR by zone), S2 Figure (crude mortality rates), S3 Table
# panels A to E (ratio of RRs, Awasi-only comparator, survey-3 person-time from 11 July 2022,
# MDA+IRS vs IRS, design-based models), S10 Table (cause of death) and the CMR dispersion ratio (S1 Text).

pacman::p_load(broom, emmeans, lme4, srvyr, survey, tidyverse)

source("code/functions.R")

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# Helpers ----

# All-age models use the optimiser chosen by an lme4::allFit() search from default starting values (not the
# bobyqa estimates): the non-singular fit with the highest log-likelihood, which also beat the model without
# displacement status. The search takes about an hour, so the chosen optimisers are fixed here.
nlminb_control <- glmerControl(optimizer = "nlminbwrap", optCtrl = list(eval.max = 1e5))
nlopt_nm_control <- glmerControl(optimizer = "nloptwrap", optCtrl = list(algorithm = "NLOPT_LN_NELDERMEAD", maxeval = 1e5))
nm_control <- glmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 1e5))

# `data_name` is the name of a data frame in the global environment: emmeans() re-evaluates the stored call.
fit_status_model <- function(formula, data_name, control) {
  eval(
    substitute(glmer(f, data = d, family = poisson, control = ctl), list(f = formula, d = as.name(data_name), ctl = control)),
    envir = globalenv()
  )
}

# Contrasts of each zone against a reference zone; mvt adjustment over the contrast family.
zone_vs_ref <- function(emm, ref = "Routine") {
  contrast(emm, method = "trt.vs.ctrl", ref = ref, adjust = "mvt") |>
    summary(infer = TRUE) |>
    as_tibble() |>
    rename(rr = ratio, lcl = asymp.LCL, ucl = asymp.UCL) |>
    mutate(zone = str_remove(as.character(contrast), str_c(" / ", ref)))
}

# Survey-weighted crude rate per 10,000 person-days, one design per survey (ids = cluster, strata = zone).
build_mort_design <- \(df) as_survey(df, ids = cluster_uid, strata = zone, weights = weight)

crude_rates <- function(data) {
  data |>
    mutate(survey_n = as.integer(as.character(survey))) |>
    group_split(survey_n) |>
    map(\(df) {
      build_mort_design(df) |>
        group_by(zone) |>
        summarise(
          deaths = unweighted(sum(outcome_died)),
          pt_sum = unweighted(sum(pt)),
          rate = survey_ratio(outcome_died * 10000, pt, vartype = "ci", na.rm = TRUE)
        ) |>
        mutate(survey_n = unique(df$survey_n))
    }) |>
    bind_rows() |>
    rename(lcl = rate_low, ucl = rate_upp)
}

# Crude RR by zone from a quasi-Poisson model on the survey design.
crude_rr <- function(data) {
  data |>
    mutate(survey_n = as.integer(as.character(survey))) |>
    group_split(survey_n) |>
    map(\(df) {
      fit <- svyglm(outcome_died ~ zone + offset(log_pt), design = build_mort_design(df), family = quasipoisson())
      tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
        filter(term != "(Intercept)") |>
        transmute(zone = str_remove(term, "^zone"), crude_rr = estimate, survey_n = unique(df$survey_n))
    }) |>
    bind_rows()
}

# Data ----

df_ind <- read_rds("data/df_ind.rds")

# Surveys 1, 3 and 5 only: the rounds with a mortality module. Person-time is used as recorded.
dat_mort <-
  df_ind |>
  filter(survey %in% c(1, 3, 5)) |>
  select(survey, cluster, index_hh, zone, area, age, sex, status, weight, died, pt) |>
  rename(age_years = age) |>
  filter(!is.na(died), !is.na(zone), pt > 0) |>
  mutate(
    survey = factor(survey),
    zone = factor(zone, levels = zone_levels),
    outcome_died = as.integer(died),
    sex = factor(sex, levels = c("Homme", "Femme")),
    status = factor(status, levels = c("Autochtone", "Deplace")),
    age_years = if_else(age_years > 100, 100, age_years),
    log_pt = log(pt),
    cluster_uid = paste(survey, cluster, sep = "_")
  )

# Area-level baseline mortality (survey 1), winsorised at the 95th percentile of the area values
# and standardised. Areas without a survey-1 value take the mean.
add_baseline <- function(data, suffix) {
  baseline <-
    data |>
    filter(survey == "1") |>
    summarise(baseline = sum(outcome_died) / sum(pt) * 10000, .by = area)
  col <- str_c("baseline_mort_", suffix)
  mean_baseline <- mean(baseline$baseline)
  p95 <- quantile(baseline$baseline, 0.95)
  data |>
    left_join(baseline |> rename(!!col := baseline), by = "area") |>
    mutate(
      !!col := replace_na(.data[[col]], mean_baseline),
      !!str_c(col, "_w") := pmin(.data[[col]], p95),
      !!str_c(col, "_z") := as.numeric(scale(.data[[str_c(col, "_w")]]))
    )
}

dat_mort <- add_baseline(dat_mort, "ac")

# U5MR uses the under-five baseline rate.
dat_mort_u5 <- dat_mort |> filter(age_years < 5) |> add_baseline("u5")

# Primary models ----

# All-age pooled CMR models: displacement status in both; the baseline-adjusted model supplies
# surveys 3 and 5, the model without the baseline covariate supplies survey 1.
m_cmr_adj <- fit_status_model(
  outcome_died ~ zone * survey + age_years + sex + status + baseline_mort_ac_z + offset(log_pt) + (1 | cluster_uid/index_hh),
  "dat_mort", nlminb_control
)
m_cmr_unadj <- fit_status_model(
  outcome_died ~ zone * survey + age_years + sex + status + offset(log_pt) + (1 | cluster_uid/index_hh),
  "dat_mort", nlminb_control
)

# Pearson dispersion of the CMR model (S1 Text 4.6); overdispersion is not modelled.
dispersion_cmr <- sum(residuals(m_cmr_adj, type = "pearson")^2) / df.residual(m_cmr_adj)
write_csv(tibble(dispersion_ratio = dispersion_cmr), "output/tables/s1_text_cmr_dispersion.csv")
message("CMR Pearson dispersion ratio: ", round(dispersion_cmr, 4))

contrasts_cmr_adj <- emmeans(m_cmr_adj, ~ zone | survey, type = "response", offset = log(10000)) |>
  zone_vs_ref() |>
  mutate(survey_n = as.integer(as.character(survey)))
contrasts_cmr_unadj <- emmeans(m_cmr_unadj, ~ zone | survey, type = "response", offset = log(10000)) |>
  zone_vs_ref() |>
  mutate(survey_n = as.integer(as.character(survey)))

# U5MR: separate single-survey models with a cluster random intercept (a household term gives a
# singular fit on the few under-five deaths). Surveys 3 and 5 adjust for baseline U5MR.
dat_mort_u5_s1 <- filter(dat_mort_u5, survey == "1")
dat_mort_u5_s3 <- filter(dat_mort_u5, survey == "3")
dat_mort_u5_s5 <- filter(dat_mort_u5, survey == "5")
m_u5_s1_unadj <- glmer(
  outcome_died ~ zone + sex + status + offset(log_pt) + (1 | cluster_uid),
  data = dat_mort_u5_s1, family = poisson, control = glmerControl(optimizer = "bobyqa")
)
m_u5_s3_adj <- glmer(
  outcome_died ~ zone + sex + status + baseline_mort_u5_z + offset(log_pt) + (1 | cluster_uid),
  data = dat_mort_u5_s3, family = poisson, control = glmerControl(optimizer = "bobyqa")
)
m_u5_s5_adj <- glmer(
  outcome_died ~ zone + sex + status + baseline_mort_u5_z + offset(log_pt) + (1 | cluster_uid),
  data = dat_mort_u5_s5, family = poisson, control = glmerControl(optimizer = "bobyqa")
)

contrasts_u5_unadj <- emmeans(m_u5_s1_unadj, ~ zone, type = "response", offset = log(10000)) |> zone_vs_ref() |> mutate(survey_n = 1L)
contrasts_u5 <- bind_rows(
  emmeans(m_u5_s3_adj, ~ zone, type = "response", offset = log(10000)) |> zone_vs_ref() |> mutate(survey_n = 3L),
  emmeans(m_u5_s5_adj, ~ zone, type = "response", offset = log(10000)) |> zone_vs_ref() |> mutate(survey_n = 5L)
)

# Table 2 ----

timepoint_lu <- survey_timing |> filter(survey_n %in% c(1, 3, 5)) |> select(survey_n, months_from_baseline)

# A mortality rate cannot be negative: the lower limit is truncated at zero.
crude_cmr <- crude_rates(dat_mort) |> mutate(lcl = pmax(lcl, 0))
crude_u5mr <- crude_rates(dat_mort_u5) |> mutate(lcl = pmax(lcl, 0))

# Survey 1 comes from the model without the baseline covariate, which is derived from survey 1 itself.
build_table2 <- function(crude, crude_rr_tbl, contrasts_adj, contrasts_unadj) {
  crude |>
    mutate(dn_fmt = sprintf("%d / %.0f", deaths, pt_sum), crude_fmt = fmt_ci(rate, lcl, ucl)) |>
    select(zone, survey_n, dn_fmt, crude_fmt) |>
    left_join(
      bind_rows(
        crude_rr_tbl |> mutate(crude_rr_ci = sprintf("%.2f", crude_rr)) |> select(zone, survey_n, crude_rr_ci),
        tibble(zone = "Routine", survey_n = c(1, 3, 5), crude_rr_ci = "Reference")
      ),
      by = c("zone", "survey_n")
    ) |>
    left_join(
      bind_rows(
        contrasts_unadj |> filter(survey_n == 1) |> mutate(adj_rr_ci = fmt_ci(rr, lcl, ucl)) |> select(zone, survey_n, adj_rr_ci, p.value),
        contrasts_adj |> filter(survey_n != 1) |> mutate(adj_rr_ci = fmt_ci(rr, lcl, ucl)) |> select(zone, survey_n, adj_rr_ci, p.value),
        tibble(zone = "Routine", survey_n = c(1, 3, 5), adj_rr_ci = "Reference", p.value = NA_real_)
      ),
      by = c("zone", "survey_n")
    ) |>
    left_join(timepoint_lu, by = "survey_n") |>
    mutate(
      zone = zone_display(zone),
      p_fmt = case_when(survey_n == 1 ~ "---", p.value < 0.001 ~ "<0.001", .default = sprintf("%.3f", p.value))
    ) |>
    arrange(survey_n, zone) |>
    select(
      Zone = zone, Survey = survey_n, `Months from baseline` = months_from_baseline, `Deaths/person-days` = dn_fmt,
      `Crude rate (95% CI)` = crude_fmt, `Crude RR` = crude_rr_ci, `RR (95% CI)` = adj_rr_ci, `p-value` = p_fmt
    )
}

table2_u5 <- build_table2(crude_u5mr, crude_rr(dat_mort_u5), contrasts_u5, contrasts_u5_unadj)
table2_allage <- build_table2(crude_cmr, crude_rr(dat_mort), contrasts_cmr_adj, contrasts_cmr_unadj)
write_csv(table2_u5, "output/tables/table2_mortality_u5.csv")
write_csv(table2_allage, "output/tables/table2_mortality_allage.csv")
print(table2_u5)
print(table2_allage)

# S2 Figure ----

fig_crude <-
  bind_rows(crude_u5mr |> mutate(outcome = "U5MR"), crude_cmr |> mutate(outcome = "CMR")) |>
  left_join(timepoint_lu, by = "survey_n") |>
  mutate(zone = factor(zone, levels = zone_display_levels), outcome = factor(outcome, levels = c("U5MR", "CMR"))) |>
  ggplot(aes(x = months_from_baseline, y = rate, colour = zone, fill = zone, group = zone)) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = zone_cols) +
  scale_fill_manual(values = zone_cols) +
  scale_x_continuous(breaks = timepoint_lu$months_from_baseline, expand = expansion(mult = 0.15)) +
  facet_grid(outcome ~ zone, scales = "free_y", labeller = labeller(zone = zone_display_labels)) +
  labs(x = "Months from baseline", y = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"), legend.position = "none")

ggsave("output/figures/s2_figure_crude_mortality.png", fig_crude, width = 9, height = 6, dpi = 300)

# S3 Table A ----

# (zone_t / Routine_t) / (zone_1 / Routine_1): the zone-by-survey interaction of the primary CMR model.
emm_joint <- emmeans(m_cmr_adj, ~ zone * survey, type = "response", offset = log(10000))
emm_joint_df <- as_tibble(summary(emm_joint))

make_cmat_did <- function(df, zone_val, survey_val) {
  v <- rep(0, nrow(df))
  v[df$zone == zone_val & as.character(df$survey) == as.character(survey_val)] <- 1
  v[df$zone == zone_val & as.character(df$survey) == "1"] <- -1
  v[df$zone == "Routine" & as.character(df$survey) == as.character(survey_val)] <- -1
  v[df$zone == "Routine" & as.character(df$survey) == "1"] <- 1
  setNames(list(v), str_c(zone_val, "_S", survey_val, "_DiD"))
}

clist_did <- c(
  map(c("3", "5"), \(s) make_cmat_did(emm_joint_df, "IRS", s)),
  map(c("3", "5"), \(s) make_cmat_did(emm_joint_df, "MDA", s))
) |> unlist(recursive = FALSE)

contrasts_did <-
  contrast(emm_joint, method = clist_did) |>
  summary(infer = TRUE) |>
  as_tibble() |>
  rename(ror = ratio, lcl = asymp.LCL, ucl = asymp.UCL) |>
  mutate(zone = str_extract(contrast, "^[A-Z]+"), survey_n = as.integer(str_extract(contrast, "(?<=_S)\\d")))

crude_cmr_fmt <-
  crude_rates(dat_mort) |>
  mutate(dn_fmt = sprintf("%d / %.0f", deaths, pt_sum), crude_fmt = fmt_ci(rate, lcl, ucl)) |>
  select(zone, survey_n, dn_fmt, crude_fmt)

s3_table_a <-
  contrasts_did |>
  mutate(ror_fmt = fmt_ci(ror, lcl, ucl), p_fmt = fmt_p(p.value)) |>
  left_join(crude_cmr_fmt, by = c("zone", "survey_n")) |>
  left_join(timepoint_lu, by = "survey_n") |>
  mutate(zone = zone_display(zone)) |>
  arrange(survey_n, zone) |>
  select(
    Zone = zone, Survey = survey_n, `Months from baseline` = months_from_baseline, `Deaths/person-days` = dn_fmt,
    `Crude rate (95% CI)` = crude_fmt, `Ratio of RRs (95% CI)` = ror_fmt, `p-value` = p_fmt
  )
write_csv(s3_table_a, "output/tables/s3_table_a_ratio_of_rrs_mortality.csv")
print(s3_table_a)

# S3 Table B ----

# Routine held fixed at Awasi, the only health area never intervened, at every survey.
dat_mort_awasi <- dat_mort |>
  filter(zone != "Routine" | area == "Awasi") |>
  mutate(zone = fct_drop(zone))

m_awasi_adj <- fit_status_model(
  outcome_died ~ zone * survey + age_years + sex + status + baseline_mort_ac_z + offset(log_pt) + (1 | cluster_uid/index_hh),
  "dat_mort_awasi", nlminb_control
)
m_awasi_unadj <- fit_status_model(
  outcome_died ~ zone * survey + age_years + sex + status + offset(log_pt) + (1 | cluster_uid/index_hh),
  "dat_mort_awasi", nlminb_control
)

contrasts_awasi_adj <- emmeans(m_awasi_adj, ~ zone | survey, type = "response", offset = log(10000)) |>
  zone_vs_ref() |>
  mutate(survey_n = as.integer(as.character(survey)))
contrasts_awasi_unadj <- emmeans(m_awasi_unadj, ~ zone | survey, type = "response", offset = log(10000)) |>
  zone_vs_ref() |>
  mutate(survey_n = as.integer(as.character(survey)))

crude_awasi_fmt <-
  crude_rates(dat_mort_awasi) |>
  mutate(dn_fmt = sprintf("%d / %.0f", deaths, pt_sum), crude_fmt = fmt_ci(rate, lcl, ucl)) |>
  select(zone, survey_n, dn_fmt, crude_fmt)

s3_table_b <-
  bind_rows(
    contrasts_awasi_unadj |> filter(survey_n == 1) |> mutate(rr_fmt = fmt_ci(rr, lcl, ucl)) |> select(zone, survey_n, rr_fmt, p.value),
    contrasts_awasi_adj |> filter(survey_n != 1) |> mutate(rr_fmt = fmt_ci(rr, lcl, ucl)) |> select(zone, survey_n, rr_fmt, p.value),
    tibble(zone = "Routine", survey_n = c(1, 3, 5), rr_fmt = "Reference", p.value = NA_real_)
  ) |>
  left_join(crude_awasi_fmt, by = c("zone", "survey_n")) |>
  left_join(timepoint_lu, by = "survey_n") |>
  mutate(
    zone = zone_display(zone),
    p_fmt = case_when(survey_n == 1 ~ "---", p.value < 0.001 ~ "<0.001", .default = sprintf("%.3f", p.value))
  ) |>
  arrange(survey_n, zone) |>
  select(
    Zone = zone, Survey = survey_n, `Months from baseline` = months_from_baseline, `Deaths/person-days` = dn_fmt,
    `Crude rate (95% CI)` = crude_fmt, `RR (95% CI)` = rr_fmt, `p-value` = p_fmt
  )
write_csv(s3_table_b, "output/tables/s3_table_b_awasi_comparator_mortality.csv")
print(s3_table_b)

# S3 Table C ----

# Person-time restarts at the later of the recall start and 11 July 2022 (MDA round 1); deaths are
# those on or after it. A single-survey fit on the full recall window is reported beside it.
window_start <- as.Date("2022-07-11")
# The baseline covariate (survey-1 rates by area, winsorised and standardised) is rebuilt on this
# data frame, so it is standardised over the survey-3 rows.
baseline_ac <- dat_mort |>
  filter(survey == "1") |>
  summarise(baseline_mort_ac = sum(outcome_died) / sum(pt) * 10000, .by = area)
baseline_ac_mean <- mean(baseline_ac$baseline_mort_ac)
baseline_ac_p95 <- quantile(baseline_ac$baseline_mort_ac, 0.95)
dat_mort3_post <-
  df_ind |>
  filter(survey == 3) |>
  select(cluster, index_hh, zone, area, age, sex, status, weight, pt_start, pt_end, date_died) |>
  mutate(
    survey = 3L,
    zone = factor(zone, levels = zone_levels),
    pt = as.numeric(as.Date(pt_end) - pmax(as.Date(pt_start), window_start)),
    died = !is.na(date_died) & as.Date(date_died) >= window_start
  ) |>
  select(survey, cluster, index_hh, zone, area, age, sex, status, weight, died, pt) |>
  rename(age_years = age) |>
  filter(!is.na(died), !is.na(zone), pt > 0) |>
  mutate(
    survey = factor(survey),
    outcome_died = as.integer(died),
    sex = factor(sex, levels = c("Homme", "Femme")),
    status = factor(status, levels = c("Autochtone", "Deplace")),
    age_years = if_else(age_years > 100, 100, age_years),
    log_pt = log(pt),
    cluster_uid = paste(survey, cluster, sep = "_")
  ) |>
  left_join(baseline_ac, by = "area") |>
  mutate(
    baseline_mort_ac = replace_na(baseline_mort_ac, baseline_ac_mean),
    baseline_mort_ac_w = pmin(baseline_mort_ac, baseline_ac_p95),
    baseline_mort_ac_z = as.numeric(scale(baseline_mort_ac_w))
  )

fit_mort3 <- \(data_name, control) fit_status_model(
  outcome_died ~ zone + age_years + sex + status + baseline_mort_ac_z + offset(log_pt) + (1 | cluster_uid/index_hh),
  data_name, control
)

dat_mort3_full <- dat_mort |> filter(survey == "3") |> mutate(survey = fct_drop(survey))
m_mort3_full <- fit_mort3("dat_mort3_full", nlopt_nm_control)
m_mort3_post <- fit_mort3("dat_mort3_post", nm_control)

contrasts_mort3_full <- emmeans(m_mort3_full, ~ zone, type = "response") |> zone_vs_ref()
contrasts_mort3_post <- emmeans(m_mort3_post, ~ zone, type = "response") |> zone_vs_ref()

crude_mort3_post <- dat_mort3_post |>
  summarise(deaths = sum(outcome_died), person_days = sum(pt), .by = zone) |>
  mutate(rate_per_10000 = deaths / person_days * 1e4)

months_s3 <- timepoint_lu$months_from_baseline[timepoint_lu$survey_n == 3]
primary_s3 <- contrasts_cmr_adj |> filter(survey_n == 3) |> transmute(zone, primary = fmt_ci(rr, lcl, ucl))
full_s3 <- contrasts_mort3_full |> transmute(zone, full = fmt_ci(rr, lcl, ucl))

s3_table_c <-
  crude_mort3_post |>
  mutate(zone = as.character(zone)) |>
  left_join(contrasts_mort3_post |> select(zone, rr, lcl, ucl, p.value), by = "zone") |>
  left_join(primary_s3, by = "zone") |>
  left_join(full_s3, by = "zone") |>
  transmute(
    Zone = zone_display(zone), Survey = 3L, `Months from baseline` = months_s3,
    `Deaths / person-days (from 11 July 2022)` = str_glue("{deaths} / {format(round(person_days), big.mark = ',')}"),
    `Crude rate per 10,000 person-days` = sprintf("%.2f", rate_per_10000),
    `Primary pooled-model adjusted RR (95% CI)` = if_else(zone == "Routine", "Reference", primary),
    `Single-survey adjusted RR, full recall window (95% CI)` = if_else(zone == "Routine", "Reference", full),
    `Adjusted RR, post-intervention person-time (95% CI)` = if_else(zone == "Routine", "Reference", fmt_ci(rr, lcl, ucl)),
    `p-value` = if_else(zone == "Routine", "---", fmt_p(p.value))
  ) |>
  arrange(Zone)
write_csv(s3_table_c, "output/tables/s3_table_c_survey3_person_time_mortality.csv")
print(s3_table_c)

# S3 Table D ----

# Reported interval and p-value are mvt-adjusted over the family of contrasts with IRS as reference
# (Routine vs IRS, MDA+IRS vs IRS); the unadjusted Wald interval is kept alongside.
mda_vs_irs <- function(emm) {
  contrast_rows <- \(adjust) {
    contrast(emm, method = "trt.vs.ctrl", ref = "IRS", adjust = adjust) |>
      summary(infer = TRUE) |>
      as_tibble() |>
      filter(str_detect(as.character(contrast), "^MDA"))
  }
  unadj <- contrast_rows("none")
  adj <- contrast_rows("mvt")
  unadj |>
    transmute(
      survey_n = if ("survey" %in% names(unadj)) as.integer(as.character(survey)) else NA_integer_,
      ratio = exp(estimate),
      lcl = exp(adj$asymp.LCL), ucl = exp(adj$asymp.UCL), p = adj$p.value,
      lcl_unadj = exp(asymp.LCL), ucl_unadj = exp(asymp.UCL), p_unadj = p.value
    )
}

s3_table_d <-
  bind_rows(
    mda_vs_irs(emmeans(m_cmr_adj, ~ zone | survey)) |> filter(survey_n %in% c(3L, 5L)) |> mutate(outcome = "CMR (all ages)"),
    mda_vs_irs(emmeans(m_u5_s3_adj, ~ zone)) |> mutate(outcome = "U5MR", survey_n = 3L),
    mda_vs_irs(emmeans(m_u5_s5_adj, ~ zone)) |> mutate(outcome = "U5MR", survey_n = 5L)
  ) |>
  left_join(survey_timing |> select(survey_n, months_from_baseline), by = "survey_n") |>
  arrange(outcome, survey_n) |>
  transmute(outcome, survey = survey_n, months = months_from_baseline, ratio, lcl, ucl, p, lcl_unadj, ucl_unadj, p_unadj)
write_csv(s3_table_d, "output/tables/s3_table_d_mda_irs_vs_irs_mortality.csv")
print(s3_table_d)

# S3 Table E ----

# The same rate ratios re-estimated by survey::svyglm (quasi-Poisson, log person-time offset,
# cluster as PSU, strata = survey x zone, survey weights), beside the mixed-model estimates.
# emmeans needs the model rows (data =) and the design degrees of freedom (df =).
u5_by_area <- dat_mort_u5 |> distinct(area, baseline_mort_u5_z)

dat_svy <-
  dat_mort |>
  select(survey, zone, cluster_uid, area, outcome_died, pt, log_pt, weight, age_years, sex, status, baseline_mort_ac_z) |>
  left_join(u5_by_area, by = "area") |>
  filter(!is.na(age_years), !is.na(sex), !is.na(status)) |>
  as.data.frame()

make_design <- function(data) {
  data$stratum <- interaction(data$survey, data$zone, drop = TRUE)
  svydesign(ids = ~cluster_uid, strata = ~stratum, weights = ~weight, data = data, nest = FALSE)
}

# Models mirror the primary objects, plus the offset. Under-five models use the under-five
# domain of the survey's design.
model_specs <- list(
  E1 = list(formula = outcome_died ~ zone * survey + age_years + sex + status + baseline_mort_ac_z + offset(log_pt), survey_n = NA_integer_, primary = m_cmr_adj),
  E2 = list(formula = outcome_died ~ zone * survey + age_years + sex + status + offset(log_pt), survey_n = NA_integer_, primary = m_cmr_unadj),
  E3 = list(formula = outcome_died ~ zone + sex + status + baseline_mort_u5_z + offset(log_pt), survey_n = 3L, primary = m_u5_s3_adj),
  E4 = list(formula = outcome_died ~ zone + sex + status + baseline_mort_u5_z + offset(log_pt), survey_n = 5L, primary = m_u5_s5_adj),
  E5 = list(formula = outcome_died ~ zone + sex + status + offset(log_pt), survey_n = 1L, primary = m_u5_s1_unadj)
)

fit_design <- function(spec) {
  if (is.na(spec$survey_n)) {
    data <- dat_svy
    design <- make_design(data)
  } else {
    dat_s <- filter(dat_svy, survey == as.character(spec$survey_n))
    data <- filter(dat_s, age_years < 5)
    design <- subset(make_design(dat_s), age_years < 5)
  }
  list(fit = svyglm(spec$formula, design = design, family = quasipoisson()), data = data)
}
fits_design <- map(model_specs, fit_design)

tidy_contrast <- function(ct) {
  s <- as_tibble(summary(ct, infer = TRUE))
  lo <- intersect(c("asymp.LCL", "lower.CL"), names(s))[1]
  hi <- intersect(c("asymp.UCL", "upper.CL"), names(s))[1]
  tibble(
    contrast = as.character(s$contrast),
    survey_n = if ("survey" %in% names(s)) as.integer(as.character(s$survey)) else NA_integer_,
    log_est = s$estimate, lcl_log = s[[lo]], ucl_log = s[[hi]], p = s$p.value
  )
}

# Zone contrasts (mvt-adjusted within the reference-zone family) and, for pooled models, the ratios of
# RRs versus survey 1.
get_contrasts <- function(fit, pooled, survey_n, data = NULL, df = NULL) {
  grid <- \(specs) if (is.null(data)) emmeans(fit, specs) else emmeans(fit, specs, data = data, df = df)
  emm <- grid(if (pooled) ~ zone | survey else ~ zone)
  zone_ct <- \(ref, keep) {
    contrast(emm, method = "trt.vs.ctrl", ref = ref, adjust = "mvt") |>
      tidy_contrast() |>
      filter(str_detect(contrast, keep)) |>
      mutate(comparison = str_c(str_extract(contrast, "^[A-Za-z]+"), " vs ", ref))
  }
  out <- bind_rows(zone_ct("Routine", "^(MDA|IRS) - Routine"), zone_ct("IRS", "^MDA - IRS")) |> select(-contrast)
  if (!pooled) out$survey_n <- survey_n
  if (pooled) {
    emm_j <- grid(~ zone * survey)
    grid_df <- as_tibble(summary(emm_j))
    clist <- c(
      map(c("3", "5"), \(s) make_cmat_did(grid_df, "IRS", s)),
      map(c("3", "5"), \(s) make_cmat_did(grid_df, "MDA", s))
    ) |> unlist(recursive = FALSE)
    did <- contrast(emm_j, method = clist) |>
      tidy_contrast() |>
      mutate(survey_n = as.integer(str_extract(contrast, "(?<=_S)\\d")), comparison = str_c("DiD ", str_extract(contrast, "^[A-Z]+"))) |>
      select(-contrast)
    out <- bind_rows(out, did)
  }
  out
}

ct_primary <- imap(model_specs, \(spec, id) {
  get_contrasts(spec$primary, pooled = is.na(spec$survey_n), survey_n = spec$survey_n) |> mutate(model = id)
}) |> bind_rows()
ct_design <- imap(model_specs, \(spec, id) {
  f <- fits_design[[id]]
  get_contrasts(f$fit, pooled = is.na(spec$survey_n), survey_n = spec$survey_n, data = f$data, df = df.residual(f$fit)) |> mutate(model = id)
}) |> bind_rows()

comparison_label <- c(
  "MDA vs Routine" = str_c(zone_display_labels[["MDA"]], " vs Routine"),
  "IRS vs Routine" = "IRS vs Routine",
  "MDA vs IRS" = str_c(zone_display_labels[["MDA"]], " vs IRS"),
  "DiD MDA" = str_c("Ratio of RRs vs survey 1: ", zone_display_labels[["MDA"]]),
  "DiD IRS" = "Ratio of RRs vs survey 1: IRS"
)

panel_rows <- bind_rows(
  tibble(outcome = "U5MR", model = "E5", survey_n = 1L, comparison = c("MDA vs Routine", "IRS vs Routine")),
  expand_grid(outcome = "U5MR", model = "E3", survey_n = 3L, comparison = c("MDA vs Routine", "IRS vs Routine", "MDA vs IRS")),
  expand_grid(outcome = "U5MR", model = "E4", survey_n = 5L, comparison = c("MDA vs Routine", "IRS vs Routine", "MDA vs IRS")),
  tibble(outcome = "CMR (all ages)", model = "E2", survey_n = 1L, comparison = c("MDA vs Routine", "IRS vs Routine")),
  expand_grid(outcome = "CMR (all ages)", model = "E1", survey_n = c(3L, 5L), comparison = names(comparison_label))
)

s3_table_e <-
  panel_rows |>
  left_join(ct_primary |> rename_with(\(x) str_c(x, "_primary"), c(log_est, lcl_log, ucl_log, p)), by = c("model", "survey_n", "comparison")) |>
  left_join(ct_design |> rename_with(\(x) str_c(x, "_design"), c(log_est, lcl_log, ucl_log, p)), by = c("model", "survey_n", "comparison")) |>
  left_join(survey_timing |> transmute(survey_n = as.integer(survey_n), months_from_baseline), by = "survey_n") |>
  mutate(
    rr_primary = exp(log_est_primary), lcl_primary = exp(lcl_log_primary), ucl_primary = exp(ucl_log_primary),
    p_primary = if_else(survey_n == 1L, NA_real_, p_primary),
    rr_design = exp(log_est_design), lcl_design = exp(lcl_log_design), ucl_design = exp(ucl_log_design),
    p_design = if_else(survey_n == 1L, NA_real_, p_design),
    ratio_design_to_primary = rr_design / rr_primary,
    note = if_else(outcome == "U5MR" & survey_n == 5L, "uninformative: 4 Routine under-five deaths; excluded from the decision rule", "")
  ) |>
  transmute(
    outcome, survey = survey_n, months_from_baseline, comparison = unname(comparison_label[comparison]),
    rr_primary, lcl_primary, ucl_primary, p_primary, rr_design, lcl_design, ucl_design, p_design, ratio_design_to_primary, note
  )
write_csv(s3_table_e, "output/tables/s3_table_e_design_based_mortality.csv", na = "---")
print(s3_table_e, n = Inf)

# S10 Table ----

# Deaths included in the mortality models, by zone, survey and age group; counts and proportions only.
cause_levels <- c(
  "Malaria and/or anaemia", "Cough/respiratory", "Diarrhoea/vomiting", "Measles",
  "Malnutrition/oedema", "Accident/violence", "Childbirth complication",
  "Other", "Unknown (don't know / not recorded)"
)

# Survey 3 spells "Ne sais pas", surveys 1 and 5 "Ne sait pas"; matched on the stem.
bin_cause <- \(x) {
  case_when(
    x == "Malaria et/ou anemie" ~ "Malaria and/or anaemia",
    str_detect(x, "^Toux") ~ "Cough/respiratory",
    str_detect(x, "^Diarhee") ~ "Diarrhoea/vomiting",
    str_detect(x, "^Rougeole") ~ "Measles",
    x %in% c("Malnutrition_Marasme", "oedeme") ~ "Malnutrition/oedema",
    x %in% c("Accident", "Violence") ~ "Accident/violence",
    str_detect(x, "^Complication") ~ "Childbirth complication",
    x == "Autres" ~ "Other",
    str_detect(x, "^Ne sai") | is.na(x) ~ "Unknown (don't know / not recorded)",
    .default = NA_character_
  ) |>
    factor(levels = cause_levels)
}

deaths <-
  df_ind |>
  filter(survey %in% c(1, 3, 5)) |>
  select(survey, zone, age, died, pt, died_cause) |>
  filter(!is.na(died), !is.na(zone), pt > 0, died) |>
  mutate(
    zone = factor(as.character(zone), levels = zone_levels),
    cause = bin_cause(died_cause),
    age_years = if_else(age > 100, 100, age)
  )
stopifnot(!anyNA(deaths$cause))

s10_table <-
  bind_rows(
    deaths |> mutate(age_group = "All ages"),
    deaths |> filter(age_years < 5) |> mutate(age_group = "Under five")
  ) |>
  mutate(age_group = factor(age_group, levels = c("Under five", "All ages"))) |>
  count(zone, survey, age_group, cause, .drop = FALSE) |>
  group_by(zone, survey, age_group) |>
  mutate(deaths_total = sum(n), pct = if_else(deaths_total > 0, 100 * n / deaths_total, NA_real_)) |>
  ungroup() |>
  filter(survey != 2L) |>
  mutate(zone = zone_display(zone)) |>
  arrange(survey, age_group, zone, cause) |>
  rename(deaths_n = n, pct_of_deaths = pct)
write_csv(s10_table, "output/tables/s10_table_cause_of_death.csv")
