# Sensitivity analyses of the malaria prevalence contrasts.
# Reproduces S2 Table panels A (ratio of ORs versus survey 1, residual survey 1 ORs), B (Awasi-only comparator),
# C (health-area random intercept), D (area-level permutation test), E (MDA+IRS versus IRS), F (E-values)
# and G (net-adjusted ORs), and S12 Table (net use by zone). Needs output/models/prevalence_model.rds from 02_prevalence.R.

pacman::p_load(lme4, emmeans, srvyr, survey, EValue, tidyverse)

source("code/functions.R")

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

# Model settings ----

glmer_control <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# Primary model ----

model <- read_rds("output/models/prevalence_model.rds")
dat_rdt <- model$dat_rdt
m4_durability <- model$m4_durability
m4_durability_adj <- model$m4_durability_adj

# The multivariate-t critical values are simulated
set.seed(1)

emm_prev <- emmeans(m4_durability_adj, ~ zone | survey, type = "response")
emm_prev_unadj <- emmeans(m4_durability, ~ zone | survey, type = "response")
contrasts_prev <- tidy_contrasts(emm_prev) |> mutate(survey_n = as.integer(as.character(survey)))

crude_prev <- prevalence_by_zone(dat_rdt)

# S2 Table A ----

# Each zone's change from its own survey 1 relative to Routine's change (the zone:survey interaction)
emm_joint <- emmeans(m4_durability_adj, ~ zone * survey, type = "response")
emm_joint_df <- as_tibble(summary(emm_joint))

did_weights <- function(zone_val, survey_val) {
  v <- rep(0, nrow(emm_joint_df))
  v[emm_joint_df$zone == zone_val & as.character(emm_joint_df$survey) == survey_val] <- 1
  v[emm_joint_df$zone == zone_val & as.character(emm_joint_df$survey) == "1"] <- -1
  v[emm_joint_df$zone == "Routine" & as.character(emm_joint_df$survey) == survey_val] <- -1
  v[emm_joint_df$zone == "Routine" & as.character(emm_joint_df$survey) == "1"] <- 1
  set_names(list(v), str_c(zone_val, "_S", survey_val, "_DiD"))
}

did_list <- c(map(as.character(2:5), \(s) did_weights("IRS", s)), map(as.character(2:5), \(s) did_weights("MDA", s))) |> unlist(recursive = FALSE)

s2_table_a <-
  contrast(emm_joint, method = did_list) |>
  summary(infer = TRUE) |>
  as_tibble() |>
  rename(ror = odds.ratio, lcl = asymp.LCL, ucl = asymp.UCL) |>
  mutate(zone = str_extract(contrast, "^[A-Z]+"), survey_n = as.integer(str_extract(contrast, "(?<=_S)\\d"))) |>
  left_join(crude_prev |> mutate(prev_pct = sprintf("%.1f%% (%.1f-%.1f%%)", 100 * prev, 100 * lcl, 100 * ucl)) |> select(zone, survey_n, prev_pct), by = c("zone", "survey_n")) |>
  left_join(season_lu, by = "survey_n") |>
  mutate(zone = zone_display(zone)) |>
  arrange(survey_n, zone) |>
  transmute(
    Zone = zone, Survey = survey_n, `Months from baseline` = months_from_baseline, Season = season,
    `Crude malaria prevalence (95% CI)` = prev_pct, `Ratio of ORs (95% CI)` = fmt_ci(ror, lcl, ucl), `p-value` = fmt_p(p.value)
  )

print(s2_table_a, n = Inf, width = Inf)
write_csv(s2_table_a, "output/tables/s2_table_a_difference_in_differences_prevalence.csv")

# Survey 1 zone contrasts that remain after adjustment for the baseline covariate: imbalance, not effects
contrast_survey1 <- function(ref) {
  tidy_contrasts(emm_prev, ref) |>
    filter(survey == "1") |>
    mutate(treated = str_remove(as.character(contrast), str_c(" / ", ref)), ref = ref)
}

residual_ors <-
  bind_rows(
    contrast_survey1("Routine") |> filter(treated %in% c("MDA", "IRS")),
    contrast_survey1("IRS") |> filter(treated == "MDA")
  ) |>
  mutate(comparison = str_c(zone_display_labels[treated], " vs ", zone_display_labels[ref])) |>
  transmute(
    survey = 1L, months_from_baseline = 0L, comparison, or, lcl, ucl, p_value = p.value,
    model = "adjusted for baseline prevalence (area level), age, sex and displacement status; multivariate-t over the two zone contrasts"
  ) |>
  arrange(factor(comparison, levels = c("MDA+IRS vs Routine", "IRS vs Routine", "MDA+IRS vs IRS")))

print(residual_ors, width = Inf)
write_csv(residual_ors, "output/tables/s2_table_a_residual_survey1_ors.csv")

# S2 Table B ----

# Routine comparator held fixed at Awasi, the one area never intervened, at every survey
dat_rdt_awasi <-
  dat_rdt |>
  filter(zone != "Routine" | area == "Awasi") |>
  mutate(zone = fct_drop(zone))

m4_awasi_adj <- glmer(
  outcome_rdt ~ zone * survey + age_years + sex + status + baseline_prev_z + (1 | cluster_uid/index_hh),
  data = dat_rdt_awasi, family = binomial, control = glmer_control
)

contrasts_awasi <-
  tidy_contrasts(emmeans(m4_awasi_adj, ~ zone | survey, type = "response")) |>
  mutate(survey_n = as.integer(as.character(survey)))

crude_prev_awasi <-
  prevalence_by_zone(dat_rdt_awasi) |>
  mutate(prev_pct = sprintf("%.1f%% (%.1f-%.1f%%)", 100 * prev, 100 * lcl, 100 * ucl)) |>
  select(zone, survey_n, n_total, prev_pct)

s2_table_b <-
  contrasts_awasi |>
  mutate(or_ci = fmt_ci(or, lcl, ucl), p_fmt = if_else(survey_n == 1, "---", fmt_p(p.value))) |>
  bind_rows(tibble(zone = "Routine", survey_n = 1:5, or_ci = "Reference", p_fmt = "---")) |>
  left_join(crude_prev_awasi, by = c("zone", "survey_n")) |>
  left_join(season_lu, by = "survey_n") |>
  mutate(zone = zone_display(zone)) |>
  arrange(survey_n, zone) |>
  select(
    Zone = zone, Survey = survey_n, `Months from baseline` = months_from_baseline, Season = season, N = n_total,
    `Crude malaria prevalence (95% CI)` = prev_pct, `OR (95% CI)` = or_ci, `p-value` = p_fmt
  )

print(s2_table_b, n = Inf, width = Inf)
write_csv(s2_table_b, "output/tables/s2_table_b_awasi_only_comparator_prevalence.csv")

# S2 Table C ----

# Primary model plus a health-area random intercept. The area term is weakly identified on 14 areas, so
# both optimisers are fitted and the non-singular fit with the higher log-likelihood is used.
area_formula <- outcome_rdt ~ zone * survey + age_years + sex + status + baseline_prev_z + (1 | area) + (1 | cluster_uid/index_hh)

area_fits <- list(
  bobyqa = glmer(area_formula, data = dat_rdt, family = binomial, control = glmer_control),
  Nelder_Mead = glmer(area_formula, data = dat_rdt, family = binomial, control = glmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 2e5)))
)
area_ok <- !map_lgl(area_fits, isSingular)
area_ll <- map_dbl(area_fits, \(f) as.numeric(logLik(f)))
area_pick <- if (any(area_ok)) names(which.max(area_ll[area_ok])) else names(which.max(area_ll))
print(c(optimiser = area_pick, area_ll))

contrasts_area <-
  tidy_contrasts(emmeans(area_fits[[area_pick]], ~ zone | survey, type = "response")) |>
  mutate(survey_n = as.integer(as.character(survey)))

s2_table_c <-
  contrasts_area |>
  filter(survey_n >= 2) |>
  left_join(contrasts_prev |> transmute(zone, survey_n, `Primary adjusted OR (95% CI)` = fmt_ci(or, lcl, ucl)), by = c("zone", "survey_n")) |>
  left_join(season_lu, by = "survey_n") |>
  transmute(
    Zone = zone_display(zone), Survey = survey_n, `Months from baseline` = months_from_baseline,
    `Primary adjusted OR (95% CI)`,
    `Adjusted OR with area random intercept (95% CI)` = fmt_ci(or, lcl, ucl), `p-value` = fmt_p(p.value)
  ) |>
  arrange(Survey, Zone)

print(s2_table_c, n = Inf, width = Inf)
write_csv(s2_table_c, "output/tables/s2_table_c_area_random_intercept_prevalence.csv")

# S2 Table D ----

# Difference in mean area-level prevalence (unweighted, percentage points) between two zones. Zone labels are
# permuted over the areas present at the survey, keeping group sizes, and every labelling is enumerated (exact p).
area_level <-
  dat_rdt |>
  summarise(n = n(), prev = 100 * mean(outcome_rdt), .by = c(survey, zone, area)) |>
  mutate(survey_n = as.integer(as.character(survey)))

perm_one <- function(df, z1, z2) {
  d <- df |> filter(zone %in% c(z1, z2))
  x <- d$prev
  is_z1 <- d$zone == z1
  n1 <- sum(is_z1)
  idx <- combn(length(x), n1)
  stats <- apply(idx, 2, \(i) mean(x[i]) - mean(x[-i]))
  obs <- mean(x[is_z1]) - mean(x[!is_z1])
  tibble(
    zone = z1, comparator = z2,
    n_areas_zone = n1, n_areas_comparator = length(x) - n1,
    n_labellings = ncol(idx),
    mean_prev_zone = mean(x[is_z1]), mean_prev_comparator = mean(x[!is_z1]),
    difference_pp = obs,
    complete_separation = max(x[is_z1]) < min(x[!is_z1]),
    p_two_sided = mean(abs(stats) >= abs(obs) - 1e-9)
  )
}

perm_results <-
  map(2:5, \(s) {
    df <- area_level |> filter(survey_n == s)
    bind_rows(perm_one(df, "MDA", "Routine"), perm_one(df, "IRS", "Routine"), perm_one(df, "MDA", "IRS")) |>
      mutate(survey_n = s, .before = 1)
  }) |>
  list_rbind()

s2_table_d <-
  perm_results |>
  left_join(season_lu, by = "survey_n") |>
  transmute(
    Survey = survey_n, `Months from baseline` = months_from_baseline,
    Zone = zone_display(zone), Comparator = zone_display(comparator),
    `Areas (zone / comparator)` = str_glue("{n_areas_zone} / {n_areas_comparator}"),
    `Mean area prevalence, % (zone / comparator)` = sprintf("%.1f / %.1f", mean_prev_zone, mean_prev_comparator),
    `Difference (pp)` = sprintf("%.1f", difference_pp),
    `Every zone area below every comparator area` = if_else(complete_separation, "Yes", "No"),
    `Distinct labellings` = n_labellings,
    `Permutation p-value` = sprintf("%.3f", p_two_sided)
  )

print(s2_table_d, n = Inf, width = Inf)
write_csv(s2_table_d, "output/tables/s2_table_d_area_permutation_prevalence.csv")

# Survey 1 precedes the interventions: the same differences without a test
survey1_row <- function(zone_name, comparator_name) {
  prev_z <- area_level |> filter(survey_n == 1, zone_display(zone) == zone_name) |> pull(prev)
  prev_c <- area_level |> filter(survey_n == 1, zone_display(zone) == comparator_name) |> pull(prev)
  tibble(
    Survey = 1L, `Months from baseline` = 0L, Zone = zone_name, Comparator = comparator_name,
    `Areas (zone / comparator)` = str_glue("{length(prev_z)} / {length(prev_c)}"),
    `Mean area prevalence, % (zone / comparator)` = sprintf("%.1f / %.1f", mean(prev_z), mean(prev_c)),
    `Difference (pp)` = sprintf("%.1f", mean(prev_z) - mean(prev_c)),
    `Every zone area below every comparator area` = if_else(max(prev_z) < min(prev_c), "Yes", "No"),
    `Distinct labellings` = NA_integer_,
    `Permutation p-value` = NA_character_
  )
}

s2_table_d_survey1 <- bind_rows(survey1_row("MDA+IRS", "Routine"), survey1_row("IRS", "Routine"), survey1_row("MDA+IRS", "IRS"))

print(s2_table_d_survey1, width = Inf)
write_csv(s2_table_d_survey1, "output/tables/s2_table_d_area_permutation_prevalence_survey1.csv", na = "NA")

# S2 Table E ----

# MDA+IRS compared with IRS: the difference associated with adding MDA where IRS was also delivered,
# not an MDA-alone effect. Survey 1 comes from the model without the baseline covariate, with no p-value.
contrast_vs_irs <- function(emm) {
  tidy_contrasts(emm, "IRS") |>
    mutate(survey_n = as.integer(as.character(survey)), contrast = str_remove(as.character(contrast), " / IRS")) |>
    filter(contrast == "MDA") |>
    select(survey_n, or, lcl, ucl, p.value)
}

s2_table_e <-
  bind_rows(
    contrast_vs_irs(emm_prev_unadj) |> filter(survey_n == 1) |> mutate(p.value = NA_real_, model = "adjusted for age, sex and displacement status; no baseline covariate (survey 1)"),
    contrast_vs_irs(emm_prev) |> filter(survey_n != 1) |> mutate(model = "adjusted for baseline prevalence, age, sex and displacement status")
  ) |>
  arrange(survey_n) |>
  left_join(survey_timing |> select(survey_n, months_from_baseline), by = "survey_n") |>
  transmute(
    survey = as.integer(survey_n), months_from_baseline, comparison = "MDA+IRS vs IRS", or, lcl, ucl, p_value = p.value, model,
    interpretation = "Difference associated with adding MDA where IRS was also delivered; not an MDA-alone effect"
  )

print(s2_table_e, width = Inf)
write_csv(s2_table_e, "output/tables/s2_table_e_mda_irs_vs_irs_prevalence.csv")

# Ratio of the two zones' ratios of ORs versus survey 1: (MDA_t / MDA_1) / (IRS_t / IRS_1), Routine cancels.
# Wald interval from the joint covariance of the two interaction coefficients, no multiplicity adjustment.
cells <- as_tibble(summary(emmeans(m4_durability_adj, ~ zone * survey))) |> mutate(zone = as.character(zone), survey = as.character(survey))

mda_irs_weights <- function(survey_t) {
  w <- rep(0, nrow(cells))
  w[cells$zone == "MDA" & cells$survey == survey_t] <- 1
  w[cells$zone == "MDA" & cells$survey == "1"] <- -1
  w[cells$zone == "IRS" & cells$survey == survey_t] <- -1
  w[cells$zone == "IRS" & cells$survey == "1"] <- 1
  w
}

s2_table_e_did <-
  contrast(emmeans(m4_durability_adj, ~ zone * survey), method = set_names(map(as.character(2:5), mda_irs_weights), 2:5), adjust = "none") |>
  summary(infer = TRUE) |>
  as_tibble() |>
  transmute(survey = as.integer(as.character(contrast)), ratio = exp(estimate), lcl = exp(asymp.LCL), ucl = exp(asymp.UCL), p = p.value) |>
  left_join(survey_timing |> select(survey = survey_n, months_from_baseline), by = "survey") |>
  select(survey, months_from_baseline, ratio, lcl, ucl, p)

print(s2_table_e_did, width = Inf)
write_csv(s2_table_e_did, "output/tables/s2_table_e_did_mda_irs_vs_irs_prevalence.csv")

# S2 Table F ----

# E-values (VanderWeele and Ding 2017) for the within-survey ORs, with survey 1 benchmark rows. The risk ratio
# is approximated by sqrt(OR) because the outcome is common (Routine prevalence 51-64%).
zone_contrasts <- function(emm, ref, adjust, ratio_col = "odds.ratio") {
  set.seed(20261104)
  contrast(emm, method = "trt.vs.ctrl", ref = ref, adjust = adjust) |>
    summary(infer = TRUE) |>
    as_tibble() |>
    rename(est = all_of(ratio_col), lcl = asymp.LCL, ucl = asymp.UCL) |>
    transmute(survey = as.integer(as.character(survey)), treated = str_remove(as.character(contrast), str_c(" / ", ref)), ref, est, lcl, ucl)
}

limit_nearer_1 <- function(est, lcl, ucl) if_else(est < 1, ucl, lcl)

evalue_rr <- function(rr) {
  r <- if_else(rr < 1, 1 / rr, rr)
  r + sqrt(r * (r - 1))
}

# Point and CI-limit E-values from the EValue package (CI limit 1 when the interval includes 1)
evalue_pkg <- function(est, lcl, ucl, type = c("OR", "RR")) {
  type <- match.arg(type)
  m <- suppressMessages(if (type == "OR") evalues.OR(est, lcl, ucl, rare = FALSE) else evalues.RR(est, lcl, ucl))
  ci <- m["E-values", c("lower", "upper")]
  tibble(point = unname(m["E-values", "point"]), ci = unname(ci[!is.na(ci)][1]))
}

or_routine <- zone_contrasts(emm_prev, "Routine", "mvt") |> filter(treated %in% c("MDA", "IRS"))
or_irs <- zone_contrasts(emm_prev, "IRS", "mvt") |> filter(treated == "MDA")
or_unadj_routine <- zone_contrasts(emm_prev, "Routine", "none")
or_unadj_irs <- zone_contrasts(emm_prev, "IRS", "none") |> filter(treated == "MDA")
or_nobase <-
  bind_rows(zone_contrasts(emm_prev_unadj, "Routine", "mvt") |> filter(treated == "MDA"), zone_contrasts(emm_prev_unadj, "IRS", "mvt") |> filter(treated == "MDA")) |>
  filter(survey == 1)

# Model-based risk ratios and reference risks (check columns)
emm_log <- regrid(emm_prev, transform = "log")
rr_model <-
  bind_rows(zone_contrasts(emm_log, "Routine", "mvt", ratio_col = "ratio"), zone_contrasts(emm_log, "IRS", "mvt", ratio_col = "ratio") |> filter(treated == "MDA")) |>
  transmute(survey, treated, ref, rr_model = est, rr_model_ci_limit = limit_nearer_1(est, lcl, ucl), rr_model_lcl = lcl, rr_model_ucl = ucl)
p0_model <-
  summary(emm_prev) |>
  as_tibble() |>
  transmute(survey = as.integer(as.character(survey)), ref = as.character(zone), p0_reference_model = prob)

adjusted_label <- "adjusted for baseline prevalence (area level), age, sex and displacement status"
nobase_label <- "adjusted for age, sex and displacement status; no baseline covariate"

ors <- bind_rows(
  bind_rows(or_routine, or_irs) |> filter(survey %in% 2:5) |> mutate(row_type = "within-survey contrast", model = adjusted_label),
  bind_rows(or_routine, or_irs) |> filter(survey == 1, treated == "MDA") |> mutate(row_type = "survey-1 benchmark", model = adjusted_label),
  or_nobase |> mutate(row_type = "survey-1 benchmark", model = nobase_label)
)

unadj_limits <-
  bind_rows(or_unadj_routine, or_unadj_irs) |>
  transmute(
    survey, treated, ref, model = adjusted_label,
    ci_limit_or_unadjusted = limit_nearer_1(est, lcl, ucl), ci_unadjusted_includes_1 = lcl <= 1 & ucl >= 1
  )

evalues <-
  ors |>
  left_join(unadj_limits, by = c("survey", "treated", "ref", "model")) |>
  mutate(
    ci_limit_or = limit_nearer_1(est, lcl, ucl),
    ci_includes_1 = lcl <= 1 & ucl >= 1,
    rr_sqrt = sqrt(est),
    rr_sqrt_ci_limit = sqrt(ci_limit_or),
    pkg = pmap(list(est, lcl, ucl), \(e, l, u) evalue_pkg(e, l, u, "OR")),
    evalue_point = map_dbl(pkg, "point"),
    evalue_ci = map_dbl(pkg, "ci"),
    evalue_ci_unadjusted = case_when(
      is.na(ci_limit_or_unadjusted) ~ NA_real_,
      ci_unadjusted_includes_1 ~ 1,
      .default = evalue_rr(sqrt(ci_limit_or_unadjusted))
    )
  ) |>
  left_join(rr_model |> mutate(model = adjusted_label), by = c("survey", "treated", "ref", "model")) |>
  left_join(p0_model |> mutate(model = adjusted_label), by = c("survey", "ref", "model")) |>
  mutate(
    pkg_model = pmap(list(rr_model, rr_model_lcl, rr_model_ucl), \(e, l, u) if (is.na(e)) tibble(point = NA_real_, ci = NA_real_) else evalue_pkg(e, l, u, "RR")),
    evalue_point_model = map_dbl(pkg_model, "point"),
    evalue_ci_model = map_dbl(pkg_model, "ci")
  )

s2_table_f <-
  evalues |>
  left_join(survey_timing |> select(survey_n, months_from_baseline), by = c("survey" = "survey_n")) |>
  mutate(
    comparison = str_c(zone_display_labels[treated], " vs ", zone_display_labels[ref]),
    # Printed in the table: MDA+IRS vs Routine at surveys 2 to 5, MDA+IRS vs IRS at surveys 4 and 5, and the survey 1 benchmarks
    printed = (row_type == "within-survey contrast" & treated == "MDA" & ((ref == "Routine" & survey %in% 2:5) | (ref == "IRS" & survey %in% 4:5))) |
      (row_type == "survey-1 benchmark" & model %in% c(adjusted_label, nobase_label)),
    method = "RR approximated as sqrt(OR) (common outcome); E-value for the point estimate and the multivariate-t 95% CI limit closer to 1"
  ) |>
  transmute(
    row_type, survey, months_from_baseline, comparison, model,
    or = est, or_lcl = lcl, or_ucl = ucl, ci_limit_or, ci_limit_or_unadjusted,
    rr_sqrt, rr_sqrt_ci_limit, evalue_point, evalue_ci, evalue_ci_unadjusted,
    p0_reference_model, rr_model, rr_model_ci_limit, evalue_point_model, evalue_ci_model,
    printed, method
  ) |>
  arrange(
    factor(row_type, levels = c("within-survey contrast", "survey-1 benchmark")),
    factor(comparison, levels = c("MDA+IRS vs Routine", "MDA+IRS vs IRS", "IRS vs Routine")),
    model != adjusted_label, survey
  )

print(s2_table_f |> filter(printed) |> select(row_type, survey, comparison, or, evalue_point, evalue_ci), n = Inf, width = Inf)
write_csv(s2_table_f, "output/tables/s2_table_f_evalues_prevalence.csv")

# S2 Table G ----

# Net use enters as (N1) the area-level survey 1 net-use share, a pre-treatment covariate, and (N2) individual
# net use at surveys 1, 3 and 5. Survey 4's net item is not used and survey 2 did not ask it.
net_any <-
  read_rds("data/df_ind.rds") |>
  filter(rdt_result %in% c("Positif", "Negatif"), !is.na(zone)) |>
  pull(net_any)
stopifnot("net use rows do not line up with the analysis data" = length(net_any) == nrow(dat_rdt))

base_net <- dat_rdt |> mutate(net_any = net_any) |> filter(survey == "1") |> summarise(baseline_net = mean(net_any, na.rm = TRUE), .by = area)
overall_net <- mean(base_net$baseline_net)
p95_net <- quantile(base_net$baseline_net, 0.95)

dat_net <-
  dat_rdt |>
  mutate(net_any = net_any, survey_n = as.integer(as.character(survey))) |>
  left_join(base_net, by = "area") |>
  mutate(
    baseline_net = replace_na(baseline_net, overall_net),
    baseline_net_w = pmin(baseline_net, p95_net),
    baseline_net_z = as.numeric(scale(baseline_net_w))
  )

# The two baseline covariates must be separable for N1 to be fitted
r_area <- dat_net |> distinct(area, baseline_prev_w, baseline_net_w) |> summarise(r = cor(baseline_prev_w, baseline_net_w)) |> pull(r)
print(c(area_level_correlation_baseline_prev_net = r_area))
stopifnot(abs(r_area) <= 0.8)

dat_s135 <- dat_net |> filter(survey_n %in% c(1L, 3L, 5L), !is.na(net_any)) |> mutate(survey = droplevels(survey))

re <- "(1 | cluster_uid/index_hh)"
f_prim <- as.formula(str_c("outcome_rdt ~ zone * survey + age_years + sex + status + baseline_prev_z + ", re))
f_n1 <- as.formula(str_c("outcome_rdt ~ zone * survey + age_years + sex + status + baseline_prev_z + baseline_net_z + ", re))
f_n2net <- as.formula(str_c("outcome_rdt ~ zone * survey + net_any * survey + age_years + sex + status + baseline_prev_z + ", re))

# Each model uses the optimiser chosen by an lme4::allFit() search from default starting values: the
# non-singular fit with the highest log-likelihood. The search takes hours, so the chosen optimisers are fixed here.
nlopt_control <- \(algorithm) glmerControl(optimizer = "nloptwrap", optCtrl = list(algorithm = algorithm, maxeval = 2e5))
nlminb_control <- glmerControl(optimizer = "nlminbwrap", optCtrl = list(eval.max = 2e5))

# The formula and data are embedded in the call so that emmeans() can re-evaluate it
fit_glmer <- \(formula, data, control) do.call(glmer, list(formula = formula, data = data, family = binomial, control = control))

fits <- list(
  prim = fit_glmer(f_prim, dat_net, nlopt_control("NLOPT_LN_BOBYQA")),
  n1 = fit_glmer(f_n1, dat_net, nlminb_control),
  n2_ref = fit_glmer(f_prim, dat_s135, nlopt_control("NLOPT_LN_NELDERMEAD")),
  n2_net = fit_glmer(f_n2net, dat_s135, nlopt_control("NLOPT_LN_BOBYQA"))
)

net_zone_contrasts <- function(fit, d) {
  contrast(emmeans(fit, ~ zone | survey, data = d), method = "trt.vs.ctrl", ref = "Routine", adjust = "mvt") |>
    summary(infer = TRUE) |>
    as_tibble() |>
    transmute(
      survey_n = as.integer(as.character(survey)), zone = str_remove(as.character(contrast), " - Routine"),
      log_or = estimate, or = exp(estimate), lcl = exp(asymp.LCL), ucl = exp(asymp.UCL), p_value = p.value
    )
}

# Adjusted prevalence difference (percentage points) at the covariate means, random effects at zero
net_prev_diff <- function(fit, d) {
  contrast(regrid(emmeans(fit, ~ zone | survey, data = d)), method = "trt.vs.ctrl", ref = "Routine", adjust = "mvt") |>
    summary(infer = TRUE) |>
    as_tibble() |>
    transmute(
      survey_n = as.integer(as.character(survey)), zone = str_remove(as.character(contrast), " - Routine"),
      pd_pp = 100 * estimate, pd_lcl = 100 * asymp.LCL, pd_ucl = 100 * asymp.UCL
    )
}

# Ratios are formed within one fitting family: N1 against the primary specification refitted here, and N2-net
# against N2-ref. The primary model's own OR is shown for reference.
prim_bobyqa <- contrasts_prev |> transmute(survey_n, zone, or_primary_bobyqa = or)
prim_selected <- net_zone_contrasts(fits$prim, dat_net) |> transmute(survey_n, zone, or_primary_selected = or)

n1_tbl <-
  net_zone_contrasts(fits$n1, dat_net) |>
  left_join(prim_bobyqa, by = c("survey_n", "zone")) |>
  left_join(prim_selected, by = c("survey_n", "zone")) |>
  transmute(
    analysis = "N1", survey_n, zone, or_primary_bobyqa, or_primary_selected,
    or_ref = NA_real_, lcl_ref = NA_real_, ucl_ref = NA_real_,
    or_adjusted = or, lcl_adjusted = lcl, ucl_adjusted = ucl, p_value_adjusted = p_value,
    ratio_n1_to_primary_selected = or / or_primary_selected, ratio_adjusted_to_ref = NA_real_,
    pd_ref_pp = NA_real_, pd_ref_lcl = NA_real_, pd_ref_ucl = NA_real_,
    pd_adjusted_pp = NA_real_, pd_adjusted_lcl = NA_real_, pd_adjusted_ucl = NA_real_
  )

n2_ref_zc <- net_zone_contrasts(fits$n2_ref, dat_s135)
n2_tbl <-
  net_zone_contrasts(fits$n2_net, dat_s135) |>
  rename(or_adjusted = or, lcl_adjusted = lcl, ucl_adjusted = ucl, p_value_adjusted = p_value) |>
  select(-log_or) |>
  left_join(n2_ref_zc |> select(survey_n, zone, or_ref = or, lcl_ref = lcl, ucl_ref = ucl), by = c("survey_n", "zone")) |>
  left_join(prim_bobyqa, by = c("survey_n", "zone")) |>
  left_join(prim_selected, by = c("survey_n", "zone")) |>
  left_join(net_prev_diff(fits$n2_ref, dat_s135) |> rename(pd_ref_pp = pd_pp, pd_ref_lcl = pd_lcl, pd_ref_ucl = pd_ucl), by = c("survey_n", "zone")) |>
  left_join(net_prev_diff(fits$n2_net, dat_s135) |> rename(pd_adjusted_pp = pd_pp, pd_adjusted_lcl = pd_lcl, pd_adjusted_ucl = pd_ucl), by = c("survey_n", "zone")) |>
  mutate(analysis = "N2", ratio_n1_to_primary_selected = NA_real_, ratio_adjusted_to_ref = or_adjusted / or_ref) |>
  select(any_of(names(n1_tbl)))

s2_table_g <-
  bind_rows(n1_tbl, n2_tbl) |>
  mutate(
    months_from_baseline = survey_timing$months_from_baseline[match(survey_n, survey_timing$survey_n)],
    zone = unname(zone_display_labels[zone]),
    ref_vs_primary_selected_within_5pct = if_else(analysis == "N2", abs(or_ref / or_primary_selected - 1) <= 0.05, NA),
    note = case_when(
      survey_n == 1L & analysis == "N2" ~ "Survey 1: or_primary_selected and or_primary_bobyqa are the adjusted primary-model survey-1 zone ORs (residual survey-1 contrast after conditioning on baseline prevalence, individual covariates and the random effects), not crude baseline contrasts.",
      survey_n == 1L & analysis == "N1" ~ "Survey 1: residual survey-1 contrast, adjusted as in the primary model.",
      .default = ""
    )
  ) |>
  relocate(months_from_baseline, .after = survey_n) |>
  arrange(analysis, survey_n, factor(zone, levels = c("MDA+IRS", "IRS"))) |>
  mutate(across(where(is.double), \(x) round(x, 4)))

print(s2_table_g |> select(analysis, survey_n, zone, or_primary_selected, or_ref, or_adjusted, ratio_n1_to_primary_selected, ratio_adjusted_to_ref), n = Inf, width = Inf)
write_csv(s2_table_g, "output/tables/s2_table_g_net_adjusted_prevalence.csv")

# S12 Table ----

# Net use and household net ownership by zone, survey-weighted, surveys 1, 3 and 5. Survey 2 did not ask the
# net-use question; the survey 4 item is not reported.
df_ind <- read_rds("data/df_ind.rds")

pop <-
  df_ind |>
  filter(survey %in% c(1L, 3L, 5L), !is.na(net_any)) |>
  mutate(
    survey_n = survey,
    intact = if_else(is.na(net_cat), NA_integer_, as.integer(net_cat == "intact")),
    u5 = !is.na(age) & age < 5
  )
pop_rdt <- dat_net |> filter(survey_n %in% c(1L, 3L, 5L), !is.na(net_any))

row_defs <- list(
  list(label = "Slept under any net last night, all current members", data = \(d) d, var = "net_any", src = "pop"),
  list(label = "Slept under an intact net last night, all current members", data = \(d) d |> filter(!is.na(intact)), var = "intact", src = "pop"),
  list(label = "Slept under any net last night, children under 5 years", data = \(d) d |> filter(u5), var = "net_any", src = "pop"),
  list(label = "Slept under any net last night, RDT-tested participants", data = \(d) d, var = "net_any", src = "rdt"),
  list(label = "Lives in a household owning at least one net", data = \(d) d |> filter(!is.na(hh_own_net)), var = "hh_own_net", src = "pop")
)

net_estimate <- function(d, var, s, z) {
  d <- d |> filter(survey_n == s, zone == z, !is.na(cluster), !is.na(weight), !is.na(.data[[var]]))
  if (nrow(d) == 0L) return(NULL)
  d$y <- d[[var]]
  est <- calculate_prop(d, key = y, value = 1)
  tibble(n = nrow(d), n_yes = sum(d$y == 1L), n_clusters = n_distinct(d$cluster), pct = est$prop, lcl = est$prop_low, ucl = est$prop_upp, deff = est$prop_deff)
}

net_long <-
  expand_grid(def_i = seq_along(row_defs), survey_n = c(1L, 3L, 5L), zone = c("Routine", "MDA", "IRS")) |>
  mutate(res = pmap(list(def_i, survey_n, zone), \(i, s, z) {
    def <- row_defs[[i]]
    net_estimate(def$data(if (def$src == "rdt") pop_rdt else pop), def$var, s, z)
  })) |>
  filter(!map_lgl(res, is.null)) |>
  unnest(res) |>
  mutate(
    measure = map_chr(def_i, \(i) row_defs[[i]]$label),
    months_from_baseline = survey_timing$months_from_baseline[match(survey_n, survey_timing$survey_n)],
    zone_label = unname(zone_display_labels[zone])
  ) |>
  arrange(def_i, survey_n, factor(zone, levels = zone_display_levels))

s12_long <-
  net_long |>
  transmute(
    survey = survey_n, months_from_baseline, zone = zone_label, measure, n, n_yes, n_clusters,
    pct = round(pct, 2), lcl = round(lcl, 2), ucl = round(ucl, 2), deff = round(deff, 2)
  )
write_csv(s12_long, "output/tables/s12_table_net_use_by_zone.csv")

# Crossover-excluded areas (no zone) per survey, for the table note
excluded_areas <-
  df_ind |>
  filter(survey %in% c(1L, 3L, 5L), is.na(zone), !is.na(area)) |>
  distinct(survey, area) |>
  summarise(areas = str_c(sort(area), collapse = ", "), .by = survey)

fmt_cell <- \(pct, lcl, ucl, n) str_c(sprintf("%.1f", pct), " (", sprintf("%.1f", lcl), "-", sprintf("%.1f", ucl), "), n = ", n)

s12_wide <-
  net_long |>
  mutate(cell = fmt_cell(pct, lcl, ucl, n)) |>
  select(def_i, measure, survey_n, months_from_baseline, zone_label, cell) |>
  pivot_wider(names_from = zone_label, values_from = cell) |>
  bind_rows(
    expand_grid(def_i = seq_along(row_defs), survey_n = c(2L, 4L)) |>
      mutate(
        measure = map_chr(def_i, \(i) row_defs[[i]]$label),
        months_from_baseline = survey_timing$months_from_baseline[match(survey_n, survey_timing$survey_n)],
        # Survey 4 has no household net-ownership item
        Routine = if_else(survey_n == 2L | (survey_n == 4L & measure == "Lives in a household owning at least one net"), "Not collected", "Not reported"),
        `MDA+IRS` = Routine, IRS = Routine
      )
  ) |>
  arrange(def_i, survey_n) |>
  left_join(excluded_areas |> rename(survey_n = survey, crossover_areas_excluded = areas), by = "survey_n") |>
  transmute(measure, survey = survey_n, months_from_baseline, Routine, `MDA+IRS`, IRS, crossover_areas_excluded = replace_na(crossover_areas_excluded, ""))

print(s12_wide, n = Inf, width = Inf)
write_csv(s12_wide, "output/tables/s12_table_net_use_by_zone_wide.csv")
