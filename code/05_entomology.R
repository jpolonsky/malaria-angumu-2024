# Entomological surveillance: Table 3 (zone x half-year indicators), Figure 5 (estimated HBR by health
# area), S3 Figure (monthly indicators by zone), the entomology numbers quoted in S1 Text and the
# S5 Table workbook (Additional file 2).

pacman::p_load(patchwork, ggrepel, glmmTMB, openxlsx, tidyverse)

source("code/functions.R")

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# Helpers ----

ento_metrics <- c("hbr", "density", "sr", "eir_annual", "hbr_person", "eir_person_annual")
ento_spz_start <- as.Date("2023-01-01")

# Health-area populations (the weights of the zone means)
ento_pop <-
  tribble(
    ~area, ~pop,
    "Abia", 16835,
    "Anzika", 13322,
    "Apala", 17843,
    "Ara", 19818,
    "Awasi", 15616,
    "Cawa", 21866,
    "Gengere", 25172,
    "Lanyi", 13560,
    "Linga", 7871,
    "Pakala", 11193,
    "Panyadong", 8810,
    "Udongo Abira", 10032,
    "Ugudo-Zii", 22166,
    "Uyandu", 14128
  )

# Month from which a Routine area leaves the zone comparison (first intervention received); Awasi, the
# last Routine area, leaves at the end of the trapping window
ento_crossover_date <-
  study_areas |>
  filter(!is.na(crossover_date)) |>
  select(area, out_from = crossover_date) |>
  bind_rows(tibble(area = "Awasi", out_from = as.Date("2024-12-01")))

# Half-year periods, each containing one survey
ento_periods <-
  tribble(
    ~period, ~p_start, ~p_end, ~survey_in_period,
    "P0", as.Date("2022-05-01"), as.Date("2022-06-30"), NA_real_,
    "P1", as.Date("2022-07-01"), as.Date("2022-12-31"), 2,
    "P2", as.Date("2023-01-01"), as.Date("2023-06-30"), 3,
    "P3", as.Date("2023-07-01"), as.Date("2023-12-31"), 4,
    "P4", as.Date("2024-01-01"), as.Date("2024-06-30"), 5,
    "P5", as.Date("2024-07-01"), as.Date("2024-11-30"), NA_real_
  )
irs_round2_end <- as.Date("2023-05-31")

# Estimation by bootstrap within health area. Health areas are fixed strata and trapping sites the
# resampling unit; a resampled site carries all its trap-nights and mosquitoes across periods.
# Each resample recomputes, per health area, zone (as treated) and period or month: estimated HBR
# (Anopheles per working trap-night), density (random houses only), sporozoite rate (positives /
# tested), annual EIR (HBR x sporozoite rate x 365) and the per-occupant HBR and EIR, then averages
# the areas of a zone with population weights. A zone-area-period cell enters only when its area
# was trapped in at least `min_months` calendar months of the period (P0 is exempt).
ento_fit <- function(tn = ento_trapnights, mq = ento_mosq, unit = c("period", "month"), cluster = c("site", "house"), min_months = 3, n_boot = 2000, seed = 20261007, pop = ento_pop) {
  unit <- match.arg(unit)
  cluster <- match.arg(cluster)
  tn <- tn |> mutate(u = as.character(.data[[unit]]), unit_id = if (cluster == "site") site_id else house_id)
  mq <- mq |> mutate(u = as.character(.data[[unit]]), unit_id = if (cluster == "site") site_id else house_id)
  keys <- c("area", "unit_id", "zone", "u")

  months_cov <-
    tn |>
    filter(!is.na(zone), !is.na(u), working, !is.na(n_anoph)) |>
    summarise(n_months = n_distinct(month), .by = c(area, zone, u))
  thin_cells <- if (unit == "period") months_cov |> filter(n_months < min_months, u != "P0") else months_cov[0, ]

  th <-
    tn |>
    filter(!is.na(zone), !is.na(u), working, !is.na(n_anoph)) |>
    summarise(
      nights = n(), anoph = sum(n_anoph), nights_r = sum(random_house), anoph_r = sum(n_anoph[random_house]),
      pnights = sum(persons, na.rm = TRUE), anoph_p = sum(n_anoph[!is.na(persons)]), .by = all_of(keys)
    )
  mh <-
    mq |>
    filter(!is.na(zone), !is.na(u)) |>
    mutate(is_tested = tested %in% TRUE, is_pos = is_tested & spz_result %in% "Positive") |>
    summarise(n_tested = sum(is_tested), n_pos = sum(is_pos), .by = all_of(keys))
  hh <-
    full_join(th, mh, by = keys) |>
    mutate(across(c(nights, anoph, nights_r, anoph_r, pnights, anoph_p, n_tested, n_pos), \(x) replace_na(x, 0))) |>
    anti_join(thin_cells, by = c("area", "zone", "u")) |>
    mutate(
      hkey = paste(area, unit_id),
      ckey = paste(area, zone, u, sep = "|"),
      zkey = paste(zone, u, sep = "|")
    )
  houses <- distinct(hh, area, hkey)
  cells <- hh |> distinct(ckey, area, zone, u, zkey) |> mutate(pop = pop$pop[match(area, pop$area)])
  zus <- cells |> distinct(zkey, zone, u) |> arrange(match(zone, zone_display_levels), u)

  hi <- match(hh$hkey, houses$hkey)
  ci <- match(hh$ckey, cells$ckey)
  vals <- as.matrix(hh[c("nights", "anoph", "nights_r", "anoph_r", "n_tested", "n_pos", "pnights", "anoph_p")])
  zi_of_cell <- match(cells$zkey, zus$zkey)
  area_house_idx <- split(seq_len(nrow(houses)), houses$area)

  cell_metrics <- function(s) {
    hbr <- if_else(s[, 1] > 0, s[, 2] / s[, 1], NA_real_)
    density <- if_else(s[, 3] > 0, s[, 4] / s[, 3], NA_real_)
    sr <- if_else(s[, 5] > 0, s[, 6] / s[, 5], NA_real_)
    hbr_person <- if_else(s[, 7] > 0, s[, 8] / s[, 7], NA_real_)
    list(hbr = hbr, density = density, sr = sr, eir_annual = hbr * sr * 365, hbr_person = hbr_person, eir_person_annual = hbr_person * sr * 365)
  }
  zone_mean <- function(x) {
    ok <- !is.na(x)
    num <- rowsum(if_else(ok, cells$pop * x, 0), zi_of_cell, reorder = TRUE)[, 1]
    den <- rowsum(if_else(ok, cells$pop, 0), zi_of_cell, reorder = TRUE)[, 1]
    if_else(den > 0, num / den, NA_real_)
  }

  # Within-zone change contrasts (P3 and P5 against P1, same season, and against P2): ratio of
  # population-weighted means over the areas present in both periods. Sporozoite rate and EIR
  # have no P1 (no usable results in 2022).
  chg_grid <- tibble()
  if (unit == "period") {
    chg_grid <-
      expand_grid(zone = zone_display_levels, p = c("P3", "P5"), ref = c("P1", "P2")) |>
      mutate(
        both = pmap(list(zone, p, ref), \(z, p, r) intersect(cells$area[cells$zone == z & cells$u == p], cells$area[cells$zone == z & cells$u == r])),
        ip = pmap(list(zone, p, both), \(z, p, a) match(paste(a, z, p, sep = "|"), cells$ckey)),
        ir = pmap(list(zone, ref, both), \(z, r, a) match(paste(a, z, r, sep = "|"), cells$ckey)),
        metric = map(ref, \(r) if (r == "P1") c("hbr", "density", "hbr_person") else ento_metrics)
      ) |>
      filter(map_lgl(both, \(a) length(a) > 0)) |>
      unnest(metric) |>
      arrange(match(metric, ento_metrics), ref, match(zone, zone_display_levels), p)
    chg_grid$routine_row <- match(paste(chg_grid$metric, chg_grid$p, chg_grid$ref, "Routine"), paste(chg_grid$metric, chg_grid$p, chg_grid$ref, chg_grid$zone))
  }
  ratio_fun <- function(x, ip, ir) {
    a <- x[ip]
    b <- x[ir]
    w <- cells$pop[ip]
    ok <- !is.na(a) & !is.na(b)
    sum(w[ok] * a[ok]) / sum(w[ok] * b[ok])
  }

  pair_defs <- tibble(a = c("MDA", "IRS", "MDA"), b = c("Routine", "Routine", "IRS"))

  one_pass <- function(mult) {
    s <- rowsum(vals * mult[hi], ci, reorder = TRUE)
    cm <- cell_metrics(s)
    est <- map(cm, zone_mean)
    out_est <- unlist(est, use.names = FALSE)
    out_chg <- numeric(0)
    out_rel <- numeric(0)
    if (nrow(chg_grid) > 0) {
      out_chg <- pmap_dbl(list(chg_grid$metric, chg_grid$ip, chg_grid$ir), \(m, ip, ir) ratio_fun(cm[[m]], ip, ir))
      out_rel <- out_chg / out_chg[chg_grid$routine_row]
    }
    out_btw <- numeric(0)
    for (m in ento_metrics) {
      e <- est[[m]]
      for (k in seq_len(nrow(pair_defs))) {
        for (uu in unique(zus$u)) {
          ia <- which(zus$zone == pair_defs$a[k] & zus$u == uu)
          ib <- which(zus$zone == pair_defs$b[k] & zus$u == uu)
          out_btw <- c(out_btw, if (length(ia) == 1 && length(ib) == 1) e[ia] / e[ib] else NA_real_)
        }
      }
    }
    c(out_est, out_chg, out_rel, out_btw)
  }

  est_key <- expand_grid(metric = ento_metrics, zus) |> transmute(type = "estimate", metric, zone, zone2 = NA_character_, u, ref = NA_character_)
  chg_key <- if (nrow(chg_grid) > 0) chg_grid |> transmute(type = "change_vs_ref", metric, zone, zone2 = NA_character_, u = p, ref) else tibble()
  rel_key <- if (nrow(chg_grid) > 0) chg_key |> mutate(type = "change_relative_to_routine") else tibble()
  btw_key <- expand_grid(metric = ento_metrics, pair_defs |> mutate(k = row_number()), u = unique(zus$u)) |> transmute(type = "between_zone_ratio", metric, zone = a, zone2 = b, u, ref = NA_character_)
  key <- bind_rows(est_key, chg_key, rel_key, btw_key)

  point <- one_pass(rep(1, nrow(houses)))
  stopifnot(length(point) == nrow(key))

  set.seed(seed)
  boot <- matrix(NA_real_, nrow = n_boot, ncol = length(point))
  for (b in seq_len(n_boot)) {
    mult <- numeric(nrow(houses))
    for (idx in area_house_idx) mult[idx] <- tabulate(sample.int(length(idx), length(idx), replace = TRUE), nbins = length(idx))
    boot[b, ] <- one_pass(mult)
  }
  # 0/0 is undefined (dropped); x/0 is +Inf and stays, so an unbounded upper tail is reported as such
  boot[is.nan(boot) | is.na(boot) | boot == -Inf] <- NA_real_
  boot_q <- \(x, p) {
    x <- x[!is.na(x)]
    if (length(x) < 0.9 * n_boot) return(NA_real_)
    quantile(x, p, names = FALSE, type = if (any(is.infinite(x))) 1 else 7)
  }

  res <-
    key |>
    mutate(
      estimate = point,
      lcl = apply(boot, 2, boot_q, p = 0.025),
      ucl = apply(boot, 2, boot_q, p = 0.975),
      n_boot_valid = colSums(!is.na(boot))
    )

  # Unweighted counts behind each zone-unit cell, from the unresampled data
  ok_tn <-
    tn |>
    filter(!is.na(zone), !is.na(u), working, !is.na(n_anoph)) |>
    mutate(ckey = paste(area, zone, u, sep = "|")) |>
    filter(ckey %in% cells$ckey)
  counts <-
    ok_tn |>
    summarise(
      areas = n_distinct(area), areas_included = str_c(sort(unique(area)), collapse = ", "),
      sites = n_distinct(site_id), houses = n_distinct(house_id), trap_nights = n(), anopheles_working_traps = sum(n_anoph),
      trap_nights_random = sum(random_house), .by = c(zone, u)
    ) |>
    left_join(hh |> summarise(mosquitoes_tested = sum(n_tested), sporozoite_positives = sum(n_pos), .by = c(zone, u)), by = c("zone", "u")) |>
    mutate(across(c(mosquitoes_tested, sporozoite_positives), \(x) if_else(x == 0 & !any(mq$tested %in% TRUE), NA_real_, x)))

  list(
    estimates = res |> filter(type == "estimate") |> select(metric, zone, u, estimate, lcl, ucl, n_boot_valid),
    change = res |> filter(type %in% c("change_vs_ref", "change_relative_to_routine")) |> select(type, metric, zone, u, ref, estimate, lcl, ucl, n_boot_valid),
    between = res |> filter(type == "between_zone_ratio") |> select(metric, zone, zone2, u, estimate, lcl, ucl, n_boot_valid),
    counts = counts
  )
}

# Zone x period (or month) estimates with their unweighted counts, one wide row per cell
fit_to_long <- function(fit, periods = paste0("P", 0:5)) {
  fit$estimates |>
    filter(u %in% periods) |>
    select(-n_boot_valid) |>
    pivot_wider(names_from = metric, values_from = c(estimate, lcl, ucl)) |>
    left_join(fit$counts |> rename(period_u = u), by = c("zone", "u" = "period_u")) |>
    rename(period = u)
}

# Exact (Clopper-Pearson) interval limits for x positives of n
exact_lcl <- \(x, n) if_else(n > 0, qbeta(0.025, x, n - x + 1), NA_real_)
exact_ucl <- \(x, n) if_else(n > 0, qbeta(0.975, x + 1, n - x), NA_real_)

# "est (lcl-ucl)" with `digits` decimals; `signif = TRUE` prints `digits` significant figures
fmt_est_ci <- \(est, lcl, ucl, scale = 1, digits = 2, signif = FALSE) {
  fmt_one <- \(x) if (signif) str_remove(formatC(signif(scale * x, digits), digits = digits, format = "fg", flag = "#"), "\\.$") else sprintf(str_c("%.", digits, "f"), scale * x)
  if_else(is.na(est), NA_character_, str_c(map_chr(est, fmt_one), " (", map_chr(lcl, fmt_one), "-", map_chr(ucl, fmt_one), ")"))
}

# Data ----

ento_trapnights <- read_rds("data/entomology_trapnights.rds")
ento_mosq <- read_rds("data/entomology_mosquitoes.rds")

n_boot_table <- 2000
n_boot_monthly <- 500

# Table 3 ----

months_between <- \(a, b) as.integer(round(as.numeric(b - a) / 30.4375))
period_lab <-
  ento_periods |>
  filter(period %in% paste0("P", 1:5)) |>
  mutate(
    months_after_irs2 = case_when(
      p_end <= irs_round2_end ~ NA_character_,
      p_start <= irs_round2_end ~ "during/at end",
      .default = str_c(months_between(irs_round2_end, p_start), "-", months_between(irs_round2_end, p_end))
    ),
    dates = str_c(format(p_start, "%b %Y"), " to ", format(p_end, "%b %Y"))
  ) |>
  select(period, dates, survey_in_period, months_after_irs2)

# Sporozoite rate: positives / tested mosquitoes pooled over the cell, exact interval. HBR, EIR and
# density: population-weighted means over the health areas of the cell, bootstrap percentile intervals.
table3_rows <- function(long) {
  long |>
    filter(period != "P0") |>
    left_join(period_lab, by = "period") |>
    mutate(zone_label = as.character(zone_display(zone))) |>
    arrange(period, match(zone, zone_display_levels)) |>
    transmute(
      period, dates, survey = survey_in_period, months_after_irs_round2 = months_after_irs2,
      zone = zone_label, health_areas = areas, areas_included, sites, houses, working_trap_nights = trap_nights, anopheles_working_traps,
      mosquitoes_tested, sporozoite_positives,
      sporozoite_rate_pooled = if_else(mosquitoes_tested > 0, sporozoite_positives / mosquitoes_tested, NA_real_),
      sporozoite_pooled_lcl = exact_lcl(sporozoite_positives, mosquitoes_tested),
      sporozoite_pooled_ucl = exact_ucl(sporozoite_positives, mosquitoes_tested),
      hbr = estimate_hbr, hbr_lcl = lcl_hbr, hbr_ucl = ucl_hbr,
      eir_annual = estimate_eir_annual, eir_lcl = lcl_eir_annual, eir_ucl = ucl_eir_annual,
      density_random_houses = estimate_density, density_lcl = lcl_density, density_ucl = ucl_density,
      hbr_per_occupant = estimate_hbr_person, hbr_per_occupant_lcl = lcl_hbr_person, hbr_per_occupant_ucl = ucl_hbr_person,
      eir_per_occupant_annual = estimate_eir_person_annual, eir_per_occupant_lcl = lcl_eir_person_annual, eir_per_occupant_ucl = ucl_eir_person_annual,
      sporozoite_rate_popweighted_sensitivity = estimate_sr, sporozoite_popweighted_lcl = lcl_sr, sporozoite_popweighted_ucl = ucl_sr
    )
}

fit_period <- ento_fit(unit = "period", n_boot = n_boot_table)
table3_long <- fit_to_long(fit_period)
table3_full <- table3_rows(table3_long)

# Main-text version: zone rows only, three measures; EIR is blank wherever the sporozoite rate is
table3 <-
  table3_full |>
  transmute(
    Period = dates,
    `Months after IRS round 2` = months_after_irs_round2,
    Zone = zone,
    `Working trap-nights` = format(working_trap_nights, big.mark = ",", trim = TRUE),
    `Anopheles in working traps` = format(anopheles_working_traps, big.mark = ",", trim = TRUE),
    `Mosquitoes tested (positive)` = if_else(mosquitoes_tested > 0, str_c(format(mosquitoes_tested, big.mark = ",", trim = TRUE), " (", format(sporozoite_positives, big.mark = ",", trim = TRUE), ")"), NA_character_),
    `Sporozoite rate % (95% CI)` = fmt_est_ci(sporozoite_rate_pooled, sporozoite_pooled_lcl, sporozoite_pooled_ucl, scale = 100, digits = 1),
    `HBR (95% CI)` = fmt_est_ci(hbr, hbr_lcl, hbr_ucl),
    `EIR per year (95% CI)` = fmt_est_ci(eir_annual, eir_lcl, eir_ucl, digits = 2, signif = TRUE)
  )
write_csv(table3, "output/tables/table3_entomology.csv")
print(table3, width = Inf)

# Change over time and between-zone ratios (bootstrap CIs), quoted in the Results and S5 Table sheet E
table3_contrasts <-
  bind_rows(
    fit_period$change |> mutate(comparison = type, zone2 = NA_character_),
    fit_period$between |> mutate(comparison = "between_zone_ratio", ref = NA_character_)
  ) |>
  filter(!(comparison == "change_relative_to_routine" & zone == "Routine")) |>
  mutate(
    zone = as.character(zone_display(zone)),
    zone2 = if_else(is.na(zone2), NA_character_, as.character(zone_display(zone2)))
  ) |>
  select(comparison, metric, zone, zone2, period = u, reference_period = ref, estimate, lcl, ucl, n_boot_valid)
write_csv(table3_contrasts, "output/tables/table3_entomology_contrasts.csv")

# MDA+IRS split into Gengere and the other four areas: same engine and seed on the subset of areas
mda_areas <- study_areas$area[study_areas$baseline_zone == "MDA"]
mda_groups <- list("Gengere" = "Gengere", "Other four areas" = setdiff(mda_areas, "Gengere"))
fit_mda_groups <-
  map(mda_groups, \(a) ento_fit(tn = ento_trapnights |> filter(area %in% a), mq = ento_mosq |> filter(area %in% a), unit = "period", n_boot = n_boot_table))

mda_split <-
  imap(fit_mda_groups, \(f, g) {
    f |>
      fit_to_long() |>
      filter(zone == "MDA") |>
      table3_rows() |>
      mutate(group = g, .after = zone)
  }) |>
  bind_rows() |>
  bind_rows(table3_full |> filter(zone == "MDA+IRS") |> mutate(group = "All MDA+IRS areas", .after = zone)) |>
  arrange(period, match(group, c("All MDA+IRS areas", "Gengere", "Other four areas")))

# P0 (about 4 weeks before IRS round 1; not a usable baseline), S5 Table sheet E
table3_p0 <- table3_long |> filter(period == "P0") |> mutate(zone = as.character(zone_display(zone)))

# Area-level HBR and density for every area; cells with fewer than three trapped months are
# kept and flagged
by_area <-
  ento_trapnights |>
  filter(!is.na(zone), period %in% paste0("P", 1:5), working, !is.na(n_anoph)) |>
  summarise(
    trap_nights = n(), trapped_months = n_distinct(month), anopheles = sum(n_anoph),
    hbr = sum(n_anoph) / n(), density = sum(n_anoph[random_house]) / sum(random_house),
    .by = c(zone, area, period)
  ) |>
  left_join(
    ento_mosq |> filter(!is.na(zone), period %in% paste0("P", 1:5)) |> summarise(mosquitoes = n(), tested = sum(tested), positives = sum(spz_result %in% "Positive"), .by = c(zone, area, period)),
    by = c("zone", "area", "period")
  ) |>
  mutate(
    sporozoite_rate = if_else(tested > 0, positives / tested, NA_real_),
    eir_annual = hbr * sporozoite_rate * 365,
    tested_share = tested / mosquitoes,
    zone = as.character(zone_display(zone)),
    in_zone_estimate = trapped_months >= 3
  ) |>
  arrange(match(zone, unname(zone_display_labels[zone_display_levels])), area, period)

# Area-level HBR ratios (P3 and P5 against P1), point estimates
area_change <-
  by_area |>
  select(zone, area, period, hbr, trapped_months) |>
  pivot_wider(names_from = period, values_from = c(hbr, trapped_months)) |>
  mutate(ratio_p3_p1 = hbr_P3 / hbr_P1, ratio_p5_p1 = hbr_P5 / hbr_P1)

# Figure 5 ----

x_limits <- c(as.Date("2022-06-01"), as.Date("2024-12-01"))
x_scale <- scale_x_date(limits = x_limits, date_breaks = "3 months", date_labels = "%b\n%Y", expand = expansion(mult = 0.01))
theme_ento <- theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank())
zone_cols_lab <- set_names(zone_cols[zone_display_levels], unname(zone_display_labels[zone_display_levels]))

# Rainfall: PNLP weekly series for the whole Angumu health zone (same value for every zone)
rain <-
  read_rds("data/climate_weekly.rds") |>
  mutate(
    jan4 = make_date(year, 1, 4),
    week_start = jan4 - (wday(jan4, week_start = 1) - 1) + 7 * (epi_week - 1)
  ) |>
  filter(week_start >= x_limits[1] - 7, week_start <= x_limits[2])

# Intervention calendar
calendar <-
  tribble(
    ~row, ~label, ~start, ~end,
    "IRS", "IRS round 1", as.Date("2022-06-27"), as.Date("2022-11-10"),
    "IRS", "IRS round 2", as.Date("2023-02-02"), as.Date("2023-05-31"),
    "MDA", "MDA rounds 1-2", as.Date("2022-07-01"), as.Date("2022-08-31"),
    "MDA", "MDA rounds 3-4", as.Date("2023-01-01"), as.Date("2023-02-28"),
    "ITN", "ITN distribution", as.Date("2023-12-01"), as.Date("2023-12-31"),
    "Survey", "Survey 2", as.Date("2022-11-01"), as.Date("2022-12-15"),
    "Survey", "Survey 3", as.Date("2023-04-01"), as.Date("2023-05-15"),
    "Survey", "Survey 4", as.Date("2023-10-01"), as.Date("2023-10-31"),
    "Survey", "Survey 5", as.Date("2024-04-01"), as.Date("2024-05-15")
  ) |>
  mutate(row = factor(row, levels = c("Survey", "ITN", "MDA", "IRS")))

# With crossover = TRUE a fifth row marks the dates on which Routine areas leave the comparison
crossover_events <-
  ento_crossover_date |>
  summarise(label = str_c(sort(area), collapse = ", "), .by = out_from) |>
  arrange(out_from)

make_cal_panel <- function(crossover = FALSE) {
  n_rows <- if (crossover) 5 else 4
  cal <- calendar |> mutate(row = factor(as.character(row), levels = c(levels(calendar$row), if (crossover) "Leaves\nRoutine")))
  p <-
    ggplot(cal) +
    geom_rect(aes(xmin = start, xmax = end, ymin = as.numeric(row) - 0.35, ymax = as.numeric(row) + 0.35), fill = "grey55") +
    geom_vline(xintercept = as.Date("2024-01-01"), linetype = "dashed", colour = "grey30") +
    annotate("text", x = as.Date("2024-01-08"), y = 4.45, label = "Routine zone = Awasi only", hjust = 0, size = 2.6, colour = "grey30")
  if (crossover) {
    p <-
      p +
      geom_point(data = crossover_events, aes(out_from, 5), shape = 25, size = 2.6, colour = zone_cols[["Routine"]], fill = zone_cols[["Routine"]]) +
      geom_text(
        # labels from January 2024 are right-aligned and shifted 8 days left to stay clear of the dashed line
        data = crossover_events |> mutate(hj = if_else(out_from >= as.Date("2024-01-01"), 1, 0.5), x_lab = out_from - if_else(hj == 1, 8, 0)),
        aes(x_lab, 5.55, label = label, hjust = hj), size = 2.6, colour = "grey30", vjust = 0
      )
  }
  p +
    scale_y_continuous(breaks = seq_len(n_rows), labels = levels(cal$row)[seq_len(n_rows)], limits = c(0.5, if (crossover) 6.2 else 4.7)) +
    x_scale + labs(y = NULL) + theme_ento +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.y = element_blank())
}
p_cal <- make_cal_panel(crossover = FALSE)
p_cal_s3 <- make_cal_panel(crossover = TRUE)

p_rain <-
  ggplot(rain, aes(week_start, precipitation)) +
  geom_col(width = 6, fill = "grey70", na.rm = TRUE) +
  x_scale + labs(y = "Rainfall\n(PNLP weekly series, mm)") + theme_ento +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# One line per health area on calendar quarters, one panel per assigned zone. Gengere is drawn
# heavier; open marks: fewer than 30 working trap-nights in the quarter; dashed segments: after the
# area received its first intervention. Pseudo-log axis keeps zero-catch quarters on the plot.
q_pts <-
  ento_trapnights |>
  filter(working, !is.na(n_anoph), date >= as.Date("2022-06-01"), date <= as.Date("2024-11-30")) |>
  mutate(q_start = floor_date(date, "quarter")) |>
  summarise(trap_nights = n(), months = n_distinct(month), hbr = sum(n_anoph) / n(), .by = c(zone_assigned, area, q_start)) |>
  filter(months >= 2) |>
  left_join(ento_crossover_date, by = "area") |>
  mutate(
    x = q_start + 45,
    post = !is.na(out_from) & x >= out_from,
    zone = factor(as.character(zone_display(zone_assigned)), levels = unname(zone_display_labels[zone_display_levels]))
  ) |>
  arrange(area, x)

q_seg <-
  q_pts |>
  group_by(area) |>
  mutate(x0 = lag(x), y0 = lag(hbr)) |>
  ungroup() |>
  filter(!is.na(x0))
q_last <- q_pts |> group_by(area) |> slice_max(x, n = 1) |> ungroup()

p_area <-
  ggplot() +
  annotate("rect", xmin = as.Date("2022-06-27"), xmax = as.Date("2022-11-10"), ymin = -Inf, ymax = Inf, fill = "grey93") +
  annotate("rect", xmin = as.Date("2023-02-02"), xmax = as.Date("2023-05-31"), ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_segment(data = q_seg, aes(x = x0, xend = x, y = y0, yend = hbr, colour = zone, linetype = post, linewidth = area == "Gengere")) +
  geom_point(data = q_pts, aes(x, hbr, colour = zone, fill = if_else(trap_nights >= 30, as.character(zone), NA_character_), size = area == "Gengere"), shape = 21, stroke = 0.7) +
  geom_text_repel(
    data = q_last, aes(x, hbr, label = area, colour = zone), size = 3, nudge_x = 30, direction = "y", hjust = 0,
    segment.size = 0.2, min.segment.length = 0, box.padding = 0.15, show.legend = FALSE
  ) +
  scale_colour_manual(values = zone_cols_lab, guide = "none") +
  scale_fill_manual(values = zone_cols_lab, na.value = "white", guide = "none") +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"), guide = "none") +
  scale_linewidth_manual(values = c(`FALSE` = 0.5, `TRUE` = 1.1), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.6, `TRUE` = 2.4), guide = "none") +
  scale_y_continuous(transform = scales::pseudo_log_trans(base = 10, sigma = 0.03), breaks = c(0, 0.1, 1, 10), labels = c("0", "0.1", "1", "10")) +
  facet_grid(zone ~ .) +
  x_scale + coord_cartesian(clip = "off") +
  labs(y = "HBR (Anopheles per working trap-night)") +
  theme_ento + theme(strip.text.y = element_text(angle = 0, hjust = 0), axis.text.x = element_blank(), axis.ticks.x = element_blank())

figure5 <-
  (p_cal / p_area / (p_rain + theme(axis.text.x = element_text(), axis.ticks.x = element_line()))) +
  plot_layout(heights = c(1, 5.5, 1)) &
  scale_x_date(limits = c(as.Date("2022-06-01"), as.Date("2025-03-01")), date_breaks = "3 months", date_labels = "%b\n%Y", expand = expansion(mult = 0.01))
ggsave("output/figures/figure5_hbr_by_health_area.png", figure5, width = 9, height = 9, dpi = 300)

# S1 Text ----

# Share of mosquitoes tested, by zone x area x month (as-treated zone on the collection date;
# denominator: all Anopheles collected, working and non-working traps)
testing <-
  ento_mosq |>
  filter(date >= ento_spz_start, !is.na(zone)) |>
  summarise(
    anopheles_collected = n(), tested_n = sum(tested), positives = sum(spz_result %in% "Positive"),
    .by = c(zone, area, month)
  ) |>
  mutate(
    tested_share = tested_n / anopheles_collected,
    denominator = "all Anopheles collected in the cell (working and non-working traps)",
    zone = as.character(zone_display(zone))
  ) |>
  rename(zone_as_treated = zone) |>
  arrange(month, zone_as_treated, area)
write_csv(testing, "output/tables/s1_text_entomology_testing_coverage.csv")

testing |>
  summarise(tested_pct = 100 * sum(tested_n) / sum(anopheles_collected), .by = zone_as_treated) |>
  print()
table3_full |>
  summarise(
    working_trap_nights = sum(working_trap_nights), anopheles = sum(anopheles_working_traps),
    tested_in_cells = sum(mosquitoes_tested, na.rm = TRUE), positive_in_cells = sum(sporozoite_positives, na.rm = TRUE)
  ) |>
  print()
ento_mosq |>
  filter(tested) |>
  summarise(tested_all = n(), positive_all = sum(spz_result %in% "Positive")) |>
  print()
# Gengere's share of MDA+IRS working trap-nights and of its Anopheles, by half-year
mda_split |>
  filter(group != "All MDA+IRS areas") |>
  select(period, group, working_trap_nights, anopheles_working_traps) |>
  pivot_wider(names_from = group, values_from = c(working_trap_nights, anopheles_working_traps)) |>
  mutate(
    gengere_trap_night_pct = 100 * working_trap_nights_Gengere / (working_trap_nights_Gengere + `working_trap_nights_Other four areas`),
    gengere_anopheles_pct = 100 * anopheles_working_traps_Gengere / (anopheles_working_traps_Gengere + `anopheles_working_traps_Other four areas`)
  ) |>
  print(width = Inf)

# S3 Figure ----

fit_month <- ento_fit(unit = "month", n_boot = n_boot_monthly)
monthly <-
  fit_month$estimates |>
  select(-n_boot_valid) |>
  pivot_wider(names_from = metric, values_from = c(estimate, lcl, ucl)) |>
  left_join(fit_month$counts |> rename(month_u = u), by = c("zone", "u" = "month_u")) |>
  mutate(month = as.Date(u), zone = as.character(zone_display(zone))) |>
  select(-u) |>
  arrange(month, zone)
write_csv(monthly, "output/tables/s3_figure_entomology_monthly.csv")

# MDA+IRS without Gengere (dashed line): the full engine on the data without Gengere
fit_month_nog <- ento_fit(tn = ento_trapnights |> filter(area != "Gengere"), mq = ento_mosq |> filter(area != "Gengere"), unit = "month", n_boot = n_boot_monthly)
monthly_nog <-
  fit_month_nog$estimates |>
  filter(zone == "MDA") |>
  select(-n_boot_valid) |>
  pivot_wider(names_from = metric, values_from = c(estimate, lcl, ucl)) |>
  left_join(fit_month_nog$counts |> rename(month_u = u), by = c("zone", "u" = "month_u")) |>
  mutate(month = as.Date(u), zone = "MDA+IRS without Gengere") |>
  select(-u) |>
  arrange(month)

# Monthly pooled sporozoite rate: positives / tested in the zone-month, exact 95% CI
sr_monthly <-
  bind_rows(monthly, monthly_nog) |>
  transmute(
    month, zone, mosquitoes_tested, sporozoite_positives,
    sr_pooled = if_else(mosquitoes_tested > 0, sporozoite_positives / mosquitoes_tested, NA_real_),
    sr_pooled_lcl = exact_lcl(sporozoite_positives, mosquitoes_tested),
    sr_pooled_ucl = exact_ucl(sporozoite_positives, mosquitoes_tested)
  ) |>
  filter(mosquitoes_tested > 0) |>
  arrange(month, match(zone, c("Routine", "MDA+IRS", "IRS", "MDA+IRS without Gengere")))

# Four series: the three zones plus MDA+IRS without Gengere (dashed). Points are nudged by a few
# days per series so the intervals do not overlap; months with fewer than 20 working trap-nights
# are drawn open.
series_levels <- c(unname(zone_display_labels["Routine"]), unname(zone_display_labels["MDA"]), unname(zone_display_labels["IRS"]), "MDA+IRS without Gengere")
monthly_plot <-
  bind_rows(monthly, monthly_nog) |>
  left_join(sr_monthly |> select(month, zone, starts_with("sr_pooled")), by = c("month", "zone")) |>
  mutate(
    zone = factor(zone, levels = series_levels),
    month_x = month + c(-6, -2, 6, 2)[as.integer(zone)]
  )
series_cols <- set_names(c(zone_cols[["Routine"]], zone_cols[["MDA"]], zone_cols[["IRS"]], zone_cols[["MDA"]]), series_levels)
series_lty <- set_names(c("solid", "solid", "solid", "22"), series_levels)
min_trap_nights <- 20

# n_tested_min: sporozoite-based panels show a zone-month only with at least this many tested
# mosquitoes; open marks below 30. trap_col: trap-night count behind the open-mark rule (random-house
# trap-nights for density). The catch-rate panels use a pseudo-log axis that keeps zero months.
panel_metric <- function(est, lcl, ucl, ylab, pseudo_log = FALSE, n_tested_min = NULL, trap_col = "trap_nights") {
  d <- monthly_plot |> mutate(v = .data[[est]], lo = .data[[lcl]], hi = .data[[ucl]], thin = .data[[trap_col]] < min_trap_nights)
  if (!is.null(n_tested_min)) d <- d |> mutate(across(c(v, lo, hi), \(x) if_else(mosquitoes_tested >= n_tested_min, x, NA_real_)), thin = mosquitoes_tested < 30)
  p <-
    ggplot(d, aes(month_x, v, colour = zone, linetype = zone)) +
    annotate("rect", xmin = as.Date("2022-06-27"), xmax = as.Date("2022-11-10"), ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.6) +
    annotate("rect", xmin = as.Date("2023-02-02"), xmax = as.Date("2023-05-31"), ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.6) +
    geom_vline(xintercept = as.Date("2024-01-01"), linetype = "dashed", colour = "grey30") +
    geom_linerange(aes(ymin = lo, ymax = hi, alpha = zone), linetype = "solid", na.rm = TRUE) +
    geom_line(na.rm = TRUE, linewidth = 0.5) +
    geom_point(aes(shape = thin), size = 1.5, fill = "white", na.rm = TRUE) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), guide = "none") +
    scale_alpha_manual(values = set_names(c(0.5, 0.5, 0.5, 0.35), series_levels), guide = "none") +
    scale_colour_manual(values = series_cols) +
    scale_linetype_manual(values = series_lty) +
    x_scale + labs(y = ylab) + theme_ento
  if (pseudo_log) {
    p <- p + scale_y_continuous(transform = scales::pseudo_log_trans(base = 10, sigma = 0.05), breaks = c(0, 0.1, 1, 10, 100), labels = c("0", "0.1", "1", "10", "100"))
  }
  p
}

no_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p_sr <- panel_metric("sr_pooled", "sr_pooled_lcl", "sr_pooled_ucl", "Sporozoite rate\n(pooled, exact CI)", n_tested_min = 10) + no_x
p_hbr <- panel_metric("estimate_hbr", "lcl_hbr", "ucl_hbr", "HBR\n(Anopheles per working\ntrap-night, \nlog scale)", pseudo_log = TRUE) + no_x
p_eir <- panel_metric("estimate_eir_annual", "lcl_eir_annual", "ucl_eir_annual", "EIR\n(infective Anopheles per\ntrap per year)", n_tested_min = 10) + no_x
p_density <- panel_metric("estimate_density", "lcl_density", "ucl_density", "Density, random houses\n(Anopheles per working\ntrap-night, \nlog scale)", pseudo_log = TRUE, trap_col = "trap_nights_random") +
  theme(axis.text.x = element_text(), axis.ticks.x = element_line())

s3_figure <-
  (p_cal_s3 / p_rain / p_sr / p_hbr / p_eir / p_density) +
  plot_layout(heights = c(1.5, 0.8, 1.2, 1.4, 1.2, 1.4), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave("output/figures/s3_figure_entomology_monthly.png", s3_figure, width = 9, height = 12, dpi = 300)

# S5 Table ----

# Sheets A to C: one table, three levels stacked. Health-area rows carry bootstrap CIs from
# ento_fit() run on that area alone (month-coverage rule off so every area-period appears).
area_fits <-
  by_area |>
  distinct(area) |>
  pull(area) |>
  set_names() |>
  map(\(a) ento_fit(tn = ento_trapnights |> filter(area == a), mq = ento_mosq |> filter(area == a), unit = "period", min_months = 1, n_boot = n_boot_table))
area_ci <-
  area_fits |>
  map(\(f) fit_to_long(f, paste0("P", 1:5))) |>
  list_rbind(names_to = "area") |>
  mutate(zone = as.character(zone_display(zone))) |>
  select(area, zone, period, starts_with("lcl_"), starts_with("ucl_"), point_hbr = estimate_hbr, point_density = estimate_density, point_eir = estimate_eir_annual, hbr_pp = estimate_hbr_person, eir_pp = estimate_eir_person_annual)

sr_cols <- \(d) d |> mutate(
  tested_share = if_else(mosquitoes_collected > 0, mosquitoes_tested / mosquitoes_collected, NA_real_),
  sr_pooled = if_else(mosquitoes_tested > 0, sporozoite_positives / mosquitoes_tested, NA_real_),
  sr_pooled_lcl = exact_lcl(sporozoite_positives, mosquitoes_tested),
  sr_pooled_ucl = exact_ucl(sporozoite_positives, mosquitoes_tested),
  sr_divergence_pp = 100 * (sr_popweighted - sr_pooled),
  sr_popweighted_to_pooled = if_else(sr_pooled > 0, sr_popweighted / sr_pooled, NA_real_)
)

s5_area <-
  by_area |>
  left_join(area_ci, by = c("zone", "area", "period")) |>
  left_join(period_lab |> select(period, dates, survey = survey_in_period), by = "period") |>
  transmute(
    level = "health area", zone, unit = area, period, dates, survey, trapped_months, in_zone_estimate,
    working_trap_nights = trap_nights, anopheles_working_traps = anopheles,
    mosquitoes_collected = replace_na(mosquitoes, 0), mosquitoes_tested = replace_na(tested, 0), sporozoite_positives = replace_na(positives, 0),
    sr_popweighted = sporozoite_rate, sr_popweighted_lcl = lcl_sr, sr_popweighted_ucl = ucl_sr,
    hbr, hbr_lcl = lcl_hbr, hbr_ucl = ucl_hbr, density, density_lcl = lcl_density, density_ucl = ucl_density,
    eir_annual, eir_lcl = lcl_eir_annual, eir_ucl = ucl_eir_annual,
    hbr_per_occupant = hbr_pp, hbr_per_occupant_lcl = lcl_hbr_person, hbr_per_occupant_ucl = ucl_hbr_person,
    eir_per_occupant_annual = eir_pp, eir_per_occupant_lcl = lcl_eir_person_annual, eir_per_occupant_ucl = ucl_eir_person_annual,
    point_hbr, point_density, point_eir
  )
# the area-level engine reproduces the descriptive point estimates
stopifnot(
  isTRUE(all.equal(s5_area$hbr, s5_area$point_hbr)), isTRUE(all.equal(s5_area$density, s5_area$point_density)),
  isTRUE(all.equal(s5_area$eir_annual, s5_area$point_eir))
)
s5_area <- s5_area |> select(-starts_with("point_")) |> sr_cols()

# Zone and subgroup rows; mosquitoes collected summed over the area-periods that enter the estimates
collected_by_group <- \(zone_label, areas, per) {
  by_area |> filter(in_zone_estimate, zone == zone_label, area %in% areas, period == per) |> summarise(n = sum(mosquitoes)) |> pull(n)
}
s5_group <-
  bind_rows(
    table3_full |> mutate(level = "zone", unit = zone),
    mda_split |> filter(group != "All MDA+IRS areas") |> mutate(level = "MDA+IRS subgroup", unit = str_c("MDA+IRS: ", group))
  ) |>
  rowwise() |>
  mutate(mosquitoes_collected = collected_by_group(zone, str_split_1(areas_included, ", "), period)) |>
  ungroup() |>
  transmute(
    level, zone, unit, period, dates, survey, trapped_months = NA_real_, in_zone_estimate = NA,
    working_trap_nights, anopheles_working_traps, mosquitoes_collected, mosquitoes_tested, sporozoite_positives,
    sr_popweighted = sporozoite_rate_popweighted_sensitivity, sr_popweighted_lcl = sporozoite_popweighted_lcl, sr_popweighted_ucl = sporozoite_popweighted_ucl,
    hbr, hbr_lcl, hbr_ucl, density = density_random_houses, density_lcl, density_ucl, eir_annual, eir_lcl, eir_ucl,
    hbr_per_occupant, hbr_per_occupant_lcl, hbr_per_occupant_ucl,
    eir_per_occupant_annual, eir_per_occupant_lcl, eir_per_occupant_ucl
  ) |>
  sr_cols()

s5_detail <-
  bind_rows(s5_area, s5_group) |>
  mutate(level = factor(level, levels = c("health area", "zone", "MDA+IRS subgroup"))) |>
  arrange(level, match(zone, unname(zone_display_labels[zone_display_levels])), unit, period) |>
  mutate(level = as.character(level)) |>
  select(
    level, zone, unit, period, dates, survey, trapped_months, in_zone_estimate,
    working_trap_nights, anopheles_working_traps, mosquitoes_collected, mosquitoes_tested, tested_share, sporozoite_positives,
    sr_pooled, sr_pooled_lcl, sr_pooled_ucl, sr_popweighted, sr_popweighted_lcl, sr_popweighted_ucl, sr_divergence_pp, sr_popweighted_to_pooled,
    hbr, hbr_lcl, hbr_ucl, density, density_lcl, density_ucl, eir_annual, eir_lcl, eir_ucl,
    hbr_per_occupant, hbr_per_occupant_lcl, hbr_per_occupant_ucl, eir_per_occupant_annual, eir_per_occupant_lcl, eir_per_occupant_ucl
  )

# Sheet D: sensitivity variants (500 resamples). Each reruns ento_fit() on a modified copy of the
# data and returns the zone x period estimates and the change contrasts.
ento_sensitivity <- function(n_boot = 500) {
  tn <- ento_trapnights
  mq <- ento_mosq

  run <- function(variant, tn_v = tn, mq_v = mq, equal = FALSE, drop_person = FALSE, ...) {
    f <- ento_fit(tn = tn_v, mq = mq_v, unit = "period", n_boot = n_boot, pop = if (equal) mutate(ento_pop, pop = 1) else ento_pop, ...)
    list(
      est = f$estimates |>
        filter(u %in% paste0("P", 1:5), !(drop_person & metric %in% c("hbr_person", "eir_person_annual"))) |>
        transmute(variant = variant, metric, zone, period = u, estimate, lcl, ucl),
      con = f$change |>
        filter(!(type == "change_relative_to_routine" & zone == "Routine"), !(drop_person & metric %in% c("hbr_person", "eir_person_annual"))) |>
        transmute(variant = variant, comparison = type, metric, zone, period = u, reference_period = ref, estimate, lcl, ucl)
    )
  }

  wins_cap <- quantile(tn$n_anoph[tn$working], 0.99, na.rm = TRUE, names = FALSE)
  vec_counts <- mq |> filter(species_group != "other") |> count(tn_index, name = "n_vec")

  # Pooled zone totals with the same month-coverage rule; point estimates only, and changes are
  # ratios of pooled totals over all areas present in each period
  pooled_cells <-
    tn |>
    filter(!is.na(zone), period %in% paste0("P", 1:5), working, !is.na(n_anoph)) |>
    mutate(n_months = n_distinct(month), .by = c(area, zone, period)) |>
    filter(n_months >= 3)
  pooled_est <-
    pooled_cells |>
    summarise(hbr = sum(n_anoph) / n(), density = sum(n_anoph[random_house]) / sum(random_house), .by = c(zone, period)) |>
    pivot_longer(c(hbr, density), names_to = "metric", values_to = "estimate")
  pooled_con <-
    pooled_est |>
    pivot_wider(names_from = period, values_from = estimate) |>
    pivot_longer(c(P3, P5), names_to = "period", values_to = "later") |>
    transmute(variant = "pooled_trap_nights", comparison = "change_vs_ref", metric, zone, period, reference_period = "P1", estimate = later / P1, lcl = NA_real_, ucl = NA_real_)
  pooled_rel <-
    pooled_con |>
    mutate(routine = estimate[zone == "Routine"][match(paste(metric, period), paste(metric, period)[zone == "Routine"])], .by = c(metric, period)) |>
    filter(zone != "Routine") |>
    transmute(variant, comparison = "change_relative_to_routine", metric, zone, period, reference_period, estimate = estimate / routine, lcl, ucl)

  runs <- list(
    run("primary"),
    run("equal_area_weights", equal = TRUE),
    run("sud_without_gengere", tn |> filter(area != "Gengere"), mq |> filter(area != "Gengere")),
    run("per_household_member", tn |> mutate(n_anoph = n_anoph / persons), drop_person = TRUE),
    run("nonworking_included", tn |> mutate(working = TRUE)),
    run("winsorised_99", tn |> mutate(n_anoph = pmin(n_anoph, wins_cap))),
    run("awasi_only_routine", tn |> filter(!(zone_assigned == "Routine" & area != "Awasi")), mq |> filter(!(zone_assigned == "Routine" & area != "Awasi"))),
    run("awasi_only_routine_without_gengere", tn |> filter(!(zone_assigned == "Routine" & area != "Awasi"), area != "Gengere"), mq |> filter(!(zone_assigned == "Routine" & area != "Awasi"), area != "Gengere")),
    run("resample_houses", cluster = "house"),
    run("all_months_cells", min_months = 1),
    run("vectors_only", tn |> left_join(vec_counts, by = "tn_index") |> mutate(n_anoph = if_else(is.na(n_anoph), NA_real_, replace_na(n_vec, 0))), mq |> filter(species_group != "other"))
  )

  list(
    estimates = bind_rows(map(runs, "est"), pooled_est |> transmute(variant = "pooled_trap_nights", metric, zone, period, estimate, lcl = NA_real_, ucl = NA_real_)) |>
      arrange(variant, metric, period, match(zone, zone_display_levels)),
    contrasts = bind_rows(map(runs, "con"), pooled_con, pooled_rel) |>
      arrange(variant, comparison, metric, period, match(zone, zone_display_levels))
  )
}
sens <- ento_sensitivity(n_boot = 500)
sens_estimates <- sens$estimates |> mutate(zone = as.character(zone_display(zone)))
sens_contrasts <- sens$contrasts |> mutate(zone = as.character(zone_display(zone)))

# Sheet D model check: negative binomial GLMM of the Anopheles catch per working trap-night, with
# random intercepts for site and house, area fixed effects and zone x period cells. Zone is a
# function of health area while the area is in the comparison, so area effects are parametrised
# relative to one reference area per zone (Awasi, Pakala, Gengere) and each change contrast is one
# multiplicative change assumed common to every area of the zone, weighted by information and not
# by population. Variants: month-of-year adjustment, and rainfall over the four complete ISO weeks
# before the trap-night's week.
rain_all <-
  read_rds("data/climate_weekly.rds") |>
  mutate(jan4 = make_date(year, 1, 4), week_start = jan4 - (wday(jan4, week_start = 1) - 1) + 7 * (epi_week - 1)) |>
  select(week_start, precipitation)

ref_area <- c(Routine = "Awasi", IRS = "Pakala", MDA = "Gengere")
d0 <- ento_trapnights |> filter(!is.na(zone), period %in% paste0("P", 1:5), working, !is.na(n_anoph))
thin <- d0 |> summarise(n_months = n_distinct(month), .by = c(area, zone, period)) |> filter(n_months < 3)
d_check <-
  d0 |>
  anti_join(thin, by = c("area", "zone", "period")) |>
  mutate(
    zp = factor(str_c(zone, ".", period)),
    area_f = factor(if_else(area %in% ref_area, "ref", area)) |> relevel("ref"),
    month_of_year = factor(month(date)),
    site_id = factor(site_id),
    house_id = factor(house_id),
    rain_4wk = map_dbl(date, \(x) {
      wk <- x - (wday(x, week_start = 1) - 1)
      sum(rain_all$precipitation[rain_all$week_start >= wk - 28 & rain_all$week_start <= wk - 7])
    }) / 100
  )

fit_one <- function(adjust) {
  f <- if (adjust == "month_of_year") n_anoph ~ 0 + zp + area_f + month_of_year + (1 | site_id) + (1 | house_id) else n_anoph ~ 0 + zp + area_f + rain_4wk + (1 | site_id) + (1 | house_id)
  m <- tryCatch(glmmTMB(f, family = nbinom2, data = d_check), error = \(e) NULL)
  list(model = m, ok = !is.null(m) && isTRUE(m$sdr$pdHess))
}

contrast_rows <- function(m, adjust) {
  b <- fixef(m)$cond
  v <- vcov(m)$cond
  nm <- \(z, p) str_c("zp", z, ".", p)
  lc <- \(w) {
    w <- w[names(w) %in% names(b)]
    est <- sum(w * b[names(w)])
    se <- sqrt(as.numeric(t(w) %*% v[names(w), names(w)] %*% w))
    tibble(estimate = exp(est), lcl = exp(est - 1.96 * se), ucl = exp(est + 1.96 * se))
  }
  chg <- \(z, p) set_names(c(1, -1), c(nm(z, p), nm(z, "P1")))
  out <- list()
  for (z in c("Routine", "MDA", "IRS")) {
    for (p in c("P3", "P5")) {
      if (all(c(nm(z, p), nm(z, "P1")) %in% names(b))) out[[length(out) + 1]] <- bind_cols(tibble(adjust = adjust, comparison = "change_vs_P1", zone = z, period = p), lc(chg(z, p)))
      if (z != "Routine" && all(c(nm(z, p), nm(z, "P1"), nm("Routine", p), nm("Routine", "P1")) %in% names(b))) {
        w <- c(chg(z, p), -chg("Routine", p))
        out[[length(out) + 1]] <- bind_cols(tibble(adjust = adjust, comparison = "change_relative_to_routine", zone = z, period = p), lc(w))
      }
    }
  }
  bind_rows(out)
}

diag_row <- function(res, adjust) {
  if (!res$ok) return(tibble(adjust = adjust, converged = FALSE, n_trap_nights = nrow(d_check), zero_share_observed = NA_real_, zero_share_simulated = NA_real_, aic = NA_real_))
  sim <- simulate(res$model, nsim = 50)
  tibble(adjust = adjust, converged = TRUE, n_trap_nights = nobs(res$model), zero_share_observed = mean(d_check$n_anoph == 0), zero_share_simulated = mean(unlist(sim) == 0), aic = AIC(res$model))
}

# Likelihood-ratio test of an area x period interaction for the P3 against P1 change, within a zone
het_row <- function(z) {
  dz <- d_check |> filter(zone == z, period %in% c("P1", "P3")) |> mutate(area = factor(area))
  both <- intersect(dz$area[dz$period == "P1"], dz$area[dz$period == "P3"])
  dz <- dz |> filter(area %in% both) |> droplevels()
  if (n_distinct(dz$area) < 2) return(tibble(zone = z, areas = n_distinct(dz$area), lrt_chisq = NA_real_, lrt_df = NA_real_, lrt_p = NA_real_))
  m0 <- glmmTMB(n_anoph ~ area + period + (1 | site_id) + (1 | house_id), family = nbinom2, data = dz)
  m1 <- glmmTMB(n_anoph ~ area * period + (1 | site_id) + (1 | house_id), family = nbinom2, data = dz)
  a <- anova(m0, m1)
  tibble(zone = z, areas = n_distinct(dz$area), lrt_chisq = a$Chisq[2], lrt_df = a$`Chi Df`[2], lrt_p = a$`Pr(>Chisq)`[2])
}
heterogeneity <- map(c("Routine", "MDA", "IRS"), \(z) tryCatch(het_row(z), error = \(e) tibble(zone = z, areas = NA_integer_, lrt_chisq = NA_real_, lrt_df = NA_real_, lrt_p = NA_real_))) |> list_rbind()

fits <- map(c("month_of_year", "rainfall_4wk") |> set_names(), fit_one)
model_check <- bind_rows(imap(fits, \(r, a) if (r$ok) contrast_rows(r$model, a) else tibble())) |> mutate(zone = unname(zone_display_labels[zone]))
model_diagnostics <- bind_rows(imap(fits, diag_row))
model_heterogeneity <- heterogeneity |> mutate(zone = as.character(zone_display(zone)))

# Workbook: Read me, A zones, B MDA+IRS split, C health areas, D sensitivity / contrasts / model
# check, E P0. Units: HBR is Anopheles per working trap-night (CDC light trap, not a human landing
# catch); EIR is infective Anopheles per trap per year. Ratios, rates and indices are given to two
# significant figures when the upper limit exceeds 10 times the lower; percentages are not.
eir_label <- "infective Anopheles per trap per year"
metric_label <- c(
  hbr = "HBR (Anopheles per working trap-night)",
  density = "Density (Anopheles per working trap-night, randomly selected houses)",
  hbr_person = "HBR per occupant-night",
  sr = "Sporozoite rate",
  eir_annual = str_c("EIR (", eir_label, ")"),
  eir_person_annual = "EIR per occupant (infective bites per year)"
)
metric_digits <- c(hbr = 2, density = 2, hbr_person = 2, sr = 1, eir_annual = 1, eir_person_annual = 1)
variant_label <- c(
  primary = "Primary (population weights, site resampling)",
  equal_area_weights = "Equal weight for every health area within a zone",
  pooled_trap_nights = "Zone totals pooled over areas (no interval)",
  sud_without_gengere = "MDA+IRS zone without Gengere",
  per_household_member = "Catches divided by household size before averaging",
  nonworking_included = "Trap-nights with a non-working trap included",
  winsorised_99 = "Catches above the 99th percentile set to it",
  awasi_only_routine = "Routine comparator fixed to Awasi",
  awasi_only_routine_without_gengere = "Routine comparator fixed to Awasi, MDA+IRS without Gengere",
  resample_houses = "Bootstrap over GPS-grouped houses instead of sites",
  all_months_cells = "Cells with fewer than three trapped months kept",
  vectors_only = "Anopheles other than gambiae/funestus s.l. removed"
)
comparison_label <- c(
  change_vs_ref = "Change vs reference period, within zone",
  change_relative_to_routine = "Change relative to Routine's change",
  change_vs_P1 = "Change vs P1, within zone",
  between_zone_ratio = "Between-zone ratio"
)
adjust_label <- c(month_of_year = "Month of year", rainfall_4wk = "Rainfall, four weeks before the trap-night")

est_ci <- \(d, est, lcl, ucl, digits) fmt_est_ci_sf(d[[est]], d[[lcl]], d[[ucl]], digits = digits)
pct_ci <- \(d, est, lcl, ucl) fmt_est_ci_sf(d[[est]], d[[lcl]], d[[ucl]], digits = 1, scale = 100, sf_rule = FALSE)
eir_col <- str_c("EIR, ", eir_label, " (95% CI)")

panel_tbl <- function(d, ratio_p5_p1 = NULL) {
  out <-
    tibble(
      Unit = if_else(!is.na(d$in_zone_estimate) & !d$in_zone_estimate, str_c(d$unit, " †"), d$unit),
      Period = d$period, Dates = d$dates, Survey = d$survey, `Trapped months` = d$trapped_months,
      `Working trap-nights` = d$working_trap_nights,
      `Anopheles caught (working traps)` = d$anopheles_working_traps,
      `Mosquitoes collected` = d$mosquitoes_collected,
      `Mosquitoes tested` = d$mosquitoes_tested,
      `Tested, %` = round(100 * d$tested_share, 1),
      Positive = d$sporozoite_positives,
      `Sporozoite rate pooled, % (95% CI)` = pct_ci(d, "sr_pooled", "sr_pooled_lcl", "sr_pooled_ucl"),
      `Sporozoite rate population-weighted, % (95% CI)` = pct_ci(d, "sr_popweighted", "sr_popweighted_lcl", "sr_popweighted_ucl"),
      `Difference, percentage points` = round(d$sr_divergence_pp, 1),
      `HBR (95% CI)` = est_ci(d, "hbr", "hbr_lcl", "hbr_ucl", 2),
      `Density (95% CI)` = est_ci(d, "density", "density_lcl", "density_ucl", 2),
      !!eir_col := est_ci(d, "eir_annual", "eir_lcl", "eir_ucl", 1),
      `HBR per occupant-night (95% CI)` = est_ci(d, "hbr_per_occupant", "hbr_per_occupant_lcl", "hbr_per_occupant_ucl", 2),
      `EIR per occupant, infective bites per year (95% CI)` = est_ci(d, "eir_per_occupant_annual", "eir_per_occupant_lcl", "eir_per_occupant_ucl", 1)
    )
  if (!is.null(ratio_p5_p1)) {
    out[["HBR ratio, P5 vs P1 (point estimate)"]] <- if_else(is.na(ratio_p5_p1), NA_character_, sprintf("%.2f", ratio_p5_p1))
  }
  out
}

panel_a <- s5_detail |> filter(level == "zone") |> panel_tbl()
panel_b <- s5_detail |> filter(level == "MDA+IRS subgroup") |> panel_tbl()

# The P5 vs P1 ratio goes on each area's P5 row
area_rows <- s5_detail |> filter(level == "health area")
area_ratio_by_row <-
  area_rows |>
  select(unit, period) |>
  left_join(area_change |> transmute(unit = area, period = "P5", ratio_p5_p1), by = c("unit", "period")) |>
  pull(ratio_p5_p1)
panel_c <- panel_tbl(area_rows, ratio_p5_p1 = area_ratio_by_row)

# Sporozoite rate is a proportion in the sensitivity estimates (shown as %, no two-significant-figure rule)
panel_d_sens <-
  sens_estimates |>
  mutate(
    est_ci = fmt_est_ci_sf(estimate, lcl, ucl, digits = unname(metric_digits[metric]), scale = if_else(metric == "sr", 100, 1), sf_rule = metric != "sr"),
    measure = if_else(metric == "sr", "Sporozoite rate, %", unname(metric_label[metric]))
  ) |>
  transmute(Variant = unname(variant_label[variant]), Measure = measure, Zone = zone, Period = period, `Estimate (95% CI)` = est_ci)

ratio_ci <- \(d) fmt_est_ci_sf(d$estimate, d$lcl, d$ucl, digits = 2)

panel_d_contrasts <-
  sens_contrasts |>
  transmute(
    Variant = unname(variant_label[variant]), Comparison = unname(comparison_label[comparison]),
    Measure = unname(metric_label[metric]), Zone = zone, Period = period, `Reference period` = reference_period,
    `Ratio (95% CI)` = ratio_ci(pick(estimate, lcl, ucl))
  )

panel_d_check <-
  model_check |>
  transmute(
    Adjustment = unname(adjust_label[adjust]),
    Comparison = unname(comparison_label[comparison]), Zone = zone, Period = period,
    `Ratio (95% CI)` = ratio_ci(pick(estimate, lcl, ucl))
  )

# Sheet E: P0 runs from 1 May to 30 June 2022, so no mosquitoes were tested and sporozoite rate and
# EIR are not estimated
p0_desc <-
  table3_p0 |>
  transmute(
    Zone = zone, Period = period, `Health areas (n)` = areas, `Health areas included` = areas_included,
    Sites = sites, Houses = houses,
    `Working trap-nights` = trap_nights, `Anopheles caught (working traps)` = anopheles_working_traps,
    `Trap-nights, randomly selected houses` = trap_nights_random,
    `Mosquitoes tested` = mosquitoes_tested, Positive = sporozoite_positives,
    `HBR (95% CI)` = fmt_est_ci_sf(estimate_hbr, lcl_hbr, ucl_hbr, 2),
    `Density (95% CI)` = fmt_est_ci_sf(estimate_density, lcl_density, ucl_density, 2),
    `HBR per occupant-night (95% CI)` = fmt_est_ci_sf(estimate_hbr_person, lcl_hbr_person, ucl_hbr_person, 2)
  )

e_metrics <- c("hbr", "density", "hbr_person")

p0_ratios <-
  table3_contrasts |>
  filter(comparison == "between_zone_ratio", period == "P0", metric %in% e_metrics) |>
  transmute(
    Comparison = str_c(zone, " vs ", zone2), Measure = unname(metric_label[metric]),
    `Ratio (95% CI)` = fmt_est_ci_sf(estimate, lcl, ucl, 2), `Bootstrap resamples` = n_boot_valid,
    metric_order = match(metric, e_metrics)
  ) |>
  arrange(metric_order, match(Comparison, c("MDA+IRS vs Routine", "IRS vs Routine", "MDA+IRS vs IRS"))) |>
  select(-metric_order)

p1_dates <- s5_detail |> filter(level == "zone", period == "P1") |> pull(dates) |> unique()
p5_dates <- s5_detail |> filter(level == "zone", period == "P5") |> pull(dates) |> unique()

change_ratios <-
  table3_contrasts |>
  filter(period == "P5", reference_period == "P1", metric %in% e_metrics, comparison %in% c("change_vs_ref", "change_relative_to_routine")) |>
  transmute(
    Comparison = unname(comparison_label[comparison]), Zone = zone, Measure = unname(metric_label[metric]),
    `Ratio, P5 vs P1 (95% CI)` = fmt_est_ci_sf(estimate, lcl, ucl, 2), `Bootstrap resamples` = n_boot_valid,
    cmp_order = match(comparison, c("change_vs_ref", "change_relative_to_routine")), metric_order = match(metric, e_metrics)
  ) |>
  arrange(cmp_order, metric_order, match(Zone, c("Routine", "MDA+IRS", "IRS"))) |>
  select(-cmp_order, -metric_order)

read_me <- c(
  "S5 Table. Entomological surveillance by zone, health area, and period, Angumu health zone, Democratic Republic of Congo, May 2022 to November 2024.",
  "Sheets A to C: A = the three zones, B = the MDA+IRS zone divided into Gengere and the other four health areas, C = each health area, by half-year period (P1 to P5); sheet C also gives each area's estimated HBR ratio for P5 against P1. Sheets D = sensitivity analyses (500 bootstrap resamples; Table 5 uses 2,000): estimates, change contrasts, and a negative-binomial model check of the HBR contrasts. Sheet E = May and June 2022 (P0), before and at the start of IRS round 1: descriptives, between-zone ratios at P0, and the change in estimated HBR, density and HBR per occupant-night from P1 to P5 (2,000 bootstrap resamples).",
  "Pooled sporozoite rate = positive mosquitoes / mosquitoes tested, with an exact 95% CI; the population-weighted rate weights each health area by its population and so up-weights health areas with few tested mosquitoes (MDA+IRS, January to June 2024: pooled 1.8%, population-weighted 7.2%).",
  str_c("HBR, estimated human biting rate: Anopheles per working trap-night in CDC light traps (not a human landing catch); per occupant-night, the catch divided by household members. EIR: estimated HBR x sporozoite rate x 365, in ", eir_label, "; per occupant, infective bites per year. Density: Anopheles per working trap-night in randomly selected houses."),
  "Ratios, rates and entomological indices are given to two significant figures, estimate and interval alike, when the upper limit is more than 10 times the lower; percentages are not.",
  "Panel C rows marked with a dagger had fewer than three trapped months and are outside the zone estimate. P0 runs from 1 May to 30 June 2022; trapping began on 31 May 2022 and IRS round 1 began on 27 June 2022.",
  "[CONFIRM (Arnaud): how mosquitoes were chosen for sporozoite testing; re-dated and held-out trap-nights; trap-night, catch, tested and positive totals are provisional.]"
)

wb <- createWorkbook()
head_style <- createStyle(textDecoration = "bold", wrapText = TRUE, valign = "bottom")
title_style <- createStyle(textDecoration = "bold")
wrap_style <- createStyle(wrapText = TRUE, valign = "top")

add_table_sheet <- function(name, d, widths = NULL) {
  addWorksheet(wb, name)
  writeData(wb, name, d, headerStyle = head_style, na.string = "")
  setColWidths(wb, name, cols = seq_len(ncol(d)), widths = widths %||% pmin(pmax(nchar(names(d)) * 0.8, 10), 32))
  freezePane(wb, name, firstRow = TRUE)
}

addWorksheet(wb, "Read me")
writeData(wb, "Read me", tibble(text = read_me), colNames = FALSE)
addStyle(wb, "Read me", wrap_style, rows = seq_along(read_me), cols = 1)
setColWidths(wb, "Read me", cols = 1, widths = 140)

add_table_sheet("A zones", panel_a)
add_table_sheet("B MDA+IRS split", panel_b)
add_table_sheet("C health areas", panel_c)
add_table_sheet("D sensitivity", panel_d_sens, widths = c(46, 46, 12, 8, 24))
add_table_sheet("D contrasts", panel_d_contrasts, widths = c(46, 38, 46, 12, 8, 12, 22))
add_table_sheet("D model check", panel_d_check, widths = c(40, 38, 12, 8, 22))

# Sheet E: three titled blocks on one sheet
sheet_e <- "E P0 May-Jun 2022"
addWorksheet(wb, sheet_e)
row <- 1
write_block <- function(title, d) {
  writeData(wb, sheet_e, title, startRow = row, startCol = 1)
  addStyle(wb, sheet_e, title_style, rows = row, cols = 1)
  writeData(wb, sheet_e, d, startRow = row + 1, startCol = 1, headerStyle = head_style, na.string = "")
  row <<- row + nrow(d) + 3
}
write_block("E1. P0 descriptives, 1 May to 30 June 2022 (no mosquitoes tested, so sporozoite rate and EIR are not estimated)", p0_desc)
write_block("E2. Between-zone ratios at P0 (2,000 bootstrap resamples)", p0_ratios)
write_block(str_c("E3. Ratio of the zone mean in P5 (", p5_dates, ") to that in P1 (", p1_dates, "), within each zone, then each intervention zone's ratio relative to Routine's ratio (2,000 bootstrap resamples)"), change_ratios)
setColWidths(wb, sheet_e, cols = 1:14, widths = c(24, 44, 16, 26, 10, 10, 14, 16, 16, 14, 10, 22, 22, 22))

saveWorkbook(wb, "output/tables/s5_table_entomology.xlsx", overwrite = TRUE)
