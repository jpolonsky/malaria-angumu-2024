# Helpers shared by the analysis scripts.

# Labels ----

# Models keep Routine as the reference level; tables and figures show this order.
zone_levels <- c("Routine", "IRS", "MDA")
zone_cols <- c(Routine = "#009E73", IRS = "#0072B2", MDA = "#D55E00")
zone_display_levels <- c("Routine", "MDA", "IRS")
zone_display_labels <- c(Routine = "Routine", MDA = "MDA+IRS", IRS = "IRS")

zone_display <- function(x) {
  factor(as.character(x), levels = zone_display_levels, labels = unname(zone_display_labels[zone_display_levels]))
}

# Health areas in the zone comparison: name in the survey data and in the GRID3 map, zone at survey 1, and the
# survey and month from which a Routine area that received an intervention leaves the comparison.
study_areas <-
  tibble::tribble(
    ~area, ~map_name, ~baseline_zone, ~crossover_survey, ~crossover_date,
    "Abia", "Abia", "Routine", 4L, "2023-08-01",
    "Ara", "Ara", "Routine", 5L, "2024-01-01",
    "Awasi", "Awasi Luga", "Routine", NA_integer_, NA_character_,
    "Cawa", "Cawa", "Routine", 4L, "2023-08-01",
    "Linga", "Linga", "Routine", 5L, "2024-01-01",
    "Anzika", "Anzika", "IRS", NA_integer_, NA_character_,
    "Pakala", "Pakala", "IRS", NA_integer_, NA_character_,
    "Panyadong", "Panyadong", "IRS", NA_integer_, NA_character_,
    "Udongo Abira", "Udongo Abira", "IRS", NA_integer_, NA_character_,
    "Apala", "Apala", "MDA", NA_integer_, NA_character_,
    "Gengere", "Gengere", "MDA", NA_integer_, NA_character_,
    "Lanyi", "Lanyi", "MDA", NA_integer_, NA_character_,
    "Ugudo-Zii", "Ugudo Zii", "MDA", NA_integer_, NA_character_,
    "Uyandu", "Uyandu", "MDA", NA_integer_, NA_character_
  ) |>
  dplyr::mutate(crossover_date = as.Date(crossover_date))

# Zone of a health area in a month: Routine areas that crossed over are excluded (NA) from their crossover month.
# Surveillance exports spell Ugudo-Zii as in the map, so both names are matched.
zone_at_date <- function(area, date) {
  i <- dplyr::coalesce(match(area, study_areas$area), match(area, study_areas$map_name))
  out_from <- study_areas$crossover_date[i]
  factor(dplyr::if_else(!is.na(out_from) & date >= out_from, NA_character_, study_areas$baseline_zone[i]), levels = zone_levels)
}

# First and last interview date of each survey; months from baseline use the midpoints.
survey_timing <-
  tibble::tribble(
    ~survey_n, ~start, ~end, ~season,
    1, "2022-03-31", "2022-04-24", "Low transmission",
    2, "2022-11-17", "2022-12-04", "High transmission",
    3, "2023-04-27", "2023-05-16", "Low transmission",
    4, "2023-10-16", "2023-10-21", "High transmission",
    5, "2024-04-18", "2024-05-04", "Low transmission"
  ) |>
  dplyr::mutate(
    start = as.Date(start),
    end = as.Date(end),
    mid = start + (end - start) / 2,
    months_from_baseline = as.integer(round(as.numeric(mid - mid[survey_n == 1]) / 30.4375))
  )

season_lu <- survey_timing |> dplyr::select(survey_n, months_from_baseline, season)

# Formatting ----

fmt_sig <- function(x, sf = 2) {
  out <- rep(NA_character_, length(x))
  ok <- !is.na(x)
  zero <- ok & x == 0
  out[zero] <- "0"
  pos <- ok & !zero
  decimals <- pmax(0, sf - 1 - floor(log10(abs(x[pos]))))
  out[pos] <- sprintf("%.*f", decimals, signif(x[pos], sf))
  out
}

# Intervals spanning more than a tenfold range are shown to two significant figures.
ci_is_wide <- function(lcl, ucl, ratio = 10) {
  !is.na(lcl) & !is.na(ucl) & ucl > 0 & (lcl <= 0 | ucl / lcl > ratio)
}

fmt_est_ci_sf <- function(est, lcl, ucl, digits = 2, sep = "–", scale = 1, sf_rule = TRUE) {
  est <- scale * est
  lcl <- scale * lcl
  ucl <- scale * ucl
  fixed <- \(x) dplyr::if_else(is.na(x), NA_character_, sprintf("%.*f", digits, x))
  wide <- ci_is_wide(lcl, ucl) & sf_rule
  e <- dplyr::if_else(wide, fmt_sig(est), fixed(est))
  l <- dplyr::if_else(wide, fmt_sig(lcl), fixed(lcl))
  u <- dplyr::if_else(wide, fmt_sig(ucl), fixed(ucl))
  dplyr::case_when(
    is.na(est) ~ NA_character_,
    is.na(lcl) | is.na(ucl) ~ e,
    .default = stringr::str_c(e, " (", l, sep, u, ")")
  )
}

fmt_ci <- \(est, lcl, ucl) sprintf("%.2f (%.2f-%.2f)", est, lcl, ucl)
fmt_p <- \(p) dplyr::if_else(p < 0.001, "<0.001", sprintf("%.3f", p))

# Survey-weighted analysis ----

# Survey-weighted percentage of `value` among `key`: single-stage cluster design, strata = zone.
calculate_prop <- function(data, key, value) {
  data |>
    dplyr::filter(!is.na(cluster), !is.na(zone)) |>
    srvyr::as_survey(ids = cluster, strata = zone, weights = weight) |>
    dplyr::summarise(prop = srvyr::survey_mean({{ key }} %in% {{ value }} * 100, vartype = "ci", deff = TRUE, na.rm = TRUE))
}

# Contrasts of each zone against `ref` within each survey, multivariate-t adjusted.
tidy_contrasts <- function(emm, ref = "Routine") {
  emmeans::contrast(emm, method = "trt.vs.ctrl", ref = ref, adjust = "mvt") |>
    summary(infer = TRUE) |>
    tibble::as_tibble() |>
    dplyr::rename(or = odds.ratio, lcl = asymp.LCL, ucl = asymp.UCL) |>
    dplyr::mutate(zone = stringr::str_remove(as.character(contrast), stringr::str_glue(" / {ref}")))
}

# Survey-weighted prevalence by zone, one design per survey (clusters, strata = zone). `dat` needs survey,
# zone, cluster_uid, weight and outcome_rdt.
prevalence_by_zone <- function(dat) {
  dat |>
    dplyr::mutate(survey_n = as.integer(as.character(survey))) |>
    dplyr::group_split(survey_n) |>
    purrr::map(\(df) {
      df |>
        srvyr::as_survey(ids = cluster_uid, strata = zone, weights = weight) |>
        dplyr::group_by(zone) |>
        dplyr::summarise(
          n_total = srvyr::unweighted(dplyr::n()),
          n_pos = srvyr::unweighted(sum(outcome_rdt)),
          prev = srvyr::survey_mean(outcome_rdt, vartype = "ci", na.rm = TRUE)
        ) |>
        dplyr::mutate(survey_n = unique(df$survey_n))
    }) |>
    purrr::list_rbind() |>
    dplyr::rename(lcl = prev_low, ucl = prev_upp)
}
