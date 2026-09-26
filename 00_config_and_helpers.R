## =============================================================================
## CAR Study -- 00_config_and_helpers.R
##
## Shared configuration, cohort variable dictionaries, person-interval
## construction, survey design objects, and estimation helpers.
##
## Conventions, per the Research Strategy:
##   - Age is the time scale in every primary model.
##   - Discrete-time hazards use a complementary log-log link with
##     log(interval length) as an offset.
##   - Coordinated driving states: 0 unrestricted, 1 restricted, 2 ceased.
##     Each cohort is coded to this common definition from its own items;
##     cohorts are combined at the estimate level, never pooled at the item
##     level.
##     NOTE this is the ASCENDING-severity coding used in all modeling. It is
##     the reverse of the HRS RwDRVST5 coding built earlier (where higher =
##     more driving) and matches the NHATS drive_state3 coding. 01_build
##     handles the conversion; nothing downstream should re-map states.
##   - Primary models are cohort-specific; pooling happens only at the
##     estimate level via random-effects meta-analysis.
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(stringr)
  library(readr); library(tibble); library(haven)
  library(survey); library(survival); library(metafor); library(splines)
})

options(survey.lonely.psu = "adjust")

ROOT <- "C:/Users/trbell/Documents/Lab/Research/CARS"
DER  <- file.path(ROOT, "_derived")
OUT  <- file.path(ROOT, "_results")
MPLUS_DIR <- file.path(ROOT, "_mplus")
walk(c(DER, OUT, MPLUS_DIR), dir.create, showWarnings = FALSE, recursive = TRUE)

SEED <- 20260924
set.seed(SEED)

COHORTS <- c("HRS", "NHATS", "ACTIVE")

## =============================================================================
## Cohort source files
## =============================================================================

PATHS <- list(
  hrs_drive   = file.path(ROOT, "HRS", "derived", "cars_hrs_driving_wide.rds"),
  hrs_rand    = file.path(ROOT, "HRS", "randhrs1992_2022v1_SPSS",
                          "randhrs1992_2022v1.sav"),
  hrs_langa   = file.path(ROOT, "HRS", "cogfinalimp_9522wide.sas7bdat"),
  hrs_tracker = file.path(ROOT, "HRS", "trk2022tr_r.sas7bdat"),
  nhats_drive = file.path(ROOT, "NHATS", "_derived", "nhats_driving_long.rds"),
  nhats_inputs= file.path(ROOT, "NHATS", "_derived",
                          "nhats_classification_inputs_raw.rds"),
  nhats_dir   = file.path(ROOT, "NHATS"),
  active_po   = file.path(ROOT, "ACTIVE", "active_person_occasion.csv"),
  hcap        = file.path(ROOT, "HRS", "HCAP", "hcap_biomarkers.sas7bdat"),
  cms_adrd    = file.path(ROOT, "CMS", "adrd_first_claim.rds")   # post-linkage
)

## =============================================================================
## HRS
##
## Driving comes from the wide file built earlier: R{w}DRVST5 / R{w}DRVST3 /
## R{w}DRVASK / R{w}DRVNEVER, keyed on HHIDPN = HHID*1000 + PN.
## RAND supplies age, proxy status, interview status, demographics, and the
## design variables. Langa-Weir supplies cogfunction{year}.
##
## Usable waves for the five-level state: 8, 9, 10, 11, 12, 14, 16
## (2006-2014, 2018, 2022). Wave 13 (2016) and wave 15 (2020) are excluded:
## the alternate-wave rotation left only 55 and 734 respondents asked.
## =============================================================================

HRS_WAVES <- tibble(
  wave = c(6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L),
  year = c(2002L, 2004L, 2006L, 2008L, 2010L, 2012L, 2014L, 2016L, 2018L,
           2020L, 2022L),
  st5_ok = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE),
  st3_ok = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)
)

HRS_VARS <- list(
  id        = "HHIDPN",
  state5    = function(w) sprintf("R%dDRVST5", w),
  state3    = function(w) sprintf("R%dDRVST3", w),
  asked     = function(w) sprintf("R%dDRVASK", w),
  never     = function(w) sprintf("R%dDRVNEVER", w),
  able      = function(w) sprintf("R%dDRVABLE", w),
  drovemo   = function(w) sprintf("R%dDRVMONTH", w),
  limit     = function(w) sprintf("R%dDRVLIM", w),
  caravail  = function(w) sprintf("R%dDRVCAR", w),
  ## RAND
  age       = function(w) sprintf("R%dAGEY_E", w),
  proxy     = function(w) sprintf("R%dPROXY", w),
  iwstat    = function(w) sprintf("R%dIWSTAT", w),
  wt_resp   = function(w) sprintf("R%dWTRESP", w),
  ## RAND time-invariant
  sex       = "RAGENDER",
  race      = "RARACEM",
  hispanic  = "RAHISPAN",
  educ_yrs  = "RAEDYRS",
  birthyear = "RABYEAR",
  deathyear = "RADYEAR",
  strata    = "RAESTRAT",
  cluster   = "RAEHSAMP",
  ## RAND time-varying predictor domains (Aim 2)
  adl       = function(w) sprintf("R%dADLA", w),
  iadl      = function(w) sprintf("R%dIADLZA", w),
  cesd      = function(w) sprintf("R%dCESD", w),
  shlt      = function(w) sprintf("R%dSHLT", w),
  conde     = function(w) sprintf("R%dCONDE", w),
  hosp      = function(w) sprintf("R%dHOSP", w),
  vigact    = function(w) sprintf("R%dVGACTX", w),
  ## Langa-Weir
  cogfun    = function(y) sprintf("cogfunction%d", y)
)

## =============================================================================
## NHATS
##
## Driving comes from the long file built earlier, keyed on spid.
## drive_state3 there is 0 = drives no avoidance, 1 = drives with avoidance,
## 2 = not driving, which already matches the coordinated ascending-severity
## coding, so it passes through unchanged.
##
## Variable stems are matched by regex rather than pasted, because the
## never-drove flag is numbered for the round the respondent's COHORT entered,
## not the file's round (R14 carries fl13dneverdrv).
## =============================================================================

NHATS_ROUNDS <- 1:14

NHATS_VARS <- list(
  id           = "spid",
  drives       = "^fl[0-9]+drives$",
  never_drove  = "^fl[0-9]+dneverdrv$",
  drive_freq   = "^dt[0-9]+oftedrive$",
  drove_year   = "^dt[0-9]+driveyr$",
  avoid_night  = "^dt[0-9]+avoidriv1$",
  avoid_alone  = "^dt[0-9]+avoidriv2$",
  avoid_hwy    = "^dt[0-9]+avoidriv3$",
  avoid_weather= "^dt[0-9]+avoidriv4$",
  resp_type    = "^is[0-9]+resptype$",
  age_band     = "^r[0-9]+d2intvrage$",
  ## dementia classification inputs (Kasper); demclas is NOT in the public
  ## files and must be built with the Technical Paper 5 addenda
  dx_dementia  = "^hc[0-9]+disescn9$",
  ad8          = "^cp[0-9]+dad8dem$",
  clock_manual = "^cg[0-9]+dclkdraw$",
  clock_machine= "^cg[0-9]+dclkdlnn$",
  clock_attempt= "^cg[0-9]+atdrwclck$",
  recall_imm   = "^cg[0-9]+dwrdimmrc$",
  recall_dly   = "^cg[0-9]+dwrddlyrc$",
  ## design
  weight       = "^w[0-9]+anfinwgt0$",
  varunit      = "^w[0-9]+varunit$",
  varstrat     = "^w[0-9]+varstrat$"
)

## HC2 battery. Items labelled "INCLUDES EVER" carry code 7 = previously
## reported; reading code 1 alone undercounts prevalence several-fold.
## Heart attack, stroke and cancer are asked fresh each round and have no
## code 7.
NHATS_HC <- tibble(
  idx       = 1:10,
  condition = c("heart_attack", "heart_disease", "hypertension", "arthritis",
                "osteoporosis", "diabetes", "lung_disease", "stroke",
                "dementia", "cancer"),
  carry_fwd = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE)
)

## Five-year age bands; midpoints are a stopgap until the Sensitive File DUA
## supplies exact date of birth.
NHATS_AGE_MID <- c(67, 72, 77, 82, 87, 92)

## =============================================================================
## ACTIVE
##
## Six occasions at years 0, 1, 2, 3, 5, 10. Control arm only (INTGRPR == 0),
## per the scope of work. The source file fans out to six rows per person.
## =============================================================================

ACTIVE_OCC <- tibble(occasion = 1:6, year_nominal = c(0, 1, 2, 3, 5, 10))

ACTIVE_VARS <- list(
  id        = "ID",
  occasion  = "occasion",
  year      = "year_nominal",
  age       = "age",
  drives    = "CURDRIV0",          # 1 = yes, 2 = no
  everdrove = "EVERDRIV",          # 2 = never driven
  avoid_night = "NIGHTDRV",        # 1 = has driven at night, 2 = avoids
  avoid_alone = "ALONDRIV",
  avoid_rain  = "RAINDRIV",
  prompted  = "LIMITDRV",          # anyone suggested limiting/stopping
  prompted2 = "LIMDRIV",
  last_yr   = "LDRIVYER",
  last_mo   = "LDRIVMON",
  miles     = "MILEDRIV",
  days      = "DAYSDRIV",
  avoid_sum = "DAVOID",
  space     = "TOTDS",
  ## person level
  arm       = "INTGRPR",           # 0 = no-contact control
  arm_full  = "INTGRP",
  age_base  = "AGEB",
  sex       = "GENDER",
  race      = "RACE_CAT",
  educ      = "EDUCLEVL",
  marital   = "MARSTAT_CAT",
  mci_ever  = "MCI_EVER_2",
  mci_base  = "MCI_B",
  mci_age   = "AGE_FIRST_MCI",
  comorb    = "COMORBIDITY_B",
  srh       = "SF36GH_B"
)

## =============================================================================
## Design variables, under the common names produced by 01_build
## =============================================================================

DESIGN <- list(
  HRS    = list(weight = "svy_wt", cluster = "svy_psu", strata = "svy_strat"),
  NHATS  = list(weight = "svy_wt", cluster = "svy_psu", strata = "svy_strat"),
  ACTIVE = list(weight = NULL,     cluster = "id",      strata = NULL)
)

## Baseline covariates Z_i
COVARS <- c("age_base", "sex", "race_eth", "educ")

## Aim 2 predictor domains, with the cohort-specific source noted. Functional
## independence excludes driving, transportation and life-space items.
DOMAINS <- tribble(
  ~domain,                   ~hrs,        ~nhats,            ~active,
  "cognition",               "cogtot",    "cog_composite",   "mci_status",
  "physical_function",       "R{w}ADLA",  "sppb",            "timed_iadl",
  "comorbidity",             "R{w}CONDE", "hc_count",        "COMORBIDITY_B",
  "hospitalization",         "R{w}HOSP",  "hosp_claims",     "hosp_claims",
  "medication_burden",       "rx_count",  "rx_count",        "rx_count",
  "functional_independence", "R{w}IADLZA","iadl_no_transp",  "iadl_no_transp",
  "social_support",          "lb_support","social_partic",   "social_partic",
  "physical_activity",       "R{w}VGACTX","phys_activity",   "phys_activity"
)

## Transitions modeled in Aim 2. Cessation is absorbing.
TRANSITIONS <- tribble(
  ~name,             ~from, ~to, ~primary,
  "unrest_to_rest",      0L,  1L,  TRUE,
  "unrest_to_cease",     0L,  2L,  TRUE,
  "rest_to_cease",       1L,  2L,  TRUE,
  "rest_to_unrest",      1L,  0L,  FALSE
)

## =============================================================================
## Utilities
## =============================================================================

read_any <- function(path) {
  if (!file.exists(path)) stop("not found: ", path)
  ext <- tolower(tools::file_ext(path))
  switch(ext,
    rds      = readRDS(path),
    sas7bdat = { cf <- sub("\\.sas7bdat$", ".sas7bcat", path)
                 if (file.exists(cf)) read_sas(path, catalog_file = cf)
                 else read_sas(path) },
    sav      = read_sav(path),
    dta      = read_dta(path),
    csv      = read_csv(path, show_col_types = FALSE),
    stop("unhandled extension: ", ext))
}

peek_names <- function(path) {
  ext <- tolower(tools::file_ext(path))
  tryCatch(tolower(names(switch(ext,
    sas7bdat = read_sas(path, n_max = 0),
    sav      = read_sav(path, n_max = 0),
    csv      = read_csv(path, n_max = 0, show_col_types = FALSE)))),
    error = function(e) character(0))
}

#' Find a column by regex, case-insensitively, returning NA if absent.
find_var <- function(nms, pattern) {
  h <- grep(pattern, nms, value = TRUE, ignore.case = TRUE)
  if (length(h)) h[1] else NA_character_
}

num <- function(x) suppressWarnings(as.numeric(zap_labels(x)))

#' HRS identifier. Pad after coercion: SAS may deliver HHID and PN as
#' numerics with leading zeros already stripped.
make_hhidpn <- function(hhid, pn) {
  as.numeric(paste0(str_pad(str_trim(as.character(hhid)), 6, pad = "0"),
                    str_pad(str_trim(as.character(pn)),   3, pad = "0")))
}

## =============================================================================
## Person-interval construction
## =============================================================================

#' Expand a coordinated person-wave file into person-intervals.
#' Input columns: id, cohort, age, state (0/1/2 ascending severity).
make_intervals <- function(d) {
  d %>%
    filter(!is.na(state), !is.na(age)) %>%
    arrange(id, age) %>%
    group_by(id) %>%
    mutate(ceased_before = cumsum(lag(state, default = 0L) == 2L)) %>%
    filter(ceased_before == 0) %>%
    mutate(state_next = lead(state),
           age_next   = lead(age),
           dt         = lead(age) - age) %>%
    ungroup() %>%
    filter(!is.na(state_next), dt > 0) %>%
    transmute(id, cohort, age_start = age, age_end = age_next, dt,
              from = state, to = state_next, log_dt = log(dt),
              across(any_of(c("svy_wt", "svy_psu", "svy_strat"))))
}

## =============================================================================
## Survey design
## =============================================================================

#' Two-level weight rescaling (Asparouhov method A) then svydesign.
make_design <- function(iv, cohort_name) {
  cfg <- DESIGN[[cohort_name]]
  if (is.null(cfg$weight)) {
    return(svydesign(ids = as.formula(paste0("~", cfg$cluster)),
                     weights = ~1, data = iv))
  }
  iv <- iv %>%
    group_by(.data[[cfg$cluster]]) %>%
    mutate(w_scaled = .data[[cfg$weight]] * n() / sum(.data[[cfg$weight]])) %>%
    ungroup()
  svydesign(ids     = as.formula(paste0("~", cfg$cluster)),
            strata  = as.formula(paste0("~", cfg$strata)),
            weights = ~w_scaled, data = iv, nest = TRUE)
}

## =============================================================================
## Discrete-time cause-specific hazard model
## =============================================================================

fit_hazard <- function(design, event, rhs, age_spec = "ns(age_start, df = 3)") {
  f <- as.formula(sprintf("%s ~ %s + %s + offset(log_dt)", event, age_spec, rhs))
  svyglm(f, design = design, family = binomial(link = "cloglog"))
}

fit_hazard_unwt <- function(dat, event, rhs, cluster = "id",
                            age_spec = "ns(age_start, df = 3)") {
  f <- as.formula(sprintf("%s ~ %s + %s + offset(log_dt)", event, age_spec, rhs))
  m <- glm(f, data = dat, family = binomial(link = "cloglog"))
  list(fit = m, vcov = sandwich::vcovCL(m, cluster = dat[[cluster]]))
}

## =============================================================================
## Competing-risks cumulative incidence from cause-specific discrete hazards
## =============================================================================

cif_from_hazards <- function(haz_list, ages) {
  H <- map_dfc(haz_list, ~approx(.x$age, .x$h, xout = ages, rule = 2)$y)
  names(H) <- names(haz_list)
  H <- as.matrix(H); H[is.na(H)] <- 0
  h_all  <- rowSums(H)
  S_prev <- c(1, head(cumprod(1 - h_all), -1))
  as_tibble(sweep(H, 1, S_prev, `*`)) %>%
    mutate(across(everything(), cumsum)) %>%
    mutate(age = ages, surv = cumprod(1 - h_all), .before = 1)
}

predict_hazard <- function(fit, newdata, ages, age_var = "age_start") {
  nd <- newdata[rep(1, length(ages)), , drop = FALSE]
  nd[[age_var]] <- ages
  nd$log_dt <- 0
  eta <- predict(fit, newdata = nd, type = "link")
  tibble(age = ages, h = 1 - exp(-exp(as.numeric(eta))))
}

rm_years_free <- function(cif_case, cif_ctrl, col, a_lo, a_hi) {
  keep <- cif_case$age >= a_lo & cif_case$age <= a_hi
  tibble(window_lo = a_lo, window_hi = a_hi,
         rm_case    = sum(1 - cif_case[[col]][keep]),
         rm_control = sum(1 - cif_ctrl[[col]][keep]),
         diff       = sum(1 - cif_case[[col]][keep]) -
                      sum(1 - cif_ctrl[[col]][keep]))
}

## =============================================================================
## Bootstrap, pooling, multiplicity, tidying
## =============================================================================

boot_persons <- function(iv, stat_fun, B = 1000) {
  ids <- unique(iv$id)
  map_dfr(seq_len(B), function(b) {
    take <- tibble(id = sample(ids, length(ids), replace = TRUE))
    as_tibble(as.list(stat_fun(
      take %>% left_join(iv, by = "id", relationship = "many-to-many"))))
  }) %>%
    summarise(across(everything(),
                     list(est = ~mean(.x, na.rm = TRUE),
                          lo  = ~quantile(.x, .025, na.rm = TRUE),
                          hi  = ~quantile(.x, .975, na.rm = TRUE))))
}

pool_estimates <- function(est_df, by = c("outcome", "term")) {
  est_df %>%
    group_by(across(all_of(by))) %>%
    group_modify(~{
      if (nrow(.x) < 2)
        return(tibble(k = nrow(.x), est = .x$est, se = .x$se,
                      lo = .x$est - 1.96*.x$se, hi = .x$est + 1.96*.x$se,
                      tau2 = NA_real_, I2 = NA_real_, Q_p = NA_real_))
      m <- rma(yi = .x$est, sei = .x$se, method = "REML")
      tibble(k = m$k, est = as.numeric(m$b), se = m$se,
             lo = m$ci.lb, hi = m$ci.ub, tau2 = m$tau2, I2 = m$I2, Q_p = m$QEp)
    }) %>% ungroup()
}

bh_within <- function(df, p_col = "p", family_cols = c("outcome", "family")) {
  df %>% group_by(across(all_of(family_cols))) %>%
    mutate(p_fdr = p.adjust(.data[[p_col]], method = "BH")) %>% ungroup()
}

tidy_svyglm <- function(fit, keep = NULL, label = NULL) {
  s <- summary(fit)$coefficients
  out <- tibble(term = rownames(s), est = s[, 1], se = s[, 2],
                z = s[, 3], p = s[, 4])
  if (!is.null(keep))  out <- out %>% filter(str_detect(term, keep))
  if (!is.null(label)) out <- out %>% mutate(!!!label)
  out %>% mutate(hr = exp(est), hr_lo = exp(est - 1.96*se),
                 hr_hi = exp(est + 1.96*se))
}

message("config loaded: ", paste(COHORTS, collapse = ", "))
