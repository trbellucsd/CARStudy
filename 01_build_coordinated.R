## =============================================================================
## CAR Study -- 01_build_coordinated.R
##
## Maps each cohort's native variables onto the common analysis schema used
## by the Aim 1-3 scripts. Consistent with the coordinated analysis framework:
## measures are operationalized study-appropriately within each cohort, and
## cohorts are compared and combined at the estimate level, not pooled into a
## single harmonized measure. Nothing downstream touches cohort-native names.
##
## Outputs written to _derived/:
##   driving_long_all_cohorts.rds   id, cohort, age, state, wave, proxy, design
##   baseline_covariates.rds        id, cohort, age_base, sex, race_eth, educ
##   predictor_domains_long.rds     id, cohort, domain, age, value
##   adrd_index_all_cohorts.rds     id, cohort, age_dx, year_dx, age_death,
##                                  age_claims_start, age_claims_end, age_last_obs
##
## State coding, ascending severity: 0 unrestricted, 1 restricted, 2 ceased.
## This REVERSES the HRS RwDRVST5 coding built earlier (higher = more driving)
## and matches NHATS drive_state3 as built.
## =============================================================================

source("00_config_and_helpers.R")

## =============================================================================
## HRS
##
## RwDRVST5: 0 never drove, 1 ceased, 2 able but no driving last month,
##           3 drives limited to nearby, 4 drives long trips
##
## Mapping to the coordinated state definition:
##   4 -> 0 (unrestricted)
##   3 -> 1 (restricted)
##   2 -> 1 (restricted; see note)
##   1 -> 2 (ceased)
##   0 -> dropped (never drove; excluded from all analyses)
##
## State 2 is folded into restricted rather than kept separate. In 2018 only
## 55.3% of people in that state had a car available, so it is roughly half
## loss of vehicle access rather than behavioral restriction, and RwDRVCAR
## does not exist from 2022 onward so the distinction cannot be carried
## forward. Sensitivity analysis below treats state 2 as unrestricted instead.
## =============================================================================

build_hrs <- function(st5_to_state = c(`4` = 0L, `3` = 1L, `2` = 1L, `1` = 2L)) {
  drv  <- read_any(PATHS$hrs_drive)
  rand <- read_any(PATHS$hrs_rand)
  names(rand) <- toupper(names(rand))

  waves <- HRS_WAVES %>% filter(st5_ok)

  long <- map_dfr(seq_len(nrow(waves)), function(i) {
    w <- waves$wave[i]; y <- waves$year[i]
    gd <- function(f) { v <- f(w); if (v %in% names(drv))  drv[[v]][match(rand$HHIDPN, drv$HHIDPN)] else NA_real_ }
    gr <- function(f) { v <- f(w); if (v %in% names(rand)) num(rand[[v]]) else NA_real_ }
    tibble(
      id     = as.numeric(rand$HHIDPN),
      cohort = "HRS",
      wave   = w, year = y,
      age    = gr(HRS_VARS$age),
      proxy  = gr(HRS_VARS$proxy),
      iwstat = gr(HRS_VARS$iwstat),
      st5    = gd(HRS_VARS$state5),
      never  = gd(HRS_VARS$never),
      caravl = gd(HRS_VARS$caravail),
      svy_wt = gr(HRS_VARS$wt_resp)
    )
  })

  never_ids <- long %>% filter(never == 1) %>% pull(id) %>% unique()

  long %>%
    filter(proxy == 0,                      # self-report only
           !id %in% never_ids,
           !is.na(st5), st5 >= 1, !is.na(age)) %>%
    mutate(state = unname(st5_to_state[as.character(st5)])) %>%
    filter(!is.na(state)) %>%
    left_join(
      tibble(id     = as.numeric(rand$HHIDPN),
             svy_psu   = num(rand[[HRS_VARS$cluster]]),
             svy_strat = num(rand[[HRS_VARS$strata]])),
      by = "id") %>%
    select(id, cohort, wave, year, age, state, proxy, caravl,
           svy_wt, svy_psu, svy_strat)
}

hrs_drv <- build_hrs()

## Sensitivity variant: state 2 treated as unrestricted
hrs_drv_alt <- build_hrs(st5_to_state = c(`4` = 0L, `3` = 1L, `2` = 0L, `1` = 2L))
saveRDS(hrs_drv_alt, file.path(DER, "hrs_driving_state2_as_unrestricted.rds"))

## ---- HRS baseline covariates -------------------------------------------
rand <- read_any(PATHS$hrs_rand); names(rand) <- toupper(names(rand))
hrs_cov <- tibble(
  id       = as.numeric(rand$HHIDPN),
  cohort   = "HRS",
  sex      = factor(num(rand[[HRS_VARS$sex]]), 1:2, c("male", "female")),
  educ     = num(rand[[HRS_VARS$educ_yrs]]),
  race_m   = num(rand[[HRS_VARS$race]]),
  hisp     = num(rand[[HRS_VARS$hispanic]]),
  birthyr  = num(rand[[HRS_VARS$birthyear]]),
  deathyr  = num(rand[[HRS_VARS$deathyear]])
) %>%
  mutate(race_eth = factor(case_when(
    hisp %in% 1:3 ~ "hispanic",
    race_m == 1   ~ "white_nh",
    race_m == 2   ~ "black_nh",
    TRUE          ~ "other_nh"), levels = c("white_nh","black_nh","hispanic","other_nh"))) %>%
  select(id, cohort, sex, educ, race_eth, birthyr, deathyr)

hrs_cov <- hrs_drv %>% group_by(id) %>% summarise(age_base = min(age), .groups = "drop") %>%
  right_join(hrs_cov, by = "id")

## ---- HRS predictor domains ---------------------------------------------
hrs_pred <- map_dfr(HRS_WAVES$wave, function(w) {
  gr <- function(f) { v <- f(w); if (v %in% names(rand)) num(rand[[v]]) else NA_real_ }
  age <- gr(HRS_VARS$age)
  tibble(id = as.numeric(rand$HHIDPN), cohort = "HRS", age = age,
         physical_function      = -gr(HRS_VARS$adl),      # oriented: higher = better
         functional_independence= -gr(HRS_VARS$iadl),
         comorbidity            =  gr(HRS_VARS$conde),
         hospitalization        =  gr(HRS_VARS$hosp),
         physical_activity      =  gr(HRS_VARS$vigact),
         social_support         =  NA_real_)              # from leave-behind module
}) %>%
  pivot_longer(-c(id, cohort, age), names_to = "domain", values_to = "value") %>%
  filter(!is.na(value), !is.na(age))

## ---- HRS cognition from Langa-Weir --------------------------------------
lw <- read_any(PATHS$hrs_langa); names(lw) <- tolower(names(lw))
hrs_cog <- map_dfr(grep("^cogfunction[0-9]{4}$", names(lw), value = TRUE), function(cc) {
  y <- as.integer(str_extract(cc, "[0-9]{4}"))
  tibble(id = make_hhidpn(lw$hhid, lw$pn), cohort = "HRS", year = y,
         cogfun = num(lw[[cc]]))
}) %>%
  filter(!is.na(cogfun)) %>%
  ## 1 normal, 2 CIND, 3 dementia -> oriented so higher = better
  mutate(domain = "cognition", value = 4 - cogfun) %>%
  left_join(hrs_drv %>% distinct(id, year, age), by = c("id", "year")) %>%
  filter(!is.na(age)) %>%
  select(id, cohort, domain, age, value)

hrs_pred <- bind_rows(hrs_pred, hrs_cog)

## =============================================================================
## NHATS
##
## drive_state3 as built: 0 drives no avoidance, 1 drives with avoidance,
## 2 not driving. Already ascending severity, so it passes through.
## Age is a five-year band until the Sensitive File DUA supplies exact DOB.
## =============================================================================

build_nhats <- function() {
  d <- read_any(PATHS$nhats_drive)
  names(d) <- tolower(names(d))

  never_ids <- d %>% filter(never_drove == 1) %>% pull(spid) %>% unique()

  ## design variables, pulled per round from the SP files
  sp_files <- list.files(PATHS$nhats_dir, pattern = "(?i)sp_file.*\\.(sas7bdat|sav)$",
                         recursive = TRUE, full.names = TRUE)
  design <- map_dfr(sp_files, function(p) {
    nms <- peek_names(p)
    wv <- find_var(nms, NHATS_VARS$weight)
    vu <- find_var(nms, NHATS_VARS$varunit)
    vs <- find_var(nms, NHATS_VARS$varstrat)
    if (is.na(wv)) return(NULL)
    r <- as.integer(str_extract(wv, "[0-9]+"))
    x <- read_any(p); names(x) <- tolower(names(x))
    tibble(spid = as.character(x$spid), round = r,
           svy_wt    = num(x[[wv]]),
           svy_psu   = if (!is.na(vu)) num(x[[vu]]) else NA_real_,
           svy_strat = if (!is.na(vs)) num(x[[vs]]) else NA_real_)
  })

  d %>%
    mutate(spid = as.character(spid)) %>%
    filter(proxy == 0, !spid %in% never_ids, !is.na(drive_state3)) %>%
    left_join(design, by = c("spid", "round")) %>%
    transmute(id = as.numeric(spid), cohort = "NHATS",
              wave = round, year = 2010L + round,
              age = NHATS_AGE_MID[as.integer(age_cat)],
              state = as.integer(drive_state3),
              proxy, svy_wt, svy_psu, svy_strat)
}

nhats_drv <- build_nhats()

## ---- NHATS comorbidity with the HC2 carry-forward -----------------------
## Items labelled INCLUDES EVER use 1 = new report, 7 = previously reported,
## 2 = no. Ever = 1 or 7. Reading code 1 alone gives hypertension as 154
## people in R14 instead of 5,073.
nh_raw <- read_any(PATHS$nhats_inputs)
nhats_comorb <- nh_raw %>%
  mutate(spid = as.character(spid)) %>%
  select(spid, round, num_range("hc", 1:10)) %>%
  pivot_longer(starts_with("hc"), names_to = "item", values_to = "val") %>%
  mutate(idx = as.integer(str_remove(item, "hc"))) %>%
  left_join(NHATS_HC, by = "idx") %>%
  mutate(ever = as.integer(val %in% c(1, 7))) %>%
  group_by(spid, round) %>%
  summarise(value = sum(ever, na.rm = TRUE), .groups = "drop") %>%
  transmute(id = as.numeric(spid), cohort = "NHATS", domain = "comorbidity",
            round, value) %>%
  left_join(nhats_drv %>% distinct(id, wave, age), by = c("id", "round" = "wave")) %>%
  filter(!is.na(age)) %>% select(-round)

nhats_cov <- nhats_drv %>%
  group_by(id) %>% summarise(age_base = min(age, na.rm = TRUE), .groups = "drop") %>%
  mutate(cohort = "NHATS")
## sex, race_eth, educ come from the NHATS demographic file; join here once
## the file path is confirmed.

## =============================================================================
## ACTIVE
##
## CURDRIV0 1 = drives, 2 = does not. Avoidance items are phrased as "have you
## driven ...", so 2 = avoids. EVERDRIV 2 = never driven. Control arm only.
## =============================================================================

build_active <- function() {
  a <- read_any(PATHS$active_po)
  V <- ACTIVE_VARS
  g <- function(v) if (v %in% names(a)) num(a[[v]]) else NA_real_

  never_ids <- a[[V$id]][which(g(V$everdrove) == 2)] %>% unique()

  a %>%
    mutate(.drv = g(V$drives), .age = g(V$age), .arm = g(V$arm),
           .nt = g(V$avoid_night), .al = g(V$avoid_alone), .rn = g(V$avoid_rain)) %>%
    filter(.arm == 0, !.data[[V$id]] %in% never_ids,
           .drv %in% c(1, 2), !is.na(.age)) %>%
    mutate(n_avoid = rowSums(cbind(.nt == 2, .al == 2, .rn == 2), na.rm = TRUE),
           state = case_when(.drv == 2 ~ 2L, n_avoid > 0 ~ 1L, TRUE ~ 0L)) %>%
    transmute(id = as.numeric(.data[[V$id]]), cohort = "ACTIVE",
              wave = .data[[V$occasion]], year = .data[[V$year]],
              age = .age, state, proxy = 0,
              svy_wt = NA_real_, svy_psu = NA_real_, svy_strat = NA_real_)
}

active_drv <- build_active()

active_raw <- read_any(PATHS$active_po)
active_cov <- active_raw %>%
  group_by(id = as.numeric(.data[[ACTIVE_VARS$id]])) %>%
  summarise(
    cohort   = "ACTIVE",
    age_base = suppressWarnings(first(na.omit(num(.data[[ACTIVE_VARS$age_base]])))),
    sex      = suppressWarnings(first(na.omit(num(.data[[ACTIVE_VARS$sex]])))),
    race     = suppressWarnings(first(na.omit(num(.data[[ACTIVE_VARS$race]])))),
    educ     = suppressWarnings(first(na.omit(num(.data[[ACTIVE_VARS$educ]])))),
    arm      = suppressWarnings(first(na.omit(num(.data[[ACTIVE_VARS$arm]])))),
    mci_ever = suppressWarnings(max(num(.data[[ACTIVE_VARS$mci_ever]]), na.rm = TRUE)),
    mci_age  = suppressWarnings(first(na.omit(num(.data[[ACTIVE_VARS$mci_age]])))),
    .groups = "drop") %>%
  filter(arm == 0) %>%
  mutate(sex = factor(sex, 1:2, c("male", "female")),
         race_eth = factor(race, 1:3, c("white_nh", "black_nh", "other_nh")),
         across(where(is.numeric), ~replace(.x, is.infinite(.x), NA_real_)))

## ACTIVE continuous measures for the complementary Aim 1 analysis
active_cont <- active_raw %>%
  transmute(id = as.numeric(.data[[ACTIVE_VARS$id]]),
            occasion = .data[[ACTIVE_VARS$occasion]],
            age = num(.data[[ACTIVE_VARS$age]]),
            arm = num(.data[[ACTIVE_VARS$arm]]),
            days_driven   = num(.data[[ACTIVE_VARS$days]]),
            driving_space = num(.data[[ACTIVE_VARS$space]]),
            avoidance     = num(.data[[ACTIVE_VARS$avoid_sum]]),
            miles         = num(.data[[ACTIVE_VARS$miles]]),
            prompted      = num(.data[[ACTIVE_VARS$prompted]])) %>%
  filter(arm == 0) %>%
  ## orient so higher = greater restriction, then standardize
  mutate(days_driven_z   = as.numeric(scale(-days_driven)),
         driving_space_z = as.numeric(scale(-driving_space)),
         avoidance_z     = as.numeric(scale(avoidance)))
saveRDS(active_cont, file.path(DER, "active_continuous_long.rds"))

## =============================================================================
## Combine
## =============================================================================

driving_long <- bind_rows(hrs_drv, nhats_drv, active_drv) %>% arrange(cohort, id, age)
saveRDS(driving_long, file.path(DER, "driving_long_all_cohorts.rds"))

covariates <- bind_rows(
  hrs_cov    %>% select(id, cohort, age_base, sex, race_eth, educ),
  nhats_cov  %>% select(id, cohort, age_base),
  active_cov %>% select(id, cohort, age_base, sex, race_eth, educ)
)
saveRDS(covariates, file.path(DER, "baseline_covariates.rds"))

predictors <- bind_rows(hrs_pred, nhats_comorb) %>% filter(!is.na(value))
saveRDS(predictors, file.path(DER, "predictor_domains_long.rds"))

## =============================================================================
## ADRD index file
##
## Primary definition: first qualifying CMS claim, subtypes combined. This is
## the first CMS-RECORDED diagnosis, not biological onset. AD-specific codes
## (G30, 331.0) and a two-claim definition are evaluated in sensitivity
## analyses.
##
## Until the linked files arrive, the survey-based classifications stand in:
## HRS Langa-Weir dementia, NHATS hc*disescn9 with carry-forward, ACTIVE
## AGE_FIRST_MCI. Swap PATHS$cms_adrd in once linkage is complete.
## =============================================================================

build_adrd_interim <- function() {
  hrs_dx <- hrs_cog %>%
    filter(value == 1) %>%                       # 4 - cogfun == 1  => dementia
    group_by(id, cohort) %>%
    summarise(age_dx = min(age), .groups = "drop")

  nh_dx <- nh_raw %>%
    mutate(spid = as.character(spid),
           pos = dx_raw %in% c(1, 7)) %>%
    filter(pos) %>%
    group_by(spid) %>% summarise(round_dx = min(round), .groups = "drop") %>%
    transmute(id = as.numeric(spid), cohort = "NHATS", round_dx) %>%
    left_join(nhats_drv %>% distinct(id, wave, age), by = c("id", "round_dx" = "wave")) %>%
    transmute(id, cohort, age_dx = age)

  ac_dx <- active_cov %>% filter(mci_ever == 1) %>%
    transmute(id, cohort, age_dx = mci_age)

  bind_rows(hrs_dx, nh_dx, ac_dx)
}

adrd_dx <- if (file.exists(PATHS$cms_adrd)) read_any(PATHS$cms_adrd) else build_adrd_interim()

last_obs <- driving_long %>% group_by(id, cohort) %>%
  summarise(age_last_obs = max(age), .groups = "drop")

death <- hrs_cov %>%
  transmute(id, cohort, age_death = ifelse(is.na(deathyr), NA_real_, deathyr - birthyr))

adrd_index <- last_obs %>%
  left_join(adrd_dx, by = c("id", "cohort")) %>%
  left_join(death,   by = c("id", "cohort")) %>%
  mutate(year_dx = NA_integer_,            # filled from CMS claim date
         ## claims observability window; replace with FFS enrollment spans
         age_claims_start = pmin(age_last_obs, 65, na.rm = TRUE),
         age_claims_end   = pmax(age_last_obs, coalesce(age_dx, 0), na.rm = TRUE))
saveRDS(adrd_index, file.path(DER, "adrd_index_all_cohorts.rds"))

## =============================================================================
## Build report
## =============================================================================

cat("\n=== coordinated driving file ===\n")
driving_long %>% group_by(cohort) %>%
  summarise(persons = n_distinct(id), obs = n(),
            waves = n_distinct(wave),
            pct_unrest = round(100*mean(state == 0), 1),
            pct_rest   = round(100*mean(state == 1), 1),
            pct_ceased = round(100*mean(state == 2), 1)) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n=== transition pairs ===\n")
make_intervals(driving_long) %>% group_by(cohort) %>%
  summarise(pairs = n(), modal_gap = as.numeric(names(sort(table(dt), decreasing = TRUE))[1]),
            to_rest = sum(from == 0 & to == 1), to_cease = sum(to == 2),
            recovery = sum(from == 1 & to == 0)) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\nbuild complete -> ", DER)
