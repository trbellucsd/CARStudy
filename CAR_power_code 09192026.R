# Aim 1: Monte Carlo power for age at FIRST REPORTED driving transition.
# Base R only. Run: Rscript CAR_Aim1_first_transition_power.R
# Table 2 fixes cohort sizes, repeated-observation counts, observed event counts,
# age slopes, and recovery percentages. Table 1 supplies the visit schedules.
# Replace make_demo_visits() with your real, visit-level analytic skeleton when ready.
# A real skeleton needs: id, cohort (HRS/NHATS/ACTIVE), age at each completed
# interview, future_adrd (0/1); retain only PRE-DIAGNOSIS visits among cases.
# Begin with current drivers age >=65. Include visits even after an observed
# driving event; the simulation generates NEW transitions at those visits.

set.seed(20260918)
B <- 500L                        # Start at 50 for a quick check; use >=1000 for grant estimates.
effects <- c(1.15, 1.25, 1.40)  # Annual first-transition hazard ratios for future ADRD.
future_adrd_fraction <- 0.25     # ILLUSTRATIVE; replace with cohort-specific observed fractions.
design_effect <- c(HRS=1.5, NHATS=1.7, ACTIVE=1.0)   # ILLUSTRATIVE; replace from survey analysis.

design <- data.frame(
  cohort=c('HRS','NHATS','ACTIVE'),
  n=c(19099L,15314L,655L),
  person_observations=c(63708L,70424L,2572L),
  persons_with_2plus=c(14154L,12361L,539L),
  transition_pairs=c(44608L,55110L,1917L),
  restriction_events=c(4400L,3663L,174L),
  cessation_events=c(2472L,2526L,58L),
  restriction_age_slope=c(0.068,0.033,0.061),
  cessation_age_slope=c(0.121,0.086,0.127),
  recovery_percent=c(7.0,7.1,5.8),
  stringsAsFactors=FALSE
)
visit_years <- list(HRS=c(seq(0,12,2),16,20), NHATS=0:13,
                    ACTIVE=c(0,1,2,3,5,10))
age_slope <- as.matrix(design[,c('restriction_age_slope','cessation_age_slope')])
dimnames(age_slope) <- list(design$cohort,c('restriction','cessation'))

assign_visit_counts <- function(n,n2,observations,max_waves) {
  # First n-n2 people have one visit; n2 have at least two. Distribute the
  # remaining Table 2 observations without exceeding Table 1's visit schedule.
  remaining <- observations-(n+n2)
  stopifnot(remaining>=0,remaining<=n2*(max_waves-2L))
  extra <- integer(n2)
  while (remaining>0) {
    draw <- as.vector(rmultinom(1,remaining,rep.int(1,n2)))
    accepted <- pmin(draw,max_waves-2L-extra)
    extra <- extra+accepted
    remaining <- remaining-sum(accepted)
  }
  sample(c(rep.int(1L,n-n2),2L+extra))
}

make_demo_visits <- function(design, future_fraction) {
  out <- vector('list',nrow(design))
  for (k in seq_len(nrow(design))) {
    z <- design[k,]; times <- visit_years[[z$cohort]]
    id <- paste0(z$cohort,'_',seq_len(z$n))
    if (z$cohort=='ACTIVE') {
      age0 <- pmin(94,pmax(65,rnorm(z$n,73.6,5.8)))
    } else {
      age0 <- runif(z$n,65,if(z$cohort=='HRS') 85 else 90)
    }
    group <- rbinom(z$n,1,future_fraction)
    # Table 2 does not report diagnosis ages. Its totals cannot determine how
    # many visits precede diagnosis; supply real prediagnosis visits for that.
    nwaves <- assign_visit_counts(z$n,z$persons_with_2plus,
                                  z$person_observations,length(times))
    rows <- lapply(seq_len(z$n), function(i) {
      a <- age0[i]+times[seq_len(nwaves[i])]
      data.frame(id=id[i],cohort=z$cohort,age=a,
                 future_adrd=group[i],stringsAsFactors=FALSE)
    })
    out[[k]] <- do.call(rbind,rows)
  }
  do.call(rbind,out)
}

make_intervals <- function(visits) {
  required <- c('id','cohort','age','future_adrd')
  stopifnot(all(required %in% names(visits)))
  visits <- visits[order(visits$cohort,visits$id,visits$age),required]
  stopifnot(!anyNA(visits),all(visits$age>=65),
            all(visits$future_adrd %in% 0:1),
            all(!duplicated(visits[c('cohort','id','age')])))
  visits$person <- match(paste(visits$cohort,visits$id),
                         unique(paste(visits$cohort,visits$id)))
  visits$wave <- ave(visits$age,visits$person,FUN=seq_along)
  byperson <- split(visits,visits$person)
  intervals <- lapply(byperson,function(v) {
    if (nrow(v)<2L) return(NULL)
    data.frame(person=v$person[-1L],cohort=v$cohort[-1L],
               future_adrd=v$future_adrd[-1L],
               wave=v$wave[-1L],start=v$age[-nrow(v)],end=v$age[-1L])
  })
  d <- do.call(rbind,intervals)
  rownames(d) <- NULL
  d$dt <- d$end-d$start
  d$age_mid <- (d$start+d$end)/2
  stopifnot(nrow(d)>0,all(d$dt>0))
  d
}

# A participant can first report restriction while still driving, cease without
# a reported restriction, or cease after restriction. Thus cessation competes
# with a FIRST REPORT of restriction. Interview loss/death is represented by
# the final visit in the supplied skeleton; it is NOT simulated independently.
simulate_transitions <- function(d, base_rates, hr_group, age_slope) {
  state <- integer(max(d$person)) # 0 unrestricted, 1 restricted, 2 ceased
  restrict_risk <- cease_risk <- restrict_event <- cease_event <- integer(nrow(d))
  for (w in sort(unique(d$wave))) {
    j <- which(d$wave==w)
    person <- d$person[j]; before <- state[person]
    r_risk <- before==0L; c_risk <- before!=2L
    restrict_risk[j] <- as.integer(r_risk)
    cease_risk[j] <- as.integer(c_risk)
    r_rate <- base_rates[d$cohort[j],'restriction'] *
      exp(age_slope[d$cohort[j],'restriction']*(d$age_mid[j]-75)) *
      hr_group['restriction']^d$future_adrd[j]
    c_rate <- base_rates[d$cohort[j],'cessation'] *
      exp(age_slope[d$cohort[j],'cessation']*(d$age_mid[j]-75)) *
      hr_group['cessation']^d$future_adrd[j]
    time_r <- rep(Inf,length(j));time_c <- rep(Inf,length(j))
    time_r[r_risk] <- rexp(sum(r_risk),rate=r_rate[r_risk])
    time_c[c_risk] <- rexp(sum(c_risk),rate=c_rate[c_risk])
    ceased <- c_risk & time_c<d$dt[j]
    # Restriction must still be observable at the END of the interview interval.
    restricted <- r_risk & time_r<d$dt[j] & !ceased
    restrict_event[j] <- as.integer(restricted)
    cease_event[j] <- as.integer(ceased)
    state[person[restricted]] <- 1L
    state[person[ceased]] <- 2L
  }
  d$restrict_risk <- restrict_risk
  d$cease_risk <- cease_risk
  d$restrict_event <- restrict_event
  d$cease_event <- cease_event
  d
}

calibrate_rates <- function(d, targets, age_slope, n_cal=20L, steps=5L) {
  person_years <- aggregate(dt~cohort,d,sum)
  rates <- matrix(NA_real_,nrow(targets),2,
                  dimnames=list(targets$cohort,c('restriction','cessation')))
  for (k in seq_len(nrow(targets))) {
    cy <- person_years$dt[match(targets$cohort[k],person_years$cohort)]
    rates[k,] <- c(targets$restriction_events[k],
                   targets$cessation_events[k])/cy
  }
  for (step in seq_len(steps)) {
    count <- matrix(0,nrow(rates),2,dimnames=dimnames(rates))
    for (b in seq_len(n_cal)) {
      sim <- simulate_transitions(d,rates,c(restriction=1,cessation=1),age_slope)
      for (k in rownames(rates)) {
        x <- sim$cohort==k
        count[k,] <- count[k,]+c(sum(sim$restrict_event[x]),
                                 sum(sim$cease_event[x]))
      }
    }
    goal <- as.matrix(targets[,c('restriction_events','cessation_events')])
    observed <- pmax(count/n_cal,1)
    rates <- rates * pmax(0.4,pmin(2.5,goal/observed))
  }
  rates
}

fit_one <- function(sim,cohort,outcome,deff) {
  risk <- if (outcome=='restriction') 'restrict_risk' else 'cease_risk'
  event <- if (outcome=='restriction') 'restrict_event' else 'cease_event'
  a <- sim[sim$cohort==cohort & sim[[risk]]==1L,]
  if (sum(a[[event]])<8L || length(unique(a$future_adrd))<2L)
    return(c(events=sum(a[[event]]),HR=NA_real_,p=NA_real_))
  a$y <- a[[event]];a$age_c <- a$age_mid-75
  fit <- tryCatch(suppressWarnings(glm(
    y ~ future_adrd + age_c + I(age_c^2) + offset(log(dt)),
    data=a,family=binomial(link='cloglog'))),error=function(e) NULL)
  if (is.null(fit)) return(c(events=sum(a$y),HR=NA_real_,p=NA_real_))
  cf <- summary(fit)$coefficients
  if (!('future_adrd' %in% rownames(cf)))
    return(c(events=sum(a$y),HR=NA_real_,p=NA_real_))
  beta <- unname(cf['future_adrd','Estimate'])
  # deff[cohort] carries a name (e.g., "HRS"). Remove it so c(p=...)
  # produces a component named exactly "p", rather than "p.HRS".
  se <- unname(cf['future_adrd','Std. Error']) *
    sqrt(unname(deff[cohort]))
  out <- setNames(c(sum(a$y),exp(beta),2*pnorm(-abs(beta/se))),
                  c('events','HR','p'))
  stopifnot(identical(names(out),c('events','HR','p')))
  out
}

run_power <- function(d,rates,effects,B,age_slope,deff) {
  results <- vector('list',length(effects)*B*2L*length(rownames(rates)))
  pos <- 0L
  for (hr in effects) {
    for (b in seq_len(B)) {
      sim <- simulate_transitions(d,rates,
                                  c(restriction=hr,cessation=hr),age_slope)
      for (outcome in c('restriction','cessation')) {
        for (cohort in rownames(rates)) {
          ans <- fit_one(sim,cohort,outcome,deff)
          pos <- pos+1L
          results[[pos]] <- data.frame(true_HR=hr,replicate=b,
                                       cohort=cohort,outcome=outcome,events=ans['events'],
                                       fitted_HR=ans['HR'],p=ans['p'])
        }
      }
    }
    cat('Finished true HR =',hr,'\n')
  }
  results <- do.call(rbind,results)
  if (any(!is.finite(results$p[results$cohort %in% c('HRS','NHATS')])))
    stop('Missing HRS/NHATS p-values: inspect fit_one() before interpreting power.')
  results$significant_correct <- with(results,
                                      !is.na(p) & p<0.05 & !is.na(fitted_HR) & fitted_HR>1)
  summary <- aggregate(cbind(events,significant_correct)~true_HR+outcome+cohort,
                       results,mean,na.action=na.pass)
  names(summary)[names(summary)=='significant_correct'] <- 'power'
  list(replicates=results,summary=summary)
}

visits <- make_demo_visits(design,future_adrd_fraction)
# To use real visit data, replace the preceding line with:
# visits <- read.csv('aim1_visit_skeleton.csv',stringsAsFactors=FALSE)
d <- make_intervals(visits)
for (k in seq_len(nrow(design))) {
  cohort <- design$cohort[k]
  stopifnot(sum(visits$cohort==cohort)==design$person_observations[k],
            length(unique(visits$id[visits$cohort==cohort]))==design$n[k],
            sum(table(visits$id[visits$cohort==cohort])>=2)==
              design$persons_with_2plus[k],
            abs(sum(d$cohort==cohort)-design$transition_pairs[k])<=1L)
}
cat('Table 2 check (HRS reports one fewer pair than observations minus persons):\n')
print(data.frame(cohort=design$cohort,table_pairs=design$transition_pairs,
                 simulated_pairs=as.integer(table(d$cohort)[design$cohort]),
                 recovery_percent=design$recovery_percent),row.names=FALSE)
rates <- calibrate_rates(d,design,age_slope)
cat('Calibrated annual hazards at age 75 (before ADRD group effect):\n')
print(round(rates,4))
calibration_counts <- replicate(30L, {
  x <- simulate_transitions(d,rates,c(restriction=1,cessation=1),age_slope)
  unlist(lapply(design$cohort,function(g) c(
    restriction=sum(x$restrict_event[x$cohort==g]),
    cessation=sum(x$cease_event[x$cohort==g]))))
})
calibration <- data.frame(
  cohort=rep(design$cohort,each=2),
  outcome=rep(c('restriction','cessation'),times=nrow(design)),
  table_events=as.vector(t(as.matrix(design[,c('restriction_events',
                                               'cessation_events')]))),
  mean_simulated_events=rowMeans(calibration_counts))
cat('Null-group calibration against Table 2 event counts:\n')
print(calibration,row.names=FALSE)
write.csv(calibration,'CAR_Aim1_transition_calibration.csv',row.names=FALSE)
ans <- run_power(d,rates,effects,B,age_slope,design_effect)
cat('First fitted results (check that p-values are numeric):\n')
print(head(ans$replicates[,c('cohort','outcome','fitted_HR','p')]))
print(ans$summary,row.names=FALSE)
write.csv(ans$summary,'CAR_Aim1_transition_power_summary.csv',row.names=FALSE)
write.csv(ans$replicates,'CAR_Aim1_transition_power_replicates.csv',row.names=FALSE)
agreement <- reshape(ans$replicates[,c('true_HR','replicate','outcome',
                                       'cohort','fitted_HR','p')],
                     idvar=c('true_HR','replicate','outcome'),timevar='cohort',direction='wide')
agreement$same_direction <- with(agreement,
                                 is.finite(fitted_HR.HRS) & is.finite(fitted_HR.NHATS) &
                                   fitted_HR.HRS>1 & fitted_HR.NHATS>1)
agreement$HRS_significant_NHATS_agrees <- with(agreement,
                                               same_direction & is.finite(p.HRS) & p.HRS<0.05)
agreement_summary <- aggregate(
  cbind(same_direction,HRS_significant_NHATS_agrees)~true_HR+outcome,
  agreement,mean)
print(agreement_summary,row.names=FALSE)
write.csv(agreement_summary,'CAR_Aim1_transition_agreement_summary.csv',
          row.names=FALSE)

# Interpret power separately for restriction and cessation. The HRS/NHATS/ACTIVE
# results are cohort-specific. These simulated values are planning estimates:
# the Table 2 recovery percentages describe reversals AFTER restriction and
# cannot inform power for the FIRST restriction event. They are shown in the
# Table 2 check but do not enter the first-transition model. Table 2 also lacks
# ages, ADRD case fraction/diagnosis timing, death dates, and design effects.
# Replace those assumptions with observed values before reporting grant power.


# ============================================================================
# Coordinated meta-analytic power
#
# Requires ans$replicates to contain:
#   true_HR, replicate, outcome, cohort, fitted_HR, p
#
# The code reconstructs the standard error of log(HR) from the fitted HR and
# two-sided p-value, then pools cohort-specific log hazard ratios within each
# simulation replicate.
# ============================================================================

required_columns <- c(
  "true_HR", "replicate", "outcome",
  "cohort", "fitted_HR", "p"
)

stopifnot(all(required_columns %in% names(ans$replicates)))

meta_input <- ans$replicates |>
  dplyr::mutate(
    beta = log(fitted_HR),

    # Recover the absolute Wald statistic from the two-sided p-value.
    p_bounded = pmin(pmax(p, 1e-15), 1 - 1e-15),
    z_abs = stats::qnorm(1 - p_bounded / 2),

    # Wald SE of log(HR). Estimates too close to zero cannot be recovered
    # reliably from rounded p-values and are excluded from pooling.
    se = dplyr::if_else(
      is.finite(beta) &
        is.finite(z_abs) &
        z_abs > 0 &
        abs(beta) > 1e-10,
      abs(beta) / z_abs,
      NA_real_
    ),
    variance = se^2
  ) |>
  dplyr::filter(
    is.finite(beta),
    is.finite(variance),
    variance > 0
  )


# ----------------------------------------------------------------------------
# Function to pool one simulated replicate
#
# Returns:
#   fixed-effect inverse-variance estimate
#   DerSimonian-Laird random-effects estimate
#   heterogeneity statistics
#
# The random-effects results should be treated as a sensitivity analysis when
# only HRS and NHATS contribute because tau-squared is imprecisely estimated
# from two cohorts.
# ----------------------------------------------------------------------------

pool_one_replicate <- function(x) {

  k <- nrow(x)

  if (k < 2) {
    return(tibble::tibble(
      k = k,
      pooled_HR_fixed = NA_real_,
      se_fixed = NA_real_,
      p_fixed = NA_real_,
      pooled_HR_random = NA_real_,
      se_random = NA_real_,
      p_random = NA_real_,
      Q = NA_real_,
      tau2 = NA_real_,
      I2 = NA_real_,
      same_direction = NA
    ))
  }

  # Fixed-effect inverse-variance model
  w_fixed <- 1 / x$variance

  beta_fixed <- sum(w_fixed * x$beta) / sum(w_fixed)
  se_fixed <- sqrt(1 / sum(w_fixed))
  z_fixed <- beta_fixed / se_fixed
  p_fixed <- 2 * stats::pnorm(-abs(z_fixed))

  # DerSimonian-Laird estimate of between-cohort heterogeneity
  Q <- sum(w_fixed * (x$beta - beta_fixed)^2)
  df_Q <- k - 1

  C_value <- sum(w_fixed) -
    sum(w_fixed^2) / sum(w_fixed)

  tau2 <- if (
    is.finite(C_value) &&
      C_value > 0
  ) {
    max(0, (Q - df_Q) / C_value)
  } else {
    0
  }

  # Random-effects model
  w_random <- 1 / (x$variance + tau2)

  beta_random <- sum(w_random * x$beta) / sum(w_random)
  se_random <- sqrt(1 / sum(w_random))
  z_random <- beta_random / se_random
  p_random <- 2 * stats::pnorm(-abs(z_random))

  I2 <- if (Q > 0) {
    max(0, 100 * (Q - df_Q) / Q)
  } else {
    0
  }

  tibble::tibble(
    k = k,
    pooled_HR_fixed = exp(beta_fixed),
    se_fixed = se_fixed,
    p_fixed = p_fixed,
    pooled_HR_random = exp(beta_random),
    se_random = se_random,
    p_random = p_random,
    Q = Q,
    tau2 = tau2,
    I2 = I2,
    same_direction = all(x$beta > 0)
  )
}


# ----------------------------------------------------------------------------
# Function to run the replicate-level meta-analysis for a selected cohort set
# ----------------------------------------------------------------------------

run_meta_power <- function(data, included_cohorts, analysis_name) {

  replicate_results <- data |>
    dplyr::filter(cohort %in% included_cohorts) |>
    dplyr::group_by(true_HR, replicate, outcome) |>
    dplyr::group_modify(~ pool_one_replicate(.x)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      analysis = analysis_name,
      significant_fixed = (
        is.finite(p_fixed) &
          p_fixed < 0.05 &
          pooled_HR_fixed > 1
      ),
      significant_random = (
        is.finite(p_random) &
          p_random < 0.05 &
          pooled_HR_random > 1
      )
    )

  summary_results <- replicate_results |>
    dplyr::group_by(analysis, true_HR, outcome) |>
    dplyr::summarise(
      completed_replicates = sum(is.finite(p_fixed)),

      mean_pooled_HR_fixed = mean(
        pooled_HR_fixed,
        na.rm = TRUE
      ),

      fixed_effect_power = mean(
        significant_fixed,
        na.rm = TRUE
      ),

      mean_pooled_HR_random = mean(
        pooled_HR_random,
        na.rm = TRUE
      ),

      random_effects_power = mean(
        significant_random,
        na.rm = TRUE
      ),

      same_direction_probability = mean(
        same_direction,
        na.rm = TRUE
      ),

      median_tau2 = median(
        tau2,
        na.rm = TRUE
      ),

      median_I2 = median(
        I2,
        na.rm = TRUE
      ),

      .groups = "drop"
    )

  list(
    replicates = replicate_results,
    summary = summary_results
  )
}


# Primary coordinated analysis: national cohorts
meta_national <- run_meta_power(
  data = meta_input,
  included_cohorts = c("HRS", "NHATS"),
  analysis_name = "HRS + NHATS"
)


# Sensitivity analysis: add the ACTIVE control arm
meta_all <- run_meta_power(
  data = meta_input,
  included_cohorts = c("HRS", "NHATS", "ACTIVE"),
  analysis_name = "HRS + NHATS + ACTIVE"
)


# Combine results
meta_power_summary <- dplyr::bind_rows(
  meta_national$summary,
  meta_all$summary
) |>
  dplyr::arrange(
    analysis,
    outcome,
    true_HR
  )

meta_power_replicates <- dplyr::bind_rows(
  meta_national$replicates,
  meta_all$replicates
)


cat("\nMeta-analytic power summary:\n")
print(
  meta_power_summary,
  n = Inf
)


# Save results
readr::write_csv(
  meta_power_summary,
  file.path(
    OUT_DIR,
    "CAR_Aim1_meta_analytic_power_summary.csv"
  )
)

readr::write_csv(
  meta_power_replicates,
  file.path(
    OUT_DIR,
    "CAR_Aim1_meta_analytic_power_replicates.csv"
  )
)


# ----------------------------------------------------------------------------
# Grant-ready values for HR = 1.15
# ----------------------------------------------------------------------------

grant_power <- meta_power_summary |>
  dplyr::filter(
    analysis == "HRS + NHATS",
    true_HR == 1.15
  ) |>
  dplyr::select(
    outcome,
    fixed_effect_power,
    random_effects_power,
    same_direction_probability,
    median_I2
  )

cat("\nGrant-ready coordinated power at HR = 1.15:\n")
print(grant_power, n = Inf)


## =============================================================================
## CAR Study -- aim1_active_standalone.R
##
## Aim 1 in ACTIVE only. Self-contained: no dependency on 00_config or on the
## coordinated build. Reads active_person_occasion.csv and produces every Aim 1
## result ACTIVE can support.
##
## Why ACTIVE runs on its own:
##   - Control arm only (INTGRPR == 0), per the UAB scope of work.
##   - No complex survey design, so models are unweighted with cluster-robust
##     standard errors rather than svyglm.
##   - Exact ages, unlike NHATS five-year bands.
##   - Occasions at years 0, 1, 2, 3, 5, 10, so intervals are unequal by
##     design and the interval-length offset matters more than elsewhere.
##   - ACTIVE alone carries graded continuous driving measures (days per week,
##     driving space, avoidance) and retrospective cessation dating, which the
##     survey cohorts do not have.
##
## Outputs (to _results/active/):
##   01_sample_flow.csv            eligibility cascade
##   02_transition_matrix.csv      observed state transitions
##   03_first_transition_events.csv event counts and person-time
##   04_hazard_models.csv          cause-specific cloglog models
##   05_cumulative_incidence.csv   competing-risks CIF by ADRD status
##   06_restricted_mean_years.csv  difference in years free, bootstrapped
##   07_continuous_trajectories.csv mixed models on graded measures
##   08_sensitivity.csv            confirmed events, drop-last, arm check
##   09_cessation_dating.csv       retrospective dating cross-check
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(stringr)
  library(readr); library(tibble); library(splines)
  library(survival); library(sandwich); library(lmtest)
  library(lme4); library(lmerTest)
})

ROOT   <- "C:/Users/trbell/Documents/Lab/Research/CARS"
INFILE <- file.path(ROOT, "ACTIVE", "active_person_occasion.csv")
OUT    <- file.path(ROOT, "_results", "active")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SEED <- 20260924
set.seed(SEED)

CONTROL_ONLY <- TRUE       # INTGRPR == 0; set FALSE only to check arm effects
N_BOOT       <- 1000
AGE_GRID     <- 70:92

num <- function(x) suppressWarnings(as.numeric(x))

## =============================================================================
## 1. Read and derive the driving state
##
## CURDRIV0  1 = currently drives, 2 = does not
## EVERDRIV  2 = never driven (person excluded entirely)
## NIGHTDRV / ALONDRIV / RAINDRIV are phrased "have you driven ...", so
##   1 = has driven in that situation, 2 = avoids it
##
## State, ascending severity:
##   0 unrestricted  drives, no avoidance endorsed
##   1 restricted    drives, one or more avoidance endorsed
##   2 ceased        not currently driving (absorbing)
## =============================================================================

raw <- read_csv(INFILE, show_col_types = FALSE)
stopifnot(all(c("ID", "occasion", "year_nominal", "age") %in% names(raw)))

flow <- tibble(step = "rows in source file", n_rows = nrow(raw),
               n_persons = n_distinct(raw$ID))

## The source file fans out to six rows per person; occasion distinguishes them
stopifnot(nrow(raw) == 6 * n_distinct(raw$ID))

d0 <- raw %>%
  transmute(
    id        = num(ID),
    occasion  = as.integer(occasion),
    year      = num(year_nominal),
    age       = num(age),
    arm       = num(INTGRPR),
    drives    = num(CURDRIV0),
    everdrove = num(EVERDRIV),
    av_night  = num(NIGHTDRV),
    av_alone  = num(ALONDRIV),
    av_rain   = num(RAINDRIV),
    prompted  = num(LIMITDRV),
    mci_ever  = num(MCI_EVER_2),
    mci_age   = num(AGE_FIRST_MCI),
    sex       = num(GENDER),
    race      = num(RACE_CAT),
    educ      = num(EDUCLEVL),
    age_base  = num(AGEB),
    comorb    = num(COMORBIDITY_B)
  )

never_ids <- d0 %>% filter(everdrove == 2) %>% pull(id) %>% unique()

d <- d0 %>%
  { if (CONTROL_ONLY) filter(., arm == 0) else . } %>%
  filter(!id %in% never_ids,
         drives %in% c(1, 2),
         !is.na(age)) %>%
  mutate(
    n_avoid = rowSums(cbind(av_night == 2, av_alone == 2, av_rain == 2), na.rm = TRUE),
    state   = case_when(drives == 2 ~ 2L, n_avoid > 0 ~ 1L, TRUE ~ 0L)
  ) %>%
  arrange(id, occasion)

flow <- bind_rows(flow,
                  tibble(step = if (CONTROL_ONLY) "control arm (INTGRPR == 0)" else "all arms",
                         n_rows = nrow(d0 %>% { if (CONTROL_ONLY) filter(., arm == 0) else . }),
                         n_persons = n_distinct(d0$id[if (CONTROL_ONLY) d0$arm == 0 else TRUE])),
                  tibble(step = "excluding never-drivers", n_rows = NA_integer_,
                         n_persons = n_distinct(d$id)),
                  tibble(step = "person-occasions with a valid driving state",
                         n_rows = nrow(d), n_persons = n_distinct(d$id)))

## =============================================================================
## 2. ADRD status
##
## Until CMS linkage is complete, incident MCI stands in for incident
## CMS-ascertained ADRD. AGE_FIRST_MCI gives an exact onset age, which no
## other cohort provides. Swap in the CMS diagnosis age when linkage lands:
## only the two lines below change.
## =============================================================================

ADRD_SOURCE <- "MCI"     # "MCI" interim, "CMS" once linked

person <- d %>%
  group_by(id) %>%
  summarise(across(c(sex, race, educ, age_base, comorb, mci_ever, mci_age),
                   ~suppressWarnings(first(na.omit(.x)))),
            age_first = min(age), age_last = max(age), n_occ = n(),
            .groups = "drop") %>%
  mutate(
    age_dx = if (ADRD_SOURCE == "MCI") mci_age else NA_real_,
    adrd   = as.integer(!is.na(age_dx)),
    sex      = factor(sex, 1:2, c("male", "female")),
    race_eth = factor(race, 1:3, c("white_nh", "black_nh", "other_nh"))
  )

cat("\nincident", ADRD_SOURCE, "cases:", sum(person$adrd),
    sprintf(" (%.1f%% of %d)\n", 100 * mean(person$adrd), nrow(person)))

## Only prediagnostic observations contribute for cases
d <- d %>%
  left_join(person %>% select(id, age_dx, adrd), by = "id") %>%
  filter(is.na(age_dx) | age < age_dx)

flow <- bind_rows(flow,
                  tibble(step = "prediagnostic observations only", n_rows = nrow(d),
                         n_persons = n_distinct(d$id)))
write_csv(flow, file.path(OUT, "01_sample_flow.csv"))

## =============================================================================
## 3. Person-intervals and the observed transition matrix
## =============================================================================

iv <- d %>%
  arrange(id, age) %>%
  group_by(id) %>%
  mutate(ceased_before = cumsum(lag(state, default = 0L) == 2L)) %>%
  filter(ceased_before == 0) %>%
  mutate(state_next = lead(state), age_next = lead(age), dt = lead(age) - age) %>%
  ungroup() %>%
  filter(!is.na(state_next), dt > 0) %>%
  transmute(id, age_start = age, age_end = age_next, dt, log_dt = log(dt),
            from = state, to = state_next, adrd, age_dx) %>%
  left_join(person %>% select(id, sex, race_eth, educ, age_base, comorb), by = "id")

tm <- table(from = iv$from, to = iv$to)
tm_pct <- round(100 * prop.table(tm, 1), 1)
cat("\nTransition matrix, row %:\n"); print(tm_pct)

as.data.frame(tm) %>%
  rename(n = Freq) %>%
  left_join(as.data.frame(tm_pct) %>% rename(pct = Freq), by = c("from", "to")) %>%
  write_csv(file.path(OUT, "02_transition_matrix.csv"))

cat("\ninterval lengths:\n"); print(table(round(iv$dt, 1)))

## =============================================================================
## 4. First-transition analytic files
##
## Restriction: entrants unrestricted at first eligible occasion; event is the
##   first subsequent report of any restriction. Cessation before restriction
##   is a competing event.
## Cessation:  entrants are current drivers; event is the first report of not
##   currently driving.
## Participants already reporting the outcome at entry are excluded.
## =============================================================================

make_first_transition <- function(iv, outcome = c("rest", "cease")) {
  outcome <- match.arg(outcome)
  entry_ok <- if (outcome == "rest") function(s) s == 0L else function(s) s < 2L
  is_ev    <- if (outcome == "rest") function(s) s == 1L else function(s) s == 2L
  
  iv %>%
    arrange(id, age_start) %>%
    group_by(id) %>%
    filter(entry_ok(first(from))) %>%
    mutate(event = as.integer(is_ev(to)),
           comp  = if (outcome == "rest") as.integer(to == 2L) else 0L,
           cum   = cumsum(event + comp)) %>%
    filter(cum == 0 | (cum == 1 & (event == 1 | comp == 1))) %>%
    ungroup()
}

ft <- list(rest  = make_first_transition(iv, "rest"),
           cease = make_first_transition(iv, "cease"))

event_summary <- imap_dfr(ft, function(x, nm) {
  tibble(outcome = nm,
         persons = n_distinct(x$id), intervals = nrow(x),
         events = sum(x$event), competing = sum(x$comp),
         person_years = round(sum(x$dt), 1),
         rate_per_1000 = round(1000 * sum(x$event) / sum(x$dt), 1),
         events_adrd = sum(x$event[x$adrd == 1]),
         events_ctrl = sum(x$event[x$adrd == 0]))
})
print(as.data.frame(event_summary), row.names = FALSE)
write_csv(event_summary, file.path(OUT, "03_first_transition_events.csv"))

## =============================================================================
## 5. Cause-specific discrete-time hazard models
##
##   log[-log(1 - p_ij)] = log(dt) + alpha(age) + beta * ADRD + gamma' X
##
## Unweighted with cluster-robust SEs at the participant level. The age spline
## uses 2 df rather than 3 because ACTIVE's event counts are modest.
## =============================================================================

AGE_SPEC <- "ns(age_start, df = 2)"
RHS      <- "adrd + age_base + sex + race_eth + educ"

fit_cs <- function(dat, event_col) {
  f <- as.formula(sprintf("%s ~ %s + %s + offset(log_dt)", event_col, AGE_SPEC, RHS))
  m <- glm(f, data = dat, family = binomial(link = "cloglog"))
  list(fit = m, ct = coeftest(m, vcov. = vcovCL(m, cluster = dat$id)))
}

haz <- imap(ft, function(x, nm) {
  if (sum(x$event) < 15) { message("too few events for ", nm); return(NULL) }
  fit_cs(x, "event")
}) %>% compact()

haz_tab <- imap_dfr(haz, function(h, nm) {
  s <- h$ct
  tibble(outcome = nm, term = rownames(s), est = s[, 1], se = s[, 2],
         z = s[, 3], p = s[, 4]) %>%
    mutate(hr = exp(est), hr_lo = exp(est - 1.96 * se), hr_hi = exp(est + 1.96 * se))
})
write_csv(haz_tab, file.path(OUT, "04_hazard_models.csv"))
cat("\nADRD hazard ratios:\n")
haz_tab %>% filter(term == "adrd") %>%
  select(outcome, hr, hr_lo, hr_hi, p) %>% as.data.frame() %>% print(row.names = FALSE)

## =============================================================================
## 6. Competing-risks cumulative incidence
##
## Cause-specific hazards are converted to CIFs on a one-year age grid:
##   h_k(a) = 1 - exp(-exp(eta_k(a)))
##   S(a)   = prod_{u <= a} (1 - sum_k h_k(u))
##   CIF_k(a) = sum_{u <= a} h_k(u) * S(u-1)
## Competing events: cessation before reported restriction, and death where
## a death age is available.
## =============================================================================

pred_h <- function(m, nd, ages) {
  nd2 <- nd[rep(1, length(ages)), , drop = FALSE]
  nd2$age_start <- ages
  nd2$log_dt <- 0
  tibble(age = ages,
         h = 1 - exp(-exp(as.numeric(predict(m, newdata = nd2, type = "link")))))
}

ref_row <- function(dat) {
  tibble(age_base = mean(dat$age_base, na.rm = TRUE),
         sex      = factor(names(sort(table(dat$sex), decreasing = TRUE))[1],
                           levels = levels(dat$sex)),
         race_eth = factor(names(sort(table(dat$race_eth), decreasing = TRUE))[1],
                           levels = levels(dat$race_eth)),
         educ     = mean(dat$educ, na.rm = TRUE))
}

cif_for <- function(dat, m_main, m_comp = NULL, ages = AGE_GRID) {
  map_dfr(c(0, 1), function(g) {
    nd <- ref_row(dat) %>% mutate(adrd = g)
    H  <- list(main = pred_h(m_main, nd, ages))
    if (!is.null(m_comp)) H$comp <- pred_h(m_comp, nd, ages)
    M <- as.matrix(as.data.frame(map(H, "h")))
    h_all  <- rowSums(M)
    S_prev <- c(1, head(cumprod(1 - h_all), -1))
    tibble(age = ages, adrd = g,
           cif_main = cumsum(M[, "main"] * S_prev),
           surv = cumprod(1 - h_all))
  })
}

cif_all <- imap_dfr(ft, function(x, nm) {
  if (is.null(haz[[nm]])) return(NULL)
  mc <- if (nm == "rest" && sum(x$comp) >= 15) fit_cs(x, "comp")$fit else NULL
  cif_for(x, haz[[nm]]$fit, mc) %>% mutate(outcome = nm)
})
write_csv(cif_all, file.path(OUT, "05_cumulative_incidence.csv"))

## =============================================================================
## 7. Restricted mean years free, with participant-level bootstrap
## =============================================================================

rmy_point <- function(x, nm) {
  if (sum(x$event) < 15) return(NA_real_)
  m <- try(fit_cs(x, "event")$fit, silent = TRUE)
  if (inherits(m, "try-error")) return(NA_real_)
  cf <- cif_for(x, m)
  lo <- max(min(x$age_start), min(AGE_GRID))
  hi <- min(max(x$age_end),   max(AGE_GRID))
  k  <- cf$age >= lo & cf$age <= hi
  sum(1 - cf$cif_main[k & cf$adrd == 1]) - sum(1 - cf$cif_main[k & cf$adrd == 0])
}

rmy <- imap_dfr(ft, function(x, nm) {
  pt  <- rmy_point(x, nm)
  ids <- unique(x$id)
  bs  <- map_dbl(seq_len(N_BOOT), function(b) {
    take <- tibble(id = sample(ids, length(ids), replace = TRUE))
    rmy_point(take %>% left_join(x, by = "id", relationship = "many-to-many"), nm)
  })
  tibble(outcome = nm, diff_years = pt,
         lo = quantile(bs, .025, na.rm = TRUE),
         hi = quantile(bs, .975, na.rm = TRUE),
         n_boot_ok = sum(!is.na(bs)))
})
print(as.data.frame(rmy), row.names = FALSE)
write_csv(rmy, file.path(OUT, "06_restricted_mean_years.csv"))

## Aalen-Johansen comparison
aj <- imap_dfr(ft, function(x, nm) {
  x2 <- x %>% mutate(status = factor(case_when(event == 1 ~ "outcome",
                                               comp  == 1 ~ "competing",
                                               TRUE ~ "censor"),
                                     levels = c("censor", "outcome", "competing")))
  f <- survfit(Surv(age_start, age_end, status) ~ adrd, data = x2, id = id)
  tibble(outcome = nm, age = f$time, cif = f$pstate[, "outcome"],
         stratum = rep(names(f$strata), f$strata))
})
write_csv(aj, file.path(OUT, "05b_aalen_johansen.csv"))

## =============================================================================
## 8. Complementary continuous measures
##
## Days driven per week, driving space, graded avoidance, standardized and
## oriented so higher values indicate greater restriction. The Age x ADRD
## interaction tests whether age-related change differs by subsequent
## diagnosis. Only prediagnostic measurements contribute for cases.
## =============================================================================

cont <- raw %>%
  transmute(id = num(ID), occasion = as.integer(occasion), age = num(age),
            arm = num(INTGRPR),
            days = num(DAYSDRIV), space = num(TOTDS), avoid = num(DAVOID),
            miles = num(MILEDRIV)) %>%
  { if (CONTROL_ONLY) filter(., arm == 0) else . } %>%
  filter(!id %in% never_ids) %>%
  left_join(person %>% select(id, adrd, age_dx, sex, race_eth, educ), by = "id") %>%
  filter(is.na(age_dx) | age < age_dx) %>%
  mutate(age_c = age - 75,
         days_z  = as.numeric(scale(-days)),
         space_z = as.numeric(scale(-space)),
         avoid_z = as.numeric(scale(avoid)),
         miles_z = as.numeric(scale(-miles)))

MEASURES <- c("days_z", "space_z", "avoid_z", "miles_z")

cont_fits <- map_dfr(MEASURES, function(y) {
  dd <- cont %>% filter(!is.na(.data[[y]]))
  if (nrow(dd) < 200) return(NULL)
  f <- as.formula(sprintf(
    "%s ~ age_c * adrd + sex + race_eth + educ + (1 + age_c | id)", y))
  m <- try(lmer(f, data = dd, REML = TRUE,
                control = lmerControl(optimizer = "bobyqa")), silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  s <- summary(m)$coefficients
  tibble(measure = y, term = rownames(s), est = s[, "Estimate"],
         se = s[, "Std. Error"], p = s[, "Pr(>|t|)"],
         n_obs = nrow(dd), n_id = n_distinct(dd$id))
})
write_csv(cont_fits, file.path(OUT, "07_continuous_trajectories.csv"))
cat("\nAge x ADRD interactions:\n")
cont_fits %>% filter(term == "age_c:adrd") %>%
  select(measure, est, se, p) %>% as.data.frame() %>% print(row.names = FALSE)

## Nonlinearity check
nonlin <- map_dfr(MEASURES, function(y) {
  dd <- cont %>% filter(!is.na(.data[[y]]))
  if (nrow(dd) < 200) return(NULL)
  m1 <- lmer(as.formula(sprintf("%s ~ age_c * adrd + (1 + age_c | id)", y)),
             data = dd, REML = FALSE)
  m2 <- lmer(as.formula(sprintf("%s ~ (age_c + I(age_c^2)) * adrd + (1 + age_c | id)", y)),
             data = dd, REML = FALSE)
  tibble(measure = y, aic_linear = AIC(m1), aic_quadratic = AIC(m2),
         lrt_p = anova(m1, m2)$`Pr(>Chisq)`[2])
})
write_csv(nonlin, file.path(OUT, "07b_nonlinearity.csv"))

## =============================================================================
## 9. Sensitivity
## =============================================================================

sens <- list()

## 9a. Require the outcome at two consecutive occasions
sens$confirmed <- map(ft, function(x) {
  x %>% group_by(id) %>%
    mutate(event = as.integer(event == 1 & lead(to, default = 9L) >= to)) %>%
    ungroup()
})

## 9b. Drop the last prediagnostic occasion
sens$drop_last <- map(ft, function(x) {
  x %>% group_by(id) %>% filter(row_number() < n()) %>% ungroup()
})

## 9c. Exclude the year 5 to 10 interval, the only gap longer than 2 years
sens$short_intervals <- map(ft, function(x) x %>% filter(dt <= 3))

sens_tab <- imap_dfr(sens, function(lst, lbl) {
  imap_dfr(lst, function(x, nm) {
    if (sum(x$event) < 15) return(NULL)
    h <- fit_cs(x, "event")
    s <- h$ct
    tibble(sensitivity = lbl, outcome = nm, term = rownames(s),
           est = s[, 1], se = s[, 2], p = s[, 4]) %>%
      filter(term == "adrd") %>%
      mutate(hr = exp(est), events = sum(x$event))
  })
})

## 9d. Arm check: does including the training arms change the ADRD estimate?
if (CONTROL_ONLY) {
  message("\nArm check skipped. Set CONTROL_ONLY <- FALSE and rerun to compare.")
} else {
  arm_check <- imap_dfr(ft, function(x, nm) {
    if (sum(x$event) < 15) return(NULL)
    m <- glm(as.formula(sprintf(
      "event ~ %s + adrd + factor(arm) + age_base + sex + race_eth + educ + offset(log_dt)",
      AGE_SPEC)), data = x, family = binomial(link = "cloglog"))
    s <- coeftest(m, vcov. = vcovCL(m, cluster = x$id))
    tibble(sensitivity = "all_arms_adjusted", outcome = nm,
           term = rownames(s), est = s[, 1], se = s[, 2], p = s[, 4]) %>%
      filter(str_detect(term, "adrd|arm")) %>% mutate(hr = exp(est))
  })
  sens_tab <- bind_rows(sens_tab, arm_check)
}
write_csv(sens_tab, file.path(OUT, "08_sensitivity.csv"))

## =============================================================================
## 10. Retrospective cessation dating
##
## LDRIVYER and LDRIVMON give years and months since last driving among
## non-drivers, so cessation can be dated rather than bracketed. Most reports
## are long-standing non-drivers (median about 6 years), so only those within
## roughly 2 years of the report plausibly date an event observed during
## follow-up. Used to check how much timing error the interval-censored
## treatment absorbs, not as a primary outcome.
## =============================================================================

if (all(c("LDRIVYER", "LDRIVMON") %in% names(raw))) {
  dating <- raw %>%
    transmute(id = num(ID), occasion = as.integer(occasion), age = num(age),
              arm = num(INTGRPR),
              yr = num(LDRIVYER), mo = num(LDRIVMON)) %>%
    { if (CONTROL_ONLY) filter(., arm == 0) else . } %>%
    filter(!is.na(yr) | !is.na(mo)) %>%
    mutate(since = coalesce(yr, 0) + coalesce(mo, 0) / 12,
           age_cease_reported = age - since) %>%
    filter(since >= 0, since < 60)
  
  ## Compare the reported cessation age with the interval bracketing it
  compare <- ft$cease %>% filter(event == 1) %>%
    select(id, age_start, age_end) %>%
    inner_join(dating %>% group_by(id) %>%
                 slice_min(since, n = 1, with_ties = FALSE) %>%
                 select(id, age_cease_reported, since), by = "id") %>%
    mutate(inside_interval = age_cease_reported >= age_start &
             age_cease_reported <= age_end,
           dist_from_midpoint = age_cease_reported - (age_start + age_end) / 2)
  
  cat("\ncessation dating cross-check:\n")
  cat("  events with a dating report:", nrow(compare), "\n")
  cat("  reported age inside the bracketing interval:",
      sum(compare$inside_interval, na.rm = TRUE),
      sprintf(" (%.1f%%)\n", 100 * mean(compare$inside_interval, na.rm = TRUE)))
  cat("  median |reported - interval midpoint|:",
      round(median(abs(compare$dist_from_midpoint), na.rm = TRUE), 2), "years\n")
  write_csv(compare, file.path(OUT, "09_cessation_dating.csv"))
} else {
  message("LDRIVYER / LDRIVMON not in the source file; dating check skipped.")
}

message("\nACTIVE Aim 1 complete -> ", OUT)


# =============================================================================
# AIM 3A: LATENT DRIVING PROFILES AND INCIDENT ADRD
# =============================================================================

# Planning landmark sample. Post-landmark ADRD events are varied because
# person-level CMS diagnosis data are not available for preliminary power.
N_HRS_LANDMARK <- 11500L
ADRD_EVENT_FRACTION_GRID <- c(0.10, 0.15, 0.20, 0.25)
HRS_DESIGN_EFFECT <- 1.5

PROFILE_PREVALENCE_GRID <- c(0.05, 0.10, 0.20)
ENTROPY_GRID <- c(0.70, 0.80, 0.90)
HR_GRID <- c(1.15, 1.25, 1.40)

# Sensitivity to residual heterogeneity is included for later coordinated
# analyses. For HRS alone, tau=0 is the relevant row.
TAU_GRID_AIM3A <- c(0, 0.05)

simulate_aim3a <- function(HR, profile_prevalence, entropy, tau,
                           events,
                           design_effect = HRS_DESIGN_EFFECT,
                           n_sim = B) {
  
  true_beta <- log(HR)
  
  # Approximate classification-error attenuation. Replace with pseudo-class
  # draws from the fitted model when posterior probabilities are available.
  observed_beta <- entropy * true_beta
  
  information <-
    (events / design_effect) *
    profile_prevalence * (1 - profile_prevalence)
  
  sampling_se <- 1 / sqrt(information)
  analysis_se <- sqrt(sampling_se^2 + tau^2)
  
  estimates <- rnorm(n_sim, mean = observed_beta, sd = sampling_se)
  z <- estimates / analysis_se
  p <- 2 * pnorm(-abs(z))
  detected <- p < ALPHA & estimates > 0
  pw <- mean(detected)
  
  tibble(
    power = pw,
    mcse = sqrt(pw * (1 - pw) / n_sim),
    attenuated_HR = exp(observed_beta),
    sampling_se = sampling_se,
    analysis_se = analysis_se
  )
}

aim3a_power <- crossing(
  HR = HR_GRID,
  profile_prevalence = PROFILE_PREVALENCE_GRID,
  entropy = ENTROPY_GRID,
  tau = TAU_GRID_AIM3A,
  adrd_event_fraction = ADRD_EVENT_FRACTION_GRID
) |>
  pmap_dfr(function(HR, profile_prevalence, entropy, tau,
                    adrd_event_fraction) {
    events <- round(N_HRS_LANDMARK * adrd_event_fraction)
    bind_cols(
      tibble(HR = HR, profile_prevalence = profile_prevalence,
             entropy = entropy, tau = tau,
             adrd_event_fraction = adrd_event_fraction),
      simulate_aim3a(HR, profile_prevalence, entropy, tau, events = events)
    )
  }) |>
  mutate(
    landmark_n = N_HRS_LANDMARK,
    adrd_events = round(landmark_n * adrd_event_fraction),
    design_effect = HRS_DESIGN_EFFECT
  ) |>
  arrange(adrd_event_fraction, profile_prevalence, entropy, tau, HR)

write_csv(aim3a_power, file.path(OUT_DIR, "Aim3A_profile_ADRD_power.csv"))

cat("\n=== AIM 3A: selected HRS scenarios ===\n")
print(
  aim3a_power |>
    filter(
      tau == 0,
      HR %in% c(1.15, 1.25),
      adrd_event_fraction == 0.20
    ) |>
    select(HR, adrd_event_fraction, adrd_events, profile_prevalence,
           entropy, power, attenuated_HR, mcse),
  n = Inf
)


# =============================================================================
# AIM 3B: HRS-HCAP BIOMARKER ASSOCIATIONS
# =============================================================================

# Plausible analytic sample sizes after intersecting driving, biomarker,
# landmark, and covariate eligibility. These remain sensitivity scenarios.
BIOMARKER_N_GRID <- c(1000L, 1500L, 2000L, 2500L, 3500L)
BIOMARKER_D_GRID <- seq(0.15, 0.50, by = 0.05)
BIOMARKER_PROFILE_PREVALENCE <- c(0.05, 0.10, 0.20)
BIOMARKER_ENTROPY <- c(0.70, 0.80, 0.90)

# Conservative family-wise alpha for three biomarkers. The planned BH-FDR
# analysis will generally be less conservative, but its exact power depends on
# correlations among A beta 42/40, p-tau181, and NfL.
ALPHA_BIOMARKER <- ALPHA / 3

normal_two_group_power <- function(n, prevalence, d, entropy, alpha) {
  attenuated_d <- entropy * d
  se <- sqrt(1 / (n * prevalence * (1 - prevalence)))
  ncp <- attenuated_d / se
  critical <- qnorm(1 - alpha / 2)
  pnorm(-critical - ncp) + 1 - pnorm(critical - ncp)
}

aim3b_power <- crossing(
  analytic_n = BIOMARKER_N_GRID,
  profile_prevalence = BIOMARKER_PROFILE_PREVALENCE,
  entropy = BIOMARKER_ENTROPY,
  standardized_difference = BIOMARKER_D_GRID
) |>
  mutate(
    expected_profile_n = analytic_n * profile_prevalence,
    power_nominal = pmap_dbl(
      list(analytic_n, profile_prevalence, standardized_difference, entropy),
      ~ normal_two_group_power(..1, ..2, ..3, ..4, ALPHA)
    ),
    power_bonferroni = pmap_dbl(
      list(analytic_n, profile_prevalence, standardized_difference, entropy),
      ~ normal_two_group_power(..1, ..2, ..3, ..4, ALPHA_BIOMARKER)
    )
  ) |>
  arrange(analytic_n, profile_prevalence, entropy, standardized_difference)

aim3b_detectable <- aim3b_power |>
  group_by(analytic_n, profile_prevalence, entropy) |>
  summarise(
    minimum_d_80_nominal = if (any(power_nominal >= 0.80))
      min(standardized_difference[power_nominal >= 0.80]) else NA_real_,
    minimum_d_80_bonferroni = if (any(power_bonferroni >= 0.80))
      min(standardized_difference[power_bonferroni >= 0.80]) else NA_real_,
    minimum_d_90_bonferroni = if (any(power_bonferroni >= 0.90))
      min(standardized_difference[power_bonferroni >= 0.90]) else NA_real_,
    .groups = "drop"
  )

write_csv(aim3b_power, file.path(OUT_DIR, "Aim3B_biomarker_power.csv"))
write_csv(aim3b_detectable, file.path(OUT_DIR, "Aim3B_detectable_effects.csv"))

cat("\n=== AIM 3B: detectable standardized differences ===\n")
print(aim3b_detectable, n = Inf)


# =============================================================================
# Compact grant-planning summary
# =============================================================================

grant_summary <- bind_rows(
  aim2_power |>
    filter(
      tau == 0,
      effect_or == 1.20,
      coverage_scenario == "expected",
      is.na(adrd_fraction) | adrd_fraction == 0.20
    ) |>
    transmute(
      aim = "Aim 2",
      scenario = paste(
        effect_type, domain, outcome, "OR=1.20",
        paste0("ADRD fraction=", if_else(is.na(adrd_fraction),
                                         "not applicable",
                                         format(adrd_fraction, nsmall = 2))),
        "coverage=expected", sep = " | "
      ),
      power = power,
      note = "HRS+NHATS inverse-variance coordinated estimate"
    ),
  aim3a_power |>
    filter(
      tau == 0,
      HR == 1.25,
      entropy == 0.80,
      adrd_event_fraction == 0.20
    ) |>
    transmute(
      aim = "Aim 3A",
      scenario = paste0("HR=1.25 | ADRD event fraction=0.20 | prevalence=",
                        profile_prevalence, " | entropy=0.80"),
      power = power,
      note = "HRS age-75 landmark planning scenario"
    ),
  aim3b_power |>
    filter(near(standardized_difference, 0.30), entropy == 0.80,
           profile_prevalence == 0.10) |>
    transmute(
      aim = "Aim 3B",
      scenario = paste0("d=0.30 | N=", analytic_n,
                        " | prevalence=0.10 | entropy=0.80"),
      power = power_bonferroni,
      note = "Conservative alpha=.05/3"
    )
)

write_csv(grant_summary, file.path(OUT_DIR, "Aims2_3_grant_power_summary.csv"))

cat("\n=== GRANT-PLANNING SUMMARY ===\n")
print(grant_summary, n = Inf)

cat("\nOutputs written to: ", OUT_DIR, "\n", sep = "")
cat("Simulations use Table 2 driving transitions and prespecified sensitivity ",
    "ranges for quantities requiring person-level CMS linkage.\n", sep = "")
