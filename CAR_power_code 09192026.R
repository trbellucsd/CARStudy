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



# =============================================================================
# CARS_Aims2_3_power_simulation.R
# Planning power for Aim 2, Aim 3A, and Aim 3B
# =============================================================================
#
# This script is deliberately separate from the data-construction pipelines.
# It uses the observed transition counts in Table 2 and explicit planning
# sensitivity ranges for quantities that Table 2 does not contain. Quantities
# requiring person-level CMS linkage are not treated as observed inputs.
#
# Aim 2:
#   Power for standardized within-person effects and their interactions with
#   incident ADRD, separately for restriction and cessation. HRS and NHATS
#   information is combined using inverse-variance meta-analysis. ADRD case
#   fraction and analytic coverage are varied in sensitivity analyses.
#
# Aim 3A:
#   Power for a latent driving profile predicting incident ADRD after the
#   age-75 landmark, allowing the post-landmark ADRD event fraction, profile
#   prevalence, classification uncertainty, and residual heterogeneity to vary.
#
# Aim 3B:
#   Power for standardized biomarker differences across plausible HRS-HCAP
#   eligibility counts, profile prevalences, and classification uncertainty.
#   A conservative alpha=.05/3 calculation is reported with nominal power.
#
# These are planning calculations. Entropy is used as an attenuation factor
# for classification error; entropy is not itself a misclassification rate.
# Final analyses should replace this approximation with pseudo-class draws.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
})

set.seed(20260920)

# ---- Output and simulation settings -----------------------------------------

OUT_DIR <- "C:/Users/trbell/Documents/Lab/Grants/CARS/power"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

B <- 10000L
ALPHA <- 0.05


# =============================================================================
# AIM 2
# =============================================================================

# Observed quantities from Table 2. ACTIVE is retained as a corroborative
# cohort and is not included in the primary HRS+NHATS meta-analytic power.
aim2_design <- tribble(
  ~cohort, ~transition_pairs, ~restriction_events, ~cessation_events, ~design_effect,
  "HRS",               44608,                4400,              2472,            1.5,
  "NHATS",             55110,                3663,              2526,            1.7,
  "ACTIVE",             1917,                 174,                58,            1.0
)

# Person-level CMS data are not available for preliminary power estimation.
# The interaction analysis therefore spans plausible incident ADRD fractions.
# The 0.20 row is the expected planning scenario used in the compact summary.
ADRD_FRACTION_GRID <- c(0.10, 0.15, 0.20, 0.25)

# Coverage is the assumed fraction of otherwise eligible transition intervals
# with a usable predictor and covariate set. These are planning values because
# the complete Aim 2 predictor files and CMS linkage are not yet available.
aim2_domains <- tribble(
  ~domain,                    ~coverage_HRS, ~coverage_NHATS, ~coverage_ACTIVE, ~within_information,
  "Cognitive function",               0.90,            0.90,             0.85,                0.70,
  "Physical function",                0.55,            0.75,             0.75,                0.65,
  "IADL function",                    0.90,            0.90,             0.85,                0.65,
  "Incident health events",           0.80,            0.90,             0.70,                0.75,
  "Medication burden",                0.55,            0.60,             0.55,                0.65,
  "Protective factors",               0.60,            0.75,             0.65,                0.60
)

# Vary the domain-specific coverage assumptions without treating currently
# unavailable CMS fields as observed missingness. Values are capped at 1.00.
coverage_scenarios <- tribble(
  ~coverage_scenario, ~coverage_multiplier,
  "conservative",                    0.80,
  "expected",                        1.00,
  "favorable",                       1.10
)

# within_information is a conservative planning multiplier for information
# retained after within-person centering, measurement error, and correlation
# among repeated observations.

MAIN_OR_GRID <- c(1.10, 1.15, 1.20, 1.25)
INTERACTION_OR_GRID <- c(1.10, 1.15, 1.20, 1.25)

# Prespecified sensitivity values for residual heterogeneity in log(OR).
# tau=0 corresponds to a common effect; tau=.05 is deliberately conservative.
TAU_GRID_AIM2 <- c(0, 0.03, 0.05)

coverage_long <- aim2_domains |>
  pivot_longer(
    starts_with("coverage_"),
    names_to = "cohort",
    names_prefix = "coverage_",
    values_to = "coverage"
  )

aim2_cells <- crossing(
  domain = aim2_domains$domain,
  cohort = aim2_design$cohort,
  outcome = c("restriction", "cessation"),
  adrd_fraction = ADRD_FRACTION_GRID,
  coverage_scenario = coverage_scenarios$coverage_scenario
) |>
  left_join(aim2_design, by = "cohort") |>
  left_join(
    coverage_long |>
      select(domain, cohort, coverage, within_information),
    by = c("domain", "cohort")
  ) |>
  left_join(coverage_scenarios, by = "coverage_scenario") |>
  mutate(
    base_coverage = coverage,
    coverage = pmin(base_coverage * coverage_multiplier, 1.00),
    events = if_else(
      outcome == "restriction",
      restriction_events,
      cessation_events
    ),
    event_probability = events / transition_pairs,
    usable_pairs = transition_pairs * coverage,
    
    # Approximate Fisher information for a standardized within-person term.
    information_main =
      (usable_pairs / design_effect) *
      event_probability * (1 - event_probability) *
      within_information,
    
    # Once the main effects are in the model, information for X_within*ADRD is
    # proportional to Var(ADRD)=p(1-p).
    information_interaction =
      information_main * adrd_fraction * (1 - adrd_fraction),
    
    se_main = 1 / sqrt(information_main),
    se_interaction = 1 / sqrt(information_interaction)
  )


# Simulate inverse-variance meta-analysis using known planning SEs. This is
# computationally light and makes Monte Carlo error negligible at B=10,000.
simulate_meta_power <- function(cell_data, effect_or, effect_type, tau,
                                n_sim = B) {
  
  stopifnot(effect_type %in% c("main", "interaction"))
  
  beta <- log(effect_or)
  se_name <- if (effect_type == "main") "se_main" else "se_interaction"
  
  d <- cell_data |>
    filter(cohort %in% c("HRS", "NHATS"))
  
  ses <- d[[se_name]]
  
  if (length(ses) != 2L) {
    stop(
      "Each Aim 2 simulation must contain exactly one HRS and one NHATS row; ",
      "received ", length(ses), ".",
      call. = FALSE
    )
  }
  if (any(!is.finite(ses))) {
    stop("Aim 2 planning standard errors must be finite.", call. = FALSE)
  }
  
  estimates <- sapply(ses, function(se) rnorm(n_sim, mean = beta, sd = se))
  weights <- 1 / (ses^2 + tau^2)
  pooled_beta <- as.vector(estimates %*% weights / sum(weights))
  pooled_se <- sqrt(1 / sum(weights))
  z <- pooled_beta / pooled_se
  p <- 2 * pnorm(-abs(z))
  detected <- p < ALPHA & pooled_beta > 0
  pw <- mean(detected)
  
  tibble(
    power = pw,
    mean_estimated_or = mean(exp(pooled_beta)),
    mean_se = pooled_se,
    mcse = sqrt(pw * (1 - pw) / n_sim)
  )
}


aim2_main <- crossing(
  domain = aim2_domains$domain,
  outcome = c("restriction", "cessation"),
  effect_or = MAIN_OR_GRID,
  tau = TAU_GRID_AIM2,
  coverage_scenario = coverage_scenarios$coverage_scenario
) |>
  mutate(effect_type = "within-person main effect") |>
  pmap_dfr(function(domain, outcome, effect_or, tau, coverage_scenario,
                    effect_type) {
    cells <- aim2_cells |>
      filter(
        .data$domain == .env$domain,
        .data$outcome == .env$outcome,
        .data$coverage_scenario == .env$coverage_scenario,
        .data$adrd_fraction == 0.20
      )
    bind_cols(
      tibble(domain = domain, outcome = outcome, effect_type = effect_type,
             effect_or = effect_or, tau = tau,
             coverage_scenario = coverage_scenario,
             adrd_fraction = NA_real_),
      simulate_meta_power(cells, effect_or, "main", tau)
    )
  })


aim2_interaction <- crossing(
  domain = aim2_domains$domain,
  outcome = c("restriction", "cessation"),
  effect_or = INTERACTION_OR_GRID,
  tau = TAU_GRID_AIM2,
  coverage_scenario = coverage_scenarios$coverage_scenario,
  adrd_fraction = ADRD_FRACTION_GRID
) |>
  mutate(effect_type = "within-person change x incident ADRD") |>
  pmap_dfr(function(domain, outcome, effect_or, tau, coverage_scenario,
                    adrd_fraction, effect_type) {
    cells <- aim2_cells |>
      filter(
        .data$domain == .env$domain,
        .data$outcome == .env$outcome,
        .data$coverage_scenario == .env$coverage_scenario,
        .data$adrd_fraction == .env$adrd_fraction
      )
    bind_cols(
      tibble(domain = domain, outcome = outcome, effect_type = effect_type,
             effect_or = effect_or, tau = tau,
             coverage_scenario = coverage_scenario,
             adrd_fraction = adrd_fraction),
      simulate_meta_power(cells, effect_or, "interaction", tau)
    )
  })

aim2_power <- bind_rows(aim2_main, aim2_interaction) |>
  arrange(effect_type, domain, outcome, coverage_scenario,
          adrd_fraction, tau, effect_or)

# Minimum OR on the evaluated grid that reaches 80% and 90% power.
aim2_thresholds <- aim2_power |>
  group_by(effect_type, domain, outcome, coverage_scenario,
           adrd_fraction, tau) |>
  summarise(
    minimum_or_80 = if (any(power >= 0.80, na.rm = TRUE))
      min(effect_or[!is.na(power) & power >= 0.80]) else NA_real_,
    minimum_or_90 = if (any(power >= 0.90, na.rm = TRUE))
      min(effect_or[!is.na(power) & power >= 0.90]) else NA_real_,
    .groups = "drop"
  )

write_csv(aim2_cells, file.path(OUT_DIR, "Aim2_planning_information.csv"))
write_csv(aim2_power, file.path(OUT_DIR, "Aim2_meta_analytic_power.csv"))
write_csv(aim2_thresholds, file.path(OUT_DIR, "Aim2_power_thresholds.csv"))

cat("\n=== AIM 2: HRS + NHATS coordinated power, tau=0 ===\n")
print(
  aim2_power |>
    filter(
      tau == 0,
      effect_or %in% c(1.15, 1.20),
      coverage_scenario == "expected",
      is.na(adrd_fraction) | adrd_fraction == 0.20
    ) |>
    select(effect_type, domain, outcome, effect_or, adrd_fraction,
           coverage_scenario, power, mcse),
  n = Inf
)


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
