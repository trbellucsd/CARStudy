# CAR Study: one-file Monte Carlo power simulations for Aims 1, 2, and 3.
# Run: source('CAR_Aims1_2_3_power_simulation.R')
# Output folders: CAR_power_results/Aim1, /Aim2, /Aim3.
# To run one aim: Sys.setenv(CAR_POWER_AIMS='2') before source().
# To check quickly: set CAR_AIM1_REPS, CAR_AIM2_REPS, CAR_AIM3_REPS
# to 5 and CAR_AIM3_DRAWS to 5. The Aim 2 default tests 48 scenarios.
# To vary planning scenarios: CAR_AIM2_FULL_GRID=true,
# CAR_AIM3_FULL_GRID=true, CAR_AIM3_OTHER_COHORTS=true.
# The three simulation programs appear in full below. No other .R script,
# user data file, or previously generated CSV is read at runtime.
# R packages required by Aim 2 and the survival package required by Aim 3
# must be installed in the R environment before running this file.
# All CMS and biomarker overlap parameters are planning assumptions.

CAR_POWER_AIMS <- trimws(strsplit(Sys.getenv('CAR_POWER_AIMS', '1,2,3'), ',')[[1L]])
if (length(CAR_POWER_AIMS) == 0L ||
    any(!CAR_POWER_AIMS %in% c('1', '2', '3'))) {
  stop('CAR_POWER_AIMS must contain one or more of 1, 2, 3.')
}
CAR_POWER_ROOT <- file.path(getwd(), 'CAR_power_results')
dir.create(CAR_POWER_ROOT, showWarnings = FALSE, recursive = TRUE)
run_power_aim <- function(number, code) {
  original_dir <- getwd()
  destination <- file.path(CAR_POWER_ROOT, paste0('Aim', number))
  dir.create(destination, showWarnings = FALSE, recursive = TRUE)
  on.exit(setwd(original_dir), add = TRUE)
  if (number == '2') {
    on.exit(if (requireNamespace('future', quietly = TRUE))
      future::plan(future::sequential), add = TRUE)
  }
  setwd(destination)
  cat('\nRunning Aim', number, 'in', destination, '\n')
  eval(substitute(code), envir = new.env(parent = globalenv()))
  invisible(destination)
}

if ('1' %in% CAR_POWER_AIMS) run_power_aim('1', {
# Aim 1 planning power: first reports before a matched ADRD index date.
# This Aim 1 section is embedded in the combined script; it uses base R.
# Pre-linkage visit and event counts come from Table 2; linkage, FFS coverage,
# diagnosis dates, baseline restriction, and mortality are explicit assumptions.
# Power tests the case-control transition coefficient. Fitted cumulative
# incidence ages are descriptive, conditional on driving at age 70 and on
# survival to the matched index. They are not separately powered contrasts.
# Set CAR_AIM1_REPS=20 for an initial run; default is 500 replicates per HR.

set.seed(20260918)
B <- as.integer(Sys.getenv('CAR_AIM1_REPS','500')) # Set to 20 for a local smoke check.
effects <- c(1.15, 1.25, 1.40)  # Case-control hazard ratios at the same age.
future_adrd_fraction <- 0.20     # ILLUSTRATIVE, among eligible linked persons.
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

# Legacy unrestricted-entry generator used only for initial event calibration.
# Final simulations below allow restriction at baseline and index matching.
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

# CURRENT AIM 1 DESIGN: incidence-density controls and shared index dates.
# All linkage, coverage, diagnosis, and mortality parameters below are planning
# assumptions. Table 2 contains unlinked driving counts, not analytic counts.
set.seed(20260925)
LINK_PROBABILITY <- 0.85
FFS_PROBABILITY <- 0.75
MATCH_AGE_TOLERANCE <- 2
CONTROLS_PER_CASE <- 2L
BASELINE_RESTRICTED_PROBABILITY <- 0.20
DEATH_HAZARD_75 <- 0.025
DEATH_AGE_SLOPE <- 0.08
CURVE_ORIGIN_AGE <- 70L
MAX_CURVE_AGE <- 90L
OUT_DIR <- 'CAR_Aim1_revised_power_results'
dir.create(OUT_DIR,showWarnings=FALSE,recursive=TRUE)
cohort_first_year <- c(HRS=2002,NHATS=2011,ACTIVE=1998)

# The original simulator always starts unrestricted. For cessation, restricted
# current drivers must also be at risk at their first eligible assessment.
# Use an illustrative baseline restricted fraction in the event generator.
simulate_transitions_current <- function(d, base_rates, hr_group, age_slope,
                                         baseline_restricted) {
  state <- as.integer(baseline_restricted)
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
    restricted <- r_risk & time_r<d$dt[j] & !ceased
    restrict_event[j] <- as.integer(restricted)
    cease_event[j] <- as.integer(ceased)
    state[person[restricted]] <- 1L
    state[person[ceased]] <- 2L
  }
  d$restrict_risk <- restrict_risk; d$cease_risk <- cease_risk
  d$restrict_event <- restrict_event;d$cease_event <- cease_event
  d
}

# The visit skeleton preserves Table 2 totals before claims eligibility.
visits <- make_demo_visits(design,0)
ordered <- visits[order(visits$cohort,visits$id,visits$age),]
people <- ordered[!duplicated(paste(ordered$cohort,ordered$id)),
                  c('cohort','id','age')]
people$person <- seq_len(nrow(people))
names(people)[names(people)=='age'] <- 'age0'
people$first_calendar <- unname(cohort_first_year[people$cohort])
people$n_visits <- as.integer(table(factor(
  match(paste(ordered$cohort,ordered$id),
        paste(people$cohort,people$id)),levels=people$person)))
people$last_calendar <- people$first_calendar +
  vapply(seq_len(nrow(people)),function(i)
    as.numeric(visit_years[[people$cohort[i]]][people$n_visits[i]]),numeric(1))
people$sex <- rbinom(nrow(people),1,.55)
people$education <- pmax(4,pmin(20,rnorm(nrow(people),13,3)))
people$race_indicator <- rbinom(nrow(people),1,.30)
stopifnot(all(people$n_visits>=1L))

# The age-dependent mortality draw determines survival to the case index.
# Because cases and controls must BOTH survive to the index, death cannot be
# observed within their retained pre-index histories. Do not interpret these
# fitted curves as competing-death cumulative incidence.
draw_people <- function(template) {
  p <- template
  h0 <- DEATH_HAZARD_75*exp(DEATH_AGE_SLOPE*(p$age0-75))
  years_to_death <- log1p(DEATH_AGE_SLOPE*rexp(nrow(p))/h0)/DEATH_AGE_SLOPE
  p$death_calendar <- p$first_calendar+years_to_death
  p$linked <- runif(nrow(p))<LINK_PROBABILITY
  p$ffs <- runif(nrow(p))<FFS_PROBABILITY
  p$claims_end <- pmin(p$death_calendar,
                       p$first_calendar+runif(nrow(p),7,22))
  p$baseline_restricted <- rbinom(nrow(p),1,
                                  BASELINE_RESTRICTED_PROBABILITY)
  p$diagnosis_calendar <- Inf
  candidates <- which(p$linked & p$ffs & p$n_visits>=2L &
    p$claims_end > p$first_calendar+3)
  chosen <- candidates[runif(length(candidates))<future_adrd_fraction]
  upper <- pmin(p$claims_end[chosen],p$last_calendar[chosen]+2)
  valid <- upper>p$first_calendar[chosen]+3
  chosen <- chosen[valid];upper <- upper[valid]
  p$diagnosis_calendar[chosen] <- p$first_calendar[chosen]+3+
    runif(length(chosen))*(upper-p$first_calendar[chosen]-3)
  p$future_adrd <- as.integer(is.finite(p$diagnosis_calendar))
  p
}

visits_before <- function(p,index) {
  years <- visit_years[[as.character(p$cohort[1L])]]
  pmin(p$n_visits,findInterval(index-p$first_calendar-1e-8,years))
}

match_windows <- function(p) {
  windows <- vector('list',sum(p$future_adrd));counter <- 0L
  for (g in design$cohort) {
    pool <- which(p$cohort==g & p$linked & p$ffs & p$n_visits>=2L)
    cases <- pool[p$future_adrd[pool]==1L]
    for (i in cases) {
      index <- p$diagnosis_calendar[i]
      if (visits_before(p[i,],index)<2L) next
      available <- pool[pool!=i & p$death_calendar[pool]>index &
        p$claims_end[pool]>=index & p$diagnosis_calendar[pool]>index]
      if (!length(available)) next
      case_age <- p$age0[i]+index-p$first_calendar[i]
      control_ages <- p$age0[available]+index-p$first_calendar[available]
      available <- available[abs(control_ages-case_age)<=MATCH_AGE_TOLERANCE &
        visits_before(p[available,],index)>=2L]
      if (!length(available)) next
      controls <- available[sample.int(length(available),
        min(CONTROLS_PER_CASE,length(available)))]
      counter <- counter+1L
      windows[[counter]] <- data.frame(
        person=c(i,controls),cohort=g,set_id=i,
        case=c(1L,rep(0L,length(controls))),index_calendar=index)
    }
  }
  if (!counter) return(NULL)
  do.call(rbind,windows[seq_len(counter)])
}

# Model with person-clustered uncertainty: controls may be used in multiple
# matched sets, and a person may later become a case.
fit_matched <- function(a,outcome,deff) {
  if (is.null(a) || !nrow(a)) return(NULL)
  at_risk <- if (outcome %in% c('restriction','direct_cessation'))
    a$restrict_risk==1L else a$cease_risk==1L
  x <- a[at_risk & a$dt>0,]
  x$y <- if (outcome=='restriction') x$restrict_event else x$cease_event
  if (sum(x$y)<8L || length(unique(x$case))<2L) return(NULL)
  x$age_c <- x$age_mid-75
  x$index_c <- x$index_calendar-2010
  model <- tryCatch(suppressWarnings(glm(
    y ~ case+age_c+I(age_c^2)+index_c+age0+sex+education+
      race_indicator+offset(log(dt)),data=x,
    family=binomial(link='cloglog'),y=TRUE)),error=function(e) NULL)
  if (is.null(model) || !isTRUE(model$converged) ||
      any(!is.finite(coef(model)))) return(NULL)
  v <- tryCatch({
    m <- model.matrix(model)
    eta <- pmax(-20,pmin(20,model$linear.predictors))
    mu <- pmin(1-1e-12,pmax(1e-12,fitted(model)))
    derivative <- exp(eta-exp(eta))
    score <- (model$y-mu)*derivative/pmax(mu*(1-mu),1e-12)
    clustered <- rowsum(m*score,x$person,reorder=FALSE)
    bread <- vcov(model)
    deff*bread%*%crossprod(clustered)%*%bread
  },error=function(e) NULL)
  if (is.null(v) || !is.finite(v['case','case']) || v['case','case']<=0)
    return(NULL)
  se <- sqrt(v['case','case']); beta <- unname(coef(model)['case'])
  list(model=model,hr=exp(beta),p=2*pnorm(-abs(beta/se)),
       events=sum(x$y),persons=length(unique(x$person)),risk=x)
}

age_curves <- function(fit,competing,outcome,g) {
  if (is.null(fit) || (outcome=='restriction' && is.null(competing)))
    return(NULL)
  x <- fit$risk
  ages <- CURVE_ORIGIN_AGE:MAX_CURVE_AGE
  ncase <- vapply(ages,function(age)
    sum(x$case==1L & x$age_mid>=age & x$age_mid<age+1),integer(1))
  ncontrol <- vapply(ages,function(age)
    sum(x$case==0L & x$age_mid>=age & x$age_mid<age+1),integer(1))
  supported <- ncase>=20L & ncontrol>=20L
  if (!supported[1L]) return(NULL)
  first_unsupported <- which(!supported)
  if (length(first_unsupported))
    ages <- ages[seq_len(first_unsupported[1L]-1L)]
  if (!length(ages)) return(NULL)
  new <- data.frame(age_c=ages+.5-75,index_c=mean(x$index_c),
    age0=mean(x$age0),sex=mean(x$sex),
    education=mean(x$education),race_indicator=mean(x$race_indicator),dt=1)
  curves <- vector('list',2L)
  for (case_value in 0:1) {
    new$case <- case_value
    lambda <- exp(predict(fit$model,newdata=new,type='link'))
    other <- if (outcome=='restriction')
      exp(predict(competing$model,newdata=new,type='link')) else
        rep(0,length(ages))
    surv <- 1; cif <- numeric(length(ages))
    for (j in seq_along(ages)) {
      total <- lambda[j]+other[j]
      event_p <- if (is.finite(total) && total>0)
        lambda[j]/total*(-expm1(-total)) else NA_real_
      cif[j] <- (if(j==1L) 0 else cif[j-1L])+surv*event_p
      surv <- surv*exp(-total)
    }
    curves[[case_value+1L]] <- data.frame(cohort=g,outcome=outcome,
      case=case_value,age=ages+1,cumulative_incidence=cif)
  }
  do.call(rbind,curves)
}

threshold_age <- function(curve,case_value,threshold) {
  if (is.null(curve)) return(NA_real_)
  age <- curve$age[curve$case==case_value &
    is.finite(curve$cumulative_incidence) &
    curve$cumulative_incidence>=threshold]
  if(length(age)) min(age) else NA_real_
}

# Calibration uses the original unlinked Table 2 driving skeleton. The new
# matched analysis subsequently applies mortality and claims selection.
d <- make_intervals(visits)
people_for_calibration <- draw_people(people)
d$future_adrd <- people_for_calibration$future_adrd[d$person]
rates <- calibrate_rates(d,design,age_slope)
# Calibrate first-event counts with baseline-restricted drivers as well.
for (step in seq_len(4L)) {
  counts <- matrix(0,nrow(design),2L,
                   dimnames=list(design$cohort,c('restriction','cessation')))
  for (b in seq_len(5L)) {
    simulated <- simulate_transitions_current(d,rates,
      c(restriction=1,cessation=1),age_slope,
      people_for_calibration$baseline_restricted)
    for (g in design$cohort) {
      part <- simulated$cohort==g
      counts[g,] <- counts[g,]+c(sum(simulated$restrict_event[part]),
                                   sum(simulated$cease_event[part]))
    }
  }
  target <- as.matrix(design[,c('restriction_events','cessation_events')])
  rates <- rates*pmax(.4,pmin(2.5,target/pmax(counts/5,1)))
}

# The original script's count check remains useful, but linked counts and
# matched controls will be smaller than these pre-linkage totals.
cat('Pre-linkage Table 2 counts (HRS has a one-pair discrepancy):\n')
print(data.frame(cohort=design$cohort,people=design$n,
                 intervals=as.integer(table(d$cohort)[design$cohort])))
print(round(rates,4))

rows <- curves <- list(); nrow_out <- ncurve <- 0L
for (effect in effects) {
  for (replicate in seq_len(B)) {
    p <- draw_people(people)
    d$future_adrd <- p$future_adrd[d$person]
    sim <- simulate_transitions_current(d,rates,
      c(restriction=effect,cessation=effect),age_slope,
      p$baseline_restricted)
    windows <- match_windows(p)
    analytic <- NULL
    if (!is.null(windows)) {
      sim$start_calendar <- people$first_calendar[sim$person]+sim$start-
        people$age0[sim$person]
      sim$end_calendar <- people$first_calendar[sim$person]+sim$end-
        people$age0[sim$person]
      analytic <- merge(windows,sim,by=c('person','cohort'))
      analytic <- analytic[analytic$end_calendar < analytic$index_calendar,]
      analytic$age0 <- p$age0[analytic$person]
      analytic$sex <- p$sex[analytic$person]
      analytic$education <- p$education[analytic$person]
      analytic$race_indicator <- p$race_indicator[analytic$person]
    }
    for (k in seq_len(nrow(design))) {
      g <- design$cohort[k]
      x <- if (!is.null(analytic))
        analytic[analytic$cohort==g,] else NULL
      main_r <- fit_matched(x,'restriction',design_effect[g])
      main_c <- fit_matched(x,'cessation',design_effect[g])
      direct <- fit_matched(x,'direct_cessation',design_effect[g])
      sets <- if (!is.null(windows))
        length(unique(windows$set_id[windows$cohort==g])) else 0L
      for (outcome in c('restriction','cessation')) {
        fit <- if (outcome=='restriction') main_r else main_c
        curve <- age_curves(fit,direct,outcome,g)
        nrow_out <- nrow_out+1L
        rows[[nrow_out]] <- data.frame(true_hr=effect,replicate=replicate,
          cohort=g,outcome=outcome,matched_sets=sets,
          fitted_hr=if(is.null(fit)) NA_real_ else fit$hr,
          p=if(is.null(fit)) NA_real_ else fit$p,
          events=if(is.null(fit)) NA_real_ else fit$events,
          persons=if(is.null(fit)) NA_real_ else fit$persons,
          threshold_case=threshold_age(curve,1L,
            if(outcome=='restriction') .10 else .05),
          threshold_control=threshold_age(curve,0L,
            if(outcome=='restriction') .10 else .05))
        if (!is.null(curve)) {
          ncurve <- ncurve+1L;curve$true_hr <- effect
          curve$replicate <- replicate;curves[[ncurve]] <- curve
        }
      }
    }
  }
  cat('Completed HR',effect,'with',B,'replicates\n')
}
replicates <- do.call(rbind,rows)
replicates$significant_correct <- is.finite(replicates$p) &
  replicates$p<.05 & is.finite(replicates$fitted_hr) &
  replicates$fitted_hr>1
groups <- split(replicates,
  list(replicates$true_hr, replicates$cohort, replicates$outcome),
  drop = TRUE)
summary_rows <- vector('list', length(groups))
for (i in seq_along(groups)) {
  x <- groups[[i]]
  valid_ages <- is.finite(x$threshold_case) &
    is.finite(x$threshold_control)
  age_difference <- x$threshold_case[valid_ages] -
    x$threshold_control[valid_ages]
  summary_rows[[i]] <- data.frame(
    true_hr = x$true_hr[1L],
    cohort = x$cohort[1L],
    outcome = x$outcome[1L],
    replicates = nrow(x),
    successful_fits = sum(is.finite(x$p)),
    power = mean(x$significant_correct),
    mean_events = if (all(!is.finite(x$events))) NA_real_ else
      mean(x$events, na.rm = TRUE),
    mean_matched_sets = mean(x$matched_sets),
    threshold_reached_both = mean(valid_ages),
    mean_threshold_age_case = if (all(!is.finite(x$threshold_case)))
      NA_real_ else mean(x$threshold_case, na.rm = TRUE),
    mean_threshold_age_control = if (all(!is.finite(x$threshold_control)))
      NA_real_ else mean(x$threshold_control, na.rm = TRUE),
    mean_case_minus_control_age = if (!length(age_difference))
      NA_real_ else mean(age_difference))
}
summary <- do.call(rbind, summary_rows)
rownames(summary) <- NULL
summary <- summary[order(summary$true_hr,summary$cohort,summary$outcome),]
print(summary,row.names=FALSE)
write.csv(replicates,file.path(OUT_DIR,'CAR_Aim1_matched_replicates.csv'),
          row.names=FALSE)
write.csv(summary,file.path(OUT_DIR,'CAR_Aim1_matched_power_summary.csv'),
          row.names=FALSE)
if (length(curves)) write.csv(do.call(rbind,curves),
  file.path(OUT_DIR,'CAR_Aim1_matched_age_curves.csv'),row.names=FALSE)
cat('Power tests the case-control transition coefficient. The 10% and 5%\n',
    'ages are conditional, descriptive estimates, not separately powered tests.\n',
    'Death governs index eligibility; it cannot occur among selected participants\n',
    'during their retained pre-index histories. Vary assumed CMS coverage,\n',
    'diagnosis timing, and baseline restriction before reporting grant power.\n')

})

if ('2' %in% CAR_POWER_AIMS) run_power_aim('2', {
# CAR Study Aim 2 three-state power simulation
#
# Primary analysis represented here:
#   unrestricted -> restricted
#   unrestricted -> ceased
#   restricted   -> ceased
#
# The fitted models use a complementary log-log link, an interval-length
# offset, a baseline predictor, a one-wave-lagged within-person
# deviation, incident ADRD status, and the lagged deviation x ADRD interaction.
# Cessation is absorbing. Return from restricted to unrestricted is not included
# because it is a secondary transition and Table 2 reports first restrictions.

required_packages <- c(
  "data.table", "sandwich", "metafor", "future", "future.apply"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(data.table)

set.seed(20260924)

# -----------------------------------------------------------------------------
# USER SETTINGS
# -----------------------------------------------------------------------------

# Use 25-50 replicates to test the script. Use at least 500 for grant results.
SIM_REPS <- as.integer(Sys.getenv('CAR_AIM2_REPS', '500'))

# Set TRUE for ADRD prevalence, coverage, and heterogeneity sensitivity ranges.
# The full grid is computationally intensive.
RUN_FULL_GRID <- identical(tolower(Sys.getenv('CAR_AIM2_FULL_GRID', 'false')), 'true')

# Number of parallel R sessions. Reduce this if memory is limited.
N_WORKERS <- max(1L, min(6L, future::availableCores() - 1L))

# Planning effects are hazard ratios per 1-SD increase in the risk-oriented
# lagged within-person predictor on the interval hazard. A hazard ratio is not
# tied to a one-year contrast; the baseline hazards are expressed per year so
# unequal interval lengths can be handled with offset(log(dt)). Protective
# measures should be reverse-oriented so HR > 1 represents lower protection or
# worsening status.
PLANNING_HR <- c(1.15, 1.20, 1.25)

# Run all Aim 2 domains, or replace this vector with one or more domain names
# while testing the code (for example, "Cognitive function").
DOMAINS_TO_RUN <- c(
  "Cognitive function", "Physical function", "Comorbidity burden",
  "Hospitalizations", "Medication burden", "Functional independence",
  "Social support", "Physical activity"
)

# The proportion of cessation events that occur directly from unrestricted
# driving is not available in Table 2. The primary value is 0.50. Re-run with
# 0.25 and 0.75 as sensitivity analyses.
DIRECT_CESSATION_FRACTION <- 0.50

# AR(1) correlation for the time-varying predictor.
WITHIN_PERSON_RHO <- 0.60

# Baseline effects used to generate realistic event dependence. They are
# adjustment terms rather than targets of the power calculation.
BETWEEN_PERSON_HR <- 1.10
ADRD_MAIN_HR <- 1.15

# -----------------------------------------------------------------------------
# TABLE 2 DESIGN INPUTS
# -----------------------------------------------------------------------------

design <- data.table(
  cohort = c("HRS", "NHATS", "ACTIVE"),
  n = c(19099L, 15314L, 655L),
  person_observations = c(63708L, 70424L, 2572L),
  persons_with_2plus = c(14154L, 12361L, 539L),
  transition_pairs = c(44608L, 55110L, 1917L),
  restriction_events = c(4400L, 3663L, 174L),
  cessation_events = c(2472L, 2526L, 58L),
  restriction_age_slope = c(0.068, 0.033, 0.061),
  cessation_age_slope = c(0.121, 0.086, 0.127),
  design_effect = c(1.5, 1.7, 1.0),
  baseline_age_mean = c(70, 77, 74),
  baseline_age_sd = c(5.5, 7.0, 6.0),
  min_age = c(65, 65, 65),
  max_age = c(95, 105, 94)
)

# Visit schedules reproduce the cohort assessment cadence. HRS has one more
# observation pair than Table 2 when pairs are calculated as observations minus
# persons; this one-pair discrepancy is retained and reported below.
visit_schedules <- list(
  HRS = c(0, 2, 4, 6, 8, 10, 12, 16, 20),
  NHATS = 0:13,
  ACTIVE = c(0, 1, 2, 3, 5, 10)
)

# Expected availability of each predictor among otherwise eligible intervals.
# These are planning assumptions, not observed linked-data retention rates.
expected_coverage <- data.table(
  domain = rep(
    c(
      "Cognitive function", "Physical function", "Comorbidity burden",
      "Hospitalizations", "Medication burden", "Functional independence",
      "Social support", "Physical activity"
    ),
    each = 3L
  ),
  cohort = rep(c("HRS", "NHATS", "ACTIVE"), times = 8L),
  expected = c(
    0.85, 0.90, 0.90,  # cognitive function
    0.40, 0.75, 0.75,  # physical function
    0.70, 0.80, 0.60,  # comorbidity burden
    0.70, 0.80, 0.60,  # hospitalizations
    0.50, 0.60, 0.50,  # medication burden / Part D
    0.80, 0.90, 0.90,  # functional independence
    0.55, 0.80, 0.90,  # social support
    0.75, 0.85, 0.00   # physical activity; no repeated ACTIVE measure
  )
)

coverage_multipliers <- data.table(
  coverage_scenario = c("conservative", "expected", "favorable"),
  multiplier = c(0.75, 1.00, 1.15)
)

# -----------------------------------------------------------------------------
# BUILD COHORT-SPECIFIC PERSON-INTERVAL SKELETONS
# -----------------------------------------------------------------------------

allocate_visit_counts <- function(
    sample_n, total_observations, n_with_2plus, max_visits, allocation_seed) {
  set.seed(allocation_seed)
  n_visits <- rep.int(1L, sample_n)
  repeated_ids <- sample.int(sample_n, n_with_2plus, replace = FALSE)
  n_visits[repeated_ids] <- 2L

  remaining <- total_observations - sum(n_visits)
  if (remaining < 0L) stop("Inconsistent observation totals.")

  while (remaining > 0L) {
    eligible <- repeated_ids[n_visits[repeated_ids] < max_visits]
    if (length(eligible) == 0L) stop("Visit schedule cannot accommodate totals.")
    selected <- sample(eligible, min(length(eligible), remaining), replace = FALSE)
    n_visits[selected] <- n_visits[selected] + 1L
    remaining <- remaining - length(selected)
  }
  n_visits
}

make_interval_skeleton <- function(design_row, schedule, skeleton_seed) {
  cohort_name <- design_row$cohort[[1L]]
  n_visits <- allocate_visit_counts(
    sample_n = design_row$n[[1L]],
    total_observations = design_row$person_observations[[1L]],
    n_with_2plus = design_row$persons_with_2plus[[1L]],
    max_visits = length(schedule),
    allocation_seed = skeleton_seed
  )

  set.seed(skeleton_seed + 10000L)
  baseline_age <- rnorm(
    design_row$n[[1L]],
    design_row$baseline_age_mean[[1L]],
    design_row$baseline_age_sd[[1L]]
  )
  baseline_age <- pmin(
    design_row$max_age[[1L]],
    pmax(design_row$min_age[[1L]], baseline_age)
  )

  visit_data <- rbindlist(lapply(seq_len(design_row$n[[1L]]), function(person_id) {
    k <- n_visits[[person_id]]
    data.table(
      cohort = cohort_name,
      id = person_id,
      visit = seq_len(k),
      study_time = schedule[seq_len(k)],
      baseline_age = baseline_age[[person_id]]
    )
  }))

  setorder(visit_data, id, visit)
  visit_data[, `:=`(
    next_time = shift(study_time, type = "lead"),
    next_visit = shift(visit, type = "lead")
  ), by = id]

  intervals <- visit_data[!is.na(next_time)]
  intervals[, `:=`(
    step = visit,
    dt = next_time - study_time,
    age_start = baseline_age + study_time,
    age_mid = baseline_age + (study_time + next_time) / 2
  )]
  intervals <- intervals[dt > 0]
  intervals[, c("next_time", "next_visit") := NULL]
  setorder(intervals, step, id)
  intervals[]
}

skeletons <- list()
for (cohort_index in seq_len(nrow(design))) {
  cohort_name <- design$cohort[[cohort_index]]
  skeletons[[cohort_name]] <- make_interval_skeleton(
    design_row = design[cohort_index],
    schedule = visit_schedules[[cohort_name]],
    skeleton_seed = 2000L + cohort_index
  )
}

table2_check <- rbindlist(lapply(names(skeletons), function(cohort_name) {
  design_row <- design[cohort == cohort_name]
  data.table(
    cohort = cohort_name,
    table_pairs = design_row$transition_pairs,
    simulated_pairs = nrow(skeletons[[cohort_name]]),
    difference = nrow(skeletons[[cohort_name]]) - design_row$transition_pairs
  )
}))

cat("Table 2 transition-pair check:\n")
print(table2_check)

# -----------------------------------------------------------------------------
# CALIBRATE BASELINE STATE-TRANSITION HAZARDS TO TABLE 2
# -----------------------------------------------------------------------------

competing_probabilities <- function(lambda_1, lambda_2, interval_length) {
  total_lambda <- lambda_1 + lambda_2
  any_event <- 1 - exp(-total_lambda * interval_length)
  p1 <- ifelse(total_lambda > 0, any_event * lambda_1 / total_lambda, 0)
  p2 <- ifelse(total_lambda > 0, any_event * lambda_2 / total_lambda, 0)
  list(p1 = p1, p2 = p2)
}

expected_transition_counts <- function(log_rates, interval_data, design_row) {
  sample_n <- design_row$n[[1L]]
  probability_u <- rep.int(1, sample_n)
  probability_r <- rep.int(0, sample_n)
  expected <- c(UR = 0, UC = 0, RC = 0)

  for (step_value in sort(unique(interval_data$step))) {
    rows <- interval_data[step == step_value]
    person_ids <- rows$id

    lambda_ur <- exp(
      log_rates[[1L]] +
        design_row$restriction_age_slope[[1L]] * (rows$age_mid - 75)
    )
    lambda_uc <- exp(
      log_rates[[2L]] +
        design_row$cessation_age_slope[[1L]] * (rows$age_mid - 75)
    )
    lambda_rc <- exp(
      log_rates[[3L]] +
        design_row$cessation_age_slope[[1L]] * (rows$age_mid - 75)
    )

    from_u <- competing_probabilities(lambda_ur, lambda_uc, rows$dt)
    p_rc <- 1 - exp(-lambda_rc * rows$dt)

    current_u <- probability_u[person_ids]
    current_r <- probability_r[person_ids]

    expected[["UR"]] <- expected[["UR"]] + sum(current_u * from_u$p1)
    expected[["UC"]] <- expected[["UC"]] + sum(current_u * from_u$p2)
    expected[["RC"]] <- expected[["RC"]] + sum(current_r * p_rc)

    probability_u[person_ids] <- current_u * (1 - from_u$p1 - from_u$p2)
    probability_r[person_ids] <- current_r * (1 - p_rc) + current_u * from_u$p1
  }
  expected
}

calibrate_cohort_rates <- function(
    interval_data, design_row, direct_cessation_fraction) {
  target <- c(
    UR = design_row$restriction_events[[1L]],
    UC = design_row$cessation_events[[1L]] * direct_cessation_fraction,
    RC = design_row$cessation_events[[1L]] * (1 - direct_cessation_fraction)
  )

  person_years <- sum(interval_data$dt)
  starting_rates <- c(
    max(target[["UR"]] / person_years, 1e-5),
    max(target[["UC"]] / person_years, 1e-5),
    max(3 * target[["RC"]] / person_years, 1e-5)
  )

  objective <- function(log_rates) {
    expected <- expected_transition_counts(log_rates, interval_data, design_row)
    sum(log((expected + 1) / (target + 1))^2)
  }

  fit <- optim(
    par = log(starting_rates),
    fn = objective,
    method = "Nelder-Mead",
    control = list(maxit = 5000, reltol = 1e-10)
  )

  if (fit$convergence != 0L) {
    warning("Calibration did not fully converge for ", design_row$cohort[[1L]])
  }

  expected <- expected_transition_counts(fit$par, interval_data, design_row)
  data.table(
    cohort = design_row$cohort[[1L]],
    transition = names(target),
    annual_hazard_at_75 = exp(fit$par),
    target_events = as.numeric(target),
    expected_events = as.numeric(expected)
  )
}

calibration <- rbindlist(lapply(names(skeletons), function(cohort_name) {
  calibrate_cohort_rates(
    interval_data = skeletons[[cohort_name]],
    design_row = design[cohort == cohort_name],
    direct_cessation_fraction = DIRECT_CESSATION_FRACTION
  )
}))

cat("\nCalibrated state-transition hazards:\n")
print(calibration)
write.csv(calibration, "CAR_Aim2_transition_calibration.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# GENERATE TIME-VARYING PREDICTORS AND INCIDENT ADRD STATUS
# -----------------------------------------------------------------------------

add_predictors <- function(
    interval_data, design_row, incident_adrd_fraction, ar_correlation,
    generation_seed) {
  set.seed(generation_seed)
  output <- copy(interval_data)
  sample_n <- design_row$n[[1L]]

  person_age <- output[, .(baseline_age = first(baseline_age)), by = id]
  n_cases <- max(1L, round(sample_n * incident_adrd_fraction))
  sampling_weight <- exp(0.06 * (person_age$baseline_age - 75))
  case_ids <- sample(
    person_age$id,
    size = n_cases,
    replace = FALSE,
    prob = sampling_weight
  )

  adrd <- as.integer(seq_len(sample_n) %in% case_ids)
  between <- rnorm(sample_n) + 0.25 * adrd
  between <- as.numeric(scale(between))

  current_w <- rnorm(sample_n)
  baseline_w <- current_w
  output[, W_lag := NA_real_]
  innovation_sd <- sqrt(1 - ar_correlation^2)

  for (step_value in sort(unique(output$step))) {
    row_index <- which(output$step == step_value)
    person_ids <- output$id[row_index]
    output$W_lag[row_index] <- current_w[person_ids]
    current_w[person_ids] <-
      ar_correlation * current_w[person_ids] +
      innovation_sd * rnorm(length(person_ids))
  }

  output[, `:=`(
    adrd = adrd[id],
    X0 = (between + baseline_w)[id]
  )]

  # Define change from the participant-specific baseline before each interval, and
  # scale it so the planning HR is interpreted per 1-SD worsening.
  output[, W_lag := W_lag - baseline_w[id]]
  within_sd <- sd(output$W_lag)
  if (!is.finite(within_sd) || within_sd <= 0) {
    stop("Unable to standardize the within-person predictor.")
  }
  output[, W_lag := W_lag / within_sd]
  output
}

# -----------------------------------------------------------------------------
# SIMULATE THREE-STATE DRIVING HISTORIES
# -----------------------------------------------------------------------------

simulate_three_state_history <- function(
    interval_data, design_row, calibrated_rates, main_log_hr,
    interaction_log_hr, simulation_seed) {
  set.seed(simulation_seed)
  output <- copy(interval_data)
  sample_n <- design_row$n[[1L]]
  current_state <- rep.int(0L, sample_n) # 0 = unrestricted, 1 = restricted, 2 = ceased

  base_rate <- setNames(
    calibrated_rates$annual_hazard_at_75,
    calibrated_rates$transition
  )

  output[, `:=`(from_state = NA_integer_, to_state = NA_integer_)]

  for (step_value in sort(unique(output$step))) {
    row_index <- which(output$step == step_value)
    rows <- output[row_index]
    person_ids <- rows$id
    starting_state <- current_state[person_ids]
    ending_state <- starting_state

    linear_effect <-
      log(BETWEEN_PERSON_HR) * rows$X0 +
      log(ADRD_MAIN_HR) * rows$adrd +
      main_log_hr * rows$W_lag +
      interaction_log_hr * rows$W_lag * rows$adrd

    unrestricted_rows <- which(starting_state == 0L)
    if (length(unrestricted_rows) > 0L) {
      lambda_ur <- base_rate[["UR"]] * exp(
        design_row$restriction_age_slope[[1L]] *
          (rows$age_mid[unrestricted_rows] - 75) +
          linear_effect[unrestricted_rows]
      )
      lambda_uc <- base_rate[["UC"]] * exp(
        design_row$cessation_age_slope[[1L]] *
          (rows$age_mid[unrestricted_rows] - 75) +
          linear_effect[unrestricted_rows]
      )
      probabilities <- competing_probabilities(
        lambda_ur, lambda_uc, rows$dt[unrestricted_rows]
      )
      draw <- runif(length(unrestricted_rows))
      restricted_event <- draw < probabilities$p1
      cessation_event <-
        draw >= probabilities$p1 &
        draw < (probabilities$p1 + probabilities$p2)
      ending_state[unrestricted_rows[restricted_event]] <- 1L
      ending_state[unrestricted_rows[cessation_event]] <- 2L
    }

    restricted_rows <- which(starting_state == 1L)
    if (length(restricted_rows) > 0L) {
      lambda_rc <- base_rate[["RC"]] * exp(
        design_row$cessation_age_slope[[1L]] *
          (rows$age_mid[restricted_rows] - 75) +
          linear_effect[restricted_rows]
      )
      cessation_probability <-
        1 - exp(-lambda_rc * rows$dt[restricted_rows])
      ending_state[restricted_rows[
        runif(length(restricted_rows)) < cessation_probability
      ]] <- 2L
    }

    output$from_state[row_index] <- starting_state
    output$to_state[row_index] <- ending_state
    current_state[person_ids] <- ending_state
  }
  output
}

# -----------------------------------------------------------------------------
# FIT THE PLANNED COMPLEMENTARY LOG-LOG TRANSITION MODELS
# -----------------------------------------------------------------------------

fit_one_transition <- function(
    simulated_data, transition_code, coverage_probability, design_effect_value,
    analysis_seed) {
  set.seed(analysis_seed)

  if (coverage_probability <= 0) return(NULL)

  if (transition_code == "UR") {
    analysis_data <- simulated_data[from_state == 0L]
    analysis_data[, event := as.integer(to_state == 1L)]
  } else if (transition_code == "UC") {
    analysis_data <- simulated_data[from_state == 0L]
    analysis_data[, event := as.integer(to_state == 2L)]
  } else if (transition_code == "RC") {
    analysis_data <- simulated_data[from_state == 1L]
    analysis_data[, event := as.integer(to_state == 2L)]
  } else {
    stop("Unknown transition code: ", transition_code)
  }

  if (nrow(analysis_data) == 0L) return(NULL)
  analysis_data <- analysis_data[runif(.N) < coverage_probability]

  n_events <- sum(analysis_data$event)
  if (
    nrow(analysis_data) < 100L || n_events < 20L ||
    (nrow(analysis_data) - n_events) < 20L ||
    uniqueN(analysis_data$id) < 50L
  ) return(NULL)

  model_fit <- tryCatch(
    suppressWarnings(
      glm(
        event ~ splines::ns(age_mid, df = 3) + X0 + W_lag + adrd +
          W_lag:adrd + offset(log(dt)),
        family = binomial(link = "cloglog"),
        data = analysis_data
      )
    ),
    error = function(e) NULL
  )

  if (is.null(model_fit) || !isTRUE(model_fit$converged)) return(NULL)

  robust_vcov <- tryCatch(
    sandwich::vcovCL(
      model_fit,
      cluster = analysis_data$id,
      type = "HC0"
    ) * design_effect_value,
    error = function(e) NULL
  )
  if (is.null(robust_vcov)) return(NULL)

  coefficient_names <- names(coef(model_fit))
  interaction_name <- coefficient_names[
    coefficient_names %in% c("W_lag:adrd", "adrd:W_lag")
  ]
  if (length(interaction_name) != 1L || !("W_lag" %in% coefficient_names)) {
    return(NULL)
  }

  extract_term <- function(term_name, reporting_name) {
    estimate <- unname(coef(model_fit)[[term_name]])
    standard_error <- sqrt(robust_vcov[term_name, term_name])
    data.table(
      term = reporting_name,
      estimate = estimate,
      standard_error = standard_error,
      fitted_hr = exp(estimate),
      p_value = 2 * pnorm(-abs(estimate / standard_error)),
      events = n_events,
      intervals = nrow(analysis_data),
      persons = uniqueN(analysis_data$id)
    )
  }

  rbind(
    extract_term("W_lag", "within_person"),
    extract_term(interaction_name, "within_person_x_adrd")
  )
}

meta_analyze <- function(cohort_results, included_cohorts) {
  usable <- cohort_results[
    cohort %in% included_cohorts &
      is.finite(estimate) & is.finite(standard_error) & standard_error > 0
  ]
  if (nrow(usable) < 2L) return(NULL)

  fixed_fit <- tryCatch(
    metafor::rma.uni(
      yi = usable$estimate,
      sei = usable$standard_error,
      method = "FE",
      test = "z"
    ),
    error = function(e) NULL
  )
  random_fit <- tryCatch(
    metafor::rma.uni(
      yi = usable$estimate,
      sei = usable$standard_error,
      method = "REML",
      test = "z"
    ),
    error = function(e) NULL
  )
  if (is.null(fixed_fit) || is.null(random_fit)) return(NULL)

  get_p <- function(fit_object) {
    if (!is.null(fit_object$pval)) as.numeric(fit_object$pval) else NA_real_
  }

  hrs_beta <- usable[cohort == "HRS", estimate]
  nhats_beta <- usable[cohort == "NHATS", estimate]
  same_direction <-
    length(hrs_beta) == 1L && length(nhats_beta) == 1L &&
    hrs_beta > 0 && nhats_beta > 0

  data.table(
    completed_cohorts = nrow(usable),
    pooled_hr_fixed = exp(as.numeric(fixed_fit$b)),
    p_fixed = get_p(fixed_fit),
    pooled_hr_random = exp(as.numeric(random_fit$b)),
    p_random = get_p(random_fit),
    tau2 = as.numeric(random_fit$tau2),
    I2 = as.numeric(random_fit$I2),
    same_direction_hrs_nhats = same_direction,
    total_events = sum(usable$events)
  )
}

# -----------------------------------------------------------------------------
# SCENARIOS
# -----------------------------------------------------------------------------

if (RUN_FULL_GRID) {
  adrd_values <- c(0.10, 0.20, 0.25)
  coverage_values <- c("conservative", "expected", "favorable")
  tau_values <- c(0, 0.025, 0.05)
} else {
  adrd_values <- 0.20
  coverage_values <- "expected"
  tau_values <- 0
}

scenario_grid <- CJ(
  effect_type = c("within_person", "within_person_x_adrd"),
  domain = DOMAINS_TO_RUN,
  coverage_scenario = coverage_values,
  adrd_fraction = adrd_values,
  tau = tau_values,
  true_hr = PLANNING_HR,
  sorted = FALSE
)
scenario_grid[, scenario_id := .I]

coverage_lookup <- expected_coverage[
  , .(
    coverage_scenario = coverage_multipliers$coverage_scenario,
    multiplier = coverage_multipliers$multiplier
  ),
  by = .(domain, cohort, expected)
]
coverage_lookup[, coverage_probability := pmin(0.98, expected * multiplier)]

run_one_replicate <- function(scenario_row, replicate_id) {
  scenario_number <- scenario_row$scenario_id[[1L]]
  replicate_seed <- 1000000L + scenario_number * 10000L + replicate_id

  set.seed(replicate_seed)
  cohort_deviation <- rnorm(nrow(design), mean = 0, sd = scenario_row$tau[[1L]])
  cohort_deviation <- cohort_deviation - mean(cohort_deviation)
  names(cohort_deviation) <- design$cohort

  cohort_fits <- list()
  fit_counter <- 1L

  for (cohort_index in seq_len(nrow(design))) {
    design_row <- design[cohort_index]
    cohort_name <- design_row$cohort[[1L]]

    if (scenario_row$effect_type[[1L]] == "within_person") {
      main_log_hr <- log(scenario_row$true_hr[[1L]]) +
        cohort_deviation[[cohort_name]]
      interaction_log_hr <- 0
      target_term <- "within_person"
    } else {
      main_log_hr <- 0
      interaction_log_hr <- log(scenario_row$true_hr[[1L]]) +
        cohort_deviation[[cohort_name]]
      target_term <- "within_person_x_adrd"
    }

    generated_data <- add_predictors(
      interval_data = skeletons[[cohort_name]],
      design_row = design_row,
      incident_adrd_fraction = scenario_row$adrd_fraction[[1L]],
      ar_correlation = WITHIN_PERSON_RHO,
      generation_seed = replicate_seed + 100L * cohort_index
    )

    simulated_data <- simulate_three_state_history(
      interval_data = generated_data,
      design_row = design_row,
      calibrated_rates = calibration[cohort == cohort_name],
      main_log_hr = main_log_hr,
      interaction_log_hr = interaction_log_hr,
      simulation_seed = replicate_seed + 1000L * cohort_index
    )

    coverage_probability <- coverage_lookup[
      domain == scenario_row$domain[[1L]] &
        cohort == cohort_name &
        coverage_scenario == scenario_row$coverage_scenario[[1L]],
      coverage_probability
    ]

    for (transition_code in c("UR", "UC", "RC")) {
      fitted <- fit_one_transition(
        simulated_data = simulated_data,
        transition_code = transition_code,
        coverage_probability = coverage_probability,
        design_effect_value = design_row$design_effect[[1L]],
        analysis_seed = replicate_seed + 10000L * cohort_index +
          match(transition_code, c("UR", "UC", "RC"))
      )
      if (!is.null(fitted)) {
        fitted[, `:=`(
          cohort = cohort_name,
          transition = transition_code,
          target_term = target_term
        )]
        cohort_fits[[fit_counter]] <- fitted
        fit_counter <- fit_counter + 1L
      }
    }
  }

  if (length(cohort_fits) == 0L) return(NULL)
  all_fits <- rbindlist(cohort_fits, fill = TRUE)
  target_fits <- all_fits[term == target_term]

  meta_rows <- list()
  meta_counter <- 1L
  for (transition_code in c("UR", "UC", "RC")) {
    transition_fits <- target_fits[transition == transition_code]
    for (analysis_name in c("HRS + NHATS", "HRS + NHATS + ACTIVE")) {
      included <- if (analysis_name == "HRS + NHATS") {
        c("HRS", "NHATS")
      } else {
        c("HRS", "NHATS", "ACTIVE")
      }
      if (
        analysis_name == "HRS + NHATS + ACTIVE" &&
        !("ACTIVE" %in% transition_fits$cohort)
      ) next
      pooled <- meta_analyze(transition_fits, included)
      if (!is.null(pooled)) {
        pooled[, `:=`(
          replicate = replicate_id,
          transition = transition_code,
          analysis = analysis_name
        )]
        meta_rows[[meta_counter]] <- pooled
        meta_counter <- meta_counter + 1L
      }
    }
  }

  list(
    cohort = target_fits[, .(
      scenario_id = scenario_number,
      replicate = replicate_id,
      cohort,
      transition,
      estimate,
      standard_error,
      fitted_hr,
      p_value,
      events,
      intervals,
      persons
    )],
    meta = if (length(meta_rows) > 0L) {
      meta_result <- rbindlist(meta_rows, fill = TRUE)
      meta_result[, scenario_id := scenario_number]
      meta_result[]
    } else NULL
  )
}

run_scenario <- function(scenario_number) {
  scenario_row <- scenario_grid[scenario_id == scenario_number]
  cat(
    "Running scenario", scenario_number, "of", nrow(scenario_grid), ":",
    scenario_row$effect_type, "|", scenario_row$domain, "| HR =",
    scenario_row$true_hr, "\n"
  )

  replicate_results <- future.apply::future_lapply(
    X = seq_len(SIM_REPS),
    FUN = function(replicate_id) {
      run_one_replicate(scenario_row, replicate_id)
    },
    future.seed = TRUE
  )

  replicate_results <- Filter(Negate(is.null), replicate_results)
  if (length(replicate_results) == 0L) {
    return(list(cohort = NULL, meta = NULL))
  }

  list(
    cohort = rbindlist(lapply(replicate_results, `[[`, "cohort"), fill = TRUE),
    meta = rbindlist(lapply(replicate_results, `[[`, "meta"), fill = TRUE)
  )
}

# -----------------------------------------------------------------------------
# RUN SIMULATIONS
# -----------------------------------------------------------------------------

future::plan(future::multisession, workers = N_WORKERS)

all_results <- vector("list", nrow(scenario_grid))
for (scenario_number in scenario_grid$scenario_id) {
  all_results[[scenario_number]] <- run_scenario(scenario_number)
}

future::plan(future::sequential)

cohort_replicates <- rbindlist(
  lapply(all_results, `[[`, "cohort"),
  fill = TRUE
)
meta_replicates <- rbindlist(
  lapply(all_results, `[[`, "meta"),
  fill = TRUE
)

cohort_replicates <- merge(
  cohort_replicates,
  scenario_grid,
  by = "scenario_id",
  all.x = TRUE
)
meta_replicates <- merge(
  meta_replicates,
  scenario_grid,
  by = "scenario_id",
  all.x = TRUE
)

# Power requires statistical significance in the hypothesized direction.
cohort_power_summary <- cohort_replicates[, .(
  completed_replicates = uniqueN(replicate),
  mean_fitted_hr = mean(fitted_hr, na.rm = TRUE),
  mean_events = mean(events, na.rm = TRUE),
  power = mean(p_value < 0.05 & estimate > 0, na.rm = TRUE)
), by = .(
  effect_type, domain, coverage_scenario, adrd_fraction, tau, true_hr,
  transition, cohort
)]

meta_power_summary <- meta_replicates[, .(
  completed_replicates = uniqueN(replicate),
  mean_pooled_hr_fixed = mean(pooled_hr_fixed, na.rm = TRUE),
  fixed_effect_power = mean(p_fixed < 0.05 & pooled_hr_fixed > 1, na.rm = TRUE),
  mean_pooled_hr_random = mean(pooled_hr_random, na.rm = TRUE),
  random_effects_power = mean(
    p_random < 0.05 & pooled_hr_random > 1,
    na.rm = TRUE
  ),
  same_direction_probability = mean(same_direction_hrs_nhats, na.rm = TRUE),
  median_tau2 = median(tau2, na.rm = TRUE),
  median_I2 = median(I2, na.rm = TRUE),
  mean_total_events = mean(total_events, na.rm = TRUE)
), by = .(
  effect_type, domain, coverage_scenario, adrd_fraction, tau, true_hr,
  transition, analysis
)]

# Minimum planning HR that reaches 80% random-effects power. na.rm=TRUE prevents
# the missing-value error that can otherwise occur when a cell has no estimates.
minimum_detectable_hr <- meta_power_summary[, {
  valid <- is.finite(random_effects_power) & is.finite(true_hr)
  reaches_80 <- valid & random_effects_power >= 0.80
  list(
    minimum_hr_for_80_percent_power = if (any(reaches_80)) {
      min(true_hr[reaches_80])
    } else {
      NA_real_
    }
  )
}, by = .(
  effect_type, domain, coverage_scenario, adrd_fraction, tau,
  transition, analysis
)]

setorder(
  meta_power_summary,
  effect_type, domain, coverage_scenario, adrd_fraction, tau,
  transition, analysis, true_hr
)

cat("\nMeta-analytic power summary:\n")
print(meta_power_summary)

cat("\nMinimum HR reaching 80% random-effects power:\n")
print(minimum_detectable_hr)

write.csv(
  cohort_replicates,
  "CAR_Aim2_three_state_cohort_replicates.csv",
  row.names = FALSE
)
write.csv(
  meta_replicates,
  "CAR_Aim2_three_state_meta_replicates.csv",
  row.names = FALSE
)
write.csv(
  cohort_power_summary,
  "CAR_Aim2_three_state_cohort_power_summary.csv",
  row.names = FALSE
)
write.csv(
  meta_power_summary,
  "CAR_Aim2_three_state_meta_power_summary.csv",
  row.names = FALSE
)
write.csv(
  minimum_detectable_hr,
  "CAR_Aim2_three_state_minimum_detectable_hr.csv",
  row.names = FALSE
)

cat("\nFiles written to: ", normalizePath(getwd()), "\n", sep = "")
cat(
  "For the grant, report random-effects power for HRS + NHATS as the primary\n",
  "coordinated estimate and HRS + NHATS + ACTIVE as a corroborative analysis.\n",
  sep = ""
)

})

if ('3' %in% CAR_POWER_AIMS) run_power_aim('3', {
# CAR Study Aim 3 planning power: age-75 profiles, subsequent ADRD, and HRS-HCAP.
# This Aim 3 section is embedded in the combined script.
# Requires the survival package. Set CAR_AIM3_REPS=10 for an initial check;
# set CAR_AIM3_REPS=500 for grant planning estimates after checking outputs.
# CMS linkage, post-landmark events, profile prevalence, class uncertainty,
# biomarker overlap, and effect sizes are assumptions, not observed data.
# This simulates inference conditional on a recoverable three-class solution.
# It does not test Mplus class enumeration, stability, or cross-cohort transport.

if (!requireNamespace('survival', quietly = TRUE)) {
  stop('Please install the survival package before running this script.')
}
set.seed(20260925)
REPS <- as.integer(Sys.getenv('CAR_AIM3_REPS', '100'))
PSEUDO_DRAWS <- as.integer(Sys.getenv('CAR_AIM3_DRAWS', '20'))
stopifnot(is.finite(REPS), REPS > 0L,
          is.finite(PSEUDO_DRAWS), PSEUDO_DRAWS >= 2L)
FULL_GRID <- identical(tolower(Sys.getenv('CAR_AIM3_FULL_GRID', 'false')), 'true')
RUN_OTHER_COHORTS <- identical(tolower(Sys.getenv('CAR_AIM3_OTHER_COHORTS', 'false')), 'true')
ALPHA <- .05
MAX_FOLLOWUP_YEARS <- 10
DEATH_ANNUAL_HAZARD <- .035
MODERATE_CLASS_HR <- 1.10
MODERATE_BIOMARKER_FRACTION <- .5
BIOMARKER_ERROR_CORRELATION <- .35
OUT_DIR <- 'CAR_Aim3_power_results'
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# HRS N=11,500 is the existing age-75 planning estimate, not linked eligibility.
# NHATS/ACTIVE Ns below are placeholders and run only when explicitly enabled.
cohorts <- data.frame(cohort = c('HRS', 'NHATS', 'ACTIVE'),
  n = c(11500L, 8000L, 350L), design_effect = c(1.5, 1.7, 1.0))
if (!RUN_OTHER_COHORTS) cohorts <- cohorts[cohorts$cohort == 'HRS', ]
if (FULL_GRID) {
  EVENT_FRACTIONS <- c(.10, .20, .25)
  PROFILE_PREVALENCES <- c(.05, .10, .20)
  ENTROPIES <- c(.70, .80, .90)
  ADRD_HRS <- c(1.15, 1.25, 1.40)
  HCAP_NS <- c(1000L, 2000L, 2500L, 3500L)
  BIOMARKER_DS <- c(.20, .30, .35, .40)
} else {
  EVENT_FRACTIONS <- .20
  PROFILE_PREVALENCES <- c(.10, .20)
  ENTROPIES <- .80
  ADRD_HRS <- c(1.25, 1.40)
  HCAP_NS <- c(2000L, 2500L)
  BIOMARKER_DS <- c(.30, .35)
}
# The third profile represents a moderate trajectory. The high-risk profile
# has the specified prevalence, and class 0 is the reference.
MODERATE_CLASS_PREVALENCE <- .20

class_probabilities <- function(high_prevalence) {
  p <- c(1 - MODERATE_CLASS_PREVALENCE - high_prevalence,
         MODERATE_CLASS_PREVALENCE, high_prevalence)
  if (any(p <= 0) || abs(sum(p) - 1) > 1e-10) stop('Invalid class prevalences.')
  p
}

entropy_value <- function(q, prior) {
  # Measurement label L is correct with probability q; otherwise it is drawn
  # from the class prevalences. Apply Bayes' rule to get P(class | L).
  post <- sapply(seq_along(prior), function(label) {
    q * as.numeric(seq_along(prior) == label) + (1 - q) * prior
  })
  post_entropy <- -colSums(post * log(pmax(post, 1e-12)))
  marginal_entropy <- -sum(prior * log(prior))
  1 - sum(prior * post_entropy) / marginal_entropy
}

label_accuracy <- function(prior, target_entropy) {
  if (target_entropy <= 0 || target_entropy >= 1)
    stop('Use a target entropy strictly between 0 and 1.')
  uniroot(function(q) entropy_value(q, prior) - target_entropy,
          interval = c(0, 1), tol = 1e-8)$root
}

one_draw <- function(probability_matrix) {
  u <- runif(nrow(probability_matrix))
  1L + rowSums(u > t(apply(probability_matrix, 1L, cumsum)))
}

simulate_class_posteriors <- function(n, prevalence, entropy) {
  prior <- class_probabilities(prevalence)
  true_class <- sample.int(3L, n, replace = TRUE, prob = prior) - 1L
  q <- label_accuracy(prior, entropy)
  # This proxy class label is an imperfect driving-history classification.
  label <- ifelse(runif(n) < q, true_class,
                  sample.int(3L, n, replace = TRUE, prob = prior) - 1L)
  posterior <- matrix(0, nrow = n, ncol = 3L)
  for (k in 0:2) posterior[, k + 1L] <-
    q * as.integer(label == k) + (1 - q) * prior[k + 1L]
  realized_entropy <- 1 - mean(-rowSums(posterior * log(posterior))) /
    (-sum(prior * log(prior)))
  list(class = true_class, posterior = posterior,
       entropy = realized_entropy, prevalence = prior)
}

combine_pseudo_estimates <- function(betas, variances, design_effect) {
  valid <- is.finite(betas) & is.finite(variances) & variances > 0
  betas <- betas[valid]; variances <- variances[valid]
  if (length(betas) < 2L) return(c(beta = NA_real_, se = NA_real_,
                                    p = NA_real_, draws = length(betas)))
  m <- length(betas)
  within <- mean(variances) * design_effect
  between <- var(betas)
  total_variance <- within + (1 + 1 / m) * between
  se <- sqrt(total_variance)
  beta <- mean(betas)
  c(beta = beta, se = se, p = 2 * pnorm(-abs(beta / se)), draws = m)
}

annual_adrd_hazard <- function(target_fraction, class_probs, high_hr) {
  class_hr <- c(1, MODERATE_CLASS_HR, high_hr)
  expected_fraction <- function(rate) {
    cause_rate <- rate * class_hr
    sum(class_probs * cause_rate / (cause_rate + DEATH_ANNUAL_HAZARD) *
      (1 - exp(-(cause_rate + DEATH_ANNUAL_HAZARD) * MAX_FOLLOWUP_YEARS)))
  }
  uniroot(function(rate) expected_fraction(rate) - target_fraction,
          interval = c(1e-8, 5), tol = 1e-9)$root
}

simulate_survival <- function(n, prevalence, entropy, high_hr,
                              target_event_fraction) {
  class_data <- simulate_class_posteriors(n, prevalence, entropy)
  true_class <- class_data$class
  sex <- rbinom(n, 1L, .55)
  race_indicator <- rbinom(n, 1L, .30)
  education <- pmax(4, pmin(20, rnorm(n, 13, 3)))
  comorbidity <- rpois(n, 2)
  base_rate <- annual_adrd_hazard(target_event_fraction,
    class_data$prevalence, high_hr)
  hr <- ifelse(true_class == 2L, high_hr,
               ifelse(true_class == 1L, MODERATE_CLASS_HR, 1))
  # Centered covariate associations allow a covariate-adjusted Cox analysis.
  rate <- base_rate * hr * exp(.10 * (sex - .55) +
    .08 * (race_indicator - .30) - .02 * (education - 13) +
    .08 * (comorbidity - 2))
  time_adrd <- rexp(n, rate = rate)
  time_death <- rexp(n, rate = DEATH_ANNUAL_HAZARD)
  observed_time <- pmin(time_adrd, time_death, MAX_FOLLOWUP_YEARS)
  event_adrd <- as.integer(time_adrd <= time_death &
    time_adrd <= MAX_FOLLOWUP_YEARS)
  event_death <- as.integer(time_death < time_adrd &
    time_death <= MAX_FOLLOWUP_YEARS)
  data <- data.frame(entry_age = 75, exit_age = 75 + observed_time,
    event_adrd = event_adrd, event_death = event_death, sex = sex,
    race_indicator = race_indicator, education = education,
    comorbidity = comorbidity, true_class = true_class)
  list(data = data, posterior = class_data$posterior,
       entropy = class_data$entropy)
}

analyze_survival <- function(simulated, design_effect) {
  d <- simulated$data; prob <- simulated$posterior
  estimates <- variances <- rep(NA_real_, PSEUDO_DRAWS)
  for (m in seq_len(PSEUDO_DRAWS)) {
    membership <- one_draw(prob) - 1L
    d$high <- as.integer(membership == 2L)
    d$moderate <- as.integer(membership == 1L)
    fit <- tryCatch(suppressWarnings(survival::coxph(
      survival::Surv(entry_age, exit_age, event_adrd) ~
        high + moderate + sex + race_indicator + education + comorbidity,
      data = d, ties = 'efron')), error = function(e) NULL)
    if (!is.null(fit) && is.finite(coef(fit)['high'])) {
      estimates[m] <- unname(coef(fit)['high'])
      variances[m] <- tryCatch(vcov(fit)['high', 'high'],
                              error = function(e) NA_real_)
    }
  }
  combined <- combine_pseudo_estimates(estimates, variances, design_effect)
  c(combined, events = sum(d$event_adrd), deaths = sum(d$event_death),
    high_n = sum(d$true_class == 2L),
    high_events = sum(d$event_adrd[d$true_class == 2L]),
    entropy = simulated$entropy)
}

simulate_biomarkers <- function(n, prevalence, entropy, effect_d) {
  class_data <- simulate_class_posteriors(n, prevalence, entropy)
  true_class <- class_data$class
  high <- as.integer(true_class == 2L)
  moderate <- as.integer(true_class == 1L)
  assay_age <- runif(n, 75, 84)
  sex <- rbinom(n, 1L, .55)
  race_indicator <- rbinom(n, 1L, .30)
  education <- pmax(4, pmin(20, rnorm(n, 13, 3)))
  bmi <- pmax(17, pmin(47, rnorm(n, 28, 5)))
  kidney <- pmax(15, pmin(120, rnorm(n, 75, 16)))
  assay_batch <- rbinom(n, 1L, .5)
  diagnosed_at_assay <- rbinom(n, 1L,
    plogis(qlogis(.20) + .45 * high + .15 * moderate))
  shared <- rnorm(n)
  residuals <- matrix(rnorm(n * 3L), nrow = n, ncol = 3L)
  rho <- BIOMARKER_ERROR_CORRELATION
  residuals <- sqrt(rho) * shared + sqrt(1 - rho) * residuals
  covariate_effect <- .12 * (assay_age - 79.5) / 4 +
    .08 * sex + .08 * race_indicator - .02 * (education - 13) +
    .03 * (bmi - 28) + .004 * (75 - kidney) + .12 * assay_batch +
    .20 * diagnosed_at_assay
  effect <- effect_d * (high + MODERATE_BIOMARKER_FRACTION * moderate)
  d <- data.frame(assay_age = assay_age, sex = sex,
    race_indicator = race_indicator, education = education, bmi = bmi,
    kidney = kidney, assay_batch = assay_batch,
    diagnosed_at_assay = diagnosed_at_assay)
  # All outcomes are standardized log-biomarkers; Aβ42/40 has a negative sign.
  d$abeta_log <- covariate_effect - effect + residuals[, 1L]
  d$ptau_log <- covariate_effect + effect + residuals[, 2L]
  d$nfl_log <- covariate_effect + effect + residuals[, 3L]
  list(data = d, posterior = class_data$posterior,
       entropy = class_data$entropy, high_n = sum(high))
}

analyze_biomarkers <- function(simulated, design_effect = 1.5) {
  d <- simulated$data; prob <- simulated$posterior
  biomarkers <- c('abeta_log', 'ptau_log', 'nfl_log')
  estimates <- variances <- matrix(NA_real_, PSEUDO_DRAWS, 3L)
  for (m in seq_len(PSEUDO_DRAWS)) {
    membership <- one_draw(prob) - 1L
    d$high <- as.integer(membership == 2L)
    d$moderate <- as.integer(membership == 1L)
    for (k in seq_along(biomarkers)) {
      formula <- reformulate(c('high', 'moderate', 'assay_age', 'sex',
        'race_indicator', 'education', 'bmi', 'kidney', 'assay_batch',
        'diagnosed_at_assay'), response = biomarkers[k])
      fit <- tryCatch(lm(formula, data = d), error = function(e) NULL)
      if (!is.null(fit) && is.finite(coef(fit)['high'])) {
        estimates[m, k] <- unname(coef(fit)['high'])
        variances[m, k] <- tryCatch(vcov(fit)['high', 'high'],
                                   error = function(e) NA_real_)
      }
    }
  }
  combined <- lapply(seq_along(biomarkers), function(k)
    combine_pseudo_estimates(estimates[, k], variances[, k], design_effect))
  p <- vapply(combined, function(x) unname(x['p']), numeric(1))
  adjusted <- rep(NA_real_, length(p))
  valid <- is.finite(p)
  if (any(valid)) adjusted[valid] <- p.adjust(p[valid], method = 'BH')
  data.frame(biomarker = biomarkers,
    estimate = vapply(combined, function(x) unname(x['beta']), numeric(1)),
    p = p, fdr_p = adjusted, high_n = simulated$high_n,
    posterior_entropy = simulated$entropy)
}

summarize_power <- function(x, grouping, significance, direction = NULL) {
  keys <- do.call(interaction,
    c(x[grouping], list(drop = TRUE, lex.order = TRUE)))
  groups <- split(x, keys)
  output <- lapply(groups, function(d) {
    passed <- is.finite(d[[significance]]) & d[[significance]] < ALPHA
    if (!is.null(direction)) passed <- passed &
      is.finite(d$estimate) & (d$estimate * direction > 0)
    row <- d[1L, grouping, drop = FALSE]
    row$replicates <- nrow(d)
    row$power <- mean(passed)
    row
  })
  do.call(rbind, output)
}

run_aim3a <- function() {
  grid <- expand.grid(cohort = cohorts$cohort,
    n = NA_integer_, prevalence = PROFILE_PREVALENCES,
    entropy = ENTROPIES, event_fraction = EVENT_FRACTIONS,
    true_hr = ADRD_HRS, stringsAsFactors = FALSE)
  grid$n <- cohorts$n[match(grid$cohort, cohorts$cohort)]
  grid$design_effect <- cohorts$design_effect[
    match(grid$cohort, cohorts$cohort)]
  rows <- vector('list', nrow(grid) * REPS)
  counter <- 0L
  for (s in seq_len(nrow(grid))) {
    z <- grid[s, ]
    for (b in seq_len(REPS)) {
      sim <- simulate_survival(z$n, z$prevalence, z$entropy,
        z$true_hr, z$event_fraction)
      fit <- analyze_survival(sim, z$design_effect)
      counter <- counter + 1L
      rows[[counter]] <- data.frame(cohort = z$cohort, n = z$n,
        prevalence = z$prevalence, target_entropy = z$entropy,
        event_fraction = z$event_fraction, true_hr = z$true_hr,
        replicate = b, estimate = unname(fit['beta']),
        p = unname(fit['p']), fitted_hr = exp(unname(fit['beta'])),
        events = unname(fit['events']), deaths = unname(fit['deaths']),
        high_n = unname(fit['high_n']),
        high_events = unname(fit['high_events']),
        entropy = unname(fit['entropy']),
        completed_draws = unname(fit['draws']))
    }
    cat('Aim 3A scenario', s, 'of', nrow(grid), 'complete\n')
  }
  replicates <- do.call(rbind, rows)
  group_names <- c('cohort', 'n', 'prevalence', 'target_entropy',
                   'event_fraction', 'true_hr')
  summary <- summarize_power(replicates, group_names, 'p', direction = 1)
  keys <- do.call(interaction,
    c(replicates[group_names], list(drop = TRUE, lex.order = TRUE)))
  summary$completed_replicates <- vapply(split(replicates, keys), function(d)
    sum(is.finite(d$p)), integer(1))
  summary$mean_events <- vapply(split(replicates, keys), function(d)
    mean(d$events), numeric(1))
  summary$mean_high_events <- vapply(split(replicates, keys), function(d)
    mean(d$high_events), numeric(1))
  write.csv(replicates, file.path(OUT_DIR, 'Aim3A_replicates.csv'),
            row.names = FALSE)
  write.csv(summary, file.path(OUT_DIR, 'Aim3A_power_summary.csv'),
            row.names = FALSE)
  print(summary, row.names = FALSE)
}

run_aim3b <- function() {
  grid <- expand.grid(n = HCAP_NS, prevalence = PROFILE_PREVALENCES,
    entropy = ENTROPIES, standardized_d = BIOMARKER_DS,
    stringsAsFactors = FALSE)
  rows <- vector('list', nrow(grid) * REPS)
  counter <- 0L
  for (s in seq_len(nrow(grid))) {
    z <- grid[s, ]
    for (b in seq_len(REPS)) {
      sim <- simulate_biomarkers(z$n, z$prevalence,
                                 z$entropy, z$standardized_d)
      fit <- analyze_biomarkers(sim)
      fit$n <- z$n; fit$prevalence <- z$prevalence
      fit$target_entropy <- z$entropy
      fit$standardized_d <- z$standardized_d
      fit$replicate <- b
      counter <- counter + 1L
      rows[[counter]] <- fit
    }
    cat('Aim 3B scenario', s, 'of', nrow(grid), 'complete\n')
  }
  replicates <- do.call(rbind, rows)
  group_names <- c('n', 'prevalence', 'target_entropy',
                   'standardized_d', 'biomarker')
  summary <- summarize_power(replicates, group_names, 'fdr_p')
  # Compute directional FDR power explicitly within each scenario.
  keys <- do.call(interaction,
    c(replicates[group_names], list(drop = TRUE, lex.order = TRUE)))
  groups <- split(replicates, keys)
  summary$power <- vapply(groups, function(d) {
    expected_sign <- if (d$biomarker[1L] == 'abeta_log') -1 else 1
    mean(is.finite(d$fdr_p) & d$fdr_p < ALPHA &
           is.finite(d$estimate) & d$estimate * expected_sign > 0)
  }, numeric(1))
  summary$mean_high_n <- vapply(groups, function(d)
    mean(d$high_n), numeric(1))
  write.csv(replicates, file.path(OUT_DIR, 'Aim3B_replicates.csv'),
            row.names = FALSE)
  write.csv(summary, file.path(OUT_DIR, 'Aim3B_FDR_power_summary.csv'),
            row.names = FALSE)
  print(summary, row.names = FALSE)
}

run_aim3a()
run_aim3b()
cat('Planning estimates are conditional on an assumed three-class solution.\n',
    'Aim 3A uses age as the Cox time scale, censors death as the competing event,\n',
    'and repeats posterior pseudo-class draws. Aim 3B uses BH-FDR across three\n',
    'biomarkers per replicate. Assay selection and actual CMS overlap are not\n',
    'known before data acquisition. Class-model recovery is not simulated.\n')

})

cat('Completed selected aims:', paste(CAR_POWER_AIMS, collapse = ', '), '\n')
