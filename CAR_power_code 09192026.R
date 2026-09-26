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
