###############################################################################
# COPD Markov model -> microsimulation: baseline technical validation
#
# Purpose:
#   1. Reproduce the deterministic logic of the simplified benchmark model.
#   2. Run the structurally equivalent patient-level microsimulation.
#   3. Compare state occupancy, discounted costs, QALYs, and strategy results.
#
# Deliberately excluded:
#   - probabilistic sensitivity analysis (parameter uncertainty)
#   - exacerbation-memory extension
#   - heemod and C++ implementations
#   - alternative cycle lengths
#
# This script uses base R only.
###############################################################################

rm(list = ls(all.names = TRUE))
options(scipen = 999)

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------
N_PATIENTS <- 100000L
SEED <- 12345L
RUN_REPEATED_VALIDATION <- FALSE
N_REPEATS <- 100L

OUTPUT_DIR <- file.path(getwd(), "validation_output")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Model parameters from Markov_model_simplified_Github_incl_checks.r
# -----------------------------------------------------------------------------
n_rows <- 241L                 # baseline row + 240 quarterly transitions
cycle_length_months <- 3
annual_discount_rate <- 0.03

# Quarterly progression probabilities
prog_smoker <- c(stage1 = 0.01403, stage2 = 0.01031, stage3 = 0.02298)
prog_former <- c(stage1 = 0.00247, stage2 = 0.00300, stage3 = 0.00686)

# Utilities
utility <- c(stage1 = 0.84, stage2 = 0.79, stage3 = 0.75, stage4 = 0.65)

# Stable-state direct costs per cycle
stable_cost <- c(stage1 = 102.64, stage2 = 184.69,
                 stage3 = 366.57, stage4 = 431.35)

# Moderate and severe exacerbation probabilities by GOLD stage
p_mod <- c(stage2 = 0.18376, stage3 = 0.26767, stage4 = 0.28057)
p_sev <- c(stage2 = 0.03854, stage3 = 0.03774, stage4 = 0.04855)

# Exacerbation disutilities
u_dec_mod <- c(stage2 = -0.01975, stage3 = -0.01875, stage4 = -0.01625)
u_dec_sev <- c(stage2 = -0.03075, stage3 = -0.02919, stage4 = -0.03502)

# Exacerbation direct costs
cost_mod <- c(stage2 = 39.86, stage3 = 52.90, stage4 = 51.52)
cost_sev <- c(stage2 = 2923.94, stage3 = 3109.05, stage4 = 3322.59)

# Strategy weights and intervention costs (per patient, paid once at baseline)
prop_former_intervention <- 0.22
prop_former_usual_care <- 0.06
intervention_cost_per_patient <- 595.92
usual_care_cost_per_patient <- 0

# -----------------------------------------------------------------------------
# Time-dependent mortality schedules
# -----------------------------------------------------------------------------
expand_schedule <- function(values) {
  stopifnot(length(values) == 10L)
  c(rep(values[1:9], each = 20L), rep(values[10], 61L))
}

mort_smoker <- data.frame(
  stage1 = expand_schedule(c(0.00122, 0.00210, 0.00302, 0.00470, 0.00805,
                             0.01358, 0.01970, 0.03400, 0.04729, 0.08840)),
  stage2 = expand_schedule(c(0.00210, 0.00360, 0.00518, 0.00804, 0.01377,
                             0.02318, 0.03354, 0.05758, 0.07969, 0.14672)),
  stage34 = expand_schedule(c(0.00497, 0.00852, 0.01225, 0.01899, 0.03238,
                              0.05417, 0.07783, 0.13138, 0.17900, 0.31397))
)

mort_former <- data.frame(
  stage1 = expand_schedule(c(0.00080, 0.00137, 0.00251, 0.00390, 0.00504,
                             0.00851, 0.01318, 0.02280, 0.03957, 0.07423)),
  stage2 = expand_schedule(c(0.00136, 0.00234, 0.00430, 0.00668, 0.00863,
                             0.01455, 0.02249, 0.03877, 0.06686, 0.12386)),
  stage34 = expand_schedule(c(0.00323, 0.00555, 0.01018, 0.01578, 0.02037,
                              0.03421, 0.05258, 0.08963, 0.15156, 0.26950))
)

# Mortality conditional on a severe exacerbation.
# The first age block is filled from row 1 onward so that the general recursion
# reproduces the benchmark model's special initialization for GOLD stage 2.
ae_mort_smoker <- data.frame(
  stage2 = expand_schedule(c(0.13819, 0.13690, 0.13552, 0.25400, 0.24967,
                             0.24244, 0.30675, 0.28906, 0.27198, 0.21479)),
  stage34 = expand_schedule(c(0.13570, 0.13261, 0.12933, 0.24568, 0.23523,
                              0.21762, 0.27345, 0.22866, 0.18392, 0.02336))
)

ae_mort_former <- data.frame(
  stage2 = expand_schedule(c(0.13883, 0.13798, 0.13629, 0.25503, 0.25356,
                             0.24908, 0.31459, 0.30298, 0.28199, 0.23529)),
  stage34 = expand_schedule(c(0.13721, 0.13520, 0.13116, 0.24813, 0.24462,
                              0.23379, 0.29281, 0.26404, 0.21031, 0.08281))
)

# -----------------------------------------------------------------------------
# Reward functions
# Baseline row is assigned zero costs and zero QALYs, matching the simplified
# benchmark model after removal of mild exacerbations.
# -----------------------------------------------------------------------------
calculate_rewards <- function(trace) {
  required <- c("stage1", "stage2", "stage3", "stage4", "death")
  stopifnot(all(required %in% names(trace)))

  # Double-bracket extraction prevents nested names such as
  # "stage1.stage1", which would make later lookups return NA.
  expected_utility <- c(
    stage1 = utility[["stage1"]],
    stage2 = utility[["stage2"]] + p_mod[["stage2"]] * u_dec_mod[["stage2"]] +
      p_sev[["stage2"]] * u_dec_sev[["stage2"]],
    stage3 = utility[["stage3"]] + p_mod[["stage3"]] * u_dec_mod[["stage3"]] +
      p_sev[["stage3"]] * u_dec_sev[["stage3"]],
    stage4 = utility[["stage4"]] + p_mod[["stage4"]] * u_dec_mod[["stage4"]] +
      p_sev[["stage4"]] * u_dec_sev[["stage4"]]
  )

  expected_cost <- c(
    stage1 = stable_cost[["stage1"]],
    stage2 = stable_cost[["stage2"]] + p_mod[["stage2"]] * cost_mod[["stage2"]] +
      p_sev[["stage2"]] * cost_sev[["stage2"]],
    stage3 = stable_cost[["stage3"]] + p_mod[["stage3"]] * cost_mod[["stage3"]] +
      p_sev[["stage3"]] * cost_sev[["stage3"]],
    stage4 = stable_cost[["stage4"]] + p_mod[["stage4"]] * cost_mod[["stage4"]] +
      p_sev[["stage4"]] * cost_sev[["stage4"]]
  )

  trace$costs <-
    trace$stage1 * expected_cost["stage1"] +
    trace$stage2 * expected_cost["stage2"] +
    trace$stage3 * expected_cost["stage3"] +
    trace$stage4 * expected_cost["stage4"]

  trace$qalys <- (
    trace$stage1 * expected_utility["stage1"] +
    trace$stage2 * expected_utility["stage2"] +
    trace$stage3 * expected_utility["stage3"] +
    trace$stage4 * expected_utility["stage4"]
  ) * cycle_length_months / 12

  trace$costs[1] <- 0
  trace$qalys[1] <- 0

  discount_factor <- 1 / (
    (1 + annual_discount_rate) ^ ((0:(nrow(trace) - 1L)) * cycle_length_months / 12)
  )

  trace$costs_discounted <- trace$costs * discount_factor
  trace$qalys_discounted <- trace$qalys * discount_factor
  trace
}

# -----------------------------------------------------------------------------
# Indirect costs from the simplified benchmark model
# These are retained for a secondary societal-perspective comparison.
# -----------------------------------------------------------------------------
calculate_indirect_costs <- function(trace) {
  costs_stage2_stable <- c(0, rep(112.50, 19), rep(376.37, 20),
                           rep(667.24, 20), rep(237.22, 20), rep(0, 161))
  costs_stage3_stable <- c(0, rep(268.06, 19), rep(896.81, 20),
                           rep(1589.90, 20), rep(565.24, 20), rep(0, 161))
  costs_stage4_stable <- c(0, rep(281.24, 19), rep(940.92, 20),
                           rep(1668.10, 20), rep(593.04, 20), rep(0, 161))

  costs_stage2_mod <- c(rep(491.94, 20), rep(446.16, 20),
                        rep(364.35, 20), rep(115.30, 20), rep(0, 161))
  costs_stage3_mod <- c(rep(550.51, 20), rep(499.27, 20),
                        rep(407.73, 20), rep(129.03, 20), rep(0, 161))
  costs_stage4_mod <- c(rep(491.94, 20), rep(446.16, 20),
                        rep(364.35, 20), rep(115.30, 20), rep(0, 161))

  costs_stage2_sev <- c(rep(367.00, 20), rep(332.85, 20),
                        rep(271.82, 20), rep(86.02, 20), rep(0, 161))
  costs_stage3_sev <- c(rep(663.73, 20), rep(601.96, 20),
                        rep(491.58, 20), rep(155.56, 20), rep(0, 161))
  costs_stage4_sev <- c(rep(648.12, 20), rep(587.79, 20),
                        rep(480.02, 20), rep(151.90, 20), rep(0, 161))

  stopifnot(length(costs_stage2_stable) == nrow(trace))

  trace$indirect_costs <-
    trace$stage2 * (costs_stage2_stable +
      p_mod["stage2"] * costs_stage2_mod +
      p_sev["stage2"] * costs_stage2_sev) +
    trace$stage3 * (costs_stage3_stable +
      p_mod["stage3"] * costs_stage3_mod +
      p_sev["stage3"] * costs_stage3_sev) +
    trace$stage4 * (costs_stage4_stable +
      p_mod["stage4"] * costs_stage4_mod +
      p_sev["stage4"] * costs_stage4_sev)

  discount_factor <- 1 / (
    (1 + annual_discount_rate) ^ ((0:(nrow(trace) - 1L)) * cycle_length_months / 12)
  )
  trace$indirect_costs_discounted <- trace$indirect_costs * discount_factor
  trace
}

# -----------------------------------------------------------------------------
# Deterministic cohort Markov model
# -----------------------------------------------------------------------------
run_markov <- function(n_patients, progression, mortality, ae_mortality) {
  trace <- data.frame(
    row = seq_len(n_rows),
    transition = 0:(n_rows - 1L),
    stage1 = numeric(n_rows),
    stage2 = numeric(n_rows),
    stage3 = numeric(n_rows),
    stage4 = numeric(n_rows),
    death = numeric(n_rows)
  )
  trace$stage1[1] <- n_patients

  for (i in seq_len(n_rows - 1L)) {
    p_death1 <- mortality$stage1[i]
    p_death2 <- mortality$stage2[i] + p_sev["stage2"] * ae_mortality$stage2[i]
    p_death3 <- mortality$stage34[i] + p_sev["stage3"] * ae_mortality$stage34[i]
    p_death4 <- mortality$stage34[i] + p_sev["stage4"] * ae_mortality$stage34[i]

    trace$stage1[i + 1L] <- trace$stage1[i] *
      (1 - progression["stage1"] - p_death1)

    trace$stage2[i + 1L] <- trace$stage2[i] *
      (1 - progression["stage2"] - p_death2) +
      trace$stage1[i] * progression["stage1"]

    trace$stage3[i + 1L] <- trace$stage3[i] *
      (1 - progression["stage3"] - p_death3) +
      trace$stage2[i] * progression["stage2"]

    trace$stage4[i + 1L] <- trace$stage4[i] * (1 - p_death4) +
      trace$stage3[i] * progression["stage3"]

    trace$death[i + 1L] <- n_patients - sum(unlist(
      trace[i + 1L, c("stage1", "stage2", "stage3", "stage4")],
      use.names = FALSE
    ))
  }

  trace <- calculate_rewards(trace)
  calculate_indirect_costs(trace)
}

# -----------------------------------------------------------------------------
# Patient-level microsimulation
# Progression and death are competing, mutually exclusive events.
# Exacerbation consequences remain expected values in the baseline model.
# -----------------------------------------------------------------------------
run_microsimulation <- function(n_patients, progression, mortality,
                                ae_mortality, seed) {
  set.seed(seed)
  state <- rep.int(1L, n_patients)

  trace <- data.frame(
    row = seq_len(n_rows),
    transition = 0:(n_rows - 1L),
    stage1 = integer(n_rows),
    stage2 = integer(n_rows),
    stage3 = integer(n_rows),
    stage4 = integer(n_rows),
    death = integer(n_rows)
  )
  trace$stage1[1] <- n_patients

  for (i in seq_len(n_rows - 1L)) {
    # Snapshot at the beginning of the transition prevents multiple moves
    # by one patient within the same cycle.
    idx1 <- which(state == 1L)
    idx2 <- which(state == 2L)
    idx3 <- which(state == 3L)
    idx4 <- which(state == 4L)

    if (length(idx1) > 0L) {
      u <- runif(length(idx1))
      move <- u < progression["stage1"]
      die <- !move & u < progression["stage1"] + mortality$stage1[i]
      state[idx1[move]] <- 2L
      state[idx1[die]] <- 5L
    }

    if (length(idx2) > 0L) {
      p_death <- mortality$stage2[i] +
        p_sev["stage2"] * ae_mortality$stage2[i]
      u <- runif(length(idx2))
      move <- u < progression["stage2"]
      die <- !move & u < progression["stage2"] + p_death
      state[idx2[move]] <- 3L
      state[idx2[die]] <- 5L
    }

    if (length(idx3) > 0L) {
      p_death <- mortality$stage34[i] +
        p_sev["stage3"] * ae_mortality$stage34[i]
      u <- runif(length(idx3))
      move <- u < progression["stage3"]
      die <- !move & u < progression["stage3"] + p_death
      state[idx3[move]] <- 4L
      state[idx3[die]] <- 5L
    }

    if (length(idx4) > 0L) {
      p_death <- mortality$stage34[i] +
        p_sev["stage4"] * ae_mortality$stage34[i]
      die <- runif(length(idx4)) < p_death
      state[idx4[die]] <- 5L
    }

    counts <- tabulate(state, nbins = 5L)
    trace[i + 1L, c("stage1", "stage2", "stage3", "stage4", "death")] <- counts
  }

  trace <- calculate_rewards(trace)
  calculate_indirect_costs(trace)
}

# -----------------------------------------------------------------------------
# Validation helpers
# -----------------------------------------------------------------------------
check_trace <- function(trace, n_patients, label, stochastic = FALSE) {
  state_cols <- c("stage1", "stage2", "stage3", "stage4", "death")
  totals <- rowSums(trace[state_cols])

  if (stochastic) {
    stopifnot(all(totals == n_patients))
  } else {
    stopifnot(max(abs(totals - n_patients)) < 1e-6)
  }
  stopifnot(all(as.matrix(trace[state_cols]) >= -1e-8))
  stopifnot(trace$costs[1] == 0, trace$qalys[1] == 0)
  stopifnot(
    all(is.finite(trace$costs)),
    all(is.finite(trace$qalys)),
    all(is.finite(trace$costs_discounted)),
    all(is.finite(trace$qalys_discounted)),
    all(is.finite(trace$indirect_costs_discounted))
  )

  invisible(TRUE)
}

compare_traces <- function(markov, microsim, label, n_patients) {
  state_cols <- c("stage1", "stage2", "stage3", "stage4", "death")
  diff_matrix <- as.matrix(microsim[state_cols]) - as.matrix(markov[state_cols])

  tvd <- 0.5 * rowSums(abs(diff_matrix)) / n_patients

  cycle_table <- data.frame(
    group = label,
    row = markov$row,
    transition = markov$transition,
    tvd = tvd,
    diff_stage1 = diff_matrix[, "stage1"],
    diff_stage2 = diff_matrix[, "stage2"],
    diff_stage3 = diff_matrix[, "stage3"],
    diff_stage4 = diff_matrix[, "stage4"],
    diff_death = diff_matrix[, "death"]
  )

  endpoint <- data.frame(
    group = label,
    metric = c("discounted_direct_costs", "discounted_indirect_costs",
               "discounted_qalys"),
    markov = c(sum(markov$costs_discounted),
               sum(markov$indirect_costs_discounted),
               sum(markov$qalys_discounted)),
    microsimulation = c(sum(microsim$costs_discounted),
                        sum(microsim$indirect_costs_discounted),
                        sum(microsim$qalys_discounted))
  )
  endpoint$difference <- endpoint$microsimulation - endpoint$markov
  endpoint$relative_difference_percent <-
    100 * endpoint$difference / endpoint$markov

  summary <- data.frame(
    group = label,
    max_tvd_percent = 100 * max(tvd),
    mean_tvd_percent = 100 * mean(tvd),
    row_of_max_tvd = which.max(tvd),
    transition_of_max_tvd = markov$transition[which.max(tvd)]
  )

  list(cycle = cycle_table, endpoint = endpoint, summary = summary)
}

subcohort_totals_per_patient <- function(trace, n_patients) {
  c(
    direct_costs = sum(trace$costs_discounted) / n_patients,
    indirect_costs = sum(trace$indirect_costs_discounted) / n_patients,
    qalys = sum(trace$qalys_discounted) / n_patients
  )
}

calculate_strategy_results <- function(smoker_trace, former_trace, n_patients) {
  smoker <- subcohort_totals_per_patient(smoker_trace, n_patients)
  former <- subcohort_totals_per_patient(former_trace, n_patients)

  weighted <- function(prop_former, intervention_cost) {
    out <- (1 - prop_former) * smoker + prop_former * former
    out["direct_costs"] <- out["direct_costs"] + intervention_cost
    out
  }

  usual_care <- weighted(prop_former_usual_care,
                         usual_care_cost_per_patient)
  intervention <- weighted(prop_former_intervention,
                           intervention_cost_per_patient)
  incremental <- intervention - usual_care

  direct_icer <- incremental["direct_costs"] / incremental["qalys"]
  societal_icer <- (
    incremental["direct_costs"] + incremental["indirect_costs"]
  ) / incremental["qalys"]

  data.frame(
    strategy = c("usual_care", "intervention", "incremental"),
    direct_costs_per_patient = c(usual_care["direct_costs"],
                                 intervention["direct_costs"],
                                 incremental["direct_costs"]),
    indirect_costs_per_patient = c(usual_care["indirect_costs"],
                                   intervention["indirect_costs"],
                                   incremental["indirect_costs"]),
    qalys_per_patient = c(usual_care["qalys"],
                          intervention["qalys"],
                          incremental["qalys"]),
    direct_icer = c(NA, NA, direct_icer),
    societal_icer = c(NA, NA, societal_icer)
  )
}

# -----------------------------------------------------------------------------
# Probability checks before model execution
# -----------------------------------------------------------------------------
check_probabilities <- function(progression, mortality, ae_mortality, label) {
  p2 <- progression["stage2"] + mortality$stage2 +
    p_sev["stage2"] * ae_mortality$stage2
  p3 <- progression["stage3"] + mortality$stage34 +
    p_sev["stage3"] * ae_mortality$stage34
  p4 <- mortality$stage34 + p_sev["stage4"] * ae_mortality$stage34

  stopifnot(
    all(progression >= 0),
    progression["stage1"] + max(mortality$stage1) <= 1,
    max(p2) <= 1,
    max(p3) <= 1,
    max(p4) <= 1
  )
  invisible(TRUE)
}

check_probabilities(prog_smoker, mort_smoker, ae_mort_smoker, "smokers")
check_probabilities(prog_former, mort_former, ae_mort_former, "former smokers")

# -----------------------------------------------------------------------------
# Main run
# -----------------------------------------------------------------------------
cat("Running deterministic Markov benchmark...\n")
markov_smoker <- run_markov(N_PATIENTS, prog_smoker,
                            mort_smoker, ae_mort_smoker)
markov_former <- run_markov(N_PATIENTS, prog_former,
                            mort_former, ae_mort_former)

cat("Running patient-level microsimulation...\n")
# Different but reproducible seeds are used for the two subcohorts.
microsim_smoker <- run_microsimulation(
  N_PATIENTS, prog_smoker, mort_smoker, ae_mort_smoker, SEED
)
microsim_former <- run_microsimulation(
  N_PATIENTS, prog_former, mort_former, ae_mort_former, SEED + 1L
)

check_trace(markov_smoker, N_PATIENTS, "Markov smokers", stochastic = FALSE)
check_trace(markov_former, N_PATIENTS, "Markov former smokers", stochastic = FALSE)
check_trace(microsim_smoker, N_PATIENTS, "Microsim smokers", stochastic = TRUE)
check_trace(microsim_former, N_PATIENTS, "Microsim former smokers", stochastic = TRUE)

comparison_smoker <- compare_traces(
  markov_smoker, microsim_smoker, "smokers", N_PATIENTS
)
comparison_former <- compare_traces(
  markov_former, microsim_former, "former_smokers", N_PATIENTS
)

validation_summary <- rbind(comparison_smoker$summary,
                            comparison_former$summary)
endpoint_comparison <- rbind(comparison_smoker$endpoint,
                             comparison_former$endpoint)

strategy_markov <- calculate_strategy_results(
  markov_smoker, markov_former, N_PATIENTS
)
strategy_markov$model <- "Markov"

strategy_microsim <- calculate_strategy_results(
  microsim_smoker, microsim_former, N_PATIENTS
)
strategy_microsim$model <- "Microsimulation"

strategy_comparison <- rbind(strategy_markov, strategy_microsim)
strategy_comparison <- strategy_comparison[
  c("model", "strategy", "direct_costs_per_patient",
    "indirect_costs_per_patient", "qalys_per_patient",
    "direct_icer", "societal_icer")
]

# Reference values reported in the simplified benchmark script.
# These are checks on the deterministic Markov implementation, not acceptance
# thresholds for the stochastic microsimulation.
reference_check <- data.frame(
  result = c("direct_costs_usual_care", "direct_costs_intervention",
             "qalys_usual_care", "qalys_intervention", "direct_icer"),
  reported_benchmark = c(11500.78, 11711.26, 12.68995, 13.23588, 385.5405),
  current_markov = c(
    strategy_markov$direct_costs_per_patient[strategy_markov$strategy == "usual_care"],
    strategy_markov$direct_costs_per_patient[strategy_markov$strategy == "intervention"],
    strategy_markov$qalys_per_patient[strategy_markov$strategy == "usual_care"],
    strategy_markov$qalys_per_patient[strategy_markov$strategy == "intervention"],
    strategy_markov$direct_icer[strategy_markov$strategy == "incremental"]
  )
)
reference_check$difference <-
  reference_check$current_markov - reference_check$reported_benchmark

# The published comments in the benchmark script are rounded. These tolerances
# allow only rounding-level deviations in the deterministic implementation.
reference_tolerance <- c(0.02, 0.02, 0.00002, 0.00002, 0.02)
if (anyNA(reference_check$current_markov) ||
    any(!is.finite(reference_check$current_markov))) {
  stop(
    "Internal validation error: deterministic benchmark outputs contain ",
    "missing or non-finite values."
  )
}
if (any(abs(reference_check$difference) > reference_tolerance)) {
  stop(
    "Deterministic Markov results do not match the reported simplified ",
    "benchmark values within rounding tolerance. Inspect reference_check."
  )
}

# -----------------------------------------------------------------------------
# Save output files
# -----------------------------------------------------------------------------
write.csv(validation_summary,
          file.path(OUTPUT_DIR, "validation_summary.csv"), row.names = FALSE)
write.csv(endpoint_comparison,
          file.path(OUTPUT_DIR, "endpoint_comparison.csv"), row.names = FALSE)
write.csv(strategy_comparison,
          file.path(OUTPUT_DIR, "strategy_comparison.csv"), row.names = FALSE)
write.csv(reference_check,
          file.path(OUTPUT_DIR, "benchmark_reference_check.csv"), row.names = FALSE)
write.csv(comparison_smoker$cycle,
          file.path(OUTPUT_DIR, "cycle_comparison_smokers.csv"), row.names = FALSE)
write.csv(comparison_former$cycle,
          file.path(OUTPUT_DIR, "cycle_comparison_former_smokers.csv"),
          row.names = FALSE)

# State-trace plots
png(file.path(OUTPUT_DIR, "state_trace_smokers.png"),
    width = 1400, height = 1000, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (state_name in c("stage1", "stage2", "stage3", "stage4")) {
  plot(markov_smoker$transition, markov_smoker[[state_name]], type = "l",
       xlab = "Quarterly transition", ylab = "Patients",
       main = paste("Smokers:", state_name))
  lines(microsim_smoker$transition, microsim_smoker[[state_name]], lty = 2)
}
dev.off()

png(file.path(OUTPUT_DIR, "state_trace_former_smokers.png"),
    width = 1400, height = 1000, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (state_name in c("stage1", "stage2", "stage3", "stage4")) {
  plot(markov_former$transition, markov_former[[state_name]], type = "l",
       xlab = "Quarterly transition", ylab = "Patients",
       main = paste("Former smokers:", state_name))
  lines(microsim_former$transition, microsim_former[[state_name]], lty = 2)
}
dev.off()

# -----------------------------------------------------------------------------
# Optional repeated-seed validation
# -----------------------------------------------------------------------------
if (RUN_REPEATED_VALIDATION) {
  repeated <- vector("list", N_REPEATS)

  markov_incremental <- strategy_markov[
    strategy_markov$strategy == "incremental", ,
    drop = FALSE
  ]

  extract_relative_difference <- function(comparison, metric_name) {
    comparison$endpoint$relative_difference_percent[
      comparison$endpoint$metric == metric_name
    ]
  }

  for (r in seq_len(N_REPEATS)) {
    seed_r <- SEED + 1000L * r

    sim_s <- run_microsimulation(
      N_PATIENTS, prog_smoker, mort_smoker, ae_mort_smoker, seed_r
    )
    sim_f <- run_microsimulation(
      N_PATIENTS, prog_former, mort_former, ae_mort_former, seed_r + 1L
    )

    cmp_s <- compare_traces(
      markov_smoker, sim_s, "smokers", N_PATIENTS
    )
    cmp_f <- compare_traces(
      markov_former, sim_f, "former_smokers", N_PATIENTS
    )

    strat <- calculate_strategy_results(sim_s, sim_f, N_PATIENTS)
    incr <- strat[strat$strategy == "incremental", , drop = FALSE]

    repeated[[r]] <- data.frame(
      repetition = r,
      seed = seed_r,

      max_tvd_smokers_percent = cmp_s$summary$max_tvd_percent,
      mean_tvd_smokers_percent = cmp_s$summary$mean_tvd_percent,
      max_tvd_former_percent = cmp_f$summary$max_tvd_percent,
      mean_tvd_former_percent = cmp_f$summary$mean_tvd_percent,

      direct_costs_rel_diff_smokers_percent =
        extract_relative_difference(cmp_s, "discounted_direct_costs"),
      qalys_rel_diff_smokers_percent =
        extract_relative_difference(cmp_s, "discounted_qalys"),
      indirect_costs_rel_diff_smokers_percent =
        extract_relative_difference(cmp_s, "discounted_indirect_costs"),

      direct_costs_rel_diff_former_percent =
        extract_relative_difference(cmp_f, "discounted_direct_costs"),
      qalys_rel_diff_former_percent =
        extract_relative_difference(cmp_f, "discounted_qalys"),
      indirect_costs_rel_diff_former_percent =
        extract_relative_difference(cmp_f, "discounted_indirect_costs"),

      incremental_direct_costs = incr$direct_costs_per_patient,
      incremental_qalys = incr$qalys_per_patient,
      direct_icer = incr$direct_icer,

      incremental_direct_costs_difference =
        incr$direct_costs_per_patient -
        markov_incremental$direct_costs_per_patient,
      incremental_qalys_difference =
        incr$qalys_per_patient -
        markov_incremental$qalys_per_patient,
      direct_icer_difference =
        incr$direct_icer - markov_incremental$direct_icer
    )
  }

  repeated_results <- do.call(rbind, repeated)

  write.csv(
    repeated_results,
    file.path(OUTPUT_DIR, "repeated_seed_validation.csv"),
    row.names = FALSE
  )

  summarise_metric <- function(x, benchmark = NA_real_) {
    has_benchmark <- is.finite(benchmark)

    data.frame(
      benchmark = if (has_benchmark) benchmark else NA_real_,
      mean = mean(x),
      standard_deviation = sd(x),
      median = median(x),
      percentile_2_5 = unname(quantile(x, 0.025)),
      percentile_97_5 = unname(quantile(x, 0.975)),
      minimum = min(x),
      maximum = max(x),
      mean_signed_difference =
        if (has_benchmark) mean(x - benchmark) else NA_real_,
      mean_absolute_difference =
        if (has_benchmark) mean(abs(x - benchmark)) else NA_real_,
      n_above_benchmark =
        if (has_benchmark) sum(x > benchmark) else NA_integer_,
      n_below_benchmark =
        if (has_benchmark) sum(x < benchmark) else NA_integer_
    )
  }

  summary_items <- list(
    max_tvd_smokers_percent = summarise_metric(
      repeated_results$max_tvd_smokers_percent
    ),
    max_tvd_former_percent = summarise_metric(
      repeated_results$max_tvd_former_percent
    ),
    direct_costs_rel_diff_smokers_percent = summarise_metric(
      repeated_results$direct_costs_rel_diff_smokers_percent,
      benchmark = 0
    ),
    qalys_rel_diff_smokers_percent = summarise_metric(
      repeated_results$qalys_rel_diff_smokers_percent,
      benchmark = 0
    ),
    direct_costs_rel_diff_former_percent = summarise_metric(
      repeated_results$direct_costs_rel_diff_former_percent,
      benchmark = 0
    ),
    qalys_rel_diff_former_percent = summarise_metric(
      repeated_results$qalys_rel_diff_former_percent,
      benchmark = 0
    ),
    incremental_direct_costs = summarise_metric(
      repeated_results$incremental_direct_costs,
      benchmark = markov_incremental$direct_costs_per_patient
    ),
    incremental_qalys = summarise_metric(
      repeated_results$incremental_qalys,
      benchmark = markov_incremental$qalys_per_patient
    ),
    direct_icer = summarise_metric(
      repeated_results$direct_icer,
      benchmark = markov_incremental$direct_icer
    )
  )

  repeated_summary <- do.call(
    rbind,
    lapply(names(summary_items), function(metric_name) {
      out <- summary_items[[metric_name]]
      out$metric <- metric_name
      out
    })
  )

  repeated_summary <- repeated_summary[
    c(
      "metric", "benchmark", "mean", "standard_deviation", "median",
      "percentile_2_5", "percentile_97_5", "minimum", "maximum",
      "mean_signed_difference", "mean_absolute_difference",
      "n_above_benchmark", "n_below_benchmark"
    )
  ]

  rownames(repeated_summary) <- NULL

  write.csv(
    repeated_summary,
    file.path(OUTPUT_DIR, "repeated_seed_summary.csv"),
    row.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# Console and text report
# -----------------------------------------------------------------------------
report <- capture.output({
  cat("\n============================================================\n")
  cat("BASELINE TECHNICAL VALIDATION\n")
  cat("============================================================\n")
  cat("Patients per subcohort:", N_PATIENTS, "\n")
  cat("Seed smokers:", SEED, "\n")
  cat("Seed former smokers:", SEED + 1L, "\n\n")

  cat("1. Deterministic benchmark reference check\n")
  print(reference_check, row.names = FALSE, digits = 10)

  cat("\n2. Cycle-specific state-distribution agreement\n")
  print(validation_summary, row.names = FALSE, digits = 6)

  cat("\n3. Discounted endpoint agreement by subcohort\n")
  print(endpoint_comparison, row.names = FALSE, digits = 10)

  cat("\n4. Strategy-level results per patient\n")
  print(strategy_comparison, row.names = FALSE, digits = 10)

  cat("\nInterpretation:\n")
  cat("- The Markov reference check should differ only through rounding.\n")
  cat("- Microsimulation endpoint differences should be small and may be\n")
  cat("  positive or negative because of first-order Monte Carlo error.\n")
  cat("- State totals are automatically checked in every row.\n")
  if (RUN_REPEATED_VALIDATION) {
    cat("\n5. Repeated-seed validation summary\n")
    cat("Number of repetitions:", N_REPEATS, "\n")
    print(repeated_summary, row.names = FALSE, digits = 10)
  } else {
    cat("- A systematic assessment across seeds is available by setting\n")
    cat("  RUN_REPEATED_VALIDATION <- TRUE near the top of the script.\n")
  }

  cat("\nOutput directory:\n", normalizePath(OUTPUT_DIR), "\n", sep = "")
})

cat(paste(report, collapse = "\n"), "\n")
writeLines(report, file.path(OUTPUT_DIR, "validation_report.txt"))

cat("\nValidation script completed without failed checks.\n")
