###############################################################################
# COPD cohort Markov model and equivalent patient-level microsimulation
# Clean baseline implementation
#
# Model lineage:
#   Original COPD model: Petra Menn
#   R implementation, simplification, and optimization: Tobias Niedermaier
#
# Purpose:
#   - Run the simplified deterministic cohort Markov model.
#   - Run the technically equivalent patient-level microsimulation.
#   - Calculate discounted direct costs, indirect costs, and QALYs.
#   - Aggregate current-smoker and former-smoker subcohorts into usual-care
#     and smoking-cessation intervention strategies.
#
# Deliberately excluded from this baseline file:
#   - parameter uncertainty / probabilistic sensitivity analysis
#   - patient-level exacerbation events
#   - exacerbation-history memory
#   - heemod and C++ implementations
#
# This file uses base R only.
###############################################################################

rm(list = ls(all.names = TRUE))
options(scipen = 999)

# =============================================================================
# 1. User settings
# =============================================================================
N_PATIENTS <- 100000L
SEED_SMOKERS <- 12345L
SEED_FORMER_SMOKERS <- 12346L

SAVE_STATE_TRACES <- TRUE
SAVE_TRACE_PLOTS <- TRUE

OUTPUT_DIR <- file.path(getwd(), "model_output_baseline")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 2. General model settings
# =============================================================================
N_ROWS <- 241L                    # baseline row + 240 quarterly transitions
CYCLE_LENGTH_MONTHS <- 3
ANNUAL_DISCOUNT_RATE <- 0.03

# The source model used a 60-year horizon represented by a baseline row and
# 240 subsequent quarterly transitions. Rewards are zero in the baseline row
# and are assigned to rows 2 to 241.

# =============================================================================
# 3. Model parameters
# =============================================================================

# Quarterly disease-progression probabilities
PROGRESSION_SMOKER <- c(
  stage1 = 0.01403,
  stage2 = 0.01031,
  stage3 = 0.02298
)

PROGRESSION_FORMER_SMOKER <- c(
  stage1 = 0.00247,
  stage2 = 0.00300,
  stage3 = 0.00686
)

# State utilities
UTILITY <- c(
  stage1 = 0.84,
  stage2 = 0.79,
  stage3 = 0.75,
  stage4 = 0.65
)

# Stable-state direct costs per quarter
STABLE_DIRECT_COST <- c(
  stage1 = 102.64,
  stage2 = 184.69,
  stage3 = 366.57,
  stage4 = 431.35
)

# Moderate and severe exacerbation probabilities per quarter
P_MODERATE_AE <- c(
  stage2 = 0.18376,
  stage3 = 0.26767,
  stage4 = 0.28057
)

P_SEVERE_AE <- c(
  stage2 = 0.03854,
  stage3 = 0.03774,
  stage4 = 0.04855
)

# Utility decrements associated with exacerbations
DISUTILITY_MODERATE_AE <- c(
  stage2 = -0.01975,
  stage3 = -0.01875,
  stage4 = -0.01625
)

DISUTILITY_SEVERE_AE <- c(
  stage2 = -0.03075,
  stage3 = -0.02919,
  stage4 = -0.03502
)

# Direct costs associated with exacerbations
DIRECT_COST_MODERATE_AE <- c(
  stage2 = 39.86,
  stage3 = 52.90,
  stage4 = 51.52
)

DIRECT_COST_SEVERE_AE <- c(
  stage2 = 2923.94,
  stage3 = 3109.05,
  stage4 = 3322.59
)

# Strategy composition and intervention costs
PROP_FORMER_SMOKER_USUAL_CARE <- 0.06
PROP_FORMER_SMOKER_INTERVENTION <- 0.22

USUAL_CARE_COST_PER_PATIENT <- 0
INTERVENTION_COST_PER_PATIENT <- 595.92

# =============================================================================
# 4. Time-dependent mortality schedules
# =============================================================================
expand_schedule <- function(values) {
  stopifnot(length(values) == 10L)
  c(rep(values[1:9], each = 20L), rep(values[10], 61L))
}

MORTALITY_SMOKER <- data.frame(
  stage1 = expand_schedule(c(
    0.00122, 0.00210, 0.00302, 0.00470, 0.00805,
    0.01358, 0.01970, 0.03400, 0.04729, 0.08840
  )),
  stage2 = expand_schedule(c(
    0.00210, 0.00360, 0.00518, 0.00804, 0.01377,
    0.02318, 0.03354, 0.05758, 0.07969, 0.14672
  )),
  stage34 = expand_schedule(c(
    0.00497, 0.00852, 0.01225, 0.01899, 0.03238,
    0.05417, 0.07783, 0.13138, 0.17900, 0.31397
  ))
)

MORTALITY_FORMER_SMOKER <- data.frame(
  stage1 = expand_schedule(c(
    0.00080, 0.00137, 0.00251, 0.00390, 0.00504,
    0.00851, 0.01318, 0.02280, 0.03957, 0.07423
  )),
  stage2 = expand_schedule(c(
    0.00136, 0.00234, 0.00430, 0.00668, 0.00863,
    0.01455, 0.02249, 0.03877, 0.06686, 0.12386
  )),
  stage34 = expand_schedule(c(
    0.00323, 0.00555, 0.01018, 0.01578, 0.02037,
    0.03421, 0.05258, 0.08963, 0.15156, 0.26950
  ))
)

# Mortality conditional on a severe exacerbation. In the baseline model,
# excess AE mortality is included as an expected probability:
# P(severe AE) * P(death | severe AE).
AE_MORTALITY_SMOKER <- data.frame(
  stage2 = expand_schedule(c(
    0.13819, 0.13690, 0.13552, 0.25400, 0.24967,
    0.24244, 0.30675, 0.28906, 0.27198, 0.21479
  )),
  stage34 = expand_schedule(c(
    0.13570, 0.13261, 0.12933, 0.24568, 0.23523,
    0.21762, 0.27345, 0.22866, 0.18392, 0.02336
  ))
)

AE_MORTALITY_FORMER_SMOKER <- data.frame(
  stage2 = expand_schedule(c(
    0.13883, 0.13798, 0.13629, 0.25503, 0.25356,
    0.24908, 0.31459, 0.30298, 0.28199, 0.23529
  )),
  stage34 = expand_schedule(c(
    0.13721, 0.13520, 0.13116, 0.24813, 0.24462,
    0.23379, 0.29281, 0.26404, 0.21031, 0.08281
  ))
)

# =============================================================================
# 5. Shared helper functions
# =============================================================================
get_discount_factor <- function(n_rows) {
  1 / (
    (1 + ANNUAL_DISCOUNT_RATE) ^
      ((0:(n_rows - 1L)) * CYCLE_LENGTH_MONTHS / 12)
  )
}

calculate_rewards <- function(trace) {
  required_columns <- c("stage1", "stage2", "stage3", "stage4", "death")
  stopifnot(all(required_columns %in% names(trace)))

  expected_utility <- c(
    stage1 = UTILITY[["stage1"]],
    stage2 = UTILITY[["stage2"]] +
      P_MODERATE_AE[["stage2"]] * DISUTILITY_MODERATE_AE[["stage2"]] +
      P_SEVERE_AE[["stage2"]] * DISUTILITY_SEVERE_AE[["stage2"]],
    stage3 = UTILITY[["stage3"]] +
      P_MODERATE_AE[["stage3"]] * DISUTILITY_MODERATE_AE[["stage3"]] +
      P_SEVERE_AE[["stage3"]] * DISUTILITY_SEVERE_AE[["stage3"]],
    stage4 = UTILITY[["stage4"]] +
      P_MODERATE_AE[["stage4"]] * DISUTILITY_MODERATE_AE[["stage4"]] +
      P_SEVERE_AE[["stage4"]] * DISUTILITY_SEVERE_AE[["stage4"]]
  )

  expected_direct_cost <- c(
    stage1 = STABLE_DIRECT_COST[["stage1"]],
    stage2 = STABLE_DIRECT_COST[["stage2"]] +
      P_MODERATE_AE[["stage2"]] * DIRECT_COST_MODERATE_AE[["stage2"]] +
      P_SEVERE_AE[["stage2"]] * DIRECT_COST_SEVERE_AE[["stage2"]],
    stage3 = STABLE_DIRECT_COST[["stage3"]] +
      P_MODERATE_AE[["stage3"]] * DIRECT_COST_MODERATE_AE[["stage3"]] +
      P_SEVERE_AE[["stage3"]] * DIRECT_COST_SEVERE_AE[["stage3"]],
    stage4 = STABLE_DIRECT_COST[["stage4"]] +
      P_MODERATE_AE[["stage4"]] * DIRECT_COST_MODERATE_AE[["stage4"]] +
      P_SEVERE_AE[["stage4"]] * DIRECT_COST_SEVERE_AE[["stage4"]]
  )

  trace$direct_costs <-
    trace$stage1 * expected_direct_cost[["stage1"]] +
    trace$stage2 * expected_direct_cost[["stage2"]] +
    trace$stage3 * expected_direct_cost[["stage3"]] +
    trace$stage4 * expected_direct_cost[["stage4"]]

  trace$qalys <- (
    trace$stage1 * expected_utility[["stage1"]] +
    trace$stage2 * expected_utility[["stage2"]] +
    trace$stage3 * expected_utility[["stage3"]] +
    trace$stage4 * expected_utility[["stage4"]]
  ) * CYCLE_LENGTH_MONTHS / 12

  # The first row is the baseline state distribution rather than a model cycle.
  trace$direct_costs[1] <- 0
  trace$qalys[1] <- 0

  discount_factor <- get_discount_factor(nrow(trace))
  trace$direct_costs_discounted <- trace$direct_costs * discount_factor
  trace$qalys_discounted <- trace$qalys * discount_factor

  trace
}

calculate_indirect_costs <- function(trace) {
  costs_stage2_stable <- c(
    0, rep(112.50, 19), rep(376.37, 20), rep(667.24, 20),
    rep(237.22, 20), rep(0, 161)
  )
  costs_stage3_stable <- c(
    0, rep(268.06, 19), rep(896.81, 20), rep(1589.90, 20),
    rep(565.24, 20), rep(0, 161)
  )
  costs_stage4_stable <- c(
    0, rep(281.24, 19), rep(940.92, 20), rep(1668.10, 20),
    rep(593.04, 20), rep(0, 161)
  )

  costs_stage2_moderate <- c(
    rep(491.94, 20), rep(446.16, 20), rep(364.35, 20),
    rep(115.30, 20), rep(0, 161)
  )
  costs_stage3_moderate <- c(
    rep(550.51, 20), rep(499.27, 20), rep(407.73, 20),
    rep(129.03, 20), rep(0, 161)
  )
  costs_stage4_moderate <- c(
    rep(491.94, 20), rep(446.16, 20), rep(364.35, 20),
    rep(115.30, 20), rep(0, 161)
  )

  costs_stage2_severe <- c(
    rep(367.00, 20), rep(332.85, 20), rep(271.82, 20),
    rep(86.02, 20), rep(0, 161)
  )
  costs_stage3_severe <- c(
    rep(663.73, 20), rep(601.96, 20), rep(491.58, 20),
    rep(155.56, 20), rep(0, 161)
  )
  costs_stage4_severe <- c(
    rep(648.12, 20), rep(587.79, 20), rep(480.02, 20),
    rep(151.90, 20), rep(0, 161)
  )

  stopifnot(
    length(costs_stage2_stable) == nrow(trace),
    length(costs_stage2_moderate) == nrow(trace),
    length(costs_stage2_severe) == nrow(trace)
  )

  trace$indirect_costs <-
    trace$stage2 * (
      costs_stage2_stable +
        P_MODERATE_AE[["stage2"]] * costs_stage2_moderate +
        P_SEVERE_AE[["stage2"]] * costs_stage2_severe
    ) +
    trace$stage3 * (
      costs_stage3_stable +
        P_MODERATE_AE[["stage3"]] * costs_stage3_moderate +
        P_SEVERE_AE[["stage3"]] * costs_stage3_severe
    ) +
    trace$stage4 * (
      costs_stage4_stable +
        P_MODERATE_AE[["stage4"]] * costs_stage4_moderate +
        P_SEVERE_AE[["stage4"]] * costs_stage4_severe
    )

  discount_factor <- get_discount_factor(nrow(trace))
  trace$indirect_costs_discounted <- trace$indirect_costs * discount_factor

  trace
}

check_probabilities <- function(progression, mortality, ae_mortality, label) {
  p_leave_stage1 <- progression[["stage1"]] + mortality$stage1
  p_leave_stage2 <- progression[["stage2"]] + mortality$stage2 +
    P_SEVERE_AE[["stage2"]] * ae_mortality$stage2
  p_leave_stage3 <- progression[["stage3"]] + mortality$stage34 +
    P_SEVERE_AE[["stage3"]] * ae_mortality$stage34
  p_leave_stage4 <- mortality$stage34 +
    P_SEVERE_AE[["stage4"]] * ae_mortality$stage34

  if (
    any(progression < 0) ||
    any(p_leave_stage1 < 0 | p_leave_stage1 > 1) ||
    any(p_leave_stage2 < 0 | p_leave_stage2 > 1) ||
    any(p_leave_stage3 < 0 | p_leave_stage3 > 1) ||
    any(p_leave_stage4 < 0 | p_leave_stage4 > 1)
  ) {
    stop("Invalid transition probabilities for ", label, ".")
  }

  invisible(TRUE)
}

check_trace <- function(trace, n_patients, label, stochastic) {
  state_columns <- c("stage1", "stage2", "stage3", "stage4", "death")
  totals <- rowSums(trace[state_columns])

  if (stochastic) {
    if (!all(totals == n_patients)) {
      stop("Population conservation failed for ", label, ".")
    }
  } else if (max(abs(totals - n_patients)) >= 1e-6) {
    stop("Population conservation failed for ", label, ".")
  }

  if (any(as.matrix(trace[state_columns]) < -1e-8)) {
    stop("Negative state occupancy detected for ", label, ".")
  }

  numeric_outputs <- c(
    "direct_costs", "direct_costs_discounted", "indirect_costs",
    "indirect_costs_discounted", "qalys", "qalys_discounted"
  )

  if (any(!is.finite(as.matrix(trace[numeric_outputs])))) {
    stop("Missing or non-finite rewards detected for ", label, ".")
  }

  if (trace$direct_costs[1] != 0 || trace$qalys[1] != 0) {
    stop("Baseline rewards are not zero for ", label, ".")
  }

  invisible(TRUE)
}

# =============================================================================
# 6. Deterministic cohort Markov model
# =============================================================================
run_markov_model <- function(n_patients, progression, mortality, ae_mortality) {
  trace <- data.frame(
    row = seq_len(N_ROWS),
    transition = 0:(N_ROWS - 1L),
    stage1 = numeric(N_ROWS),
    stage2 = numeric(N_ROWS),
    stage3 = numeric(N_ROWS),
    stage4 = numeric(N_ROWS),
    death = numeric(N_ROWS)
  )

  trace$stage1[1] <- n_patients

  for (i in seq_len(N_ROWS - 1L)) {
    p_death_stage1 <- mortality$stage1[i]
    p_death_stage2 <- mortality$stage2[i] +
      P_SEVERE_AE[["stage2"]] * ae_mortality$stage2[i]
    p_death_stage3 <- mortality$stage34[i] +
      P_SEVERE_AE[["stage3"]] * ae_mortality$stage34[i]
    p_death_stage4 <- mortality$stage34[i] +
      P_SEVERE_AE[["stage4"]] * ae_mortality$stage34[i]

    trace$stage1[i + 1L] <- trace$stage1[i] * (
      1 - progression[["stage1"]] - p_death_stage1
    )

    trace$stage2[i + 1L] <-
      trace$stage2[i] * (
        1 - progression[["stage2"]] - p_death_stage2
      ) +
      trace$stage1[i] * progression[["stage1"]]

    trace$stage3[i + 1L] <-
      trace$stage3[i] * (
        1 - progression[["stage3"]] - p_death_stage3
      ) +
      trace$stage2[i] * progression[["stage2"]]

    trace$stage4[i + 1L] <-
      trace$stage4[i] * (1 - p_death_stage4) +
      trace$stage3[i] * progression[["stage3"]]

    trace$death[i + 1L] <- n_patients - sum(unlist(
      trace[i + 1L, c("stage1", "stage2", "stage3", "stage4")],
      use.names = FALSE
    ))
  }

  trace <- calculate_rewards(trace)
  trace <- calculate_indirect_costs(trace)
  trace
}

# =============================================================================
# 7. Patient-level microsimulation
# =============================================================================
run_microsimulation <- function(
    n_patients, progression, mortality, ae_mortality, seed) {
  set.seed(seed)

  # State coding: 1 = GOLD 1, 2 = GOLD 2, 3 = GOLD 3,
  # 4 = GOLD 4, 5 = death.
  state <- rep.int(1L, n_patients)

  trace <- data.frame(
    row = seq_len(N_ROWS),
    transition = 0:(N_ROWS - 1L),
    stage1 = integer(N_ROWS),
    stage2 = integer(N_ROWS),
    stage3 = integer(N_ROWS),
    stage4 = integer(N_ROWS),
    death = integer(N_ROWS)
  )

  trace$stage1[1] <- n_patients

  for (i in seq_len(N_ROWS - 1L)) {
    # Snapshot indices prevent an individual from moving more than once during
    # a single quarterly transition.
    index_stage1 <- which(state == 1L)
    index_stage2 <- which(state == 2L)
    index_stage3 <- which(state == 3L)
    index_stage4 <- which(state == 4L)

    if (length(index_stage1) > 0L) {
      random_number <- runif(length(index_stage1))
      progress <- random_number < progression[["stage1"]]
      die <- !progress & random_number < (
        progression[["stage1"]] + mortality$stage1[i]
      )

      state[index_stage1[progress]] <- 2L
      state[index_stage1[die]] <- 5L
    }

    if (length(index_stage2) > 0L) {
      p_death <- mortality$stage2[i] +
        P_SEVERE_AE[["stage2"]] * ae_mortality$stage2[i]

      random_number <- runif(length(index_stage2))
      progress <- random_number < progression[["stage2"]]
      die <- !progress & random_number < (
        progression[["stage2"]] + p_death
      )

      state[index_stage2[progress]] <- 3L
      state[index_stage2[die]] <- 5L
    }

    if (length(index_stage3) > 0L) {
      p_death <- mortality$stage34[i] +
        P_SEVERE_AE[["stage3"]] * ae_mortality$stage34[i]

      random_number <- runif(length(index_stage3))
      progress <- random_number < progression[["stage3"]]
      die <- !progress & random_number < (
        progression[["stage3"]] + p_death
      )

      state[index_stage3[progress]] <- 4L
      state[index_stage3[die]] <- 5L
    }

    if (length(index_stage4) > 0L) {
      p_death <- mortality$stage34[i] +
        P_SEVERE_AE[["stage4"]] * ae_mortality$stage34[i]

      die <- runif(length(index_stage4)) < p_death
      state[index_stage4[die]] <- 5L
    }

    state_counts <- tabulate(state, nbins = 5L)
    trace[i + 1L, c("stage1", "stage2", "stage3", "stage4", "death")] <-
      state_counts
  }

  trace <- calculate_rewards(trace)
  trace <- calculate_indirect_costs(trace)
  trace
}

# =============================================================================
# 8. Strategy aggregation and economic interpretation
# =============================================================================
subcohort_results_per_patient <- function(trace, n_patients) {
  c(
    direct_costs = sum(trace$direct_costs_discounted) / n_patients,
    indirect_costs = sum(trace$indirect_costs_discounted) / n_patients,
    societal_costs = (
      sum(trace$direct_costs_discounted) +
        sum(trace$indirect_costs_discounted)
    ) / n_patients,
    qalys = sum(trace$qalys_discounted) / n_patients
  )
}

evaluate_increment <- function(delta_costs, delta_qalys) {
  tolerance <- 1e-12

  if (abs(delta_qalys) <= tolerance) {
    return(list(status = "No material QALY difference", icer = NA_real_))
  }

  if (delta_qalys > 0 && delta_costs < 0) {
    return(list(status = "Dominant", icer = NA_real_))
  }

  if (delta_qalys < 0 && delta_costs > 0) {
    return(list(status = "Dominated", icer = NA_real_))
  }

  status <- if (delta_qalys > 0 && delta_costs > 0) {
    "More effective and more costly"
  } else {
    "Less effective and less costly"
  }

  list(status = status, icer = delta_costs / delta_qalys)
}

calculate_strategy_results <- function(
    smoker_trace, former_smoker_trace, n_patients, model_name) {
  smoker <- subcohort_results_per_patient(smoker_trace, n_patients)
  former_smoker <- subcohort_results_per_patient(
    former_smoker_trace, n_patients
  )

  mix_subcohorts <- function(prop_former_smoker, additional_direct_cost) {
    result <-
      (1 - prop_former_smoker) * smoker +
      prop_former_smoker * former_smoker

    result[["direct_costs"]] <-
      result[["direct_costs"]] + additional_direct_cost
    result[["societal_costs"]] <-
      result[["direct_costs"]] + result[["indirect_costs"]]

    result
  }

  usual_care <- mix_subcohorts(
    PROP_FORMER_SMOKER_USUAL_CARE,
    USUAL_CARE_COST_PER_PATIENT
  )

  intervention <- mix_subcohorts(
    PROP_FORMER_SMOKER_INTERVENTION,
    INTERVENTION_COST_PER_PATIENT
  )

  incremental <- intervention - usual_care

  direct_evaluation <- evaluate_increment(
    incremental[["direct_costs"]], incremental[["qalys"]]
  )
  societal_evaluation <- evaluate_increment(
    incremental[["societal_costs"]], incremental[["qalys"]]
  )

  output <- data.frame(
    model = model_name,
    strategy = c("usual_care", "intervention", "incremental"),
    direct_costs_per_patient = c(
      usual_care[["direct_costs"]],
      intervention[["direct_costs"]],
      incremental[["direct_costs"]]
    ),
    indirect_costs_per_patient = c(
      usual_care[["indirect_costs"]],
      intervention[["indirect_costs"]],
      incremental[["indirect_costs"]]
    ),
    societal_costs_per_patient = c(
      usual_care[["societal_costs"]],
      intervention[["societal_costs"]],
      incremental[["societal_costs"]]
    ),
    qalys_per_patient = c(
      usual_care[["qalys"]],
      intervention[["qalys"]],
      incremental[["qalys"]]
    ),
    direct_result = c(NA_character_, NA_character_, direct_evaluation$status),
    direct_icer = c(NA_real_, NA_real_, direct_evaluation$icer),
    societal_result = c(
      NA_character_, NA_character_, societal_evaluation$status
    ),
    societal_icer = c(NA_real_, NA_real_, societal_evaluation$icer),
    stringsAsFactors = FALSE
  )

  output
}

calculate_subcohort_table <- function(
    smoker_trace, former_smoker_trace, n_patients, model_name) {
  smoker <- subcohort_results_per_patient(smoker_trace, n_patients)
  former_smoker <- subcohort_results_per_patient(
    former_smoker_trace, n_patients
  )

  data.frame(
    model = model_name,
    subcohort = c("current_smokers", "former_smokers"),
    direct_costs_per_patient = c(
      smoker[["direct_costs"]], former_smoker[["direct_costs"]]
    ),
    indirect_costs_per_patient = c(
      smoker[["indirect_costs"]], former_smoker[["indirect_costs"]]
    ),
    societal_costs_per_patient = c(
      smoker[["societal_costs"]], former_smoker[["societal_costs"]]
    ),
    qalys_per_patient = c(smoker[["qalys"]], former_smoker[["qalys"]]),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 9. Run models
# =============================================================================
check_probabilities(
  PROGRESSION_SMOKER,
  MORTALITY_SMOKER,
  AE_MORTALITY_SMOKER,
  "current smokers"
)
check_probabilities(
  PROGRESSION_FORMER_SMOKER,
  MORTALITY_FORMER_SMOKER,
  AE_MORTALITY_FORMER_SMOKER,
  "former smokers"
)

cat("Running deterministic cohort Markov model...\n")
markov_smokers <- run_markov_model(
  N_PATIENTS,
  PROGRESSION_SMOKER,
  MORTALITY_SMOKER,
  AE_MORTALITY_SMOKER
)
markov_former_smokers <- run_markov_model(
  N_PATIENTS,
  PROGRESSION_FORMER_SMOKER,
  MORTALITY_FORMER_SMOKER,
  AE_MORTALITY_FORMER_SMOKER
)

cat("Running patient-level microsimulation...\n")
microsimulation_smokers <- run_microsimulation(
  N_PATIENTS,
  PROGRESSION_SMOKER,
  MORTALITY_SMOKER,
  AE_MORTALITY_SMOKER,
  SEED_SMOKERS
)
microsimulation_former_smokers <- run_microsimulation(
  N_PATIENTS,
  PROGRESSION_FORMER_SMOKER,
  MORTALITY_FORMER_SMOKER,
  AE_MORTALITY_FORMER_SMOKER,
  SEED_FORMER_SMOKERS
)

check_trace(
  markov_smokers, N_PATIENTS, "Markov model: current smokers", FALSE
)
check_trace(
  markov_former_smokers, N_PATIENTS,
  "Markov model: former smokers", FALSE
)
check_trace(
  microsimulation_smokers, N_PATIENTS,
  "Microsimulation: current smokers", TRUE
)
check_trace(
  microsimulation_former_smokers, N_PATIENTS,
  "Microsimulation: former smokers", TRUE
)

# =============================================================================
# 10. Calculate results
# =============================================================================
subcohort_results <- rbind(
  calculate_subcohort_table(
    markov_smokers, markov_former_smokers, N_PATIENTS, "Markov"
  ),
  calculate_subcohort_table(
    microsimulation_smokers, microsimulation_former_smokers,
    N_PATIENTS, "Microsimulation"
  )
)

strategy_results <- rbind(
  calculate_strategy_results(
    markov_smokers, markov_former_smokers, N_PATIENTS, "Markov"
  ),
  calculate_strategy_results(
    microsimulation_smokers, microsimulation_former_smokers,
    N_PATIENTS, "Microsimulation"
  )
)

# One-run comparison between the deterministic expectation and the selected
# microsimulation seed. Repeated-seed validation remains archived separately.
markov_strategy <- strategy_results[
  strategy_results$model == "Markov",
  ,
  drop = FALSE
]
microsimulation_strategy <- strategy_results[
  strategy_results$model == "Microsimulation",
  ,
  drop = FALSE
]

model_comparison <- merge(
  markov_strategy,
  microsimulation_strategy,
  by = "strategy",
  suffixes = c("_markov", "_microsimulation"),
  sort = FALSE
)

model_comparison$direct_cost_difference <-
  model_comparison$direct_costs_per_patient_microsimulation -
  model_comparison$direct_costs_per_patient_markov
model_comparison$qaly_difference <-
  model_comparison$qalys_per_patient_microsimulation -
  model_comparison$qalys_per_patient_markov

# =============================================================================
# 11. Regression checks against the frozen deterministic benchmark
# =============================================================================
EXPECTED_BENCHMARK <- c(
  direct_costs_usual_care = 11500.78,
  direct_costs_intervention = 11711.26,
  qalys_usual_care = 12.68995,
  qalys_intervention = 13.23588,
  direct_icer = 385.5405
)

markov_incremental <- strategy_results[
  strategy_results$model == "Markov" &
    strategy_results$strategy == "incremental",
  ,
  drop = FALSE
]
markov_usual_care <- strategy_results[
  strategy_results$model == "Markov" &
    strategy_results$strategy == "usual_care",
  ,
  drop = FALSE
]
markov_intervention <- strategy_results[
  strategy_results$model == "Markov" &
    strategy_results$strategy == "intervention",
  ,
  drop = FALSE
]

CURRENT_BENCHMARK <- c(
  direct_costs_usual_care = markov_usual_care$direct_costs_per_patient,
  direct_costs_intervention = markov_intervention$direct_costs_per_patient,
  qalys_usual_care = markov_usual_care$qalys_per_patient,
  qalys_intervention = markov_intervention$qalys_per_patient,
  direct_icer = markov_incremental$direct_icer
)

BENCHMARK_TOLERANCE <- c(
  direct_costs_usual_care = 0.02,
  direct_costs_intervention = 0.02,
  qalys_usual_care = 0.00002,
  qalys_intervention = 0.00002,
  direct_icer = 0.02
)

benchmark_check <- data.frame(
  result = names(EXPECTED_BENCHMARK),
  expected = unname(EXPECTED_BENCHMARK),
  current = unname(CURRENT_BENCHMARK),
  difference = unname(CURRENT_BENCHMARK - EXPECTED_BENCHMARK),
  tolerance = unname(BENCHMARK_TOLERANCE),
  stringsAsFactors = FALSE
)
benchmark_check$passed <-
  abs(benchmark_check$difference) <= benchmark_check$tolerance

if (any(!benchmark_check$passed)) {
  stop(
    "Regression check failed: the deterministic model no longer reproduces ",
    "the frozen simplified benchmark. Inspect benchmark_check.csv."
  )
}

# =============================================================================
# 12. Save outputs
# =============================================================================
write.csv(
  subcohort_results,
  file.path(OUTPUT_DIR, "subcohort_results.csv"),
  row.names = FALSE
)
write.csv(
  strategy_results,
  file.path(OUTPUT_DIR, "strategy_results.csv"),
  row.names = FALSE
)
write.csv(
  model_comparison,
  file.path(OUTPUT_DIR, "markov_vs_microsimulation_one_run.csv"),
  row.names = FALSE
)
write.csv(
  benchmark_check,
  file.path(OUTPUT_DIR, "benchmark_check.csv"),
  row.names = FALSE
)

if (SAVE_STATE_TRACES) {
  write.csv(
    markov_smokers,
    file.path(OUTPUT_DIR, "trace_markov_current_smokers.csv"),
    row.names = FALSE
  )
  write.csv(
    markov_former_smokers,
    file.path(OUTPUT_DIR, "trace_markov_former_smokers.csv"),
    row.names = FALSE
  )
  write.csv(
    microsimulation_smokers,
    file.path(OUTPUT_DIR, "trace_microsimulation_current_smokers.csv"),
    row.names = FALSE
  )
  write.csv(
    microsimulation_former_smokers,
    file.path(OUTPUT_DIR, "trace_microsimulation_former_smokers.csv"),
    row.names = FALSE
  )
}

plot_state_traces <- function(markov_trace, microsimulation_trace, title_prefix,
                              filename) {
  png(filename, width = 1400, height = 1000, res = 150)
  on.exit(dev.off(), add = TRUE)

  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  for (state_name in c("stage1", "stage2", "stage3", "stage4")) {
    plot(
      markov_trace$transition,
      markov_trace[[state_name]],
      type = "l",
      xlab = "Quarterly transition",
      ylab = "Patients",
      main = paste(title_prefix, state_name)
    )
    lines(
      microsimulation_trace$transition,
      microsimulation_trace[[state_name]],
      lty = 2
    )
  }
}

if (SAVE_TRACE_PLOTS) {
  plot_state_traces(
    markov_smokers,
    microsimulation_smokers,
    "Current smokers:",
    file.path(OUTPUT_DIR, "state_trace_current_smokers.png")
  )
  plot_state_traces(
    markov_former_smokers,
    microsimulation_former_smokers,
    "Former smokers:",
    file.path(OUTPUT_DIR, "state_trace_former_smokers.png")
  )
}

model_results <- list(
  settings = list(
    n_patients_per_subcohort = N_PATIENTS,
    seed_current_smokers = SEED_SMOKERS,
    seed_former_smokers = SEED_FORMER_SMOKERS,
    n_rows = N_ROWS,
    quarterly_transitions = N_ROWS - 1L,
    cycle_length_months = CYCLE_LENGTH_MONTHS,
    annual_discount_rate = ANNUAL_DISCOUNT_RATE
  ),
  markov = list(
    current_smokers = markov_smokers,
    former_smokers = markov_former_smokers
  ),
  microsimulation = list(
    current_smokers = microsimulation_smokers,
    former_smokers = microsimulation_former_smokers
  ),
  subcohort_results = subcohort_results,
  strategy_results = strategy_results,
  model_comparison = model_comparison,
  benchmark_check = benchmark_check
)

saveRDS(
  model_results,
  file.path(OUTPUT_DIR, "model_results_baseline.rds")
)

# =============================================================================
# 13. Human-readable report
# =============================================================================
report_lines <- capture.output({
  cat("============================================================\n")
  cat("COPD BASELINE MODEL RESULTS\n")
  cat("============================================================\n")
  cat("Patients per smoking-status subcohort:", N_PATIENTS, "\n")
  cat("Quarterly transitions:", N_ROWS - 1L, "\n")
  cat("Time horizon:", (N_ROWS - 1L) * CYCLE_LENGTH_MONTHS / 12,
      "years\n")
  cat("Discount rate:", 100 * ANNUAL_DISCOUNT_RATE, "% per year\n\n")

  cat("1. Regression check against frozen Markov benchmark\n")
  print(benchmark_check, row.names = FALSE, digits = 10)

  cat("\n2. Results by smoking-status subcohort\n")
  print(subcohort_results, row.names = FALSE, digits = 10)

  cat("\n3. Strategy-level results\n")
  print(strategy_results, row.names = FALSE, digits = 10)

  cat("\n4. Interpretation\n")
  cat(
    "- Direct-cost perspective: report the numeric ICER when the intervention\n",
    "  is both more effective and more costly.\n",
    sep = ""
  )
  cat(
    "- Societal perspective: report 'Dominant' when total societal costs are\n",
    "  lower and QALYs are higher; do not report a negative ICER.\n",
    sep = ""
  )
  cat(
    "- The selected microsimulation seed is one stochastic realization. The\n",
    "  archived repeated-seed validation establishes agreement in expectation.\n",
    sep = ""
  )

  cat("\nOutput directory:\n", normalizePath(OUTPUT_DIR), "\n", sep = "")
})

writeLines(
  report_lines,
  con = file.path(OUTPUT_DIR, "model_run_report.txt")
)
cat(paste(report_lines, collapse = "\n"), "\n")

###############################################################################
# End of file
###############################################################################
