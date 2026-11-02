###############################################################################
# COPD patient-level microsimulation: acute-exacerbation memory extension
#
# This file is a separate extension of the frozen, technically validated
# baseline model. It does not alter the baseline implementation.
#
# Required file in the current working directory:
#   COPD_Markov_microsimulation_main_baseline.R
#
# The extension compares:
#   1. The frozen expected-value baseline microsimulation.
#   2. An event-level AE microsimulation without memory.
#   3. The same event-level AE microsimulation with patient-level memory.
#
# The memory scenario is illustrative and assumption-based, not empirically
# calibrated. Base R only.
###############################################################################

# =============================================================================
# 1. Load the frozen baseline model
# =============================================================================
BASELINE_FILENAME <- "COPD_Markov_microsimulation_main_baseline.R"

if (!file.exists(BASELINE_FILENAME)) {
  stop(
    "Required baseline file not found in the current working directory: ",
    BASELINE_FILENAME
  )
}

# The baseline file intentionally clears the global environment and then runs
# the validated benchmark and baseline microsimulation. Therefore, extension
# settings are defined only after sourcing it.
source(BASELINE_FILENAME)

# =============================================================================
# 2. Memory-extension settings
# =============================================================================
OUTPUT_DIR_MEMORY <- file.path(getwd(), "model_output_AE_memory")
dir.create(OUTPUT_DIR_MEMORY, showWarnings = FALSE, recursive = TRUE)

# Use the same patient numbers and subcohort-specific seeds as the frozen
# baseline. The event-level no-memory and memory models use common random
# numbers, which reduces random noise in their direct comparison.
EVENT_SEED_SMOKERS <- SEED_SMOKERS + 20000L
EVENT_SEED_FORMER_SMOKERS <- SEED_FORMER_SMOKERS + 20000L

# Illustrative memory assumption from the manuscript draft:
# after any moderate or severe AE, the combined probability of a subsequent AE
# is increased by at most 50%. The excess risk applies in the next three
# quarterly cycles and decays linearly:
#   first subsequent cycle:  +50.0%
#   second subsequent cycle: +33.3%
#   third subsequent cycle:  +16.7%
#   thereafter:               +0.0%
MEMORY_MAX_RELATIVE_RISK_INCREASE <- 0.50
MEMORY_DECAY_WEIGHTS <- c(1, 2 / 3, 1 / 3)
MEMORY_DURATION_CYCLES <- length(MEMORY_DECAY_WEIGHTS)

SAVE_MEMORY_TRACES <- TRUE
SAVE_MEMORY_PLOTS <- TRUE

# =============================================================================
# 3. Helper functions for event-level rewards
# =============================================================================
get_indirect_cost_schedules <- function() {
  list(
    stable = list(
      stage2 = c(
        0, rep(112.50, 19), rep(376.37, 20), rep(667.24, 20),
        rep(237.22, 20), rep(0, 161)
      ),
      stage3 = c(
        0, rep(268.06, 19), rep(896.81, 20), rep(1589.90, 20),
        rep(565.24, 20), rep(0, 161)
      ),
      stage4 = c(
        0, rep(281.24, 19), rep(940.92, 20), rep(1668.10, 20),
        rep(593.04, 20), rep(0, 161)
      )
    ),
    moderate = list(
      stage2 = c(
        rep(491.94, 20), rep(446.16, 20), rep(364.35, 20),
        rep(115.30, 20), rep(0, 161)
      ),
      stage3 = c(
        rep(550.51, 20), rep(499.27, 20), rep(407.73, 20),
        rep(129.03, 20), rep(0, 161)
      ),
      stage4 = c(
        rep(491.94, 20), rep(446.16, 20), rep(364.35, 20),
        rep(115.30, 20), rep(0, 161)
      )
    ),
    severe = list(
      stage2 = c(
        rep(367.00, 20), rep(332.85, 20), rep(271.82, 20),
        rep(86.02, 20), rep(0, 161)
      ),
      stage3 = c(
        rep(663.73, 20), rep(601.96, 20), rep(491.58, 20),
        rep(155.56, 20), rep(0, 161)
      ),
      stage4 = c(
        rep(648.12, 20), rep(587.79, 20), rep(480.02, 20),
        rep(151.90, 20), rep(0, 161)
      )
    )
  )
}

INDIRECT_COST_SCHEDULES <- get_indirect_cost_schedules()

stopifnot(
  all(vapply(
    INDIRECT_COST_SCHEDULES$stable,
    length,
    integer(1)
  ) == N_ROWS),
  all(vapply(
    INDIRECT_COST_SCHEDULES$moderate,
    length,
    integer(1)
  ) == N_ROWS),
  all(vapply(
    INDIRECT_COST_SCHEDULES$severe,
    length,
    integer(1)
  ) == N_ROWS)
)

get_memory_multiplier <- function(
    patient_indices, current_row, last_ae_row, include_memory) {
  multiplier <- rep.int(1, length(patient_indices))

  if (!include_memory || length(patient_indices) == 0L) {
    return(multiplier)
  }

  previous_ae_row <- last_ae_row[patient_indices]
  lag <- current_row - previous_ae_row

  active <- !is.na(previous_ae_row) &
    lag >= 1L & lag <= MEMORY_DURATION_CYCLES

  if (any(active)) {
    multiplier[active] <- 1 +
      MEMORY_MAX_RELATIVE_RISK_INCREASE *
      MEMORY_DECAY_WEIGHTS[lag[active]]
  }

  multiplier
}

cap_competing_ae_probabilities <- function(p_moderate, p_severe) {
  total <- p_moderate + p_severe
  above_one <- total > 1

  if (any(above_one)) {
    scaling <- (1 - .Machine$double.eps) / total[above_one]
    p_moderate[above_one] <- p_moderate[above_one] * scaling
    p_severe[above_one] <- p_severe[above_one] * scaling
  }

  list(moderate = p_moderate, severe = p_severe)
}

calculate_event_level_rewards <- function(
    trace, row_index, state_counts, moderate_counts, severe_counts) {
  stage_names <- c("stage1", "stage2", "stage3", "stage4")

  stable_direct_cost <- sum(
    state_counts[stage_names] * STABLE_DIRECT_COST[stage_names]
  )

  ae_direct_cost <-
    moderate_counts[["stage2"]] * DIRECT_COST_MODERATE_AE[["stage2"]] +
    moderate_counts[["stage3"]] * DIRECT_COST_MODERATE_AE[["stage3"]] +
    moderate_counts[["stage4"]] * DIRECT_COST_MODERATE_AE[["stage4"]] +
    severe_counts[["stage2"]] * DIRECT_COST_SEVERE_AE[["stage2"]] +
    severe_counts[["stage3"]] * DIRECT_COST_SEVERE_AE[["stage3"]] +
    severe_counts[["stage4"]] * DIRECT_COST_SEVERE_AE[["stage4"]]

  stable_utility <- sum(state_counts[stage_names] * UTILITY[stage_names])

  ae_disutility <-
    moderate_counts[["stage2"]] * DISUTILITY_MODERATE_AE[["stage2"]] +
    moderate_counts[["stage3"]] * DISUTILITY_MODERATE_AE[["stage3"]] +
    moderate_counts[["stage4"]] * DISUTILITY_MODERATE_AE[["stage4"]] +
    severe_counts[["stage2"]] * DISUTILITY_SEVERE_AE[["stage2"]] +
    severe_counts[["stage3"]] * DISUTILITY_SEVERE_AE[["stage3"]] +
    severe_counts[["stage4"]] * DISUTILITY_SEVERE_AE[["stage4"]]

  stable_indirect_cost <-
    state_counts[["stage2"]] *
      INDIRECT_COST_SCHEDULES$stable$stage2[row_index] +
    state_counts[["stage3"]] *
      INDIRECT_COST_SCHEDULES$stable$stage3[row_index] +
    state_counts[["stage4"]] *
      INDIRECT_COST_SCHEDULES$stable$stage4[row_index]

  ae_indirect_cost <-
    moderate_counts[["stage2"]] *
      INDIRECT_COST_SCHEDULES$moderate$stage2[row_index] +
    moderate_counts[["stage3"]] *
      INDIRECT_COST_SCHEDULES$moderate$stage3[row_index] +
    moderate_counts[["stage4"]] *
      INDIRECT_COST_SCHEDULES$moderate$stage4[row_index] +
    severe_counts[["stage2"]] *
      INDIRECT_COST_SCHEDULES$severe$stage2[row_index] +
    severe_counts[["stage3"]] *
      INDIRECT_COST_SCHEDULES$severe$stage3[row_index] +
    severe_counts[["stage4"]] *
      INDIRECT_COST_SCHEDULES$severe$stage4[row_index]

  trace$direct_costs[row_index] <- stable_direct_cost + ae_direct_cost
  trace$indirect_costs[row_index] <-
    stable_indirect_cost + ae_indirect_cost
  trace$qalys[row_index] <-
    (stable_utility + ae_disutility) * CYCLE_LENGTH_MONTHS / 12

  trace
}

# =============================================================================
# 4. Event-level microsimulation with optional AE memory
# =============================================================================
run_event_level_microsimulation <- function(
    n_patients,
    progression,
    mortality,
    ae_mortality,
    seed,
    include_memory = FALSE) {
  set.seed(seed)

  # State coding: 1 = GOLD 1, 2 = GOLD 2, 3 = GOLD 3,
  # 4 = GOLD 4, 5 = death.
  state <- rep.int(1L, n_patients)

  # A severe AE sampled for the current trace row affects mortality during the
  # subsequent transition. This matches the reward/transition timing of the
  # frozen source model: events and rewards belong to the current trace row.
  severe_ae_current <- rep.int(FALSE, n_patients)

  # Individual history variables.
  last_ae_row <- rep.int(NA_integer_, n_patients)
  cumulative_moderate_ae <- integer(n_patients)
  cumulative_severe_ae <- integer(n_patients)

  trace <- data.frame(
    row = seq_len(N_ROWS),
    transition = 0:(N_ROWS - 1L),
    stage1 = integer(N_ROWS),
    stage2 = integer(N_ROWS),
    stage3 = integer(N_ROWS),
    stage4 = integer(N_ROWS),
    death = integer(N_ROWS),
    moderate_ae_stage2 = integer(N_ROWS),
    moderate_ae_stage3 = integer(N_ROWS),
    moderate_ae_stage4 = integer(N_ROWS),
    severe_ae_stage2 = integer(N_ROWS),
    severe_ae_stage3 = integer(N_ROWS),
    severe_ae_stage4 = integer(N_ROWS),
    total_moderate_ae = integer(N_ROWS),
    total_severe_ae = integer(N_ROWS),
    active_memory_patients = integer(N_ROWS),
    mean_risk_multiplier = numeric(N_ROWS),
    direct_costs = numeric(N_ROWS),
    indirect_costs = numeric(N_ROWS),
    qalys = numeric(N_ROWS)
  )

  trace$stage1[1] <- n_patients
  trace$mean_risk_multiplier[1] <- 1

  for (i in seq_len(N_ROWS - 1L)) {
    # Fixed-size patient-level random-number vectors preserve common random
    # numbers between the no-memory and memory variants, even after their state
    # distributions begin to diverge.
    random_transition <- runif(n_patients)
    random_ae <- runif(n_patients)

    index_stage1 <- which(state == 1L)
    index_stage2 <- which(state == 2L)
    index_stage3 <- which(state == 3L)
    index_stage4 <- which(state == 4L)

    if (length(index_stage1) > 0L) {
      u <- random_transition[index_stage1]
      progress <- u < progression[["stage1"]]
      die <- !progress & u < (
        progression[["stage1"]] + mortality$stage1[i]
      )

      state[index_stage1[progress]] <- 2L
      state[index_stage1[die]] <- 5L
    }

    if (length(index_stage2) > 0L) {
      p_death <- mortality$stage2[i] +
        as.numeric(severe_ae_current[index_stage2]) *
        ae_mortality$stage2[i]

      if (any(progression[["stage2"]] + p_death > 1)) {
        stop("Invalid individual transition probability in GOLD 2.")
      }

      u <- random_transition[index_stage2]
      progress <- u < progression[["stage2"]]
      die <- !progress & u < (progression[["stage2"]] + p_death)

      state[index_stage2[progress]] <- 3L
      state[index_stage2[die]] <- 5L
    }

    if (length(index_stage3) > 0L) {
      p_death <- mortality$stage34[i] +
        as.numeric(severe_ae_current[index_stage3]) *
        ae_mortality$stage34[i]

      if (any(progression[["stage3"]] + p_death > 1)) {
        stop("Invalid individual transition probability in GOLD 3.")
      }

      u <- random_transition[index_stage3]
      progress <- u < progression[["stage3"]]
      die <- !progress & u < (progression[["stage3"]] + p_death)

      state[index_stage3[progress]] <- 4L
      state[index_stage3[die]] <- 5L
    }

    if (length(index_stage4) > 0L) {
      p_death <- mortality$stage34[i] +
        as.numeric(severe_ae_current[index_stage4]) *
        ae_mortality$stage34[i]

      if (any(p_death > 1)) {
        stop("Invalid individual transition probability in GOLD 4.")
      }

      die <- random_transition[index_stage4] < p_death
      state[index_stage4[die]] <- 5L
    }

    current_row <- i + 1L
    state_count_vector <- tabulate(state, nbins = 5L)
    names(state_count_vector) <- c(
      "stage1", "stage2", "stage3", "stage4", "death"
    )

    trace[current_row, c(
      "stage1", "stage2", "stage3", "stage4", "death"
    )] <- state_count_vector

    # Sample at most one AE per person and quarter. Moderate and severe AEs are
    # mutually exclusive categorical outcomes. Their marginal probabilities
    # equal the probabilities in the frozen expected-value model when memory is
    # disabled.
    moderate_ae_next <- rep.int(FALSE, n_patients)
    severe_ae_next <- rep.int(FALSE, n_patients)

    multiplier_sum <- 0
    multiplier_n <- 0L
    active_memory_n <- 0L

    for (stage_number in 2:4) {
      stage_name <- paste0("stage", stage_number)
      patient_indices <- which(state == stage_number)

      if (length(patient_indices) == 0L) {
        next
      }

      multiplier <- get_memory_multiplier(
        patient_indices,
        current_row,
        last_ae_row,
        include_memory
      )

      p_moderate <-
        P_MODERATE_AE[[stage_name]] * multiplier
      p_severe <-
        P_SEVERE_AE[[stage_name]] * multiplier

      capped <- cap_competing_ae_probabilities(p_moderate, p_severe)
      p_moderate <- capped$moderate
      p_severe <- capped$severe

      u_ae <- random_ae[patient_indices]
      severe_event <- u_ae < p_severe
      moderate_event <- !severe_event &
        u_ae < (p_severe + p_moderate)

      severe_ae_next[patient_indices[severe_event]] <- TRUE
      moderate_ae_next[patient_indices[moderate_event]] <- TRUE

      trace[current_row, paste0("moderate_ae_", stage_name)] <-
        sum(moderate_event)
      trace[current_row, paste0("severe_ae_", stage_name)] <-
        sum(severe_event)

      multiplier_sum <- multiplier_sum + sum(multiplier)
      multiplier_n <- multiplier_n + length(multiplier)
      active_memory_n <- active_memory_n + sum(multiplier > 1)
    }

    any_ae_next <- moderate_ae_next | severe_ae_next
    last_ae_row[any_ae_next] <- current_row
    cumulative_moderate_ae <- cumulative_moderate_ae +
      as.integer(moderate_ae_next)
    cumulative_severe_ae <- cumulative_severe_ae +
      as.integer(severe_ae_next)

    trace$total_moderate_ae[current_row] <- sum(moderate_ae_next)
    trace$total_severe_ae[current_row] <- sum(severe_ae_next)
    trace$active_memory_patients[current_row] <- active_memory_n
    trace$mean_risk_multiplier[current_row] <- if (multiplier_n > 0L) {
      multiplier_sum / multiplier_n
    } else {
      1
    }

    moderate_counts <- c(
      stage2 = trace$moderate_ae_stage2[current_row],
      stage3 = trace$moderate_ae_stage3[current_row],
      stage4 = trace$moderate_ae_stage4[current_row]
    )
    severe_counts <- c(
      stage2 = trace$severe_ae_stage2[current_row],
      stage3 = trace$severe_ae_stage3[current_row],
      stage4 = trace$severe_ae_stage4[current_row]
    )

    trace <- calculate_event_level_rewards(
      trace,
      current_row,
      state_count_vector,
      moderate_counts,
      severe_counts
    )

    severe_ae_current <- severe_ae_next
  }

  # Baseline row contains no cycle reward or event.
  trace$direct_costs[1] <- 0
  trace$indirect_costs[1] <- 0
  trace$qalys[1] <- 0

  discount_factor <- get_discount_factor(N_ROWS)
  trace$direct_costs_discounted <- trace$direct_costs * discount_factor
  trace$indirect_costs_discounted <- trace$indirect_costs * discount_factor
  trace$qalys_discounted <- trace$qalys * discount_factor

  patient_total_ae <- cumulative_moderate_ae + cumulative_severe_ae

  history_summary <- data.frame(
    include_memory = include_memory,
    total_moderate_ae = sum(cumulative_moderate_ae),
    total_severe_ae = sum(cumulative_severe_ae),
    total_all_ae = sum(patient_total_ae),
    patients_with_any_ae = sum(patient_total_ae > 0L),
    proportion_with_any_ae = mean(patient_total_ae > 0L),
    mean_ae_per_patient = mean(patient_total_ae),
    median_ae_per_patient = median(patient_total_ae),
    maximum_ae_per_patient = max(patient_total_ae),
    stringsAsFactors = FALSE
  )

  check_trace(
    trace,
    n_patients,
    if (include_memory) {
      "Event-level microsimulation with AE memory"
    } else {
      "Event-level microsimulation without AE memory"
    },
    TRUE
  )

  list(trace = trace, history_summary = history_summary)
}

# =============================================================================
# 5. Internal assumption and probability checks
# =============================================================================
stopifnot(
  MEMORY_MAX_RELATIVE_RISK_INCREASE >= 0,
  length(MEMORY_DECAY_WEIGHTS) >= 1L,
  all(MEMORY_DECAY_WEIGHTS >= 0),
  all(diff(MEMORY_DECAY_WEIGHTS) <= 0)
)

maximum_multiplier <- 1 +
  MEMORY_MAX_RELATIVE_RISK_INCREASE * max(MEMORY_DECAY_WEIGHTS)
maximum_total_ae_probability <- max(
  (P_MODERATE_AE + P_SEVERE_AE) * maximum_multiplier
)

if (maximum_total_ae_probability >= 1) {
  warning(
    "The maximum combined AE probability reaches or exceeds 1. ",
    "Probabilities will be proportionally capped."
  )
}

# =============================================================================
# 6. Run event-level models
# =============================================================================
cat("\nRunning event-level microsimulation without AE memory...\n")
event_no_memory_smokers <- run_event_level_microsimulation(
  N_PATIENTS,
  PROGRESSION_SMOKER,
  MORTALITY_SMOKER,
  AE_MORTALITY_SMOKER,
  EVENT_SEED_SMOKERS,
  include_memory = FALSE
)
event_no_memory_former <- run_event_level_microsimulation(
  N_PATIENTS,
  PROGRESSION_FORMER_SMOKER,
  MORTALITY_FORMER_SMOKER,
  AE_MORTALITY_FORMER_SMOKER,
  EVENT_SEED_FORMER_SMOKERS,
  include_memory = FALSE
)

cat("Running event-level microsimulation with AE memory...\n")
event_memory_smokers <- run_event_level_microsimulation(
  N_PATIENTS,
  PROGRESSION_SMOKER,
  MORTALITY_SMOKER,
  AE_MORTALITY_SMOKER,
  EVENT_SEED_SMOKERS,
  include_memory = TRUE
)
event_memory_former <- run_event_level_microsimulation(
  N_PATIENTS,
  PROGRESSION_FORMER_SMOKER,
  MORTALITY_FORMER_SMOKER,
  AE_MORTALITY_FORMER_SMOKER,
  EVENT_SEED_FORMER_SMOKERS,
  include_memory = TRUE
)

# =============================================================================
# 7. Results and comparisons
# =============================================================================
expected_value_strategy_results <- strategy_results[
  strategy_results$model == "Microsimulation",
  ,
  drop = FALSE
]
expected_value_strategy_results$model <- "Expected-value baseline"

event_no_memory_strategy_results <- calculate_strategy_results(
  event_no_memory_smokers$trace,
  event_no_memory_former$trace,
  N_PATIENTS,
  "Event-level no memory"
)

event_memory_strategy_results <- calculate_strategy_results(
  event_memory_smokers$trace,
  event_memory_former$trace,
  N_PATIENTS,
  "Event-level with memory"
)

memory_strategy_results <- rbind(
  expected_value_strategy_results,
  event_no_memory_strategy_results,
  event_memory_strategy_results
)

memory_subcohort_results <- rbind(
  calculate_subcohort_table(
    event_no_memory_smokers$trace,
    event_no_memory_former$trace,
    N_PATIENTS,
    "Event-level no memory"
  ),
  calculate_subcohort_table(
    event_memory_smokers$trace,
    event_memory_former$trace,
    N_PATIENTS,
    "Event-level with memory"
  )
)

compare_model_specifications <- function(reference, alternative) {
  merged <- merge(
    reference,
    alternative,
    by = "strategy",
    suffixes = c("_no_memory", "_memory"),
    sort = FALSE
  )

  merged$direct_cost_difference <-
    merged$direct_costs_per_patient_memory -
    merged$direct_costs_per_patient_no_memory
  merged$indirect_cost_difference <-
    merged$indirect_costs_per_patient_memory -
    merged$indirect_costs_per_patient_no_memory
  merged$societal_cost_difference <-
    merged$societal_costs_per_patient_memory -
    merged$societal_costs_per_patient_no_memory
  merged$qaly_difference <-
    merged$qalys_per_patient_memory -
    merged$qalys_per_patient_no_memory

  merged$direct_cost_relative_difference_percent <-
    100 * merged$direct_cost_difference /
    merged$direct_costs_per_patient_no_memory
  merged$qaly_relative_difference_percent <-
    100 * merged$qaly_difference /
    merged$qalys_per_patient_no_memory

  merged
}

memory_effect_comparison <- compare_model_specifications(
  event_no_memory_strategy_results,
  event_memory_strategy_results
)

compare_event_to_expected <- function(event_results, expected_results) {
  merged <- merge(
    expected_results,
    event_results,
    by = "strategy",
    suffixes = c("_expected", "_event"),
    sort = FALSE
  )

  merged$direct_cost_difference <-
    merged$direct_costs_per_patient_event -
    merged$direct_costs_per_patient_expected
  merged$qaly_difference <-
    merged$qalys_per_patient_event -
    merged$qalys_per_patient_expected
  merged$direct_icer_difference <-
    merged$direct_icer_event - merged$direct_icer_expected

  merged
}

event_no_memory_vs_expected <- compare_event_to_expected(
  event_no_memory_strategy_results,
  expected_value_strategy_results
)

history_summary <- rbind(
  cbind(
    subcohort = "current_smokers",
    event_no_memory_smokers$history_summary
  ),
  cbind(
    subcohort = "former_smokers",
    event_no_memory_former$history_summary
  ),
  cbind(
    subcohort = "current_smokers",
    event_memory_smokers$history_summary
  ),
  cbind(
    subcohort = "former_smokers",
    event_memory_former$history_summary
  )
)
history_summary$model <- ifelse(
  history_summary$include_memory,
  "Event-level with memory",
  "Event-level no memory"
)
history_summary <- history_summary[c(
  "model", "subcohort", "total_moderate_ae", "total_severe_ae",
  "total_all_ae", "patients_with_any_ae", "proportion_with_any_ae",
  "mean_ae_per_patient", "median_ae_per_patient",
  "maximum_ae_per_patient"
)]

# =============================================================================
# 8. Save outputs
# =============================================================================
write.csv(
  memory_strategy_results,
  file.path(OUTPUT_DIR_MEMORY, "memory_strategy_results.csv"),
  row.names = FALSE
)
write.csv(
  memory_subcohort_results,
  file.path(OUTPUT_DIR_MEMORY, "memory_subcohort_results.csv"),
  row.names = FALSE
)
write.csv(
  memory_effect_comparison,
  file.path(OUTPUT_DIR_MEMORY, "memory_effect_comparison.csv"),
  row.names = FALSE
)
write.csv(
  event_no_memory_vs_expected,
  file.path(OUTPUT_DIR_MEMORY, "event_no_memory_vs_expected_baseline.csv"),
  row.names = FALSE
)
write.csv(
  history_summary,
  file.path(OUTPUT_DIR_MEMORY, "AE_history_summary.csv"),
  row.names = FALSE
)

if (SAVE_MEMORY_TRACES) {
  write.csv(
    event_no_memory_smokers$trace,
    file.path(OUTPUT_DIR_MEMORY, "trace_event_no_memory_smokers.csv"),
    row.names = FALSE
  )
  write.csv(
    event_memory_smokers$trace,
    file.path(OUTPUT_DIR_MEMORY, "trace_event_memory_smokers.csv"),
    row.names = FALSE
  )
  write.csv(
    event_no_memory_former$trace,
    file.path(OUTPUT_DIR_MEMORY, "trace_event_no_memory_former_smokers.csv"),
    row.names = FALSE
  )
  write.csv(
    event_memory_former$trace,
    file.path(OUTPUT_DIR_MEMORY, "trace_event_memory_former_smokers.csv"),
    row.names = FALSE
  )
}

plot_event_comparison <- function(
    no_memory_trace, memory_trace, title_text, filename) {
  png(filename, width = 1400, height = 1000, res = 150)
  on.exit(dev.off(), add = TRUE)

  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  plot(
    no_memory_trace$transition,
    cumsum(no_memory_trace$total_moderate_ae),
    type = "l",
    xlab = "Quarterly transition",
    ylab = "Cumulative moderate AEs",
    main = paste(title_text, "moderate AEs")
  )
  lines(
    memory_trace$transition,
    cumsum(memory_trace$total_moderate_ae),
    lty = 2
  )

  plot(
    no_memory_trace$transition,
    cumsum(no_memory_trace$total_severe_ae),
    type = "l",
    xlab = "Quarterly transition",
    ylab = "Cumulative severe AEs",
    main = paste(title_text, "severe AEs")
  )
  lines(
    memory_trace$transition,
    cumsum(memory_trace$total_severe_ae),
    lty = 2
  )

  plot(
    no_memory_trace$transition,
    cumsum(no_memory_trace$direct_costs_discounted),
    type = "l",
    xlab = "Quarterly transition",
    ylab = "Cumulative discounted direct costs",
    main = paste(title_text, "direct costs")
  )
  lines(
    memory_trace$transition,
    cumsum(memory_trace$direct_costs_discounted),
    lty = 2
  )

  plot(
    no_memory_trace$transition,
    cumsum(no_memory_trace$qalys_discounted),
    type = "l",
    xlab = "Quarterly transition",
    ylab = "Cumulative discounted QALYs",
    main = paste(title_text, "QALYs")
  )
  lines(
    memory_trace$transition,
    cumsum(memory_trace$qalys_discounted),
    lty = 2
  )
}

if (SAVE_MEMORY_PLOTS) {
  plot_event_comparison(
    event_no_memory_smokers$trace,
    event_memory_smokers$trace,
    "Current smokers:",
    file.path(OUTPUT_DIR_MEMORY, "memory_comparison_current_smokers.png")
  )
  plot_event_comparison(
    event_no_memory_former$trace,
    event_memory_former$trace,
    "Former smokers:",
    file.path(OUTPUT_DIR_MEMORY, "memory_comparison_former_smokers.png")
  )
}

memory_results <- list(
  settings = list(
    n_patients_per_subcohort = N_PATIENTS,
    seed_current_smokers = EVENT_SEED_SMOKERS,
    seed_former_smokers = EVENT_SEED_FORMER_SMOKERS,
    max_relative_risk_increase = MEMORY_MAX_RELATIVE_RISK_INCREASE,
    decay_weights = MEMORY_DECAY_WEIGHTS,
    memory_duration_cycles = MEMORY_DURATION_CYCLES,
    maximum_total_ae_probability = maximum_total_ae_probability
  ),
  event_no_memory = list(
    current_smokers = event_no_memory_smokers,
    former_smokers = event_no_memory_former
  ),
  event_memory = list(
    current_smokers = event_memory_smokers,
    former_smokers = event_memory_former
  ),
  strategy_results = memory_strategy_results,
  subcohort_results = memory_subcohort_results,
  memory_effect_comparison = memory_effect_comparison,
  event_no_memory_vs_expected = event_no_memory_vs_expected,
  history_summary = history_summary
)

saveRDS(
  memory_results,
  file.path(OUTPUT_DIR_MEMORY, "model_results_AE_memory.rds")
)

# =============================================================================
# 9. Human-readable report
# =============================================================================
report_lines_memory <- capture.output({
  cat("============================================================\n")
  cat("COPD ACUTE-EXACERBATION MEMORY EXTENSION\n")
  cat("============================================================\n")
  cat("Patients per smoking-status subcohort:", N_PATIENTS, "\n")
  cat("Quarterly transitions:", N_ROWS - 1L, "\n")
  cat("Event-level AE sampling: at most one moderate or severe AE", "\n")
  cat("per patient and quarter.\n\n")

  cat("1. Illustrative memory assumptions\n")
  cat(
    "Maximum relative increase in subsequent AE probability:",
    100 * MEMORY_MAX_RELATIVE_RISK_INCREASE,
    "%\n"
  )
  cat(
    "Excess-risk decay weights:",
    paste(round(MEMORY_DECAY_WEIGHTS, 6), collapse = ", "),
    "\n"
  )
  cat("Memory duration:", MEMORY_DURATION_CYCLES, "quarterly cycles\n")
  cat(
    "Maximum combined moderate/severe AE probability:",
    round(maximum_total_ae_probability, 6),
    "\n"
  )
  cat(
    "These are proof-of-concept assumptions and are not empirically",
    "calibrated.\n"
  )

  cat("\n2. Strategy-level results\n")
  print(memory_strategy_results, row.names = FALSE, digits = 10)

  cat("\n3. Effect of adding memory to the event-level model\n")
  print(memory_effect_comparison, row.names = FALSE, digits = 10)

  cat("\n4. Event-level no-memory model versus expected-value baseline\n")
  print(event_no_memory_vs_expected, row.names = FALSE, digits = 10)

  cat("\n5. Exacerbation-history summary\n")
  print(history_summary, row.names = FALSE, digits = 10)

  cat("\nInterpretation notes\n")
  cat(
    "- The event-level no-memory model should approximate the frozen",
    "expected-value baseline, but a single run also contains AE event",
    "sampling error.\n"
  )
  cat(
    "- The no-memory and memory event models use common random numbers.",
    "Their difference therefore isolates the memory assumption more",
    "efficiently than unrelated simulation seeds.\n"
  )
  cat(
    "- Do not calculate an ICER between no-memory and memory model",
    "specifications. Report changes in costs, QALYs, and event counts.\n"
  )

  cat(
    "\nOutput directory:\n",
    normalizePath(OUTPUT_DIR_MEMORY),
    "\n",
    sep = ""
  )
})

writeLines(
  report_lines_memory,
  con = file.path(OUTPUT_DIR_MEMORY, "memory_extension_report.txt")
)
cat(paste(report_lines_memory, collapse = "\n"), "\n")

###############################################################################
# End of file
###############################################################################
