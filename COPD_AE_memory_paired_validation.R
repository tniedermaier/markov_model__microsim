###############################################################################
# Paired repeated-seed validation of the COPD AE-memory extension
#
# Required files in the current working directory:
#   COPD_Markov_microsimulation_main_baseline.R
#   COPD_microsimulation_AE_memory_extension.R
#
# This script:
#   1. Checks whether the event-level model without memory reproduces the
#      deterministic Markov benchmark in expectation.
#   2. Quantifies the paired effect of adding AE memory across repeated seeds.
#
# The no-memory and memory variants use identical random-number seeds within
# each repetition (common random numbers). Base R only.
###############################################################################

MEMORY_EXTENSION_FILENAME <- "COPD_microsimulation_AE_memory_extension.R"

if (!file.exists(MEMORY_EXTENSION_FILENAME)) {
  stop(
    "Required memory-extension file not found in the current working directory: ",
    MEMORY_EXTENSION_FILENAME
  )
}

# This also sources and runs the frozen baseline model.
source(MEMORY_EXTENSION_FILENAME)

# =============================================================================
# 1. Settings
# =============================================================================
N_MEMORY_REPEATS <- 20L
PAIRED_SEED_BASE <- 700000L

OUTPUT_DIR_MEMORY_VALIDATION <- file.path(
  getwd(),
  "model_output_AE_memory_validation"
)
dir.create(
  OUTPUT_DIR_MEMORY_VALIDATION,
  showWarnings = FALSE,
  recursive = TRUE
)

stopifnot(
  N_MEMORY_REPEATS >= 2L,
  N_PATIENTS >= 1L
)

# =============================================================================
# 2. Helper functions
# =============================================================================
get_strategy_row <- function(results, strategy_name) {
  out <- results[
    results$strategy == strategy_name,
    ,
    drop = FALSE
  ]

  if (nrow(out) != 1L) {
    stop("Could not uniquely identify strategy row: ", strategy_name)
  }

  out
}

get_markov_reference <- function() {
  reference <- strategy_results[
    strategy_results$model == "Markov",
    ,
    drop = FALSE
  ]

  if (nrow(reference) != 3L) {
    stop("The frozen Markov strategy reference could not be identified.")
  }

  reference
}

summarise_numeric_metric <- function(x) {
  data.frame(
    mean = mean(x),
    standard_deviation = sd(x),
    standard_error = sd(x) / sqrt(length(x)),
    median = median(x),
    percentile_2_5 = unname(quantile(x, 0.025)),
    percentile_97_5 = unname(quantile(x, 0.975)),
    minimum = min(x),
    maximum = max(x),
    n_positive = sum(x > 0),
    n_negative = sum(x < 0),
    n_zero = sum(x == 0),
    stringsAsFactors = FALSE
  )
}

make_summary_table <- function(results, metric_names) {
  rows <- lapply(metric_names, function(metric_name) {
    summary_row <- summarise_numeric_metric(results[[metric_name]])
    summary_row$metric <- metric_name
    summary_row
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  out[c(
    "metric", "mean", "standard_deviation", "standard_error", "median",
    "percentile_2_5", "percentile_97_5", "minimum", "maximum",
    "n_positive", "n_negative", "n_zero"
  )]
}

# =============================================================================
# 3. Frozen deterministic reference
# =============================================================================
markov_reference <- get_markov_reference()
markov_uc <- get_strategy_row(markov_reference, "usual_care")
markov_int <- get_strategy_row(markov_reference, "intervention")
markov_incr <- get_strategy_row(markov_reference, "incremental")

# =============================================================================
# 4. Repeated paired simulations
# =============================================================================
paired_results <- vector("list", N_MEMORY_REPEATS)

cat(
  "\nRunning", N_MEMORY_REPEATS,
  "paired AE-memory validation repetitions...\n"
)

for (r in seq_len(N_MEMORY_REPEATS)) {
  seed_smokers <- PAIRED_SEED_BASE + 1000L * r
  seed_former <- seed_smokers + 1L

  no_memory_smokers <- run_event_level_microsimulation(
    N_PATIENTS,
    PROGRESSION_SMOKER,
    MORTALITY_SMOKER,
    AE_MORTALITY_SMOKER,
    seed_smokers,
    include_memory = FALSE
  )
  no_memory_former <- run_event_level_microsimulation(
    N_PATIENTS,
    PROGRESSION_FORMER_SMOKER,
    MORTALITY_FORMER_SMOKER,
    AE_MORTALITY_FORMER_SMOKER,
    seed_former,
    include_memory = FALSE
  )

  memory_smokers <- run_event_level_microsimulation(
    N_PATIENTS,
    PROGRESSION_SMOKER,
    MORTALITY_SMOKER,
    AE_MORTALITY_SMOKER,
    seed_smokers,
    include_memory = TRUE
  )
  memory_former <- run_event_level_microsimulation(
    N_PATIENTS,
    PROGRESSION_FORMER_SMOKER,
    MORTALITY_FORMER_SMOKER,
    AE_MORTALITY_FORMER_SMOKER,
    seed_former,
    include_memory = TRUE
  )

  no_memory_strategy <- calculate_strategy_results(
    no_memory_smokers$trace,
    no_memory_former$trace,
    N_PATIENTS,
    "Event-level no memory"
  )
  memory_strategy <- calculate_strategy_results(
    memory_smokers$trace,
    memory_former$trace,
    N_PATIENTS,
    "Event-level with memory"
  )

  no_uc <- get_strategy_row(no_memory_strategy, "usual_care")
  no_int <- get_strategy_row(no_memory_strategy, "intervention")
  no_incr <- get_strategy_row(no_memory_strategy, "incremental")

  mem_uc <- get_strategy_row(memory_strategy, "usual_care")
  mem_int <- get_strategy_row(memory_strategy, "intervention")
  mem_incr <- get_strategy_row(memory_strategy, "incremental")

  history_no_smokers <- no_memory_smokers$history_summary
  history_no_former <- no_memory_former$history_summary
  history_mem_smokers <- memory_smokers$history_summary
  history_mem_former <- memory_former$history_summary

  paired_results[[r]] <- data.frame(
    repetition = r,
    seed_smokers = seed_smokers,
    seed_former = seed_former,

    # No-memory event model minus deterministic Markov reference
    no_memory_uc_direct_cost_deviation =
      no_uc$direct_costs_per_patient -
      markov_uc$direct_costs_per_patient,
    no_memory_int_direct_cost_deviation =
      no_int$direct_costs_per_patient -
      markov_int$direct_costs_per_patient,
    no_memory_incremental_direct_cost_deviation =
      no_incr$direct_costs_per_patient -
      markov_incr$direct_costs_per_patient,

    no_memory_uc_qaly_deviation =
      no_uc$qalys_per_patient -
      markov_uc$qalys_per_patient,
    no_memory_int_qaly_deviation =
      no_int$qalys_per_patient -
      markov_int$qalys_per_patient,
    no_memory_incremental_qaly_deviation =
      no_incr$qalys_per_patient -
      markov_incr$qalys_per_patient,

    no_memory_direct_icer_deviation =
      no_incr$direct_icer - markov_incr$direct_icer,

    # Paired memory effect: memory minus no memory
    memory_effect_uc_direct_cost =
      mem_uc$direct_costs_per_patient -
      no_uc$direct_costs_per_patient,
    memory_effect_int_direct_cost =
      mem_int$direct_costs_per_patient -
      no_int$direct_costs_per_patient,
    memory_effect_incremental_direct_cost =
      mem_incr$direct_costs_per_patient -
      no_incr$direct_costs_per_patient,

    memory_effect_uc_indirect_cost =
      mem_uc$indirect_costs_per_patient -
      no_uc$indirect_costs_per_patient,
    memory_effect_int_indirect_cost =
      mem_int$indirect_costs_per_patient -
      no_int$indirect_costs_per_patient,

    memory_effect_uc_societal_cost =
      mem_uc$societal_costs_per_patient -
      no_uc$societal_costs_per_patient,
    memory_effect_int_societal_cost =
      mem_int$societal_costs_per_patient -
      no_int$societal_costs_per_patient,
    memory_effect_incremental_societal_cost =
      mem_incr$societal_costs_per_patient -
      no_incr$societal_costs_per_patient,

    memory_effect_uc_qaly =
      mem_uc$qalys_per_patient -
      no_uc$qalys_per_patient,
    memory_effect_int_qaly =
      mem_int$qalys_per_patient -
      no_int$qalys_per_patient,
    memory_effect_incremental_qaly =
      mem_incr$qalys_per_patient -
      no_incr$qalys_per_patient,

    memory_effect_direct_icer =
      mem_incr$direct_icer - no_incr$direct_icer,

    # Event-count effects by smoking-status subcohort
    memory_effect_total_ae_smokers =
      history_mem_smokers$total_all_ae -
      history_no_smokers$total_all_ae,
    memory_effect_total_ae_former =
      history_mem_former$total_all_ae -
      history_no_former$total_all_ae,

    memory_effect_severe_ae_smokers =
      history_mem_smokers$total_severe_ae -
      history_no_smokers$total_severe_ae,
    memory_effect_severe_ae_former =
      history_mem_former$total_severe_ae -
      history_no_former$total_severe_ae,

    memory_effect_total_ae_smokers_percent =
      100 * (
        history_mem_smokers$total_all_ae /
          history_no_smokers$total_all_ae - 1
      ),
    memory_effect_total_ae_former_percent =
      100 * (
        history_mem_former$total_all_ae /
          history_no_former$total_all_ae - 1
      ),

    difference_in_proportion_with_any_ae_smokers =
      history_mem_smokers$proportion_with_any_ae -
      history_no_smokers$proportion_with_any_ae,
    difference_in_proportion_with_any_ae_former =
      history_mem_former$proportion_with_any_ae -
      history_no_former$proportion_with_any_ae,

    stringsAsFactors = FALSE
  )

  cat(
    "Completed repetition", r, "of", N_MEMORY_REPEATS,
    "- paired memory effect on UC QALYs:",
    sprintf("%.5f", paired_results[[r]]$memory_effect_uc_qaly),
    "\n"
  )
}

paired_results <- do.call(rbind, paired_results)
rownames(paired_results) <- NULL

# =============================================================================
# 5. Summaries
# =============================================================================
validation_metric_names <- c(
  "no_memory_uc_direct_cost_deviation",
  "no_memory_int_direct_cost_deviation",
  "no_memory_incremental_direct_cost_deviation",
  "no_memory_uc_qaly_deviation",
  "no_memory_int_qaly_deviation",
  "no_memory_incremental_qaly_deviation",
  "no_memory_direct_icer_deviation"
)

memory_effect_metric_names <- c(
  "memory_effect_uc_direct_cost",
  "memory_effect_int_direct_cost",
  "memory_effect_incremental_direct_cost",
  "memory_effect_uc_indirect_cost",
  "memory_effect_int_indirect_cost",
  "memory_effect_uc_societal_cost",
  "memory_effect_int_societal_cost",
  "memory_effect_incremental_societal_cost",
  "memory_effect_uc_qaly",
  "memory_effect_int_qaly",
  "memory_effect_incremental_qaly",
  "memory_effect_direct_icer",
  "memory_effect_total_ae_smokers",
  "memory_effect_total_ae_former",
  "memory_effect_severe_ae_smokers",
  "memory_effect_severe_ae_former",
  "memory_effect_total_ae_smokers_percent",
  "memory_effect_total_ae_former_percent",
  "difference_in_proportion_with_any_ae_smokers",
  "difference_in_proportion_with_any_ae_former"
)

no_memory_validation_summary <- make_summary_table(
  paired_results,
  validation_metric_names
)

paired_memory_effect_summary <- make_summary_table(
  paired_results,
  memory_effect_metric_names
)

# =============================================================================
# 6. Save outputs
# =============================================================================
write.csv(
  paired_results,
  file.path(
    OUTPUT_DIR_MEMORY_VALIDATION,
    "paired_memory_validation.csv"
  ),
  row.names = FALSE
)

write.csv(
  no_memory_validation_summary,
  file.path(
    OUTPUT_DIR_MEMORY_VALIDATION,
    "event_no_memory_validation_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  paired_memory_effect_summary,
  file.path(
    OUTPUT_DIR_MEMORY_VALIDATION,
    "paired_memory_effect_summary.csv"
  ),
  row.names = FALSE
)

report_file <- file.path(
  OUTPUT_DIR_MEMORY_VALIDATION,
  "paired_memory_validation_report.txt"
)

capture.output({
  cat("============================================================\n")
  cat("PAIRED REPEATED-SEED VALIDATION OF AE MEMORY EXTENSION\n")
  cat("============================================================\n")
  cat("Patients per smoking-status subcohort:", N_PATIENTS, "\n")
  cat("Number of paired repetitions:", N_MEMORY_REPEATS, "\n")
  cat(
    "The no-memory and memory variants use identical seeds within each\n",
    "repetition. Differences therefore estimate the memory effect with\n",
    "reduced Monte Carlo noise.\n\n",
    sep = ""
  )

  cat("1. Event-level no-memory model minus deterministic Markov benchmark\n")
  print(
    no_memory_validation_summary,
    row.names = FALSE,
    digits = 10
  )

  cat("\n2. Paired effect of adding memory (memory minus no memory)\n")
  print(
    paired_memory_effect_summary,
    row.names = FALSE,
    digits = 10
  )

  cat("\nInterpretation notes\n")
  cat(
    "- Mean no-memory deviations should be close to zero. This checks that\n",
    "  explicit AE event sampling preserves the expected-value benchmark.\n",
    sep = ""
  )
  cat(
    "- A stable memory effect should have the same direction in nearly all\n",
    "  paired repetitions.\n",
    sep = ""
  )
  cat(
    "- The proportion with any AE is expected to be identical by construction:\n",
    "  memory can affect only recurrent events after a first AE has occurred.\n",
    sep = ""
  )
  cat(
    "- The implemented memory is recency-based and non-cumulative: the most\n",
    "  recent AE resets the three-cycle decay profile; effects do not stack.\n",
    sep = ""
  )

  cat(
    "\nOutput directory:\n",
    normalizePath(OUTPUT_DIR_MEMORY_VALIDATION),
    "\n",
    sep = ""
  )
}, file = report_file)

cat("\nPaired AE-memory validation completed.\n")
cat("Output directory:", normalizePath(OUTPUT_DIR_MEMORY_VALIDATION), "\n")
cat("Main report:", normalizePath(report_file), "\n")
