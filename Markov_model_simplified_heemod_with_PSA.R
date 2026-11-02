# Corrected heemod implementation of the simplified COPD Markov model
# Goal: reproduce Markov_model_simplified_Github_incl_checks.r as closely as possible.
# Deterministic model plus corrected probabilistic sensitivity analysis (PSA).

rm(list = ls(all = TRUE))
options(scipen = 999)

library(heemod)

# Mixture parameters used again in the post-hoc intervention/UC calculation
eff_int <- 0.22
eff_uc  <- 0.06
c_int   <- 595.92
c_uc    <- 0

# -----------------------------------------------------------------------------
# 1. Parameters
# -----------------------------------------------------------------------------

param <- define_parameters(
  # Progression probabilities per 3-month cycle
  prog1_smo   = 0.01403,
  prog2_smo   = 0.01031,
  prog3_smo   = 0.02298,
  prog1_fosmo = 0.00247,
  prog2_fosmo = 0.00300,
  prog3_fosmo = 0.00686,

  # Exacerbation probabilities and utility decrements
  p_mod2 = 0.18376,
  d_mod2 = -0.01975,
  p_sev2 = 0.03854,
  d_sev2 = -0.03075,

  p_mod3 = 0.26767,
  d_mod3 = -0.01875,
  p_sev3 = 0.03774,
  d_sev3 = -0.02919,

  p_mod4 = 0.28057,
  d_mod4 = -0.01625,
  p_sev4 = 0.04855,
  d_sev4 = -0.03502,

  # State utilities
  u_1 = 0.84,
  u_2 = 0.79,
  u_3 = 0.75,
  u_4 = 0.65,

  # Direct costs: these are used AS-IS per model cycle in the Base-R reference.
  # Therefore they are NOT multiplied by cycle_length / 12.
  c_1    = 102.64,
  c_2    = 184.69,
  c_mod2 = 39.86,
  c_sev2 = 2923.94,
  c_3    = 366.57,
  c_mod3 = 52.90,
  c_sev3 = 3109.05,
  c_4    = 431.35,
  c_mod4 = 51.52,
  c_sev4 = 3322.59,

  # Intervention / usual-care mixture parameters
  c_int   = c_int,
  c_uc    = c_uc,
  eff_int = eff_int,
  eff_uc  = eff_uc,

  # Cycle length and discount rate
  cycle_length = 3,
  disc_e       = 0.03,

  # Exact quarterly discount factor used in the Base-R reference:
  # Cycle 1 -> exponent 0.25, cycle 2 -> 0.50, ..., cycle 240 -> 60.
  disc_factor = 1 / ((1 + disc_e) ^ (model_time * cycle_length / 12)),

  # ---------------------------------------------------------------------------
  # Time-dependent mortality: smokers
  # Transition parameters use blocks 1:20, 21:40, ... exactly as in Base R.
  # ---------------------------------------------------------------------------
  mort1_smo = ifelse(model_time <= 20, 0.00122,
               ifelse(model_time <= 40, 0.00210,
               ifelse(model_time <= 60, 0.00302,
               ifelse(model_time <= 80, 0.00470,
               ifelse(model_time <= 100, 0.00805,
               ifelse(model_time <= 120, 0.01358,
               ifelse(model_time <= 140, 0.01970,
               ifelse(model_time <= 160, 0.03400,
               ifelse(model_time <= 180, 0.04729, 0.08840))))))))),

  mort2_smo = ifelse(model_time <= 20, 0.00210,
               ifelse(model_time <= 40, 0.00360,
               ifelse(model_time <= 60, 0.00518,
               ifelse(model_time <= 80, 0.00804,
               ifelse(model_time <= 100, 0.01377,
               ifelse(model_time <= 120, 0.02318,
               ifelse(model_time <= 140, 0.03354,
               ifelse(model_time <= 160, 0.05758,
               ifelse(model_time <= 180, 0.07969, 0.14672))))))))),

  mort3_smo = ifelse(model_time <= 20, 0.00497,
               ifelse(model_time <= 40, 0.00852,
               ifelse(model_time <= 60, 0.01225,
               ifelse(model_time <= 80, 0.01899,
               ifelse(model_time <= 100, 0.03238,
               ifelse(model_time <= 120, 0.05417,
               ifelse(model_time <= 140, 0.07783,
               ifelse(model_time <= 160, 0.13138,
               ifelse(model_time <= 180, 0.17900, 0.31397))))))))),

  ae_mort2_smo = ifelse(model_time <= 20, 0.13819,
                  ifelse(model_time <= 40, 0.13690,
                  ifelse(model_time <= 60, 0.13552,
                  ifelse(model_time <= 80, 0.25400,
                  ifelse(model_time <= 100, 0.24967,
                  ifelse(model_time <= 120, 0.24244,
                  ifelse(model_time <= 140, 0.30675,
                  ifelse(model_time <= 160, 0.28906,
                  ifelse(model_time <= 180, 0.27198, 0.21479))))))))),

  ae_mort3_smo = ifelse(model_time <= 20, 0.13570,
                  ifelse(model_time <= 40, 0.13261,
                  ifelse(model_time <= 60, 0.12933,
                  ifelse(model_time <= 80, 0.24568,
                  ifelse(model_time <= 100, 0.23523,
                  ifelse(model_time <= 120, 0.21762,
                  ifelse(model_time <= 140, 0.27345,
                  ifelse(model_time <= 160, 0.22866,
                  ifelse(model_time <= 180, 0.18392, 0.02336))))))))),

  # ---------------------------------------------------------------------------
  # Time-dependent mortality: former smokers
  # ---------------------------------------------------------------------------
  mort1_fosmo = ifelse(model_time <= 20, 0.00080,
                 ifelse(model_time <= 40, 0.00137,
                 ifelse(model_time <= 60, 0.00251,
                 ifelse(model_time <= 80, 0.00390,
                 ifelse(model_time <= 100, 0.00504,
                 ifelse(model_time <= 120, 0.00851,
                 ifelse(model_time <= 140, 0.01318,
                 ifelse(model_time <= 160, 0.02280,
                 ifelse(model_time <= 180, 0.03957, 0.07423))))))))),

  mort2_fosmo = ifelse(model_time <= 20, 0.00136,
                 ifelse(model_time <= 40, 0.00234,
                 ifelse(model_time <= 60, 0.00430,
                 ifelse(model_time <= 80, 0.00668,
                 ifelse(model_time <= 100, 0.00863,
                 ifelse(model_time <= 120, 0.01455,
                 ifelse(model_time <= 140, 0.02249,
                 ifelse(model_time <= 160, 0.03877,
                 ifelse(model_time <= 180, 0.06686, 0.12386))))))))),

  mort3_fosmo = ifelse(model_time <= 20, 0.00323,
                 ifelse(model_time <= 40, 0.00555,
                 ifelse(model_time <= 60, 0.01018,
                 ifelse(model_time <= 80, 0.01578,
                 ifelse(model_time <= 100, 0.02037,
                 ifelse(model_time <= 120, 0.03421,
                 ifelse(model_time <= 140, 0.05258,
                 ifelse(model_time <= 160, 0.08963,
                 ifelse(model_time <= 180, 0.15156, 0.26950))))))))),

  ae_mort2_fosmo = ifelse(model_time <= 20, 0.13883,
                    ifelse(model_time <= 40, 0.13798,
                    ifelse(model_time <= 60, 0.13629,
                    ifelse(model_time <= 80, 0.25503,
                    ifelse(model_time <= 100, 0.25356,
                    ifelse(model_time <= 120, 0.24908,
                    ifelse(model_time <= 140, 0.31459,
                    ifelse(model_time <= 160, 0.30298,
                    ifelse(model_time <= 180, 0.28199, 0.23529))))))))),

  ae_mort3_fosmo = ifelse(model_time <= 20, 0.13721,
                    ifelse(model_time <= 40, 0.13520,
                    ifelse(model_time <= 60, 0.13116,
                    ifelse(model_time <= 80, 0.24813,
                    ifelse(model_time <= 100, 0.24462,
                    ifelse(model_time <= 120, 0.23379,
                    ifelse(model_time <= 140, 0.29281,
                    ifelse(model_time <= 160, 0.26404,
                    ifelse(model_time <= 180, 0.21031, 0.08281))))))))),

  # ---------------------------------------------------------------------------
  # Indirect costs.
  # IMPORTANT: The Base-R vectors contain an explicit value for t = 0.
  # With method = "end", heemod model_time = 1 corresponds to Base-R row 2.
  # Therefore the cut-points are 19, 39, 59, 79 (NOT 20, 40, 60, 80).
  # ---------------------------------------------------------------------------
  ind_stable2 = ifelse(model_time <= 19, 112.50,
                 ifelse(model_time <= 39, 376.37,
                 ifelse(model_time <= 59, 667.24,
                 ifelse(model_time <= 79, 237.22, 0)))),
  ind_stable3 = ifelse(model_time <= 19, 268.06,
                 ifelse(model_time <= 39, 896.81,
                 ifelse(model_time <= 59, 1589.90,
                 ifelse(model_time <= 79, 565.24, 0)))),
  ind_stable4 = ifelse(model_time <= 19, 281.24,
                 ifelse(model_time <= 39, 940.92,
                 ifelse(model_time <= 59, 1668.10,
                 ifelse(model_time <= 79, 593.04, 0)))),

  ind_mod2 = ifelse(model_time <= 19, 491.94,
              ifelse(model_time <= 39, 446.16,
              ifelse(model_time <= 59, 364.35,
              ifelse(model_time <= 79, 115.30, 0)))),
  ind_mod3 = ifelse(model_time <= 19, 550.51,
              ifelse(model_time <= 39, 499.27,
              ifelse(model_time <= 59, 407.73,
              ifelse(model_time <= 79, 129.03, 0)))),
  ind_mod4 = ifelse(model_time <= 19, 491.94,
              ifelse(model_time <= 39, 446.16,
              ifelse(model_time <= 59, 364.35,
              ifelse(model_time <= 79, 115.30, 0)))),

  ind_sev2 = ifelse(model_time <= 19, 367.00,
              ifelse(model_time <= 39, 332.85,
              ifelse(model_time <= 59, 271.82,
              ifelse(model_time <= 79, 86.02, 0)))),
  ind_sev3 = ifelse(model_time <= 19, 663.73,
              ifelse(model_time <= 39, 601.96,
              ifelse(model_time <= 59, 491.58,
              ifelse(model_time <= 79, 155.56, 0)))),
  ind_sev4 = ifelse(model_time <= 19, 648.12,
              ifelse(model_time <= 39, 587.79,
              ifelse(model_time <= 59, 480.02,
              ifelse(model_time <= 79, 151.90, 0))))
)

# -----------------------------------------------------------------------------
# 2. Transition matrices
# -----------------------------------------------------------------------------

mat_smo <- define_transition(
  1 - prog1_smo - mort1_smo,                         prog1_smo, 0, 0, mort1_smo,
  0, 1 - prog2_smo - mort2_smo - p_sev2*ae_mort2_smo, prog2_smo, 0, mort2_smo + p_sev2*ae_mort2_smo,
  0, 0, 1 - prog3_smo - mort3_smo - p_sev3*ae_mort3_smo, prog3_smo, mort3_smo + p_sev3*ae_mort3_smo,
  0, 0, 0, 1 - mort3_smo - p_sev4*ae_mort3_smo, mort3_smo + p_sev4*ae_mort3_smo,
  0, 0, 0, 0, 1,
  state_names = c("stage1", "stage2", "stage3", "stage4", "death")
)

mat_fosmo <- define_transition(
  1 - prog1_fosmo - mort1_fosmo,                         prog1_fosmo, 0, 0, mort1_fosmo,
  0, 1 - prog2_fosmo - mort2_fosmo - p_sev2*ae_mort2_fosmo, prog2_fosmo, 0, mort2_fosmo + p_sev2*ae_mort2_fosmo,
  0, 0, 1 - prog3_fosmo - mort3_fosmo - p_sev3*ae_mort3_fosmo, prog3_fosmo, mort3_fosmo + p_sev3*ae_mort3_fosmo,
  0, 0, 0, 1 - mort3_fosmo - p_sev4*ae_mort3_fosmo, mort3_fosmo + p_sev4*ae_mort3_fosmo,
  0, 0, 0, 0, 1,
  state_names = c("stage1", "stage2", "stage3", "stage4", "death")
)

# -----------------------------------------------------------------------------
# 3. State values
# -----------------------------------------------------------------------------
# QALYs are annual utilities multiplied by 3/12.
# Direct and indirect costs are NOT multiplied by 3/12, because the Base-R
# reference uses the supplied cost values directly in every 3-month cycle.
#
# LYG note: the Base-R reference manually sets both t=0 and t=240 to zero.
# t=0 is naturally omitted by method="end". The model_time==240 rule below
# reproduces the second zero exactly. Remove it if you later decide that the
# last-cycle zero was only an implementation artefact in the reference model.

state_stage1 <- define_state(
  qaly             = u_1 * cycle_length/12,
  qaly_disc        = u_1 * cycle_length/12 * disc_factor,
  lyg              = ifelse(model_time == 240, 0, cycle_length/12),
  lyg_disc         = ifelse(model_time == 240, 0, cycle_length/12) * disc_factor,
  direct_cost      = c_1,
  direct_cost_disc = c_1 * disc_factor,
  indirect_cost      = 0,
  indirect_cost_disc = 0,
  total_cost_disc    = c_1 * disc_factor
)

state_stage2 <- define_state(
  qaly      = (u_2 + p_mod2*d_mod2 + p_sev2*d_sev2) * cycle_length/12,
  qaly_disc = (u_2 + p_mod2*d_mod2 + p_sev2*d_sev2) * cycle_length/12 * disc_factor,
  lyg       = ifelse(model_time == 240, 0, cycle_length/12),
  lyg_disc  = ifelse(model_time == 240, 0, cycle_length/12) * disc_factor,
  direct_cost      = c_2 + p_mod2*c_mod2 + p_sev2*c_sev2,
  direct_cost_disc = (c_2 + p_mod2*c_mod2 + p_sev2*c_sev2) * disc_factor,
  indirect_cost      = ind_stable2 + p_mod2*ind_mod2 + p_sev2*ind_sev2,
  indirect_cost_disc = (ind_stable2 + p_mod2*ind_mod2 + p_sev2*ind_sev2) * disc_factor,
  total_cost_disc = (c_2 + p_mod2*c_mod2 + p_sev2*c_sev2 +
                     ind_stable2 + p_mod2*ind_mod2 + p_sev2*ind_sev2) * disc_factor
)

state_stage3 <- define_state(
  qaly      = (u_3 + p_mod3*d_mod3 + p_sev3*d_sev3) * cycle_length/12,
  qaly_disc = (u_3 + p_mod3*d_mod3 + p_sev3*d_sev3) * cycle_length/12 * disc_factor,
  lyg       = ifelse(model_time == 240, 0, cycle_length/12),
  lyg_disc  = ifelse(model_time == 240, 0, cycle_length/12) * disc_factor,
  direct_cost      = c_3 + p_mod3*c_mod3 + p_sev3*c_sev3,
  direct_cost_disc = (c_3 + p_mod3*c_mod3 + p_sev3*c_sev3) * disc_factor,
  indirect_cost      = ind_stable3 + p_mod3*ind_mod3 + p_sev3*ind_sev3,
  indirect_cost_disc = (ind_stable3 + p_mod3*ind_mod3 + p_sev3*ind_sev3) * disc_factor,
  total_cost_disc = (c_3 + p_mod3*c_mod3 + p_sev3*c_sev3 +
                     ind_stable3 + p_mod3*ind_mod3 + p_sev3*ind_sev3) * disc_factor
)

state_stage4 <- define_state(
  qaly      = (u_4 + p_mod4*d_mod4 + p_sev4*d_sev4) * cycle_length/12,
  qaly_disc = (u_4 + p_mod4*d_mod4 + p_sev4*d_sev4) * cycle_length/12 * disc_factor,
  lyg       = ifelse(model_time == 240, 0, cycle_length/12),
  lyg_disc  = ifelse(model_time == 240, 0, cycle_length/12) * disc_factor,
  direct_cost      = c_4 + p_mod4*c_mod4 + p_sev4*c_sev4,
  direct_cost_disc = (c_4 + p_mod4*c_mod4 + p_sev4*c_sev4) * disc_factor,
  indirect_cost      = ind_stable4 + p_mod4*ind_mod4 + p_sev4*ind_sev4,
  indirect_cost_disc = (ind_stable4 + p_mod4*ind_mod4 + p_sev4*ind_sev4) * disc_factor,
  total_cost_disc = (c_4 + p_mod4*c_mod4 + p_sev4*c_sev4 +
                     ind_stable4 + p_mod4*ind_mod4 + p_sev4*ind_sev4) * disc_factor
)

state_death <- define_state(
  qaly = 0,
  qaly_disc = 0,
  lyg = 0,
  lyg_disc = 0,
  direct_cost = 0,
  direct_cost_disc = 0,
  indirect_cost = 0,
  indirect_cost_disc = 0,
  total_cost_disc = 0
)

# -----------------------------------------------------------------------------
# 4. Strategies and model run
# -----------------------------------------------------------------------------

strategy_smo <- define_strategy(
  transition = mat_smo,
  stage1 = state_stage1,
  stage2 = state_stage2,
  stage3 = state_stage3,
  stage4 = state_stage4,
  death  = state_death
)

strategy_fosmo <- define_strategy(
  transition = mat_fosmo,
  stage1 = state_stage1,
  stage2 = state_stage2,
  stage3 = state_stage3,
  stage4 = state_stage4,
  death  = state_death
)

# The Base-R data frames contain t=0 plus 240 transitions (241 rows total).
# Therefore heemod must compute 240 Markov cycles, using END-of-cycle counts.
res <- run_model(
  smo   = strategy_smo,
  fosmo = strategy_fosmo,
  parameters = param,
  init = c(stage1 = 100000, stage2 = 0, stage3 = 0, stage4 = 0, death = 0),
  cycles = 240,
  method = "end",
  cost   = direct_cost_disc,
  effect = qaly_disc
)

print(summary(res))

# -----------------------------------------------------------------------------
# 5. Reconstruct a 241-row state trace (t = 0,...,240) for direct comparison
#    with dat_smo / dat_fosmo from the Base-R reference.
# -----------------------------------------------------------------------------

make_trace <- function(model_result, strategy_name) {
  x <- subset(get_counts(model_result), .strategy_names == strategy_name)
  m <- xtabs(count ~ model_time + state_names, data = x)
  m <- m[, c("stage1", "stage2", "stage3", "stage4", "death"), drop = FALSE]

  out <- rbind(
    c(stage1 = 100000, stage2 = 0, stage3 = 0, stage4 = 0, death = 0),
    m
  )
  rownames(out) <- 0:240
  out
}

trace_smo   <- make_trace(res, "smo")
trace_fosmo <- make_trace(res, "fosmo")

# Example plots
matplot(0:240, trace_smo, type = "l", lty = 1,
        xlab = "3-month cycle", ylab = "Participants", main = "Smokers")
legend("right", legend = colnames(trace_smo), lty = 1, bty = "n")

matplot(0:240, trace_fosmo, type = "l", lty = 1,
        xlab = "3-month cycle", ylab = "Participants", main = "Former smokers")
legend("right", legend = colnames(trace_fosmo), lty = 1, bty = "n")

# -----------------------------------------------------------------------------
# 6. Extract totals and compare with the Base-R reference
# -----------------------------------------------------------------------------

values_long <- get_values(res)

get_total <- function(strategy_name, value_name) {
  z <- subset(values_long,
              .strategy_names == strategy_name & value_names == value_name)
  sum(z$value)
}

results <- data.frame(
  outcome = c("LYG", "QALY", "direct_cost", "LYG_discounted",
              "QALY_discounted", "direct_cost_discounted",
              "indirect_cost", "indirect_cost_discounted"),
  smo = c(
    get_total("smo", "lyg"),
    get_total("smo", "qaly"),
    get_total("smo", "direct_cost"),
    get_total("smo", "lyg_disc"),
    get_total("smo", "qaly_disc"),
    get_total("smo", "direct_cost_disc"),
    get_total("smo", "indirect_cost"),
    get_total("smo", "indirect_cost_disc")
  ),
  fosmo = c(
    get_total("fosmo", "lyg"),
    get_total("fosmo", "qaly"),
    get_total("fosmo", "direct_cost"),
    get_total("fosmo", "lyg_disc"),
    get_total("fosmo", "qaly_disc"),
    get_total("fosmo", "direct_cost_disc"),
    get_total("fosmo", "indirect_cost"),
    get_total("fosmo", "indirect_cost_disc")
  )
)

print(results, digits = 12)

# Reference values obtained by executing the deterministic Base-R equations.
# Small floating-point differences are acceptable.
reference <- data.frame(
  outcome = results$outcome,
  smo = c(
    2186346.178464,
    1775640.995481,
    1757475693.190276,
    1531263.688566,
    1248522.796181,
    1164532502.082168,
    1338300535.003893,
    949410662.279107
  ),
  fosmo = c(
    3020695.648851,
    2520181.512542,
    1494553273.078705,
    1903962.912134,
    1589728.066062,
    923630943.568449,
    278945552.885145,
    196146845.466009
  )
)

comparison <- transform(
  results,
  ref_smo = reference$smo,
  ref_fosmo = reference$fosmo,
  diff_smo = smo - reference$smo,
  diff_fosmo = fosmo - reference$fosmo
)
print(comparison, digits = 12)

# -----------------------------------------------------------------------------
# 7. Reproduce the final intervention/usual-care calculations from Base R
# -----------------------------------------------------------------------------

n_cohort <- 100000

smo_qaly_disc   <- get_total("smo", "qaly_disc") / n_cohort
fosmo_qaly_disc <- get_total("fosmo", "qaly_disc") / n_cohort

smo_lyg_disc   <- get_total("smo", "lyg_disc") / n_cohort
fosmo_lyg_disc <- get_total("fosmo", "lyg_disc") / n_cohort

smo_dir_disc   <- get_total("smo", "direct_cost_disc") / n_cohort
fosmo_dir_disc <- get_total("fosmo", "direct_cost_disc") / n_cohort

smo_ind_disc   <- get_total("smo", "indirect_cost_disc") / n_cohort
fosmo_ind_disc <- get_total("fosmo", "indirect_cost_disc") / n_cohort

# Direct-cost perspective
cost_dir_uc <- smo_dir_disc * (1 - eff_uc) + fosmo_dir_disc * eff_uc + c_uc
cost_dir_int <- smo_dir_disc * (1 - eff_int) + fosmo_dir_disc * eff_int + c_int

qaly_uc <- smo_qaly_disc * (1 - eff_uc) + fosmo_qaly_disc * eff_uc
qaly_int <- smo_qaly_disc * (1 - eff_int) + fosmo_qaly_disc * eff_int

lyg_uc <- smo_lyg_disc * (1 - eff_uc) + fosmo_lyg_disc * eff_uc
lyg_int <- smo_lyg_disc * (1 - eff_int) + fosmo_lyg_disc * eff_int

icer_direct <- (cost_dir_int - cost_dir_uc) / (qaly_int - qaly_uc)

# Societal perspective: direct + indirect costs
cost_soc_uc <- (smo_dir_disc + smo_ind_disc) * (1 - eff_uc) +
               (fosmo_dir_disc + fosmo_ind_disc) * eff_uc + c_uc
cost_soc_int <- (smo_dir_disc + smo_ind_disc) * (1 - eff_int) +
                (fosmo_dir_disc + fosmo_ind_disc) * eff_int + c_int

icer_societal <- (cost_soc_int - cost_soc_uc) / (qaly_int - qaly_uc)

final_results <- data.frame(
  quantity = c("Direct costs UC", "Direct costs intervention",
               "QALY UC", "QALY intervention",
               "LYG UC", "LYG intervention",
               "ICER direct costs",
               "Societal costs UC", "Societal costs intervention",
               "ICER societal"),
  value = c(cost_dir_uc, cost_dir_int,
            qaly_uc, qaly_int,
            lyg_uc, lyg_int,
            icer_direct,
            cost_soc_uc, cost_soc_int,
            icer_societal)
)
print(final_results, digits = 12)

# Expected approximate final results from the Base-R reference:
# Direct costs UC            11500.7841
# Direct costs intervention  11711.2616
# QALY UC                       12.6899511
# QALY intervention             13.2358796
# LYG UC                        15.5362564
# LYG intervention              16.1325752
# ICER direct costs            385.5405 EUR/QALY
# Societal costs UC          20542.9324
# Societal costs intervention 19548.1878
# ICER societal             -1822.1154 EUR/QALY

# -----------------------------------------------------------------------------
# 8. Probabilistic sensitivity analysis (PSA)
# -----------------------------------------------------------------------------
# Methodological choices:
# 1) One value is sampled for each uncertain parameter per PSA iteration and is
#    then held constant over all 240 cycles of that iteration.
# 2) The SAME sampled values for common parameters (e.g. utilities, AE risks,
#    direct costs) are used for the smoker and former-smoker submodels.
# 3) Smoker-specific and former-smoker-specific progression probabilities are
#    separate uncertain parameters and are sampled independently.
# 4) Intervention and usual care are NOT the two heemod strategies above.
#    They are mixtures of the smoker and former-smoker submodels. Therefore the
#    final intervention-vs-UC PSA is constructed post hoc from the paired
#    subgroup results, using the same sampled eff_int, eff_uc and c_int.
# 5) No uncertainty distributions were specified in the reference code for
#    background mortality, AE excess mortality, utility decrements, discount
#    rate, or the unit values used for indirect costs. These remain fixed.
#    The indirect-cost totals are nevertheless uncertain because p_mod/p_sev
#    and the state distributions are uncertain.
# 6) No correlation matrix was specified in the reference PSA. Therefore all
#    DISTINCT resampled parameters are assumed independent. Parameters shared
#    by both submodels are shared by definition, not resampled twice.

# Helper functions: convert a mean and standard error into beta shape
# parameters. These are the same formulae used in the Base-R reference.
beta_shape1 <- function(mu, se) {
  mu^2 * (1 - mu) / se^2 - mu
}

beta_shape2 <- function(mu, se) {
  mu * (1 - mu)^2 / se^2 - (1 - mu)
}

# Define all parameter distributions in ONE PSA object. This is important:
# run_psa() then creates one joint parameter draw per iteration and evaluates
# both submodels with that same draw.
psa_def <- define_psa(
  # Progression: smokers
  prog1_smo ~ beta(
    shape1 = beta_shape1(0.01403, 0.00053),
    shape2 = beta_shape2(0.01403, 0.00053)
  ),
  prog2_smo ~ beta(
    shape1 = beta_shape1(0.01031, 0.00034),
    shape2 = beta_shape2(0.01031, 0.00034)
  ),
  prog3_smo ~ beta(
    shape1 = beta_shape1(0.02298, 0.00162),
    shape2 = beta_shape2(0.02298, 0.00162)
  ),

  # Progression: former smokers
  prog1_fosmo ~ beta(
    shape1 = beta_shape1(0.00247, 0.00036),
    shape2 = beta_shape2(0.00247, 0.00036)
  ),
  prog2_fosmo ~ beta(
    shape1 = beta_shape1(0.00300, 0.00024),
    shape2 = beta_shape2(0.00300, 0.00024)
  ),
  prog3_fosmo ~ beta(
    shape1 = beta_shape1(0.00686, 0.00036),
    shape2 = beta_shape2(0.00686, 0.00036)
  ),

  # Exacerbation probabilities
  p_mod2 ~ beta(
    shape1 = beta_shape1(0.18376, 0.01090),
    shape2 = beta_shape2(0.18376, 0.01090)
  ),
  p_sev2 ~ beta(
    shape1 = beta_shape1(0.03854, 0.00246),
    shape2 = beta_shape2(0.03854, 0.00246)
  ),
  p_mod3 ~ beta(
    shape1 = beta_shape1(0.26767, 0.01412),
    shape2 = beta_shape2(0.26767, 0.01412)
  ),
  p_sev3 ~ beta(
    shape1 = beta_shape1(0.03774, 0.002250),
    shape2 = beta_shape2(0.03774, 0.002250)
  ),
  p_mod4 ~ beta(
    shape1 = beta_shape1(0.28057, 0.009560),
    shape2 = beta_shape2(0.28057, 0.009560)
  ),
  p_sev4 ~ beta(
    shape1 = beta_shape1(0.04855, 0.001890),
    shape2 = beta_shape2(0.04855, 0.001890)
  ),

  # Utilities
  u_1 ~ beta(
    shape1 = beta_shape1(0.84, 0.0300),
    shape2 = beta_shape2(0.84, 0.0300)
  ),
  u_2 ~ beta(
    shape1 = beta_shape1(0.79, 0.0100),
    shape2 = beta_shape2(0.79, 0.0100)
  ),
  u_3 ~ beta(
    shape1 = beta_shape1(0.75, 0.0100),
    shape2 = beta_shape2(0.75, 0.0100)
  ),
  u_4 ~ beta(
    shape1 = beta_shape1(0.65, 0.0300),
    shape2 = beta_shape2(0.65, 0.0300)
  ),

  # Direct costs
  c_1    ~ gamma(mean = 102.64,  sd = 20.00),
  c_2    ~ gamma(mean = 184.69,  sd = 25.00),
  c_3    ~ gamma(mean = 366.57,  sd = 30.00),
  c_4    ~ gamma(mean = 431.35,  sd = 35.00),
  c_mod2 ~ gamma(mean = 39.86,   sd = 10.00),
  c_mod3 ~ gamma(mean = 52.90,   sd = 10.00),
  c_mod4 ~ gamma(mean = 51.52,   sd = 10.00),
  c_sev2 ~ gamma(mean = 2923.94, sd = 500.00),
  c_sev3 ~ gamma(mean = 3109.05, sd = 500.00),
  c_sev4 ~ gamma(mean = 3322.59, sd = 500.00),

  # Intervention cost and smoking-cessation proportions
  c_int ~ gamma(mean = 595.92, sd = 50.00),
  eff_int ~ beta(
    shape1 = beta_shape1(0.22, 0.07),
    shape2 = beta_shape2(0.22, 0.07)
  ),
  eff_uc ~ beta(
    shape1 = beta_shape1(0.06, 0.02),
    shape2 = beta_shape2(0.06, 0.02)
  )
)

# PSA settings. Increase N_PSA for final analyses if desired.
N_PSA    <- 1000L
PSA_SEED <- 1L

# Optional heemod parallelisation. The current heemod vignette recommends
# using no more than about four cores because gains beyond that can be small.
USE_PARALLEL <- FALSE
N_CORES <- parallel::detectCores(logical = TRUE)
if (is.na(N_CORES)) N_CORES <- 1L
N_CORES <- max(1L, min(4L, as.integer(N_CORES)))

run_joint_psa <- function() {
  if (USE_PARALLEL && N_CORES > 1L) {
    use_cluster(N_CORES)
    on.exit(close_cluster(), add = TRUE)
  }

  set.seed(PSA_SEED)
  run_psa(
    model = res,
    psa   = psa_def,
    N     = N_PSA,
    keep  = FALSE
  )
}

psa_joint <- run_joint_psa()

# Important: summary(psa_joint) and plot(psa_joint, ...) compare the two
# SUBMODELS 'smo' and 'fosmo'. They do NOT compare the actual decision
# strategies intervention vs usual care. The latter are constructed below.

# -----------------------------------------------------------------------------
# 9. Pair smoker/former-smoker outcomes within each PSA iteration
# -----------------------------------------------------------------------------

psa_tab <- psa_joint$psa

required_cols <- c(
  ".index", ".strategy_names",
  "qaly_disc", "lyg_disc",
  "direct_cost_disc", "indirect_cost_disc",
  "eff_int", "eff_uc", "c_int"
)

missing_cols <- setdiff(required_cols, names(psa_tab))
if (length(missing_cols) > 0L) {
  stop(
    "Expected columns missing from psa_joint$psa: ",
    paste(missing_cols, collapse = ", ")
  )
}

psa_smo <- psa_tab[psa_tab$.strategy_names == "smo", , drop = FALSE]
psa_fosmo <- psa_tab[psa_tab$.strategy_names == "fosmo", , drop = FALSE]

psa_smo <- psa_smo[order(psa_smo$.index), , drop = FALSE]
psa_fosmo <- psa_fosmo[order(psa_fosmo$.index), , drop = FALSE]

if (nrow(psa_smo) != N_PSA || nrow(psa_fosmo) != N_PSA) {
  stop("Unexpected number of PSA rows for smoker/former-smoker submodels.")
}
if (!identical(psa_smo$.index, psa_fosmo$.index)) {
  stop("PSA iteration indices do not match between submodels.")
}

# Check that common sampled parameters are genuinely paired between both
# submodels. run_psa() should use the same parameter table for every strategy.
shared_draws <- c(
  "p_mod2", "p_sev2", "p_mod3", "p_sev3", "p_mod4", "p_sev4",
  "u_1", "u_2", "u_3", "u_4",
  "c_1", "c_2", "c_3", "c_4",
  "c_mod2", "c_mod3", "c_mod4", "c_sev2", "c_sev3", "c_sev4",
  "c_int", "eff_int", "eff_uc"
)
shared_draws <- intersect(shared_draws, names(psa_smo))

same_draw <- vapply(
  shared_draws,
  function(v) identical(psa_smo[[v]], psa_fosmo[[v]]),
  logical(1)
)
if (!all(same_draw)) {
  stop(
    "Common PSA draws differ between submodels for: ",
    paste(names(same_draw)[!same_draw], collapse = ", ")
  )
}

# Convert cohort totals to per-patient quantities before the intervention/UC
# mixture calculation, exactly as in the deterministic Base-R reference.
psa_decision <- data.frame(
  sim = psa_smo$.index,

  eff_int = psa_smo$eff_int,
  eff_uc  = psa_smo$eff_uc,
  c_int   = psa_smo$c_int,

  qaly_smo   = psa_smo$qaly_disc / n_cohort,
  qaly_fosmo = psa_fosmo$qaly_disc / n_cohort,

  lyg_smo   = psa_smo$lyg_disc / n_cohort,
  lyg_fosmo = psa_fosmo$lyg_disc / n_cohort,

  dir_cost_smo   = psa_smo$direct_cost_disc / n_cohort,
  dir_cost_fosmo = psa_fosmo$direct_cost_disc / n_cohort,

  ind_cost_smo   = psa_smo$indirect_cost_disc / n_cohort,
  ind_cost_fosmo = psa_fosmo$indirect_cost_disc / n_cohort
)

# -----------------------------------------------------------------------------
# 10. Actual decision strategies: intervention vs usual care
# -----------------------------------------------------------------------------

# Effects
psa_decision$qaly_int <- with(
  psa_decision,
  qaly_smo * (1 - eff_int) + qaly_fosmo * eff_int
)
psa_decision$qaly_uc <- with(
  psa_decision,
  qaly_smo * (1 - eff_uc) + qaly_fosmo * eff_uc
)

psa_decision$lyg_int <- with(
  psa_decision,
  lyg_smo * (1 - eff_int) + lyg_fosmo * eff_int
)
psa_decision$lyg_uc <- with(
  psa_decision,
  lyg_smo * (1 - eff_uc) + lyg_fosmo * eff_uc
)

# Direct-cost perspective
psa_decision$cost_dir_int <- with(
  psa_decision,
  dir_cost_smo * (1 - eff_int) + dir_cost_fosmo * eff_int + c_int
)
psa_decision$cost_dir_uc <- with(
  psa_decision,
  dir_cost_smo * (1 - eff_uc) + dir_cost_fosmo * eff_uc + c_uc
)

# Societal perspective: direct + indirect costs
psa_decision$cost_soc_int <- with(
  psa_decision,
  (dir_cost_smo + ind_cost_smo) * (1 - eff_int) +
    (dir_cost_fosmo + ind_cost_fosmo) * eff_int + c_int
)
psa_decision$cost_soc_uc <- with(
  psa_decision,
  (dir_cost_smo + ind_cost_smo) * (1 - eff_uc) +
    (dir_cost_fosmo + ind_cost_fosmo) * eff_uc + c_uc
)

# Incremental results: intervention minus usual care
psa_decision$delta_qaly <- psa_decision$qaly_int - psa_decision$qaly_uc
psa_decision$delta_lyg  <- psa_decision$lyg_int  - psa_decision$lyg_uc
psa_decision$delta_cost_dir <- psa_decision$cost_dir_int - psa_decision$cost_dir_uc
psa_decision$delta_cost_soc <- psa_decision$cost_soc_int - psa_decision$cost_soc_uc

# ICERs are retained for descriptive use, but CE probabilities below are
# calculated from incremental net monetary benefit (INMB), not by comparing
# ICERs with a threshold. INMB remains valid in all quadrants of the CE plane.
eps_effect <- 1e-12
psa_decision$icer_direct <- ifelse(
  abs(psa_decision$delta_qaly) > eps_effect,
  psa_decision$delta_cost_dir / psa_decision$delta_qaly,
  NA_real_
)
psa_decision$icer_societal <- ifelse(
  abs(psa_decision$delta_qaly) > eps_effect,
  psa_decision$delta_cost_soc / psa_decision$delta_qaly,
  NA_real_
)

# -----------------------------------------------------------------------------
# 11. PSA summaries
# -----------------------------------------------------------------------------

summarise_vector <- function(x) {
  c(
    mean   = mean(x, na.rm = TRUE),
    sd     = sd(x, na.rm = TRUE),
    q025   = unname(quantile(x, 0.025, na.rm = TRUE)),
    median = median(x, na.rm = TRUE),
    q975   = unname(quantile(x, 0.975, na.rm = TRUE))
  )
}

summary_vars <- c(
  "qaly_uc", "qaly_int", "delta_qaly",
  "lyg_uc", "lyg_int", "delta_lyg",
  "cost_dir_uc", "cost_dir_int", "delta_cost_dir",
  "cost_soc_uc", "cost_soc_int", "delta_cost_soc"
)

psa_summary <- do.call(
  rbind,
  lapply(summary_vars, function(v) {
    z <- summarise_vector(psa_decision[[v]])
    data.frame(
      outcome = v,
      mean    = unname(z["mean"]),
      sd      = unname(z["sd"]),
      q025    = unname(z["q025"]),
      median  = unname(z["median"]),
      q975    = unname(z["q975"]),
      row.names = NULL,
      check.names = FALSE
    )
  })
)

cat("\nPSA summary (per patient; 95% simulation interval):\n")
print(psa_summary, digits = 8, row.names = FALSE)

# Useful decision probabilities that do not require choosing a WTP threshold.
psa_decision_probabilities <- data.frame(
  quantity = c(
    "Intervention more effective",
    "Intervention direct-cost saving",
    "Intervention societal-cost saving",
    "Intervention dominant - direct costs",
    "Intervention dominant - societal costs",
    "Intervention dominated - direct costs",
    "Intervention dominated - societal costs"
  ),
  probability = c(
    mean(psa_decision$delta_qaly > 0),
    mean(psa_decision$delta_cost_dir < 0),
    mean(psa_decision$delta_cost_soc < 0),
    mean(psa_decision$delta_qaly > 0 & psa_decision$delta_cost_dir < 0),
    mean(psa_decision$delta_qaly > 0 & psa_decision$delta_cost_soc < 0),
    mean(psa_decision$delta_qaly < 0 & psa_decision$delta_cost_dir > 0),
    mean(psa_decision$delta_qaly < 0 & psa_decision$delta_cost_soc > 0)
  )
)

cat("\nPSA decision probabilities:\n")
print(psa_decision_probabilities, digits = 5, row.names = FALSE)

# -----------------------------------------------------------------------------
# 12. CEAC and EVPI from incremental net monetary benefit
# -----------------------------------------------------------------------------

# WTP is in EUR per QALY. Adjust the grid to the decision context as needed.
WTP_GRID <- seq(0, 50000, by = 250)

prob_ce <- function(wtp, delta_effect, delta_cost) {
  mean(wtp * delta_effect - delta_cost > 0, na.rm = TRUE)
}

evpi_two_strategy <- function(wtp, qaly_int, qaly_uc, cost_int, cost_uc) {
  nmb_int <- wtp * qaly_int - cost_int
  nmb_uc  <- wtp * qaly_uc  - cost_uc

  mean(pmax(nmb_int, nmb_uc), na.rm = TRUE) -
    max(mean(nmb_int, na.rm = TRUE), mean(nmb_uc, na.rm = TRUE))
}

ceac <- data.frame(
  wtp = WTP_GRID,
  prob_ce_direct = vapply(
    WTP_GRID,
    prob_ce,
    numeric(1),
    delta_effect = psa_decision$delta_qaly,
    delta_cost   = psa_decision$delta_cost_dir
  ),
  prob_ce_societal = vapply(
    WTP_GRID,
    prob_ce,
    numeric(1),
    delta_effect = psa_decision$delta_qaly,
    delta_cost   = psa_decision$delta_cost_soc
  ),
  evpi_direct = vapply(
    WTP_GRID,
    evpi_two_strategy,
    numeric(1),
    qaly_int = psa_decision$qaly_int,
    qaly_uc  = psa_decision$qaly_uc,
    cost_int = psa_decision$cost_dir_int,
    cost_uc  = psa_decision$cost_dir_uc
  ),
  evpi_societal = vapply(
    WTP_GRID,
    evpi_two_strategy,
    numeric(1),
    qaly_int = psa_decision$qaly_int,
    qaly_uc  = psa_decision$qaly_uc,
    cost_int = psa_decision$cost_soc_int,
    cost_uc  = psa_decision$cost_soc_uc
  )
)

# Selected WTP values for quick inspection
WTP_REPORT <- c(0, 1000, 5000, 10000, 20000, 30000, 50000)
ceac_report <- ceac[match(WTP_REPORT, ceac$wtp), , drop = FALSE]
cat("\nCEAC / EVPI at selected WTP values:\n")
print(ceac_report, digits = 6, row.names = FALSE)

# -----------------------------------------------------------------------------
# 13. Plots
# -----------------------------------------------------------------------------

# Cost-effectiveness plane: direct-cost perspective
plot(
  psa_decision$delta_qaly,
  psa_decision$delta_cost_dir,
  pch = 16,
  cex = 0.55,
  xlab = "Incremental QALYs (intervention - usual care)",
  ylab = "Incremental direct costs (EUR)",
  main = "PSA: cost-effectiveness plane - direct costs"
)
abline(h = 0, v = 0, lty = 2)
points(
  qaly_int - qaly_uc,
  cost_dir_int - cost_dir_uc,
  pch = 4,
  cex = 1.4,
  lwd = 2
)

# Cost-effectiveness plane: societal perspective
plot(
  psa_decision$delta_qaly,
  psa_decision$delta_cost_soc,
  pch = 16,
  cex = 0.55,
  xlab = "Incremental QALYs (intervention - usual care)",
  ylab = "Incremental societal costs (EUR)",
  main = "PSA: cost-effectiveness plane - societal costs"
)
abline(h = 0, v = 0, lty = 2)
points(
  qaly_int - qaly_uc,
  cost_soc_int - cost_soc_uc,
  pch = 4,
  cex = 1.4,
  lwd = 2
)

# CEACs
matplot(
  ceac$wtp,
  ceac[, c("prob_ce_direct", "prob_ce_societal")],
  type = "l",
  lty = 1,
  lwd = 2,
  ylim = c(0, 1),
  xlab = "Willingness-to-pay (EUR/QALY)",
  ylab = "Probability intervention is cost-effective",
  main = "Cost-effectiveness acceptability curves"
)
legend(
  "bottomright",
  legend = c("Direct-cost perspective", "Societal perspective"),
  lty = 1,
  lwd = 2,
  bty = "n"
)

# EVPI per patient
matplot(
  ceac$wtp,
  ceac[, c("evpi_direct", "evpi_societal")],
  type = "l",
  lty = 1,
  lwd = 2,
  xlab = "Willingness-to-pay (EUR/QALY)",
  ylab = "EVPI per patient (EUR)",
  main = "Expected value of perfect information"
)
legend(
  "topright",
  legend = c("Direct-cost perspective", "Societal perspective"),
  lty = 1,
  lwd = 2,
  bty = "n"
)

# -----------------------------------------------------------------------------
# 14. Optional exports
# -----------------------------------------------------------------------------
# Uncomment if desired.
# write.csv(psa_decision, "COPD_PSA_draws.csv", row.names = FALSE)
# write.csv(psa_summary, "COPD_PSA_summary.csv", row.names = FALSE)
# write.csv(ceac, "COPD_PSA_CEAC_EVPI.csv", row.names = FALSE)

