# markov_model__microsim

This repository contains various files related to the COPD Markov model developed by Petra Menn and Rolf Holle in 2008/2009. It contains the following files:
- Markov_model_Menn_Benchmark_Github.r: The original model implemented in Base R. The top half of the codes implements the deterministic model, the bottom half the probabilistic model.
- Markov_model_Menn_Benchmark_Github_incl_checks.r: Basically the same, but with additional consistency checks.
- Markov_model_simplified_Github_incl_checks.r: The Base R implementation of the simplified Markov model where two states (LTx: Lung transplantation, and LVR: Lung volume reduction) were removed, in addition to mild exacerbations as events.
- Markov_model_simplified_heemod_with_PSA.R: The heemod implementation of the simplified Markov model.
- COPD_Markov_microsimulation_main_baseline.R: The simplified Markov model and the corresponding microsimulation model.
- COPD_microsimulation_AE_memory_extension.R: A proof-of-concept memory extension of the microsimulation model where acute axecerbations in one cycle still elevate the risk in subsequent cycles.

The remaining files show the output of technical model validations. (Check if the choice of the random seed is influential etc.)
