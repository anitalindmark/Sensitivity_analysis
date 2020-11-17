# Sensitivity_analysis
Code for simulations and empirical example for the paper "Sensitivity analysis to unobserved confounding in mediation analysis allowing for effect modification, censoring and truncation".

The **simulations** folder contains the code for the simulation study. It consists of:

- The subfolder **functions for censoring and truncation**, which contains functions to perform the modified maximum likelihood estimation of regression parameters under censoring or truncation ("Modified_ML_censored_outcome.R" and "Modified_ML_truncated_outcome.R") as well as implementations of the log-likelihoods and gradients ("loglik_and_gradient_censored_outcome.R" and "loglik_and_gradient_truncated_outcome.R").
- The subfolder **mediator-outcome** which contains the code for the simulations with mediator-outcome confounding, one code file for each Scenario a-e.
- The subfolder **exposure-mediator** which contains the code for the simulations with exposure-mediator confounding, one code file for each Scenario a, d and e.
- The subfolder **exposure-outcome** which contains the code for the simulations with exposure-outcome confounding, one code file for each Scenario a, d and e.

The file "empirical_example.R" contains the code for the empirical example based on the UPBdata available through the R package medflex: https://cran.r-project.org/package=medflex 
