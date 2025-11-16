# Data analysis in the `SurvNPSA` paper

This directory contains the R code to reproduce the sensitivity analysis for the effect of Elective Neck Dissection (END) on mortality in the paper "Nonparametric Sensitivity Analysis for Unobserved Confounding with Survival Outcomes".

------------------------------------------------------------------------

## Directory structure

-   `END_data_analysis.R`\
    Main script for generating Fig 3 and all other analyses reported in Section 6 and Web Appendix D of the paper.

-   `END_estimate.R`\
    Scripts used to compute the cross-fitted one-step estimators of the observed components with nuisance functions estimated using SuperLearner and survSuperLearner.

-   `END_senspar_cluster.R`\
    Scripts used to compute the sensitivity parameters using benchmarking based on `d` observed covariates. utils/END_senspar_rst.R computes the product of the sensitivity parameters.

-   `utils/`\
    Additional helper functions. The most up-to-date versions of these functions are provided in the SurvNPSA R package.

------------------------------------------------------------------------

To reproduce Section 6 and Web Appendix D, run the following csteps in sequence:

1.  Compute the estimators.
    -   Run `END_estimate.R`
2.  Compute the product of sensitivity parameters for each `d`.
    -   To run a single repetition with j = 1 for the leave-8-out observed confounding (d=8), use `Rscript END_senspar_cluster.R 1 8 1`

        -   Since this step only uses observed confounding as a reference to interpret the plausibility of the robustness values, and need to be repeated 100 times for each `d`, we recommend use simple nuisance estimators such as `SL.gam` and `SL.coxph` for simplicity.

    -   We recommend submit the task to a Slurm cluster using `sbatch —array 1-100%100 submit_sim.sh 8 1`
3.  Load the results in Step 1 and 2. Process and visualize the results.
    -   Effect bounds estimates and inference.

    -   Robustness values and interpretation.
