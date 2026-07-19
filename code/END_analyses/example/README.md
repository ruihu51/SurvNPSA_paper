# RHC example analysis

This folder contains a small example analysis using public RHC data. The purpose of this example is to show that the code workflow used in the END analysis can be run locally with example data. This example is not intended to reproduce the END application results in the paper.

The original END analytic dataset and fitted END analysis object are not included in the Supplementary Materials because they are based on data subject to a confidentiality/data use agreement. The scripts in the parent folder show how the END analysis objects and results were generated. The scripts in this folder use 1000 observations from the public RHC data and `fit.times = 1:30` so that the full workflow can be run locally.

The example data file was constructed from the publicly available RHC data from the Vanderbilt Biostatistics data repository (https://hbiostat.org/data/repo/rhc). These data were used in Connors et al. (1996) to study the effectiveness of right heart catheterization in the initial care of critically ill patients. The example file contains `time`, `event`, `treat`, and 11 baseline covariates. Here, `treat` indicates receipt of RHC, and `event` indicates death.

------------------------------------------------------------------------

## Directory structure

-   `RHC_estimate.R`\
    Script used to compute the cross-fitted one-step estimators of the observed components.

-   `RHC_senspar_cluster.R`\
    Script used to compute the sensitivity parameters using benchmarking based on `d` observed covariates.

-   `utils/RHC_senspar_rst.R`\
    Script used to combine the benchmarking repetition and compute the product of the sensitivity parameters.

-   `RHC_data_analysis.R`\
    Main script for loading the fitted RHC example analysis object and benchmarking summaries, then processing and visualizing the results.

-   `data/rhc_example.RData`\
    Example RHC data with 1000 observations.

-   `utils/`\
    Example script for combining the benchmarking repetition. Other helper functions are sourced from the parent `../utils/` folder.

------------------------------------------------------------------------

The paths in these scripts are relative to this directory. The full local example can be run as follows:

1.  Compute the estimators.
    -   Run `Rscript RHC_estimate.R`
    -   This step loads `data/rhc_example.RData` if it is already present. Otherwise, it downloads the public RHC data using `Hmisc::getHdata()` and constructs `data/rhc_example.RData`.
    -   This step computes `outputRHC/result.RHC.RData`, which is used by `RHC_data_analysis.R`.
2.  Compute the product of sensitivity parameters.
    -   To run one repetition with j = 1 for the leave-5-out observed confounding (d=5), use `Rscript RHC_senspar_cluster.R 1 5 1`
    -   When `SL.version = 1`, the propensity score library is `SL.gam` and `SL.mean`, and the survival and censoring libraries are `survSL.km` and `survSL.coxph`.
3.  Combine the benchmarking repetition.
    -   Run `Rscript utils/RHC_senspar_rst.R`
    -   This step computes `data/senspar.df.RHC.cluster.RData`, which is used by `RHC_data_analysis.R`.
    -   The intermediate combined files are saved in `data/senspar/`.
4.  Load the fitted RHC example analysis object and benchmarking summaries. Process and visualize the results.
    -   Run `Rscript RHC_data_analysis.R`
    -   This step generates example effect bounds, robustness values, and a figure in `figures/`.
