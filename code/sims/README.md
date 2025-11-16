# Numerical experiments in the `SurvNPSA` paper

This directory contains the R code to reproduce the numerical experiments for the paper "Nonparametric Sensitivity Analysis for Unobserved Confounding with Survival Outcomes".

------------------------------------------------------------------------

## Directory structure

-   `sim_main.R`\
    Main script for running the simulations.

-   `sim_results.R`\
    Scripts used to generate `Fig 1` and `Fig 2` in Section 5 of the paper.

-   `compute_data/`\
    Intermediate objects such as seed data and true target parameters.

-   `utils/`\
    Helper functions used in the simulations (estimation, inference, computing true target parameters, etc.). For the most up-to-date versions of these functions, see the `SurvNPSA` R package.

------------------------------------------------------------------------

To reproduce the simulation results, run the following code in sequence:

1.  Run the simulation
    -   To run a single experiment (repetition j = 501, sample size n = 1000), use `Rscript sim_main.R 501 1000`

    -   You can also submit experiments to a Slurm cluster using `sbatch —array 501-1500%100 submit_sim.sh 1000`
2.  Process and visualize the results
    -   Simulation outputs from the previous step will be saved in a directory named `output.paper/`

    -   In `sim_results.R,` you can integrate all simulation outputs into a single data frame and generate the figures used in the paper.
