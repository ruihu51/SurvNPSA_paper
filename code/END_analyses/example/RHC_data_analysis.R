rm(list=ls())

library(ggplot2)
library(dplyr)

# This script gives a local example using public RHC data. It follows the same
# workflow as the END analysis scripts, but uses 1000 observations and
# fit.times = 1:30 so that the full example can be run locally.
# The numerical results from this example are not used in the paper.

# Load the fitted RHC example analysis object.
load("outputRHC/result.RHC.RData")
result <- result.RHC

source("../utils/npsa_summary.R")
source("../utils/npsa_utils.R")

set.seed(102025)
plot.times <- 1:30

# effect bounds estimates and inference
# under no unobserved confounding
bounds.df <- .report.bounds(plot.times, result, rmst = FALSE)
bounds.df$bounds.df %>%
    filter(times==30) %>%
    select(ptwise.trans.lower, theta.obs, ptwise.trans.upper) %>%
    mutate(across(everything(), ~ round(.x, 3)))
# Example output at 30 days under no unobserved confounding:
# ptwise.trans.lower theta.obs ptwise.trans.upper
#             -0.058     0.009              0.076

# Load the precomputed benchmarking summaries from the example.
# In the END analysis, this object is computed from many repetitions.
# Here, one leave-5-out repetition is used so that the workflow runs quickly.
load("data/senspar.df.RHC.cluster.RData")
senspar.df <- senspar

# under one level of observed confounding used as an example benchmark
bounds.df.sens <- .report.bounds(plot.times, result,
                                 sens.df.mean = senspar.df$sens.df.mean,
                                 num_drop = c(5),
                                 pct_drop = NULL,
                                 n_var = 11,
                                 rmst = FALSE)

bounds.df$bounds.df <- rbind(bounds.df$bounds.df, bounds.df.sens$bounds.df)

############################################################################
# Figure: RHC example estimated survival difference, effect bounds and inference
# under one example level of observed confounding
###########################################################################
# This follows the plotting structure used for the END application.
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

fig.RHC.bounds <- bounds.df$bounds.df %>%
    mutate(
        setting = case_when(
            d == 0  ~ "No unobserved confounding",
            d == 5  ~ "Example leave-5-out confounding"
        ),
        setting = factor(setting, levels = c(
            "No unobserved confounding","Example leave-5-out confounding"
        ))
    ) %>%
    ggplot(aes(x = times)) +
    geom_step(aes(y = theta.obs, linetype = "Observed Effect", color = "Observed Effect")) +
    geom_step(aes(y = effect.lower, linetype = "Effect Bounds", color = "Effect Bounds")) +
    geom_step(aes(y = effect.upper, linetype = "Effect Bounds", color = "Effect Bounds")) +
    geom_step(aes(y = ptwise.trans.lower, linetype = "Pointwise CI", color = "Pointwise CI")) +
    geom_step(aes(y = ptwise.trans.upper, linetype = "Pointwise CI", color = "Pointwise CI")) +
    geom_step(aes(y = uniform.trans.lower, linetype = "Uniform Bands", color = "Uniform Bands")) +
    geom_step(aes(y = uniform.trans.upper, linetype = "Uniform Bands", color = "Uniform Bands")) +
    geom_hline(yintercept = 0, color = "grey") +
    scale_color_manual(values = c(
        "Observed Effect" = "black",
        "Effect Bounds"   = "red",
        "Pointwise CI"    = "blue",
        "Uniform Bands"   = "brown"
    )) +
    scale_linetype_manual(values = c(
        "Observed Effect" = "solid",
        "Effect Bounds"   = "dashed",
        "Pointwise CI"    = "dotdash",
        "Uniform Bands"   = "longdash"
    )) +
    labs(linetype = "Type", color = "Type", x = "Days since ICU admission",
         y = "Survival difference (RHC - No RHC)") +
    theme_bw() +
    theme(
        legend.position = "bottom",
        text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.width = unit(0.8, "cm"),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "grey")
    ) +
    facet_wrap(~setting, scales = "free_y")

ggsave(file="figures/Fig_RHC_example_bounds.png", width = 260,
       height = 160, dpi=300, units="mm", limitsize = FALSE, fig.RHC.bounds)

####################
# Robustness Values
###################
# RV and MIRV at t = 30.
# In this small example, the pointwise confidence interval at t = 30 covers
# zero; the RV output is included to show the computation workflow.
set.seed(4262021)
rv.times <- c(30)
out <- list(result = result.RHC, senspar.df = senspar.df, bounds.df = bounds.df)
out$res.RV <- .report.RV(rv.times, result, rho=1, unif = TRUE, q.01 = 0, q.99 = 30)
summary(out$res.RV)
# Robustness Value Report
# ------------------------
# t0 theta    RV MIRV conf.level lower.b rho
# 30     0 0.009    0       0.95      NA   1
# The p-value under no unobserved confounding is 0.7758.

### leave-d-out
# Example benchmarking value at t = 30.
senspar$sens.df.mean %>%
    filter(t==30) %>%
    mutate(sens.par = round(sens.par, 4))
# d t  sens.par
# 5 30   0.0021
