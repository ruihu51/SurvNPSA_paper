rm(list=ls())

library(cubature)
library(ggplot2)
library(dplyr)
library(R.utils)

load("outputEND/result.END.RData")
result <- result.END

source("../SurvNPSA/R/npsa_summary.R")
source("../SurvNPSA/R/npsa_utils.R")

set.seed(102025)
plot.times <- 1:60

# effect bounds estimates and inference
# under no unobserved confounding
bounds.df <- .report.bounds(plot.times, result, rmst = TRUE)
bounds.df$bounds.df %>%
    filter(times==12) %>%
    select(ptwise.trans.lower, theta.obs, ptwise.trans.upper) %>%
    mutate(across(everything(), ~ round(.x, 3)))
# paper - survival difference at one year post-diagnosis
# ptwise.trans.lower theta.obs ptwise.trans.upper
# 1              0.013     0.055              0.097

load("../data/app_rst/senspar.df.END.cluster.RData")
senspar.df <- senspar

# under different levels of unobserved confounding
bounds.df.sens <- .report.bounds(plot.times, result,
                                 sens.df.mean = senspar.df$sens.df.mean,
                                 num_drop = c(3, 8, 13),
                                 pct_drop = NULL,
                                 n_var = 16,
                                 rmst = TRUE,
                                 sens.rmst.df.mean = senspar.df$sens.rmst.df.mean)

bounds.df$bounds.df <- rbind(bounds.df$bounds.df, bounds.df.sens$bounds.df)

############################################################################
# Figure 3 (END) Estimated survival difference, effect bounds and inference
# under different levels of unobserved confounding
###########################################################################
library(scales)

fig1.END.bounds <- bounds.df$bounds.df %>%
    mutate(
        setting = case_when(
            d == 0  ~ "No unobserved confounding",
            d == 3  ~ "Weak confounding",
            d == 8  ~ "Moderate confounding",
            d == 13 ~ "Strong confounding"
        ),
        setting = factor(setting, levels = c(
            "No unobserved confounding","Weak confounding",
            "Moderate confounding","Strong confounding"
        ))
    ) %>%
    ggplot(aes(x = times/12)) +
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
    scale_y_continuous(
        breaks = breaks_width(0.1),
        minor_breaks = NULL
    ) +
    labs(linetype = "Type", color = "Type", x = "Years post-diagnosis",
         y = "Survival difference (END - No END)") +
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

ggsave(file="../figures/Fig1_END_bounds.eps", width = 260,
       height = 240, dpi=300, units="mm", device=cairo_ps, limitsize = FALSE, fig1.END.bounds)


####################
# Robustness Values
###################
set.seed(4262021)
rv.times <- c(12, 60)
out <- list(result = result.END, senspar.df = NULL, bounds.df = bounds.df)
# RV and MIRV at t=12 (1 year) and t=60 (5 year)
out$res.RV <- .report.RV(rv.times, result, rho=1, unif = TRUE, q.01 = 0, q.99 = 60)
summary(out$res.RV)
# Robustness Value Report
# ------------------------
#     t0 theta    RV  MIRV conf.level lower.b rho
# 12     0 0.032 0.008       0.95    TRUE   1
# 60     0 0.002 0.000       0.95      NA   1

# interpret RV
(sp.point <- 0.032*0.032/(1-0.032)) # 0.001057851
(sp.l.pw <- 0.008*0.008/(1-0.008)) # 6.451613e-05

### leave-one-out
senspar$sens.df %>%
    mutate(sens.par = C.Y.sq * C.A.sq) %>%
    filter(t==12 & d==1) %>%
    mutate(sig.point=sens.par>sp.point,
           sig.l.pw=sens.par>sp.l.pw)

### leave-d-out
senspar$sens.df.mean %>%
    filter(t==12) %>%
    mutate(sig.point=sens.par>sp.point,
           sig.l.pw=sens.par>sp.l.pw)

# d = 2 and 3 for single point
# d = 1 already explain away for CI

### leave-8-out
mean(senspar$sens.df %>%
         mutate(sens.par = C.Y.sq * C.A.sq) %>%
         filter(t==12 & d==8) %>%
         mutate(value = sens.par <= sp.point) %>%
         pull(value))

set.seed(4262021)
# URV over the 1st year
out$res.RV <- .report.RV(rv.times, result, rho=1, unif = TRUE, q.01 = 0, q.99 = 12)
# The p-value under no unobserved confounding is: 0.0149
out$res.RV$unif.RV
# 0.005575041
unif.RV <- round(out$res.RV$unif.RV,3)
(sp.unif <- unif.RV*unif.RV/(1-unif.RV)) # 3.62173e-05

###########################
# Supporting Information
##########################

###################################
# Fig D1 (leave-d-out confounding)
##################################
sens.df.q <- senspar$sens.df %>%
    mutate(sens.par = C.Y.sq * C.A.sq) %>%
    group_by(d, t) %>%
    summarise(
        Q25 = quantile(sens.par, probs = 0.25),
        Median = quantile(sens.par, probs = 0.50),
        Q75 = quantile(sens.par, probs = 0.75),
        Mean = mean(sens.par),
        .groups = 'drop'
    )

library(RColorBrewer)
data4plot <- tidyr::pivot_longer(sens.df.q, cols = c(Q25, Median, Q75, Mean),
                                 names_to = "Quantile", values_to = "Value")

new_rows <- expand.grid(d = c(1,3,8,13), t = 0, Quantile = unique(data4plot$Quantile)) %>%
    mutate(Value = 0)

library(ggthemes)

pal_d <- c("#5F0F40", "#9A031E", "#0F4C5C")
pal_q <- c("#FFBE0B", "#FF006E", "#3A86FF")

p1 <- data4plot %>%
    bind_rows(new_rows) %>%
    filter(Quantile == "Mean", d %in% c(3, 8, 13), t <= 60) %>%
    mutate(d = factor(d, levels = c(3, 8, 13))) %>%   # lock order to map colors/linetypes
    ggplot(aes(x = t/12, y = Value, group = d, color = d, linetype = d)) +
    geom_line(size = 0.9) +
    scale_color_manual(values = pal_d, name = "d") +
    scale_linetype_manual(values = c("solid", "dashed", "dotdash"), name = "d") +
    coord_cartesian(ylim = c(0, 0.04)) +
    labs(
        title = "Average leave-d-out confounding",
        x = "Years post-diagnosis",
        y = "Confounding level"
    ) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 12),
        panel.grid.minor = element_blank()
    )

p2 <- data4plot %>%
    bind_rows(new_rows) %>%
    filter(d == 8, t <= 60, Quantile != "Mean") %>%
    mutate(Quantile = factor(Quantile, levels = c("Q25", "Median", "Q75"))) %>%
    ggplot(aes(x = t/12, y = Value, group = Quantile, color = Quantile, linetype = Quantile)) +
    geom_line(size = 0.9) +
    scale_color_manual(values = pal_q, name = "Quantiles") +
    scale_linetype_manual(values = c("solid", "dashed", "dotdash"), name = "Quantiles") +
    coord_cartesian(ylim = c(0, 0.04)) +
    labs(
        title = "Quantiles of leave-8-out confounding",
        x = "Years post-diagnosis",
        y = "Confounding level"
    ) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 12),
        panel.grid.minor = element_blank()
    )

library(patchwork)
fig.senspar.plot.2 <- p1 + p2 + plot_layout(ncol = 2)

ggsave(file="../figures/Fig2_END_senspar.eps", width = 260,
       height = 120, dpi=300, units="mm", device=cairo_ps, limitsize = FALSE, fig.senspar.plot.2)

#################
# empirical rho
################
# go to senspar_result_END.R

# rho=0.52
#################
set.seed(4262021)
out$res.RV <- .report.RV(rv.times, result, rho=0.52, unif = TRUE, q.01 = 0, q.99 = 12)
summary(out$res.RV)
# Robustness Value Report
# ------------------------
#     t0 theta    RV  MIRV conf.level lower.b  rho
# 12     0 0.061 0.015       0.95    TRUE 0.52
# 60     0 0.004 0.000       0.95      NA 0.52

(sp.l.pw <- 0.015*0.015/(1-0.015)) # 0.0002284264

### leave-one-out
senspar$sens.df %>%
    mutate(sens.par = C.Y.sq * C.A.sq) %>%
    filter(t==12 & d==1) %>%
    mutate(sig.point=sens.par>sp.point,
           sig.l.pw=sens.par>sp.l.pw)


load("../data/END/prepped_data.Rdata")
W$neck.dissection <- NULL

set.seed(4262021)

time <- Y
event <- Delta
treat <- A.neck.dissection
confounders <- W
names(confounders)

var_list <- list(
    1,
    c(2, 3, 4, 5),
    6,
    c(7,8),
    c(9,10,11),
    12, 13, 14, 15,
    c(16,17),
    18,
    c(19,20,21,22,23),
    c(24,25),
    c(26,27),
    c(28,29),
    30
)

var_list[c(1,2,6)]
names(confounders)[1] #"age"
names(confounders)[2] #"surgery"
names(confounders)[12] #"high.t.stage"




