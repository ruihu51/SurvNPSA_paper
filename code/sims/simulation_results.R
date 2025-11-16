library(patchwork)
library(tidyr)
library(dplyr)
library(ggplot2)

# load simulation outputs
df_all <- data.frame()
folder_path <- "../clusters/simulation_paper_1105/output.paper/res."
for (n in c(500, 1000, 2500, 5000)){
    for (j in 501:1500){
        if (j%%200==0) cat(n, j, "\n")
        file_name <- paste0(folder_path, n, ".",  j, ".RData")
        tryCatch({
            load(file_name)
            df_all <- rbind(df_all, res_df)
        }, error = function(e) {
            message("An error occurred: ", e$message)
            message("An error occurred: ", file_name)
        })

    }
}
sim.s3.v3 <- df_all

# load true values
load("true.value.s3.final.RData")
# true bound for the true effect
(theta.l.true <- true.value$theta.obs.true.1 -
        sqrt(true.value$g.gs.sq.true*true.value$a.as.sq.true))
(theta.u.true <- true.value$theta.obs.true.1 +
        sqrt(true.value$g.gs.sq.true*true.value$a.as.sq.true))

# true_df for simulation result calculation
true_df <- data.frame(times = true.value$times)
true_df$theta.l.true <- theta.l.true
true_df$theta.u.true <- theta.u.true
true_df$psi.true <- true.value$e1.gs.sq.true
true_df$theta.obs.true <- true.value$theta.obs.true.1
true_df$tau.true <- true.value$as.sq.true

# combine simulation results and true parameters
sim.s3.v3$times <- as.numeric(as.character(sim.s3.v3$times))
true_df$times <- as.numeric(as.character(true_df$times))
df <- left_join(sim.s3.v3, true_df, by = "times")


check_time <- c(0.5,1,1.5,2)
##########################
# Figures in Section 5
#########################

#########################################
# Figure 1 (Bias and MSE) in Section 5
########################################
summary.2.plot <- df %>%
    mutate(bias.l = effect.lower - theta.l.true,
           bias.u = effect.upper - theta.u.true) %>%
    group_by(times, n) %>%
    summarise(mean.l=mean(bias.l, na.rm=TRUE),
              sd.l=sd(bias.l),
              mean.u=mean(bias.u, na.rm=TRUE),
              sd.u=sd(bias.u),
              cnt=n()) %>%
    mutate(ll.l=mean.l-1.96*sd.l/sqrt(cnt),
           ul.l=mean.l+1.96*sd.l/sqrt(cnt),
           ll.u=mean.u-1.96*sd.u/sqrt(cnt),
           ul.u=mean.u+1.96*sd.u/sqrt(cnt),
           ll.l.rootn=sqrt(n)*mean.l-1.96*sqrt(n)*sd.l/sqrt(cnt),
           ul.l.rootn=sqrt(n)*mean.l+1.96*sqrt(n)*sd.l/sqrt(cnt),
           ll.u.rootn=sqrt(n)*mean.u-1.96*sqrt(n)*sd.u/sqrt(cnt),
           ul.u.rootn=sqrt(n)*mean.u+1.96*sqrt(n)*sd.u/sqrt(cnt))

p2 <- summary.2.plot %>%
    filter(times %in% check_time) %>%
    ggplot(aes(x=factor(n))) +
    geom_line(aes(y = sqrt(n)*mean.l, color = "Lower", linetype = "Lower"), group = 1) +
    geom_point(aes(y = sqrt(n)*mean.l, color = "Lower"), size = 0.5) +
    geom_errorbar(aes(ymin = ll.l.rootn, ymax = ul.l.rootn, color = "Lower", linetype = "Lower"),
                  width = 0.5, size = 0.5) +
    geom_line(aes(y = sqrt(n)*mean.u, color = "Upper", linetype = "Upper"), group = 1) +
    geom_point(aes(y = sqrt(n)*mean.u, color = "Upper"), size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
    geom_errorbar(aes(ymin = ll.u.rootn, ymax = ul.u.rootn, color = "Upper", linetype = "Upper"),
                  width = 0.5, size = 0.5) +
    scale_color_manual(values = c("Lower" = "firebrick3", "Upper" = "deepskyblue4"),
                       labels = c("Lower", "Upper"),
                       name = "Estimators") +
    scale_linetype_manual(values = c("Lower" = "solid", "Upper" = "dashed"),
                          labels = c("Lower", "Upper"),
                          name = "Estimators") +
    facet_grid(cols = vars(times),
               labeller = labeller(times = function(x) paste0("t=", x))) +
    labs(x = "Sample Size (n)", y = expression(sqrt(n)%*%Bias)) +
    theme_bw() +
    ylim(-0.15, 0.2) +
    theme(legend.position = "none",
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12),
          legend.text = element_text(size = 12))

summary.3.plot <- df %>%
    mutate(bias.l.sq = (effect.lower - theta.l.true)^2,
           bias.u.sq = (effect.upper - theta.u.true)^2) %>%
    group_by(times, n) %>%
    summarise(mean.l=mean(bias.l.sq, na.rm=TRUE),
              sd.l=sd(bias.l.sq),
              mean.u=mean(bias.u.sq, na.rm=TRUE),
              sd.u=sd(bias.u.sq),
              cnt=n()) %>%
    mutate(ll.l=mean.l-1.96*sd.l/sqrt(cnt),
           ul.l=mean.l+1.96*sd.l/sqrt(cnt),
           ll.u=mean.u-1.96*sd.u/sqrt(cnt),
           ul.u=mean.u+1.96*sd.u/sqrt(cnt),
           ll.l.rootn=n*mean.l-1.96*n*sd.l/sqrt(cnt),
           ul.l.rootn=n*mean.l+1.96*n*sd.l/sqrt(cnt),
           ll.u.rootn=n*mean.u-1.96*n*sd.u/sqrt(cnt),
           ul.u.rootn=n*mean.u+1.96*n*sd.u/sqrt(cnt))

p3 <- summary.3.plot %>%
    filter(times %in% check_time) %>%
    ggplot(aes(x=factor(n))) +
    geom_line(aes(y = n*mean.l, color = "Lower", linetype = "Lower"), group = 1) +
    geom_point(aes(y = n*mean.l, color = "Lower"), size = 0.5) +
    geom_errorbar(aes(ymin = n*ll.l, ymax = n*ul.l, color = "Lower", linetype = "Lower"),
                  width = 0.5, size = 0.5) +
    geom_line(aes(y = n*mean.u, color = "Upper", linetype = "Upper"), group = 1) +
    geom_point(aes(y = n*mean.u, color = "Upper"), size = 0.5) +
    geom_errorbar(aes(ymin = n*ll.u, ymax = n*ul.u, color = "Upper", linetype = "Upper"),
                  width = 0.5, size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
    scale_color_manual(values = c("Lower" = "firebrick3", "Upper" = "deepskyblue4"),
                       labels = c("Lower", "Upper"),
                       name = "Estimators") +
    scale_linetype_manual(values = c("Lower" = "solid", "Upper" = "dashed"),
                          labels = c("Lower", "Upper"),
                          name = "Estimators") +
    facet_grid(cols = vars(times),
               labeller = labeller(times = function(x) paste0("t=", x))) +
    labs(x = "Sample Size (n)", y = expression(n%*%MSE)) +
    theme_bw() +
    theme(legend.position = "none",
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12),
          legend.text = element_text(size = 12))

fig.plot.bias <- p2/p3 + theme(legend.position = "bottom")

ggsave(file="../figures/Fig_sim_bias.eps", width = 260,
       height = 160, dpi=300, units="mm", device=cairo_ps, limitsize = FALSE, fig.plot.bias)


##############################################
# Figure 2 (Empirical coverage) in Section 5
#############################################

#########################
# Confidence intervals coverage
########################
summary.pw.1.plot <- df %>%
    mutate(l.coverage = ptwise.bounds.lower < theta.l.true,
           u.coverage = ptwise.bounds.upper > theta.u.true,
           b.coverage = (ptwise.bounds.lower < theta.l.true) & (ptwise.bounds.upper > theta.u.true),
           l.t.coverage = ptwise.trans.lower < theta.l.true,
           u.t.coverage = ptwise.trans.upper > theta.u.true,) %>%
    group_by(times, n) %>%
    summarise(mean.l = mean(l.coverage),
              mean.u = mean(u.coverage),
              mean.b = mean(b.coverage),
              sd.l=sd(l.coverage),
              sd.u=sd(u.coverage),
              sd.b=sd(b.coverage),
              cnt=n()) %>%
    mutate(ll.l=mean.l-1.96*sd.l/sqrt(cnt),
           ul.l=mean.l+1.96*sd.l/sqrt(cnt),
           ll.u=mean.u-1.96*sd.u/sqrt(cnt),
           ul.u=mean.u+1.96*sd.u/sqrt(cnt),
           ll.b=mean.b-1.96*sd.b/sqrt(cnt),
           ul.b=mean.b+1.96*sd.b/sqrt(cnt))

summary.pw.2.plot <- df %>%
    mutate(l.coverage = ptwise.bounds.lower < theta.l.true,
           u.coverage = ptwise.bounds.upper > theta.u.true,

           l.t.coverage = ptwise.trans.lower < theta.l.true,
           u.t.coverage = ptwise.trans.upper > theta.u.true,
           b.t.coverage = (ptwise.trans.lower < theta.l.true) & (ptwise.trans.upper > theta.u.true),) %>%
    group_by(times, n) %>%
    summarise(mean.l = mean(l.t.coverage),
              mean.u = mean(u.t.coverage),
              mean.b = mean(b.t.coverage),
              sd.l=sd(l.t.coverage),
              sd.u=sd(u.t.coverage),
              sd.b=sd(b.t.coverage),
              cnt=n()) %>%
    mutate(ll.l=mean.l-1.96*sd.l/sqrt(cnt),
           ul.l=mean.l+1.96*sd.l/sqrt(cnt),
           ll.u=mean.u-1.96*sd.u/sqrt(cnt),
           ul.u=mean.u+1.96*sd.u/sqrt(cnt),
           ll.b=mean.b-1.96*sd.b/sqrt(cnt),
           ul.b=mean.b+1.96*sd.b/sqrt(cnt))

summary.pw.plot <- bind_rows(summary.pw.1.plot %>%
                                 select(times, n, mean.b, ll.b, ul.b) %>%
                                 mutate(type="Pointwise"),
                             summary.pw.2.plot %>%
                                 select(times, n, mean.b, ll.b, ul.b) %>%
                                 mutate(type="Transformed pointwise"))

fig.plot.pwCI <- summary.pw.plot %>%
    filter(times %in% check_time) %>%
    ggplot(aes(x=factor(n))) +
    geom_line(aes(y = mean.b, color = type, linetype = type, group = interaction(times, type))) +
    geom_point(aes(y = mean.b, color = type), size = 0.5) +
    geom_errorbar(aes(ymin = ll.b, ymax = ul.b, color = type, linetype = type),
                  width = 0.5, size = 0.5) +
    geom_hline(yintercept = 0.95, linetype = "dotted", color = "black", size = 0.5) +
    facet_grid(cols = vars(times),
               labeller = labeller(times = function(x) paste0("t=", x))) +
    labs(x = "Sample Size (n)", y = "Empirical Coverages",
         color = "Type", linetype = "Type") +
    scale_color_manual(values = c("Pointwise" = "#008080", "Transformed pointwise" = "#ff7f0e")) +
    scale_linetype_manual(values = c("Pointwise" = "solid", "Transformed pointwise" = "dashed")) +
    theme_bw() +
    ylim(0.85, 1) +
    theme(legend.position = "bottom",
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12),
          legend.text = element_text(size = 12))

#########################
# Uniform bands coverage
########################
summary.u.1 <- df %>%
    mutate(l.unif.coverage = uniform.bounds.lower < theta.l.true,
           u.unif.coverage = uniform.bounds.upper > theta.u.true,
           b.unif.coverage = (uniform.bounds.lower < theta.l.true) & (uniform.bounds.upper > theta.u.true),
           l.unif.t.coverage = uniform.trans.lower < theta.l.true,
           u.unif.t.coverage = uniform.trans.upper > theta.u.true,
           b.unif.t.coverage = (uniform.trans.lower < theta.l.true) & (uniform.trans.upper > theta.u.true)) %>%
    group_by(n, j) %>%
    summarise(l.coverage = ifelse(mean(l.unif.coverage)==1, 1, 0),
              u.coverage = ifelse(mean(u.unif.coverage)==1, 1, 0),
              b.coverage = ifelse(mean(b.unif.coverage)==1, 1, 0),
              l.t.coverage = ifelse(mean(l.unif.t.coverage)==1, 1, 0),
              u.t.coverage = ifelse(mean(u.unif.t.coverage)==1, 1, 0),
              b.t.coverage = ifelse(mean(b.unif.t.coverage)==1, 1, 0))

summary.u.3.plot <- summary.u.1 %>%
    group_by(n) %>%
    summarise(mean.b = mean(b.coverage),
              sd.b=sd(b.coverage),
              cnt=n()) %>%
    mutate(ll.b=mean.b-1.96*sd.b/sqrt(cnt),
           ul.b=mean.b+1.96*sd.b/sqrt(cnt))

summary.u.4.plot <- summary.u.1 %>%
    group_by(n) %>%
    summarise(mean.b= mean(b.t.coverage),
              sd.b=sd(b.t.coverage),
              cnt=n()) %>%
    mutate(ll.b=mean.b-1.96*sd.b/sqrt(cnt),
           ul.b=mean.b+1.96*sd.b/sqrt(cnt)
    )

summary.u.plot <- bind_rows(summary.u.3.plot %>%
                                select(n, mean.b, ll.b, ul.b) %>%
                                mutate(type="Uniform"),
                            summary.u.4.plot %>%
                                select(n, mean.b, ll.b, ul.b) %>%
                                mutate(type="Transformed uniform"))

fig.plot.unifCI <- summary.u.plot %>%
    ggplot(aes(x=factor(n))) +
    geom_line(aes(y = mean.b, color = type, linetype = type, group = type)) +
    geom_point(aes(y = mean.b, color = type), size = 0.5) +
    geom_errorbar(aes(ymin = ll.b, ymax = ul.b, color = type, linetype = type),
                  width = 0.5, size = 0.5) +
    geom_hline(yintercept = 0.95, linetype = "dotted", color = "black", size = 0.5) +
    labs(x = "Sample Size (n)", y = "Empirical Coverages",
         color = "Type", linetype = "Type") +
    scale_color_manual(values = c("Uniform" = "#7851A9", "Transformed uniform" = "#5D432C")) +
    scale_linetype_manual(values = c("Uniform" = "solid", "Transformed uniform" = "dashed")) +
    theme_bw() +
    ylim(0.85, 1) +
    theme(legend.position = "bottom",
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12),
          legend.text = element_text(size = 12))


(fig.plot.pwCI / fig.plot.unifCI) +
    plot_layout(heights = c(1, 1))

layout <- "
AAAAAA
#BBBB#
"
fig.plot.CI <- (fig.plot.pwCI / fig.plot.unifCI) +
    plot_layout(design = layout)

ggsave(file="../figures/Fig_sim_CI.eps", width = 260,
       height = 160, dpi=300, units="mm", device=cairo_ps, limitsize = FALSE, fig.plot.CI)
