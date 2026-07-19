rm(list=ls())

# Prepare the public RHC example data.
if (file.exists("data/rhc_example.RData")) {
    load("data/rhc_example.RData")
} else {
    if (!requireNamespace("Hmisc", quietly = TRUE)) {
        stop("Please install the Hmisc package to download the public RHC data.")
    }

    Hmisc::getHdata(rhc)

    vars <- c("age", "sex", "race", "meanbp1", "hrt1", "resp1",
              "temp1", "wblc1", "alb1", "crea1", "pafi1")

    time <- ifelse(rhc$death == "Yes",
                   rhc$dthdte - rhc$sadmdte,
                   rhc$lstctdte - rhc$sadmdte)
    event <- as.numeric(rhc$death == "Yes")
    treat <- as.numeric(rhc$swang1 == "RHC")

    keep <- is.finite(time) & !is.na(event) & !is.na(treat) &
        stats::complete.cases(rhc[, vars])

    set.seed(20260719)
    idx <- sort(sample(which(keep), 1000))

    rhc.example <- data.frame(time = as.numeric(time[idx]),
                              event = event[idx],
                              treat = treat[idx],
                              rhc[idx, vars],
                              check.names = FALSE)
    rhc.example$sex <- factor(rhc.example$sex)
    rhc.example$race <- factor(rhc.example$race)

    dir.create("data", showWarnings = FALSE, recursive = TRUE)
    save(rhc.example, file = "data/rhc_example.RData")
}
# Example data summary: n = 1000, RHC treated = 391, deaths = 669.

source("../utils/estimate_nuisances.R")
source("../utils/estimate_obs_components.R")

g.SL.library <- list(c("SL.gam", "All"),
                     c("SL.mean", "All"))

surv.SL.library <- list(c("survSL.km", "All"),
                        c("survSL.coxph", "All"))

set.seed(4262021)

time <- rhc.example$time
event <- rhc.example$event
treat <- rhc.example$treat
confounders <- rhc.example[, !(names(rhc.example) %in% c("time", "event", "treat"))]
fit.times <- 1:30

cat("start estimating nuisances:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
result.RHC <- .get.nuisances.est(time = time,
                                 event = event,
                                 treat = treat,
                                 confounders = confounders,
                                 fit.times = fit.times,
                                 nuisance.options = list(prop.SL.library = g.SL.library,
                                                         cens.SL.library = surv.SL.library,
                                                         event.SL.library = surv.SL.library,
                                                         cross.fit = TRUE,
                                                         V = 5,
                                                         survSL.control=list(initWeightAlg = "survSL.coxph"),
                                                         survSL.cvControl = list(V = 5),
                                                         save.nuis.fits = FALSE),
                                 verbose = TRUE)

folder_path <- "outputRHC/result.RHC"

dir.create("outputRHC", showWarnings = FALSE, recursive = TRUE)

file_name <- paste0(folder_path, ".RData")

save(result.RHC, file=file_name)

cat("start estimating target:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
result.RHC <- .get.obs.comps(time=time, event=event, treat=treat, result=result.RHC,
                             psi.type = "hybrid",
                             verbose = TRUE)

save(result.RHC, file=file_name)

cat("Complete!")
