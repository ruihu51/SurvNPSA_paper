rm(list=ls())

load("data/prepped_data.Rdata")

source("estimate_nuisances.R")
source("estimate_obs_components.R")

W$neck.dissection <- NULL

screen.marginal.05 <- function(Y, X, family, obsWeights, id, ...) {
    p.vals <- apply(X, 2, function(col) summary(glm(Y ~ col, family=family))$coef[2,4])
    whichVariable <- p.vals <= .05
    if(sum(whichVariable) == 0) {
        whichVariable <- FALSE
        whichVariable[order(p.vals)[1:2]] <- TRUE
    }
    whichVariable
}
screen.marginal.10 <- function(Y, X, family, obsWeights, id, ...) {
    p.vals <- apply(X, 2, function(col) summary(glm(Y ~ col, family=family))$coef[2,4])
    whichVariable <- p.vals <= .10
    if(sum(whichVariable) == 0) {
        whichVariable <- FALSE
        whichVariable[order(p.vals)[1:2]] <- TRUE
    }
    whichVariable
}

g.SL.library <- c(lapply(c("SL.glm", "SL.step"), function(alg) c(alg, c("All", "screen.marginal.05", "screen.marginal.10"))),
                  lapply(c("SL.ranger", "SL.gam", "SL.earth", "SL.xgboost"), function(alg) c(alg, c("screen.marginal.05", "screen.marginal.10"))),
                  list(c("SL.mean", "All")))

survSL.pchSL1 <- function(...) {
    survSL.pchSL(breaks = 1, SL.library = g.SL.library, ...)
}

survSL.pchSL2 <- function(...) {
    survSL.pchSL(breaks = 2, SL.library = g.SL.library, ...)
}

survSL.pchSL3 <- function(...) {
    survSL.pchSL(breaks = 3, SL.library = g.SL.library, ...)
}

survSL.pchSL4 <- function(...) {
    survSL.pchSL(breaks = 4, SL.library = g.SL.library, ...)
}

survSL.pchSL5 <- function(...) {
    survSL.pchSL(breaks = 5, SL.library = g.SL.library, ...)
}

survSL.gam.cts9 <- function(time, event, X, newX, new.times, ...) {
    survSL.gam(time = time,
               event = event,
               X = X,
               newX = newX,
               new.times = new.times,
               cts.num = 9,
               ...)
}

surv.SL.library <- c(lapply(c("survSL.coxph", "survSL.loglogreg", "survSL.expreg", "survSL.weibreg"), function(alg) c(alg, c("All", "survscreen.marg", "survscreen.glmnet"))),
                     lapply(c("survSL.gam.cts9", "survSL.rfsrc"), function(alg) c(alg, c("survscreen.marg", "survscreen.glmnet"))),
                     lapply(c("survSL.km"), function(alg) c(alg, "All")))

set.seed(4262021)

time <- Y
event <- Delta
treat <- A.neck.dissection
confounders <- W
fit.times <- 1:120

cat("start estimating nuisances:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
result.END <- .get.nuisances.est(time = time,
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
                                                         save.nuis.fits = TRUE),
                                 verbose = TRUE)

folder_path <- "outputEND/result.END"

file_name <- paste0(folder_path, ".RData")

save(result.END, file=file_name)

cat("start estimating target:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
result.END <- .get.obs.comps(time=time, event=event, treat=treat, result=result.END,
                             psi.type = "hybrid",
                             verbose = TRUE)

save(result.END, file=file_name)

cat("Complete!")
