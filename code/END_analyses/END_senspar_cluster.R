rm(list=ls())

args = commandArgs(TRUE)
task.id = as.numeric(args[1])

j <- task.id
d <- as.numeric(args[2])
SL.version <- as.numeric(args[3])

cat("d =", d, " j =", j, " Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

source("utils/estimate_nuisances.R")
source("utils/estimate_obs_components.R")

library(cubature)

.get.S.hat.int.vals <- function(t, S.hat.obs, tol, result.sim.drop){
    # step function
    theta.obs.func <- stepfun(c(0, result.sim.drop$nuisance$eval.times[-length(result.sim.drop$nuisance$eval.times)]),
                              c(1, S.hat.obs), right = FALSE)
    # integration
    hcubature(theta.obs.func, lowerLimit = c(0),
              upperLimit = c(t), tol=tol)$integral
}

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

if (SL.version==1){
    g.SL.library <- list(c("SL.gam", "All"),
                         c("SL.mean", "All"))

    surv.SL.library <- list(c("survSL.km", "All"),
                            c("survSL.coxph", "All"))
}


set.seed(072825)

load("data/drop.index.END.RData")
drop.index <- drop.index.list[[d]]

# Load the END analytic dataset. This file is not included in the
# Supplementary Materials because the END data are subject to a
# confidentiality/data use agreement. The loaded object should contain
# Y, Delta, A.neck.dissection, and W.
load("path/to/END_analytic_data.RData")
W$neck.dissection <- NULL

time <- Y
event <- Delta
treat <- A.neck.dissection
confounders <- W
fit.times <- 1:120

load("outputEND/result.END.RData")

theta.obs = result.END$obs.comps.df$theta.obs
psi = result.END$obs.comps.df$psi
tau = result.END$tau
S.hat.obs = result.END$nuisance$event.pred
g.hat.obs = result.END$nuisance$prop.pred
eval.times <- result.END$nuisance$eval.times

fit.times.rmst = result.END$fit.times.rmst
gamma = result.END$gamma.est
max_gap = 1
tol = 0.01

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
n_var <- length(var_list)

if (is.matrix(drop.index)){
    confounders.drop <- confounders[,-unlist(var_list[drop.index[,j]])]
} else {
    confounders.drop <- confounders[,-(var_list[[drop.index[j]]])]
}

cat(names(confounders)[!(names(confounders) %in% names(confounders.drop))], "\n")

cat("start estimating nuisances:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
result.sim.drop <- .get.nuisances.est(time = time,
                                      event = event,
                                      treat = treat,
                                      confounders = confounders.drop,
                                      fit.times = fit.times,
                                      nuisance.options = list(prop.SL.library = g.SL.library,
                                                              cens.SL.library = surv.SL.library,
                                                              event.SL.library = surv.SL.library,
                                                              cross.fit = TRUE,
                                                              V = 5,
                                                              eval.times = eval.times,
                                                              survSL.control=list(initWeightAlg = "survSL.coxph"),
                                                              survSL.cvControl = list(V = 5),
                                                              save.nuis.fits = FALSE),
                                      verbose = TRUE)
cat("start estimating target:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
result.sim.drop <- .get.obs.comps(time=time, event=event, treat=treat,
                                  result=result.sim.drop,
                                  psi.type = "hybrid",
                                  verbose = TRUE)

eval.idx <- sapply(result.sim.drop$fit.times, function(ft) {
    which.min(ifelse(result.sim.drop$nuisance$eval.times >= ft,
                     result.sim.drop$nuisance$eval.times, Inf))
})

fit.idx <- sapply(result.sim.drop$fit.times, function(ft) {
    which.min(ifelse(fit.times >= ft,
                     fit.times, Inf))
})

S.hat.obs.drop <- result.sim.drop$nuisance$event.pred

g.hat.obs.drop <- result.sim.drop$nuisance$prop.pred

alpha.drop <- (treat - g.hat.obs.drop)/(g.hat.obs.drop*(1 - g.hat.obs.drop))
alpha.obs <- (treat - g.hat.obs)/(g.hat.obs*(1 - g.hat.obs))

theta.change.matrix <- result.sim.drop$obs.comps.df$theta.obs - theta.obs

V.g.matrix.psi <- colMeans((S.hat.obs - S.hat.obs.drop)^2)
V.a.vector <- mean((alpha.obs - alpha.drop)^2)
tau.drop <- ifelse(result.sim.drop$tau>0, result.sim.drop$tau, mean(alpha.drop^2))


Gain.out.matrix = pmax(0, V.g.matrix.psi[eval.idx] / psi[fit.idx]) # 1*t
Gain.trt.vector = pmax(0, V.a.vector / tau.drop) # 1*1

folder_path <- "output.senspar/Fold"

dir.create("output.senspar", showWarnings = FALSE, recursive = TRUE)

# save
file_name.theta.change <- paste0(folder_path, d, ".", j, ".theta.change.RData")
file_name.V.g.psi <- paste0(folder_path, d, ".", j, ".V.g.psi.RData")
file_name.V.a <- paste0(folder_path, d, ".", j,".V.a.RData")
file_name.V.h.gamma <- paste0(folder_path, d, ".", j, ".V.h.gamma.RData")
file_name.Gain.out <- paste0(folder_path, d, ".", j,".Gain.out.RData")
file_name.Gain.trt <- paste0(folder_path, d, ".", j,".Gain.trt.RData")
file_name.Gain.out.phi <- paste0(folder_path, d, ".", j,".Gain.out.phi.RData")

save(theta.change.matrix, file=file_name.theta.change)
V.g.matrix.psi <- V.g.matrix.psi[eval.idx]
save(V.g.matrix.psi, file=file_name.V.g.psi)
rm(V.g.matrix.psi)
save(V.a.vector, file=file_name.V.a)

save(Gain.out.matrix, file=file_name.Gain.out)
save(Gain.trt.vector, file=file_name.Gain.trt)

cat("RMST", "\n")
h.hat.obs <- t(apply(S.hat.obs, 1, function(row) {
    sapply(fit.times.rmst, function(t) .get.S.hat.int.vals(t, row, tol = tol, result.sim.drop = result.sim.drop))
}))
h.hat.obs.drop <- t(apply(S.hat.obs.drop, 1, function(row) {
    sapply(fit.times.rmst, function(t) .get.S.hat.int.vals(t, row, tol = tol, result.sim.drop= result.sim.drop))
}))

V.h.matrix.gamma <- colMeans((h.hat.obs - h.hat.obs.drop)^2)

Gain.out.phi.matrix = pmax(0, V.h.matrix.gamma / gamma) # 1*t


save(V.h.matrix.gamma, file=file_name.V.h.gamma)
save(Gain.out.phi.matrix, file=file_name.Gain.out.phi)

cat("Complete!")
