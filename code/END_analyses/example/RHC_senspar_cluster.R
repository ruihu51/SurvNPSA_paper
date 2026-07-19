rm(list=ls())

args = commandArgs(TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript RHC_senspar_cluster.R <j> <d> <SL.version>")
}
task.id = as.numeric(args[1])

j <- task.id
d <- as.numeric(args[2])
SL.version <- as.numeric(args[3])

cat("d =", d, " j =", j, " Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

source("../utils/estimate_nuisances.R")
source("../utils/estimate_obs_components.R")

g.SL.library <- list(c("SL.gam", "All"),
                     c("SL.mean", "All"))

surv.SL.library <- list(c("survSL.km", "All"),
                        c("survSL.coxph", "All"))

if (SL.version==1){
    g.SL.library <- list(c("SL.gam", "All"),
                         c("SL.mean", "All"))

    surv.SL.library <- list(c("survSL.km", "All"),
                            c("survSL.coxph", "All"))
}

set.seed(072825)

load("data/drop.index.RHC.RData")
drop.index <- drop.index.list[[d]]

load("data/rhc_example.RData")

time <- rhc.example$time
event <- rhc.example$event
treat <- rhc.example$treat
confounders <- rhc.example[, !(names(rhc.example) %in% c("time", "event", "treat"))]
fit.times <- 1:30

load("outputRHC/result.RHC.RData")

theta.obs = result.RHC$obs.comps.df$theta.obs
psi = result.RHC$obs.comps.df$psi
tau = result.RHC$tau
S.hat.obs = result.RHC$nuisance$event.pred
g.hat.obs = result.RHC$nuisance$prop.pred
eval.times <- result.RHC$nuisance$eval.times

var_list <- as.list(1:ncol(confounders))
n_var <- length(var_list)

if (is.matrix(drop.index)){
    confounders.drop <- confounders[,-unlist(var_list[drop.index[,j]])]
} else {
    confounders.drop <- confounders[,-(var_list[[drop.index[j]]])]
}

cat(names(confounders)[!(names(confounders) %in% names(confounders.drop))], "\n")
# For Rscript RHC_senspar_cluster.R 1 5 1, the dropped covariates are:
# age sex race meanbp1 hrt1

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
file_name.Gain.out <- paste0(folder_path, d, ".", j,".Gain.out.RData")
file_name.Gain.trt <- paste0(folder_path, d, ".", j,".Gain.trt.RData")

save(theta.change.matrix, file=file_name.theta.change)
V.g.matrix.psi <- V.g.matrix.psi[eval.idx]
save(V.g.matrix.psi, file=file_name.V.g.psi)
rm(V.g.matrix.psi)
save(V.a.vector, file=file_name.V.a)

save(Gain.out.matrix, file=file_name.Gain.out)
save(Gain.trt.vector, file=file_name.Gain.trt)

# The example run saves one leave-5-out benchmarking repetition in output.senspar/.
cat("Complete!")
