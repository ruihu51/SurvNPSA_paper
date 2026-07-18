rm(list=ls())

args = commandArgs(TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript sim_main.R <j> <n> [gam|superlearner]")
}
task.id = as.numeric(args[1])

j <- task.id
n <- as.numeric(args[2])
method <- ifelse(length(args) >= 3, tolower(args[3]), "gam")
if (!(method %in% c("gam", "superlearner"))) {
    stop("method must be either 'gam' or 'superlearner'.")
}

load("compute_data/seed.data.paper.RData")
seed <- as.integer(seed.data[seed.data$n == n & seed.data$j == j, "seed"])
if (length(seed) != 1 || is.na(seed)) {
    stop("No matching seed found. Use j in 501,...,1500 and n in 500, 1000, 2500, or 5000.")
}

set.seed(seed)
cat(n, j, seed, '\n')
cat("method:", method, "\n")

library(dplyr)

source("utils/estimate_nuisances.R")
source("utils/estimate_obs_components.R")
source("utils/npsa_utils.R")

expit <- function(x) 1 / (1 + exp(-x))

# generate data setting
alpha_w <- c(0.2, -0.1)
beta_w <- c(0.1, 0.2)
zeta_w <- c(0.3, -0.1)

alpha_u <- c(0.55, 0.5)
beta_u <- c(-0.5, -1.75)

load("compute_data/true.value.s3.final.RData")
sens.out.true <-  true.value$g.gs.sq.true/true.value$e1.gs.sq.true
sens.trt.true <- true.value$a.as.sq.true/true.value$as.sq.true

senspar <- data.frame(time = true.value$times)
senspar$sens.out.true <- sens.out.true

pi_uw <- function(u, w, alpha_u, alpha_w){
    return(expit(0.2 - w%*%alpha_w - u%*%alpha_u))
}

lambda_uw <- function(a, u, w, beta_u, beta_w){
    if (is.matrix(u)){
        return(exp(0.15 - 0.25*a - cbind(sqrt(w[,1]), w[,2])%*%beta_w
                   - cbind(sqrt(u[,1]), exp(-3+u[,2]/2))%*%beta_u))
    } else {
        return(exp(0.15 - 0.25*a - cbind(sqrt(w[1]), w[2])%*%beta_w
                   - cbind(sqrt(u[1]), exp(-3+u[2]/2))%*%beta_u))
    }
}

lambda_cens_w <- function(a, w, zeta_w){
    return(exp(-0.5 - 0.15*a - w%*%zeta_w))
}

sim.data.surv <- function(n) {
    U1 <- runif(n, 0, 1)
    U2 <- runif(n, -2, 2)
    W1 <- rbeta(n, 2*U1, 1)
    W2 <- runif(n, 0, 1)
    U <- cbind(U1, U2)
    W <- cbind(W1, W2)

    piuw <- pi_uw(u=U, w=W, alpha_u, alpha_w)
    A <- rbinom(n, 1, prob=piuw)

    lambda0 <- lambda_uw(0, U, W, beta_u, beta_w)
    lambda1 <- lambda_uw(1, U, W, beta_u, beta_w)
    T0 <- rexp(n, lambda0)
    T1 <- rexp(n, lambda1)

    lambda0_cens <- lambda_cens_w(0, W, zeta_w)
    lambda1_cens <- lambda_cens_w(1, W, zeta_w)
    C0 <- rexp(n, lambda0_cens)
    C1 <- rexp(n, lambda1_cens)

    T <- A * T1 + (1-A) * T0
    C <- A * C1 + (1-A) * C0

    Y <- pmin(T, C)
    D <- as.numeric(T<=C)

    list(W=W, U=U, A=A, Y=Y, D=D, T0=T0, T1=T1, T=T, lambda0=lambda0, lambda1=lambda1,
         C0=C0, C1=C1, C=C, lambda0_cens=lambda0_cens, lambda1_cens=lambda1_cens,
         piuw=piuw)
}

sim.data <- sim.data.surv(n)
time  <- sim.data$Y
event <- sim.data$D
treat <- sim.data$A
confounders <- sim.data$W
fit.times <- c(seq(0.1,0.9,0.1),1,1.5,2)

if (method == "gam") {
    # generalized additive logistic regression model
    g.SL.library <- list(c("SL.gam.custom", "All"))
    # generalized additive Cox regression model
    surv.SL.library <- list(c("survSL.gam.custom", "All"))
} else {
    SL.xgboost.new <- function(Y, X, newX, family, obsWeights, id,
                               nrounds = 100,
                               max_depth = 3,
                               learning_rate = 0.1,
                               subsample = 1,
                               colsample_bytree = 1,
                               min_child_weight = 1,
                               nthread = 1,
                               ...) {
        requireNamespace("xgboost")

        X_mat <- data.matrix(X)
        newX_mat <- data.matrix(newX)

        dtrain <- xgboost::xgb.DMatrix(
            data = X_mat,
            label = Y,
            weight = obsWeights
        )

        params <- list(
            objective = if (family$family == "binomial") "binary:logistic" else "reg:squarederror",
            max_depth = max_depth,
            learning_rate = learning_rate,
            subsample = subsample,
            colsample_bytree = colsample_bytree,
            min_child_weight = min_child_weight,
            nthread = nthread
        )

        fit <- xgboost::xgb.train(
            params = params,
            data = dtrain,
            nrounds = nrounds,
            verbose = 0
        )

        pred <- predict(fit, newdata = newX_mat)

        fit_obj <- list(object = fit)
        class(fit_obj) <- "SL.xgboost.new"

        out <- list(pred = pred, fit = fit_obj)
        class(out$fit) <- "SL.xgboost.new"
        out
    }

    predict.SL.xgboost.new <- function(object, newdata, ...) {
        newX_mat <- data.matrix(newdata)
        predict(object$object, newdata = newX_mat)
    }

    # Propensity score SuperLearner library used for the additional simulation reported in Web Appendix E.
    g.SL.library <- c(lapply(
        c("SL.glm", "SL.step", "SL.ranger", "SL.gam.custom",
          "SL.earth", "SL.xgboost.new", "SL.mean"),
        function(alg) c(alg, "All")
    ))

    # Survival SuperLearner library used for the additional simulation reported in Web Appendix E.
    survSL.gam.cts9 <- function(time, event, X, newX, new.times, ...) {
        survSL.gam(time = time,
                   event = event,
                   X = X,
                   newX = newX,
                   new.times = new.times,
                   cts.num = 9,
                   ...)
    }

    surv.SL.library <- c(lapply(
        c("survSL.coxph", "survSL.loglogreg", "survSL.expreg", "survSL.weibreg",
          "survSL.gam.cts9", "survSL.rfsrc", "survSL.km"),
        function(alg) c(alg, "All")
    ))
}

nuisance.options = list(prop.SL.library = g.SL.library,
                        cens.SL.library = surv.SL.library,
                        event.SL.library = surv.SL.library,
                        cross.fit = TRUE,
                        V = 5,
                        survSL.control=list(initWeightAlg = "survSL.coxph"),
                        survSL.cvControl = list(V = 5),
                        save.nuis.fits = FALSE)

# estimate nuisances functions
cat("start estimating nuisances:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
result.sim <- .get.nuisances.est(time = time,
                                 event = event,
                                 treat = treat,
                                 confounders = confounders,
                                 fit.times = fit.times,
                                 nuisance.options = nuisance.options,
                                 verbose = FALSE)

# estimate observed components
cat("start estimating target:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
result.sim <- .get.obs.comps(time=time, event=event, treat=treat, result=result.sim,
                             psi.type = "hybrid", tau.type = "hybrid",
                             verbose = FALSE)

join_df <- left_join(result.sim$obs.comps.df, senspar, by="time")
sens.out.true.input <- join_df$sens.out.true

# estimate causal effect bounds
effect.bounds <-.get.effect.bounds(fit.times = result.sim$fit.times,
                                   theta.obs = result.sim$obs.comps.df$theta.obs,
                                   psi = result.sim$obs.comps.df$psi,
                                   tau = result.sim$tau,
                                   sens.out = sens.out.true.input,
                                   sens.trt = sens.trt.true,
                                   rho = 1)

# obtain confidence intervals and bands
bounds.conf.int <- .bounds.confints(effect.bounds,
                                    psi = result.sim$obs.comps.df$psi, tau = result.sim$tau,
                                    IF.vals.theta.obs = result.sim$IF.vals.theta.obs,
                                    IF.vals.psi = result.sim$IF.vals.psi,
                                    IF.vals.tau = result.sim$IF.vals.tau,
                                    rho = 1,
                                    band.end.pts=c(0,Inf),
                                    conf.level=.95,
                                    boot=10000)

# combine results
t_len <- length(bounds.conf.int$times)
cols_names <- lapply(bounds.conf.int, function(x) if(length(x) == t_len) x else NULL)
cols_names <- cols_names[!sapply(cols_names, is.null)]

res_df <- data.frame(cols_names)

res_df$theta.obs <- result.sim$obs.comps.df$theta.obs
res_df$surv.1 <- result.sim$obs.comps.df$surv.1
res_df$surv.0 <- result.sim$obs.comps.df$surv.0
res_df$psi <- result.sim$obs.comps.df$psi
res_df$tau <- result.sim$tau

res_df$n <- n
res_df$j <- j
res_df$seed <- seed

row.names(res_df) <- NULL

# save outputs
output.dir <- ifelse(method == "gam", "output.paper", "output.paper.superlearner")
if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)
file_name <- file.path(output.dir, paste0("res.", n, ".", j, ".RData"))
save(res_df, file=file_name)

cat("Complete!")
