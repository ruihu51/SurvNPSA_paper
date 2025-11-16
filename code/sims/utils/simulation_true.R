#####################
# These scripts compute the true parameter values.
# The results are saved in `true.value.s3.final.RData`.
# Note that calculating theta.obs.true, g.gs.sq.true, and e1.gs.sq.true
# is computationally intensive; we recommend running these steps
# on the cluster in parallel for each time point.
###################

rm(list=ls())

expit <- function(x) 1 / (1 + exp(-x))

library(cubature)

true.value <- list()

alpha_w <- c(0.2, -0.1)
beta_w <- c(0.1, 0.2)
zeta_w <- c(0.3, -0.1)

alpha_u <- c(0.55, 0.5)
beta_u <- c(-0.5, -1.75)

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

s_t_uw <- function(t, a, u, w, beta_u, beta_w){
    return(exp(-lambda_uw(a, u, w, beta_u, beta_w)*t))
}

lambda_cens_w <- function(a, w, zeta_w){
    return(exp(-0.5 - 0.15*a - w%*%zeta_w))
}

G_t_w <- function(t, a, w, zeta_w){
    return(exp(-lambda_cens_w(a, w, zeta_w)*t))
}

pi_w <- function(w){
    return(p_aw(a=1,w)/p_w(w))
}

p_uw <- function(u, w){
    if (is.matrix(u)) {
        return(dbeta(w[,1], 2*u[,1], 1)*0.25)
    } else {
        return(dbeta(w[1], 2*u[1], 1)*0.25)
    }
}

p_w <- function(w){
    f <- function(u) dbeta(w[1], 2*u, 1)
    return(hcubature(f, lowerLimit = c(0), upperLimit = c(1), tol = 1e-04)$integral)
}

p_auw <- function(a, u, w){
    return(a*pi_uw(u, w, alpha_u, alpha_w)*p_uw(u,w)+(1-a)*(1-pi_uw(u, w, alpha_u, alpha_w))*p_uw(u,w))
}

p_aw <- function(a, w){
    f <- function(u) p_auw(a, u, w)
    return(hcubature(f, lowerLimit = c(0,-2), upperLimit = c(1,2), tol = 1e-04)$integral)
}

p_u_aw <- function(a, u, w){
    return(p_auw(a,u,w)/p_aw(a,w))
}

p_a <- function(a){
    f <- function(x) p_auw(a, x[1:2], x[3:4])
    return(hcubature(f, lowerLimit = c(0,-2,0,0),
                     upperLimit = c(1,2,1,1), tol = 1e-04)$integral)
}

p_uw_a <- function(a, u, w){
    return(p_auw(a, u, w)/p_a(a))
}

p_w_a <- function(a, w){
    return(p_aw(a, w)/p_a(a))
}

s_t_w <- function(t, a, w, beta_u, beta_w){
    f <- function(u) s_t_uw(t, a, u, w, beta_u, beta_w)*p_u_aw(a, u, w)
    return(hcubature(f, lowerLimit = c(0,-2), upperLimit = c(1,2), tol = 1e-04)$integral)
}

alpha <- function(a, u, w, alpha_u, alpha_w){
    return(a/pi_uw(u, w, alpha_u, alpha_w) - (1-a)/(1-pi_uw(u, w, alpha_u, alpha_w)))
}

alpha_s <- function(a, w){
    return(a/pi_w(w) - (1-a)/(1-pi_w(w)))
}

fit.times <- c(seq(0.1,0.9,by=0.1),1,1.5,2)
true.value$times <- fit.times

# theta.true
get.theta.true <- function(t){
    f <- function(x) (s_t_uw(t, a=1, x[1:2], x[3:4], beta_u, beta_w)
                      - s_t_uw(t, a=0, x[1:2], x[3:4], beta_u, beta_w))*p_uw(x[1:2], x[3:4])
    return(hcubature(f, lowerLimit = c(0,-2,0,0),
                     upperLimit = c(1,2,1,1), tol = 1e-04)$integral)
}

theta.true <- numeric(length(fit.times))
for (i in 1:length(fit.times)){
    # cat(i, "\n")
    theta.true[i] <- get.theta.true(fit.times[i])
}
true.value$theta.true <- theta.true

# theta.obs.true - computationally intensive; recommend running on the cluster for better performance
get.theta.obs.true <- function(t){
    f <- function(w) (s_t_w(t, a=1, w, beta_u, beta_w)
                      - s_t_w(t, a=0, w, beta_u, beta_w))*p_w(w)
    return(hcubature(f, lowerLimit = c(0,0),
                     upperLimit = c(1,1), tol = 1e-04)$integral)
}

theta.obs.true <- numeric(length(fit.times))
for (i in 1:length(fit.times)){
    # cat(i, "\n")
    theta.obs.true[i] <- get.theta.obs.true(fit.times[i])
}
true.value$theta.obs.true.1 <- theta.obs.true

# E(g-gs)^2 - computationally intensive; recommend running on the cluster for better performance
get.g.gs.sq.true <- function(t){
    f1 <- function(x) (s_t_uw(t, a=1, x[1:2], x[3:4], beta_u, beta_w) -
                           s_t_w(t, a=1, x[3:4], beta_u, beta_w))^2*p_uw_a(a=1, x[1:2], x[3:4])
    f0 <- function(x) (s_t_uw(t, a=0, x[1:2], x[3:4], beta_u, beta_w) -
                           s_t_w(t, a=0, x[3:4], beta_u, beta_w))^2*p_uw_a(a=0, x[1:2], x[3:4])
    return(hcubature(f1, lowerLimit = c(0,-2,0,0),
                     upperLimit = c(1,2,1,1), tol = 1e-03)$integral*p_a(1) +
               hcubature(f0, lowerLimit = c(0,-2,0,0),
                         upperLimit = c(1,2,1,1), tol = 1e-03)$integral*p_a(0))
}

g.gs.sq.true <- numeric(length(fit.times))
for (i in 1:length(fit.times)){
    # cat(i, "\n")
    g.gs.sq.true[i] <- get.g.gs.sq.true(fit.times[i])
}
true.value$g.gs.sq.true <- g.gs.sq.true

# E(alpha-alpha_s)^2
f1 <- function(x) (alpha(a=1, x[1:2], x[3:4], alpha_u, alpha_w) -
                       alpha_s(a=1, x[3:4]))^2*p_uw_a(a=1, x[1:2], x[3:4])
f0 <- function(x) (alpha(a=0, x[1:2], x[3:4], alpha_u, alpha_w) -
                       alpha_s(a=0, x[3:4]))^2*p_uw_a(a=0, x[1:2], x[3:4])
(a.as.sq.true <- hcubature(f1, lowerLimit = c(0,-2,0,0),
                           upperLimit = c(1,2,1,1), tol = 1e-04)$integral*p_a(1) +
        hcubature(f0, lowerLimit = c(0,-2,0,0),
                  upperLimit = c(1,2,1,1), tol = 1e-04)$integral*p_a(0))
true.value$a.as.sq.true <- a.as.sq.true

# E alpha_s^2 = tau
f1 <- function(w) alpha_s(a=1, w, alpha_u, alpha_w)^2*p_w_a(a=1, w)
f0 <- function(w) alpha_s(a=0, w, alpha_u, alpha_w)^2*p_w_a(a=0, w)
(as.sq.true <- hcubature(f1, lowerLimit = c(0,0),
                         upperLimit = c(1,1), tol = 1e-04)$integral*p_a(1) +
        hcubature(f0, lowerLimit = c(0,0),
                  upperLimit = c(1,1), tol = 1e-04)$integral*p_a(0))
true.value$as.sq.true <- as.sq.true

# E{1(T>t)-S(t|A,W)}^2 = E[S(t|A,W)*{1-S(t|A,W)}] = psi- computationally intensive; recommend running on the cluster for better performance
get.1.gs.sq.true <- function(t){
    f1 <- function(w) s_t_w(t, a=1, w, beta_u, beta_w)*(1-s_t_w(t, a=1, w,
                                                                beta_u, beta_w))*p_w_a(a=1, w)
    f0 <- function(w) s_t_w(t, a=0, w, beta_u, beta_w)*(1-s_t_w(t, a=0, w,
                                                                beta_u, beta_w))*p_w_a(a=0, w)
    return(hcubature(f1, lowerLimit = c(0,0),
                     upperLimit = c(1,1), tol = 1e-04)$integral*p_a(1) +
               hcubature(f0, lowerLimit = c(0,0),
                         upperLimit = c(1,1), tol = 1e-04)$integral*p_a(0))
}

e1.gs.sq.true <- numeric(length(fit.times))
for (i in 1:length(fit.times)){
    # cat(i, "\n")
    e1.gs.sq.true[i] <- get.1.gs.sq.true(fit.times[i])
}
true.value$e1.gs.sq.true <- e1.gs.sq.true

# average censoring rate
get.censor.rate <- function(t){
    f1 <- function(w) (1-G_t_w(t, a=1, w, zeta_w))*p_w_a(a=1, w)
    f0 <- function(w) (1-G_t_w(t, a=0, w, zeta_w))*p_w_a(a=0, w)
    return(hcubature(f1, lowerLimit = c(0,0),
                     upperLimit = c(1,1), tol = 1e-03)$integral*p_a(1) +
               hcubature(f0, lowerLimit = c(0,0),
                         upperLimit = c(1,1), tol = 1e-03)$integral*p_a(0))
}
true.value$censor.rate.true <- sapply(c(0.5,1,1.5,2), get.censor.rate)






