V.g.matrix.psi.list <- list()
theta.change.matrix.list <- list()
Gain.out.matrix.list <- list()
V.g.matrix.list <- list()
Gain.trt.vector.list <- list()
V.a.vector.list <- list()

load("../data/drop.index.END.RData")

folder_path <- 'output.senspar/Fold'

t_len <- 120

for (v in c(1,2,3,8,13)){
    if (v %in% c(1)){
        J <- 16
    } else {
        J <- 100
    }
    V.g.matrix.psi.all <- matrix(NA, nrow=J, ncol=t_len)
    for (j in 1:J) {
        # if (j %% 10 == 0) cat(v, j, "\n")
        file_name.V.g.psi <- paste0(folder_path, v, ".", j, ".V.g.psi.RData")
        load_result <- try(load(file_name.V.g.psi), silent = TRUE)
        if (inherits(load_result, "try-error")) next  # skip to next j
        V.g.matrix.psi.all[j, ] <- V.g.matrix.psi
    }
    V.g.matrix.psi.all <- V.g.matrix.psi.all[complete.cases(V.g.matrix.psi.all), ]
    cat(v, nrow(V.g.matrix.psi.all), "\n")
    V.g.matrix.psi.list[[v]] <- V.g.matrix.psi.all
}

for (v in c(1,2,3,8,13)){
    if (v %in% c(1)){
        J <- 16
    } else {
        J <- 100
    }
    theta.change.matrix.all <- matrix(NA, nrow=J, ncol=t_len)
    for (j in 1:J){
        if (j%%10==0) cat(v, j, "\n")
        file_name.theta.change <- paste0(folder_path, v, ".", j, ".theta.change.RData")
        load_result <- try(load(file_name.theta.change), silent = TRUE)
        if (inherits(load_result, "try-error")) next  # skip to next j
        theta.change.matrix.all[j,] <- theta.change.matrix
    }
    theta.change.matrix.all <- theta.change.matrix.all[complete.cases(theta.change.matrix.all), ]
    theta.change.matrix.list[[v]] <- theta.change.matrix.all
}

for (v in c(1,2,3,8,13)){
    if (v %in% c(1)){
        J <- 16
    } else {
        J <- 100
    }
    Gain.out.matrix.all <- matrix(NA, nrow=J, ncol=t_len)
    for (j in 1:J){
        if (j%%10==0) cat(v, j, "\n")
        file_name.Gain.out <- paste0(folder_path, v, ".", j,".Gain.out.RData")
        load_result <- try(load(file_name.Gain.out), silent = TRUE)
        if (inherits(load_result, "try-error")) next  # skip to next j
        Gain.out.matrix.all[j,] <- Gain.out.matrix
    }
    Gain.out.matrix.all <- Gain.out.matrix.all[complete.cases(Gain.out.matrix.all), ]
    Gain.out.matrix.list[[v]] <- Gain.out.matrix.all
}

for (v in c(1,2,3,8,13)){
    if (v %in% c(1)){
        J <- 16
    } else {
        J <- 100
    }
    V.g.matrix.all <- matrix(NA, nrow=J, ncol=t_len)
    for (j in 1:J){
        if (j%%10==0) cat(v, j, "\n")
        file_name.V.g <- paste0(folder_path, v, ".", j,".V.g.psi.RData")
        load_result <- try(load(file_name.V.g), silent = TRUE)
        if (inherits(load_result, "try-error")) next  # skip to next j
        V.g.matrix.all[j,] <- V.g.matrix.psi
    }
    V.g.matrix.all <- V.g.matrix.all[complete.cases(V.g.matrix.all), ]
    V.g.matrix.list[[v]] <- V.g.matrix.all
}

for (v in c(1,2,3,8,13)){
    if (v %in% c(1)){
        J <- 16
    } else {
        J <- 100
    }
    V.a.vector.all <- numeric(J)
    for (j in 1:J){
        if (j%%10==0) cat(v, j, "\n")
        file_name.V.a <- paste0(folder_path, v, ".", j,".V.a.RData")
        load_result <- try(load(file_name.V.a), silent = TRUE)
        if (inherits(load_result, "try-error")) next  # skip to next j
        V.a.vector.all[j] <- V.a.vector
    }
    V.a.vector.all <- V.a.vector.all[!is.na(V.a.vector.all)]
    V.a.vector.list[[v]] <- V.a.vector.all
}

for (v in c(1,2,3,8,13)){
    if (v %in% c(1)){
        J <- 16
    } else {
        J <- 100
    }
    Gain.trt.vector.all <- numeric(J)
    for (j in 1:J){
        if (j%%10==0) cat(v, j, "\n")
        file_name.Gain.trt <- paste0(folder_path, v, ".", j,".Gain.trt.RData")
        load_result <- try(load(file_name.Gain.trt), silent = TRUE)
        if (inherits(load_result, "try-error")) next  # skip to next j
        Gain.trt.vector.all[j] <- Gain.trt.vector
    }
    Gain.trt.vector.all <- Gain.trt.vector.all[!is.na(Gain.trt.vector.all)]
    Gain.trt.vector.list[[v]] <- Gain.trt.vector.all
}

V.h.matrix.gamma.list <- list()
Gain.out.phi.matrix.list <- list()

for (v in c(1,2,3,8,13)){
    if (v %in% c(1)){
        J <- 16
    } else {
        J <- 100
    }
    V.h.matrix.gamma.all <- matrix(NA, nrow=J, ncol=2)
    for (j in 1:J){
        if (j%%10==0) cat(v, j, "\n")
        file_name.V.h.gamma <- paste0(folder_path, v, ".", j, ".V.h.gamma.RData")
        load_result <- try(load(file_name.V.h.gamma), silent = TRUE)
        if (inherits(load_result, "try-error")) next  # skip to next j
        V.h.matrix.gamma.all[j,] <- V.h.matrix.gamma
    }
    V.h.matrix.gamma.all <- V.h.matrix.gamma.all[complete.cases(V.h.matrix.gamma.all), ]
    V.h.matrix.gamma.list[[v]] <- V.h.matrix.gamma.all
}

for (v in c(1,2,3,8,13)){
    if (v %in% c(1)){
        J <- 16
    } else {
        J <- 100
    }
    Gain.out.phi.matrix.all <- matrix(NA, nrow=J, ncol=2)
    for (j in 1:J){
        if (j%%10==0) cat(v, j, "\n")
        file_name.Gain.out.phi <- paste0(folder_path, v, ".", j,".Gain.out.phi.RData")
        load_result <- try(load(file_name.Gain.out.phi), silent = TRUE)
        if (inherits(load_result, "try-error")) next  # skip to next j
        Gain.out.phi.matrix.all[j,] <- Gain.out.phi.matrix
    }
    Gain.out.phi.matrix.all <- Gain.out.phi.matrix.all[complete.cases(Gain.out.phi.matrix.all), ]
    Gain.out.phi.matrix.list[[v]] <- Gain.out.phi.matrix.all
}

###########################################
# Estimate s_{c,T}(t)*s_{c,A}/(1-s_{c,A})
#########################################
Gain.out.df <- data.frame(value = numeric(), t = numeric(), j = integer(), d = integer())

for(i in c(1,2,3,8,13)) {

    cur_mat <- Gain.out.matrix.list[[i]]

    df <- data.frame(C.Y.sq = as.vector(cur_mat),
                     t = rep(1:ncol(cur_mat), each = nrow(cur_mat)),
                     j = rep(1:nrow(cur_mat), ncol(cur_mat)),
                     d = rep(i, length(cur_mat)))

    Gain.out.df <- rbind(Gain.out.df, df)
}

Gain.trt.df <- do.call(rbind, lapply(c(1,2,3,8,13), function(i) {
    data.frame(C.A.sq = Gain.trt.vector.list[[i]] / (1 - Gain.trt.vector.list[[i]]),
               j = 1:length(Gain.trt.vector.list[[i]]),
               d = as.factor(i))
}))

sens.df <- merge(Gain.out.df, Gain.trt.df, by = c("j", "d"))

save(sens.df, file = "../data/senspar/sens.df.RData")

sens.df.mean <- sens.df %>%
    mutate(sens.par = C.Y.sq * C.A.sq) %>%
    group_by(d, t) %>%
    summarize(sens.par = mean(sens.par)) %>%
    ungroup()

save(sens.df.mean, file = "../data/senspar/sens.df.mean.RData")

## RMST
# load("../data/senspar/V.h.matrix.gamma.list.RData")
# load("../data/senspar/Gain.out.phi.matrix.list.RData")
Gain.out.phi.df <- data.frame(value = numeric(), t = numeric(), j = integer(), d = integer())
for(i in c(1,2,3,8,13)) {

    cur_mat <- Gain.out.phi.matrix.list[[i]]

    df <- data.frame(C.Y.sq = as.vector(cur_mat),
                     t = rep(1:ncol(cur_mat), each = nrow(cur_mat)),
                     j = rep(1:nrow(cur_mat), ncol(cur_mat)),
                     d = rep(i, length(cur_mat)))

    Gain.out.phi.df <- rbind(Gain.out.phi.df, df)
}

sens.rmst.df <- merge(Gain.out.phi.df, Gain.trt.df, by = c("j", "d"))

save(sens.rmst.df, file = "../data/senspar/sens.rmst.df.RData")

sens.rmst.df.mean <- sens.rmst.df %>%
    mutate(sens.par = C.Y.sq * C.A.sq) %>%
    group_by(d, t) %>%
    summarize(sens.par = mean(sens.par)) %>%
    ungroup()

save(sens.rmst.df.mean, file = "../data/senspar/sens.rmst.df.mean.RData")

senspar <- list()
senspar$sens.df <- sens.df
senspar$sens.df.mean <- sens.df.mean
senspar$sens.rmst.df <- sens.rmst.df
senspar$sens.rmst.df.mean <- sens.rmst.df.mean
save(senspar, file = "../data/app_rst/senspar.df.END.cluster.RData")

###################
# Empirical rho
###################
valid.matrix <- matrix(NA, nrow=8, ncol=120)
cor.t.list <- list()
for (v in c(1,2,3,8)){
    if (v %in% c(1)){
        J <- 16
    } else {
        J <- 100
    }
    cor.t <- matrix(NA, nrow=J, ncol=120)
    for (t in 1:120){
        V.g <- V.g.matrix.psi.list[[v]][,t]
        V.a <- V.a.vector.list[[v]]
        valid <- V.g > 0 & V.a > 0
        theta.change <- theta.change.matrix.list[[v]]
        cor.t[valid,t] <- (abs(theta.change[valid,t])/sqrt(V.g[valid]*V.a[valid]))*sign(theta.change[valid,t])
        valid.matrix[v,t] <- mean(valid)
    }
    cor.t.list[[v]] <- cor.t
}

round(range(colMeans(abs(cor.t.list[[8]][,1:12]))),2) # [0.34, 0.56]



