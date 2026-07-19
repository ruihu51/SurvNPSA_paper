library(dplyr)

V.g.matrix.psi.list <- list()
theta.change.matrix.list <- list()
Gain.out.matrix.list <- list()
V.g.matrix.list <- list()
Gain.trt.vector.list <- list()
V.a.vector.list <- list()

folder_path <- 'output.senspar/Fold'

t_len <- 30

for (v in c(5)){
    J <- 1
    V.g.matrix.psi.all <- matrix(NA, nrow=J, ncol=t_len)
    for (j in 1:J) {
        file_name.V.g.psi <- paste0(folder_path, v, ".", j, ".V.g.psi.RData")
        load_result <- try(load(file_name.V.g.psi), silent = TRUE)
        if (inherits(load_result, "try-error")) next
        V.g.matrix.psi.all[j, ] <- V.g.matrix.psi
    }
    V.g.matrix.psi.all <- V.g.matrix.psi.all[complete.cases(V.g.matrix.psi.all), , drop = FALSE]
    cat(v, nrow(V.g.matrix.psi.all), "\n")
    V.g.matrix.psi.list[[v]] <- V.g.matrix.psi.all
}

for (v in c(5)){
    J <- 1
    theta.change.matrix.all <- matrix(NA, nrow=J, ncol=t_len)
    for (j in 1:J){
        file_name.theta.change <- paste0(folder_path, v, ".", j, ".theta.change.RData")
        load_result <- try(load(file_name.theta.change), silent = TRUE)
        if (inherits(load_result, "try-error")) next
        theta.change.matrix.all[j,] <- theta.change.matrix
    }
    theta.change.matrix.all <- theta.change.matrix.all[complete.cases(theta.change.matrix.all), , drop = FALSE]
    theta.change.matrix.list[[v]] <- theta.change.matrix.all
}

for (v in c(5)){
    J <- 1
    Gain.out.matrix.all <- matrix(NA, nrow=J, ncol=t_len)
    for (j in 1:J){
        file_name.Gain.out <- paste0(folder_path, v, ".", j,".Gain.out.RData")
        load_result <- try(load(file_name.Gain.out), silent = TRUE)
        if (inherits(load_result, "try-error")) next
        Gain.out.matrix.all[j,] <- Gain.out.matrix
    }
    Gain.out.matrix.all <- Gain.out.matrix.all[complete.cases(Gain.out.matrix.all), , drop = FALSE]
    Gain.out.matrix.list[[v]] <- Gain.out.matrix.all
}

for (v in c(5)){
    J <- 1
    V.g.matrix.all <- matrix(NA, nrow=J, ncol=t_len)
    for (j in 1:J){
        file_name.V.g <- paste0(folder_path, v, ".", j,".V.g.psi.RData")
        load_result <- try(load(file_name.V.g), silent = TRUE)
        if (inherits(load_result, "try-error")) next
        V.g.matrix.all[j,] <- V.g.matrix.psi
    }
    V.g.matrix.all <- V.g.matrix.all[complete.cases(V.g.matrix.all), , drop = FALSE]
    V.g.matrix.list[[v]] <- V.g.matrix.all
}

for (v in c(5)){
    J <- 1
    V.a.vector.all <- numeric(J)
    for (j in 1:J){
        file_name.V.a <- paste0(folder_path, v, ".", j,".V.a.RData")
        load_result <- try(load(file_name.V.a), silent = TRUE)
        if (inherits(load_result, "try-error")) next
        V.a.vector.all[j] <- V.a.vector
    }
    V.a.vector.all <- V.a.vector.all[!is.na(V.a.vector.all)]
    V.a.vector.list[[v]] <- V.a.vector.all
}

for (v in c(5)){
    J <- 1
    Gain.trt.vector.all <- numeric(J)
    for (j in 1:J){
        file_name.Gain.trt <- paste0(folder_path, v, ".", j,".Gain.trt.RData")
        load_result <- try(load(file_name.Gain.trt), silent = TRUE)
        if (inherits(load_result, "try-error")) next
        Gain.trt.vector.all[j] <- Gain.trt.vector
    }
    Gain.trt.vector.all <- Gain.trt.vector.all[!is.na(Gain.trt.vector.all)]
    Gain.trt.vector.list[[v]] <- Gain.trt.vector.all
}

###########################################
# Estimate s_{c,T}(t)*s_{c,A}/(1-s_{c,A})
#########################################
Gain.out.df <- data.frame(value = numeric(), t = numeric(), j = integer(), d = integer())

for(i in c(5)) {

    cur_mat <- Gain.out.matrix.list[[i]]

    df <- data.frame(C.Y.sq = as.vector(cur_mat),
                     t = rep(1:ncol(cur_mat), each = nrow(cur_mat)),
                     j = rep(1:nrow(cur_mat), ncol(cur_mat)),
                     d = rep(i, length(cur_mat)))

    Gain.out.df <- rbind(Gain.out.df, df)
}

Gain.trt.df <- do.call(rbind, lapply(c(5), function(i) {
    data.frame(C.A.sq = Gain.trt.vector.list[[i]] / (1 - Gain.trt.vector.list[[i]]),
               j = 1:length(Gain.trt.vector.list[[i]]),
               d = as.factor(i))
}))

sens.df <- merge(Gain.out.df, Gain.trt.df, by = c("j", "d"))

dir.create("data/senspar", showWarnings = FALSE, recursive = TRUE)

save(sens.df, file = "data/senspar/sens.df.RData")

sens.df.mean <- sens.df %>%
    mutate(sens.par = C.Y.sq * C.A.sq) %>%
    group_by(d, t) %>%
    summarize(sens.par = mean(sens.par), .groups = "drop") %>%
    ungroup()

save(sens.df.mean, file = "data/senspar/sens.df.mean.RData")

senspar <- list()
senspar$sens.df <- sens.df
senspar$sens.df.mean <- sens.df.mean
save(senspar, file = "data/senspar.df.RHC.cluster.RData")
# For the supplied example repetition, at t = 30 and d = 5,
# round(senspar$sens.df.mean$sens.par, 4) equals 0.0021.
