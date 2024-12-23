source("./gom.R")
library(gtools)
library(sirt)
library(parallel)
#library(MASS)

library(ggplot2)
library(latex2exp)
library(ggpubr)
library(dplyr)


C <- 3 # the number of possible outcomes
K <- 3
ns <- c(200, 1000, 2000, 3000, 4000, 5000)
n_rep <- 100
block <- 10 # number of blocks
rho <- 0.5 # autoregressive parameter

#### generate data ####
set.seed(234)

R_list <- T0_list <- P0_list <- list()
for(ni in 1:length(ns)) {
  N <- ns[ni]
  J <- N / 5
  R_list[[ni]] <- list()
  T0_list[[ni]] <- list()
  P0_list[[ni]] <- list()
  
  for(rep in 1:n_rep) {
    cat("ni=", ns[ni], "rep =", rep, "\n")
    N <- ns[ni]
    J <- floor(N/5)
    P0 <- rdirichlet(N, rep(1, K))
    P0[1:K, ] <- diag(1, K)
    T0 <- matrix(rbeta(J*K, 0.2, 0.2), nr=J, nc=K)
    T0_list[[ni]][[rep]] <- T0
    P0_list[[ni]][[rep]] <- P0
    
    R0 <- P0 %*% t(T0)
    # obtain quantiles
    D <- matrix(NA, nrow=N, ncol=J)
    for(j in 1:J) {
      D[, j] <- sapply(R0[, j], qnorm)
    }
    # generate locally dependent data
    R <- matrix(NA, nrow=N, ncol=J)
    for(i in 1:N) {
      Y <- mvrnorm(J/block, mu=rep(0, block), Sigma=ar1_cor(block, rho))
      Y <- as.vector(t(Y))
      R[i, ] <- as.numeric(Y < D[i, ])
    }
    R_list[[ni]][[rep]] <- R
  }
}

# save(R_list, T0_list, P0_list, file=paste0("./data/sim_loc_rho=",rho,".RData"))


#### GoM ####
set.seed(1234)

r <- 10
q <- 0.4
e <- 0.2
eps <- 0

load(paste0("./data/sim_loc_rho=", rho, ".RData"))
rm(P0_list, T0_list)
res <- list()
for(ni in 1:length(ns)) {
  cat("N =", ns[ni], "\n")
  N <- ns[ni]
  J <- floor(N/5)
  
  res[[ni]] <- mclapply(R_list[[ni]], function(R) {
    t1 <- Sys.time()
    res <- svds(R, k=K)
    U <- res$u; V <- res$v; d <- res$d
    
    res <- gomSVD(U, V, d, J, r=r, q=q, e=e, eps=eps)
    t2 <- Sys.time()
    res$t <- t2-t1
    return(res)
  })
}

# save(res, file=paste0("./data/sim_loc_res_rho=", rho, ".RData"))


#### results ####
perms <- list()
perms[[1]] <- perm_mat(1:3)
perms[[2]] <- perm_mat(c(1, 3, 2))
perms[[3]] <- perm_mat(c(2, 1, 3))
perms[[4]] <- perm_mat(c(2, 3, 1))
perms[[5]] <- perm_mat(c(3, 1, 2))
perms[[6]] <- perm_mat(c(3, 2, 1))

load(paste0("./data/sim_loc_rho=", rho, ".RData"))
load(paste0("./data/sim_loc_res_rho=", rho, ".RData"))
rm(R_list)

err_T = err_P = err_max_T = err_max_P <- matrix(NA, n_rep, length(ns))
for(ni in 1:length(ns)) {
  for(rep in 1:n_rep) {
    P0 <- P0_list[[ni]][[rep]]
    T0 <- T0_list[[ni]][[rep]]
    
    P_hat <- res[[ni]][[rep]]$P_hat
    T_hat <- res[[ni]][[rep]]$T_hat
    
    err_T_try <- rep(NA, 6)
    for(try in 1:6) {
      perm <- perms[[try]]
      err_T_try[try] <- mean(abs(T_hat - T0 %*% perm))
    }
    perm <- perms[[which.min(err_T_try)]]
    err_T[rep, ni] <- mean(abs(T_hat - T0 %*% perm))
    err_P[rep, ni] <- mean(abs(P_hat - P0 %*% perm))
    err_max_T[rep, ni] <- max(abs(T_hat - T0 %*% perm))
    err_max_P[rep, ni] <- max(apply(P_hat - P0 %*% perm, 1, function(x) sqrt(sum(x^2))))
  }
} 

tm <- theme(axis.text=element_text(size=12), axis.title=element_text(size=12), 
            plot.title = element_text(size=15), legend.title = element_blank())

err_T_dat <- data.frame(err=as.vector(err_max_T), rep=rep(1:n_rep, length(ns)), n=as.factor(rep(ns, each=n_rep)))
err_P_dat <- data.frame(err=as.vector(err_max_P), rep=rep(1:n_rep, length(ns)), n=as.factor(rep(ns, each=n_rep)))
p1 <- ggplot(data=err_P_dat, aes(x=n, y=err)) + 
  geom_boxplot(outlier.size=0.3, fatten = 1.5) + 
  xlab("N") + ylab("") + 
  ggtitle(TeX(paste0('2-to-infinity error for $\\hat{\\Pi}$, $\\rho=', rho, '$'))) + tm
p2 <- ggplot(data=err_T_dat, aes(x=n, y=err)) + 
  geom_boxplot(outlier.size=0.3, fatten = 1.5) + 
  xlab("N") + ylab("") + 
  ggtitle(TeX(paste0('Max absolute error for $\\hat{\\Theta}$, $\\rho=', rho, '$'))) + tm


p_err <- ggarrange(p1, p2, ncol = 2, nrow = 1)
p_err
## this figure corresponds to Figure 3
# ggsave("./figures/sim_loc_err_rho=0.5.png", plot=p_err, width=10, height=5)

