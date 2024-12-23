source("./gom.R")
source("./gibbs_util.R")
library(gtools)
library(sirt)
library(mixedMem)

library(foreach)
library(parallel)

library(ggplot2)
library(ggpubr)
library(latex2exp)


C <- 3 # the number of possible outcomes
K <- 3 # the number of extreme profiles
ns <- c(200, 1000, 2000, 3000, 4000, 5000) # sample sizes
n_rep <- 100

# permutations of the extreme profiles
perms <- list()
perms[[1]] <- perm_mat(1:3)
perms[[2]] <- perm_mat(c(1, 3, 2))
perms[[3]] <- perm_mat(c(2, 1, 3))
perms[[4]] <- perm_mat(c(2, 3, 1))
perms[[5]] <- perm_mat(c(3, 1, 2))
perms[[6]] <- perm_mat(c(3, 2, 1))


#### simulate data ####
set.seed(1234)

for(ni in 1:length(ns)) {
  N <- ns[ni]
  J <- floor(N/5)
  print(N)
  
  myCluster <- makeCluster(detectCores()-1, type = "PSOCK")
  data <- foreach(i_rep = 1:n_rep, .packages=c('gtools')) %dopar% {
    P0 <- rdirichlet(N, rep(1, K))
    P0[1:K, ] <- diag(1, K)
    
    T0 <- matrix(NA, nc=K, nr=J*C)
    for(j in 1:J) {
      for (k in 1:K) {
        T0[((j-1)*C+1):(j*C), k] <- rdirichlet(1, rep(0.2, 3))
      }
    }
    
    R = matrix(NA, N, C*J)
    for(i in 1:N) {
      for(j in 1:J) {
        R[i, ((j-1)*C+1):(j*C)] <- rmultinom(1, 1, matrix(P0[i, ], nr=1) %*% t(T0[((j-1)*C+1):(j*C), ]))
      }
    }
    
    return(list(R=R, T0=T0, P0=P0))
  }
  
  #save(data, file=paste0("./data/sim_N=", N, "_dirichlet.RData"))
  stopCluster(myCluster)
}


#### GoM ####
# pruning and thresholding parameters
r <- 10
q <- 0.4
e <- 0.2
eps <- 0

for(ni in 1:length(ns)) {
  cat("N =", ns[ni], "\n")
  N <- ns[ni]
  J <- floor(N/5)
  load(paste0("./data/sim_N=", N, "_dirichlet.RData"))
  
  myCluster <- makeCluster(detectCores()-1, type = "PSOCK")
  results <- foreach(i_rep=1:n_rep, .packages=c('sirt', 'gtools', 'RSpectra')) %dopar% {
    t1 <- Sys.time()
    res <- svd(data[[i_rep]][[1]], nu=K, nv=K)
    U <- res$u[, 1:K]; V <- res$v[, 1:K]; d <- res$d[1:K]
    
    res <- gomSVD(U, V, d, J, r=r, q=q, e=e, eps=eps)
    
    T_hat <- res$T_hat
    C <- nrow(V) / J
    for(j in 1:J) {
      T_hat[((j-1)*C+1):(j*C), ] <- apply(T_hat[((j-1)*C+1):(j*C), ], 2, function(x) x/sum(x))
    }
    res$T_hat <- T_hat
    
    t2 <- Sys.time()
    res$t <- t2-t1
    return(res)
  }
  
  stopCluster(myCluster)
  
  #save(results, file=paste0("./data/sim_res_gom_dirichlet_N=", N, ".RData"))
}


#### VI with mixedMem ####
set.seed(123)
alpha <- rep(0.5, K)

res_em <- list()
for(ni in 1:3) {
  N <- ns[ni]
  J <- floor(N/5)
  cat("N =", N, "J =", J, "\n")
  load(paste0("./data/sim_N=", N, "_dirichlet.RData"))
  
  Rj <- rep(1, J)
  Nijr <- array(1, dim = c(N, J, 1))
  Vj <- rep(C, J)
  dist <- rep("multinomial", J)
  theta <- array(0, dim = c(J, K, C))
  for (j in 1:J) theta[j, , ] <- rdirichlet(K, rep(0.3, C))
  
  myCluster <- makeCluster(detectCores()-1, type = "PSOCK")
  res_em[[ni]] <- foreach(i_rep=1:n_rep, .packages=c('mixedMem')) %dopar% {
    # R matrix to array
    R <- data[[i_rep]][[1]]
    obs <- array(NA, dim=c(N, J, 1, 1))
    for(i in 1:N) {
      for(j in 1:J) {
        obs[i, j, , ] <- which(R[i, ((j-1)*C+1):(j*C)] == 1) - 1
      }
    }
    
    test_model <- mixedMemModel(N, J, Rj, Nijr, K, Vj, alpha, theta, dist=dist, obs=obs)
    t1 <- Sys.time()
    out <- mmVarFit(test_model)
    t2 <- Sys.time(); t <- t2-t1
    out$t <- t
    return(out)
  }
  stopCluster(myCluster)
  
  #save(res_em, file=paste0("./results/sim_res_em_dirichlet.RData"))
}



#### Gibbs ####
alphas <- rep(1, K)
betas <- rep(1, C)
n_burnin <- 5000
n_sample <- 2000

set.seed(123)

res_gibbs_list <- list()
for(ni in 1:3) {
  N <- ns[ni]
  J <- N / 5
  print(N)
  load(paste0("./data/sim_N=", N, "_dirichlet.RData"))
  
  myCluster <- makeCluster(detectCores()-1, type = "PSOCK", outfile="./gibbs_pol.log")
  registerDoParallel(myCluster)
  
  T1 <- Sys.time()
  foreach(i_rep = 1:n_rep, .packages=c('stats', 'gtools')) %dopar% {
    print(i_rep)
    
    R <- data[[i_rep]][[1]]
    R_bin <- array(0, c(N, J, C))
    for(i in 1:N) {
      for(j in 1:J) {
        R_bin[i, j, ] <- R[i, ((j-1)*C+1):(j*C)]
      }
    }
    
    t1 <- Sys.time()
    theta <- array(NA, dim=c(J, K, C))
    for (j in 1:J) {
      for(k in 1:K) {
        theta[j, k, ] <- rdirichlet(1, betas)
      }
    }
    
    pi <- rdirichlet(N, alphas)
    
    z <- matrix(NA, N, J)
    for(i in 1:N) z[i, ] <- sample(1:K, J, TRUE, pi[i, ])
    
    res_gibbs <- gibbs_gom_pol(theta, pi, z, R_bin, alphas, betas, n_burnin, n_sample)
    t2 <- Sys.time()
    print(t2 - t1)
    res_gibbs$t <- t2 - t1
    
    # save(res_gibbs, file=paste0("./data/gibbs_results/res_", ni, "_", i_rep, "_dirichlet.RData"))
  }
  T2 <- Sys.time()
  print(T2-T1)
  
  stopCluster(myCluster)
}


#### results VI ####
# normalized phi is the estimated membership
err_T_em <- matrix(NA, n_rep, 3)
err_P_em <- matrix(NA, n_rep, 3)
t_em_res <- matrix(NA, n_rep, 3)

for(ni in 1:3) {
  print(ni)
  N <- ns[ni]
  J <- floor(N/5)
  load(paste0("./data/sim_res_em_dirichlet_", N, ".RData"))
  load(paste0("./data/sim_N=", N, "_dirichlet.RData"))
  
  for(rep in 1:n_rep) {
    P0 <- data[[rep]][[3]]
    T0 <- data[[rep]][[2]]
    
    T_hat0 <- res_em[[rep]]$theta
    print(res_em[[rep]]$t)
    t_em_res[rep, ni] <- res_em[[rep]]$t
    T_hat <- matrix(NA, nr=nrow(T0), nc=ncol(T0))
    for(j in 1:J) {
      T_hat[((j-1)*C+1):(j*C), ] <- t(T_hat0[j, , ])
    }
    
    P_hat <- res_em[[rep]]$phi
    P_hat <- t(apply(P_hat, 1, function(x) x / sum(x)))
    
    err_T_try <- rep(NA, 6)
    for(try in 1:6) {
      perm <- perms[[try]]
      err_T_try[try] <- mean(abs(T_hat - T0 %*% perm))
    }
    perm <- perms[[which.min(err_T_try)]]
    err_T_em[rep, ni] <- mean(abs(T_hat - T0 %*% perm))
    err_P_em[rep, ni] <- mean(abs(P_hat - P0 %*% perm))
  }
}

t_em_res[, 2:3] <- t_em_res[, 2:3]*60 # minutes to seconds

## the following values correspond to the second row of Table 1
print(apply(t_em_res, 2, mean))
## the following values correspond to the second row of Table 2
print(apply(err_T_em, 2, mean)) # Theta
print(apply(err_P_em, 2, mean)) # Pi


#### results GoM ####
err_T = err_P = err_max_T = err_max_P <- matrix(NA, n_rep, length(ns))
t_res <- matrix(NA, n_rep, length(ns))

for(ni in 1:length(ns)) {
  N <- ns[ni]
  J <- floor(N/5)
  cat("N=", N, "\n")
  load(paste0("./data/sim_N=", N, "_dirichlet.RData"))
  load(paste0("./data/sim_res_gom_dirichlet_N=", N, ".RData"))
  
  for(i_rep in 1:n_rep) {
    P0 <- data[[i_rep]][[3]]
    T0 <- data[[i_rep]][[2]]
    
    P_hat <- results[[i_rep]]$P_hat
    T_hat <- results[[i_rep]]$T_hat
    print(results[[i_rep]]$t)
    t_res[i_rep, ni] <- as.numeric(results[[i_rep]]$t)
    
    err_T_try <- rep(NA, 6)
    for(try in 1:6) {
      perm <- perms[[try]]
      err_T_try[try] <- max(abs(T_hat - T0 %*% perm))
    }
    perm <- perms[[which.min(err_T_try)]]
    err_T[i_rep, ni] <- mean(abs(T_hat - T0 %*% perm))
    err_P[i_rep, ni] <- mean(abs(P_hat - P0 %*% perm))
    
    err_max_T[i_rep, ni] <- max(abs(T_hat - T0 %*% perm))
    err_max_P[i_rep, ni] <- max(apply(P_hat - P0 %*% perm, 1, function(x) sqrt(sum(x^2))))
  }
}

t_res[, 5:6] <- t_res[, 5:6]*60 # minutes to seconds
## the following values correspond to the first row of Table 1
print(apply(t_res[, 1:3], 2, mean))
## the following values correspond to the first row of Table 2
print(apply(err_T[, 1:3], 2, mean)) # Theta
print(apply(err_P[, 1:3], 2, mean)) # Pi


# plots
err_T_dat <- data.frame(err=as.vector(err_max_T), 
                        rep=rep(1:n_rep, length(ns)), 
                        n=as.factor(rep(ns, each=n_rep)))
err_P_dat <- data.frame(err=as.vector(err_max_P), 
                        rep=rep(1:n_rep, length(ns)), 
                        n=as.factor(rep(ns, each=n_rep)))

p1 <- ggplot(data=err_P_dat, aes(x=n, y=err)) + 
  geom_boxplot(outlier.size=0.3, fatten = 1.5) + 
  xlab("N") + ylab("") + 
  ggtitle(TeX('2-to-infinity error for $\\Pi$')) + 
  theme(legend.title = element_blank())

p2 <- ggplot(data=err_T_dat, aes(x=n, y=err)) + 
  geom_boxplot(outlier.size=0.3, fatten = 1.5) + 
  xlab("N") + ylab("") + 
  ggtitle(TeX('Max absolute error for $\\Theta$')) + 
  theme(legend.title = element_blank())

p_err <- ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = T)
p_err

## this figure corresponds to Figure 1
# ggsave("./figures/flatten_gom_err_dirichlet.png", plot=p_err, width=10, height=4.5)


#### results Gibbs ####
err_T_gibbs = err_P_gibbs <- matrix(NA, n_rep, 3)
t_res_gibbs <- matrix(NA, n_rep, 3)

for(ni in 1:3) {
  N <- ns[ni]
  J <- N / 5
  load(paste0("./data/sim_N=", N, "_dirichlet.RData"))
  
  for(i_rep in 1:n_rep) {
    load(paste0("./data/gibbs_results/res_", ni, "_", i_rep, "_dirichlet.RData"))
    P0 <- data[[i_rep]][[3]]
    T0 <- data[[i_rep]][[2]]
    
    # Gibbs
    print(res_gibbs$t)
    pi_hat <- res_gibbs$pi_est
    theta_hat <- res_gibbs$theta_est
    t_res_gibbs[i_rep, ni] <- res_gibbs$t
    
    err_P_try <- rep(NA, 6)
    for(try in 1:6) {
      perm <- perms[[try]]
      err_P_try[try] <- mean(abs(pi_hat - P0 %*% perm))
    }
    perm <- perms[[which.min(err_P_try)]]
    
    theta0 <- array(NA, dim=c(J, K, C))
    for (j in 1:J) {
      for (k in 1:K) {
        theta0[j, k, ] <- (T0 %*% perm)[((j-1)*C+1):(j*C), k]
      }
    }
    
    err_T_gibbs[i_rep, ni] <- mean(abs(theta_hat - theta0))
    err_P_gibbs[i_rep, ni] <- mean(abs(pi_hat - P0 %*% perm))
  }
}

t_res_gibbs[, 1] <- t_res_gibbs[, 1]*60 # minutes to seconds
t_res_gibbs[, 2] <- t_res_gibbs[, 2]*60*60 # hours to seconds
t_res_gibbs[, 3] <- t_res_gibbs[, 3]*60*60 # hours to seconds

## the following values correspond to the third row of Table 1
print(apply(t_res_gibbs, 2, function(x) mean(x, na.rm=T)))
## the following values correspond to the third row of Table 2
print(apply(err_T_gibbs, 2, function(x) mean(x, na.rm=T)))
print(apply(err_P_gibbs, 2, function(x) mean(x, na.rm=T)))

