source("./gomSVD.R")
library(gtools)
library(foreach)
library(doParallel)
library(parallel)
library(MASS)

K <- 3
ns <- c(200, 400, 600, 800, 1000)
n_rep <- 100

#### generate data ####
set.seed(123)

R_list <- T0_list <- P0_list <- list()
for(ni in 1:length(ns)) {
  N <- ns[ni]
  J <- N * 5
  R_list[[ni]] <- list()
  T0_list[[ni]] <- list()
  P0_list[[ni]] <- list()
  cat("ni=", ns[ni], "\n")
  
  for(rep in 1:n_rep) {
    P0 <- rdirichlet(N, rep(1, K))
    P0[1:K, ] <- diag(1, K)
    T0 <- matrix(rgamma(J*K, 1, 2), nr=J, nc=K)
    T0_list[[ni]][[rep]] <- T0
    P0_list[[ni]][[rep]] <- P0
    
    R0 <- P0 %*% t(T0)
    R <- matrix(NA, nrow=N, ncol=J)
    for(i in 1:N) {
      for(j in 1:J) {
        R[i, j] <- rpois(1, R0[i,j])
      }
    }
    R_list[[ni]][[rep]] <- R
  }
}

save(R_list, T0_list, P0_list, file="./data/sim_pois_J=5N.RData")


#### GoM ####
r <- 10
q <- 0.4
e <- 0.2
eps <- 0

set.seed(321)
load("./data/sim_pois_J=5N.RData")

results <- list()
for(ni in 1:length(ns)) {
  N <- ns[ni]
  J <- N * 5
  cat("ni=", ni, "\n")
  
  myCluster <- makeCluster(detectCores()-1, type = "PSOCK")
  registerDoParallel(myCluster)
  results[[ni]] <- foreach(i_rep=1:n_rep, .packages=c('sirt', 'gtools', 'RSpectra')) %dopar% {
    t1 <- Sys.time()
    res <- svd(R_list[[ni]][[i_rep]], nu=K, nv=K)
    U <- res$u[, 1:K]; V <- res$v[, 1:K]; d <- res$d[1:K]
    
    res <- gomSVD_Pois(U, V, d, J, r=r, q=q, e=e, eps=eps)
    
    t2 <- Sys.time()
    res$t <- t2-t1
    return(res)
  }
  stopCluster(myCluster)
  
  save(results, file="./data/sim_res_gom_Pois_J=5N.RData")
}


#### NMF ####
set.seed(321)
load("./data/sim_pois_J=5N.RData")

results <- list()
for(ni in 1:length(ns)) {
  N <- ns[ni]
  J <- N * 5
  cat("ni=", ni, "\n")
  
  myCluster <- makeCluster(detectCores()-1)
  registerDoParallel(myCluster)
  results[[ni]] <- foreach(i_rep=1:n_rep, .packages=c('NMF')) %dopar% {
    t1 <- Sys.time()
    
    test <- nmf(R_list[[ni]][[i_rep]], K, 'lee') 
    R0_hat <- fitted(test)
    P_hat <- basis(test)
    T_hat <- t(coef(test))
    
    res <- list()
    res$R0_hat <- R0_hat
    res$P_hat <- P_hat
    res$T_hat <- T_hat
    
    t2 <- Sys.time()
    res$t <- t2-t1
    return(res)
  }
  stopCluster(myCluster)
  
  save(results, file="./data/sim_res_nmf_Pois_J=5N.RData")
}


#### res ####
library(ggplot2)
library(latex2exp)
library(ggpubr)
library(dplyr)

load("./data/sim_pois_J=5N.RData")
results_gom <- get(load("./data/sim_res_gom_Pois_J=5N.RData"))
results_nmf <- get(load("./data/sim_res_nmf_Pois_J=5N.RData"))

t_res <- array(NA, dim=c(n_rep, 4, 2),
               dimnames=list(rep=1:n_rep, N=ns[1:4], Method=c('Proposed', 'NMF')))
err_res <- array(NA, dim=c(n_rep, 4, 2),
                 dimnames=list(rep=1:n_rep, N=ns[1:4], Method=c('Proposed', 'NMF')))

for(ni in 1:4) {
  N <- ns[ni]
  J <- N * 5
  print(ni)
  
  for(rep in 1:n_rep) {
    T0 <- T0_list[[ni]][[rep]]
    P0 <- P0_list[[ni]][[rep]]
    R0 <- P0 %*% t(T0)
    
    # Proposed
    res <- results_gom[[ni]][[rep]]
    P_hat <- res$P_hat
    T_hat <- res$T_hat
    R0_hat <- P_hat %*% t(T_hat)
    #print(res$t)
    t_res[rep, ni, 'Proposed'] <- as.numeric(res$t)
    err_res[rep, ni, 'Proposed'] <- sqrt(sum((R0 - R0_hat)^2)) / sqrt(N*J)
    
    # NMF
    res <- results_nmf[[ni]][[rep]]
    #print(res$t)
    t <- as.numeric(res$t)
    if ((ni > 1) & (t < 5)) t <- t * 60 
    t_res[rep, ni, 'NMF'] <- t
    err_res[rep, ni, 'NMF'] <- sqrt(sum((R0 - res$R0_hat)^2)) / sqrt(N*J)
  }
}

err_dat2 <- as.data.frame.table(err_res)
err_dat2$Method <- factor(err_dat2$Method, levels=c('NMF', 'Proposed'))

p3 <- ggplot(data=err_dat2, aes(x=N, y=Freq, fill=Method)) + 
  geom_boxplot(outlier.size=0.3, fatten = 1) + 
  xlab("N (J=5N)") + ylab("") + 
  ggtitle(TeX('Scaled Frobenius error of $\\hat{R}$')) + 
  theme(legend.title = element_blank()) + 
  scale_fill_manual(values=c("#388EC2", "#E94F42"))

t_dat2 <- as.data.frame.table(t_res)
t_dat2$Method <- factor(t_dat2$Method, levels=c('NMF', 'Proposed'))

#p4 <- ggplot(data=t_dat2, aes(x=N, y=Freq, fill=Method)) + 
#  geom_boxplot(outlier.size=0.3, fatten = 1) + 
#  xlab("N (J=5N)") + ylab("") + 
#  ggtitle(TeX('Computation time in seconds')) + 
#  theme(legend.title = element_blank()) + 
#  scale_fill_manual(values=c("#999999", "#E69F00"))
df <- t_dat2 %>%
  group_by(N, Method) %>%
  summarise(
    q25 = quantile(Freq, 0.25, na.rm = TRUE),
    q75 = quantile(Freq, 0.75, na.rm = TRUE),
    t = median(Freq)
  )

df$Method <- factor(df$Method, levels = c('NMF', 'Proposed'))

p4 <- ggplot(df, aes(N, t, color = Method, lty = Method)) +
  geom_line(aes(group = Method), position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.3, 
                position = position_dodge(0.3)) +
  scale_linetype_manual(values = c("longdash", "solid")) +
  scale_color_manual(values = c("#388EC2", "#E94F42")) +
  ylab("") + 
  xlab("N (J=5N)") + 
  ggtitle(TeX('Computation time in seconds'))

ggarrange(p1, p3, p2, p4, nrow=1, ncol=4, common.legend=T)
ggsave("./figures/sim_pois_compare.png", width=12, height=3.5)
