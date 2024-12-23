gibbs_sample_pol <- function(theta, pi, z, R_bin, alphas, betas) {
  # theta is J*K*C
  # pi is N*K
  # z is N*J
  # R_bin is N*J*C
  # alphas is of length K
  # betas is of length C
  
  N <- nrow(pi)
  J <- dim(theta)[1]
  K <- dim(theta)[2]
  C <- dim(theta)[3]
  
  # update theta
  for(j in 1:J) {
    for(k in 1:K) {
      z_bin <- as.numeric(z[, j] == k)
      betas_updated <- matrix(z_bin, nr=1) %*% R_bin[, j, ] + betas
      theta[j, k, ] <- rdirichlet(1, betas_updated)
    }
  }
  
  # update pi
  for(i in 1:N) {
    pi[i, ] <- rdirichlet(1, sapply(1:K, function(k) sum(z[i, ] == k)) + alphas)
  }
  
  # update z
  for(i in 1:N) {
    for(j in 1:J) {
      weights <- theta[j, , which(R_bin[i, j, ] == 1)] * pi[i, ]
      weights <- weights / sum(weights)
      z[i, j] <- sample(1:K, 1, T, weights)
    }
  }
  
  return(list(theta=theta, pi=pi, z=z))
}


gibbs_gom_pol <- function(theta, pi, z, R_bin, alphas, betas, n_burnin, n_sample) {
  # theta is J*K*C
  # pi is N*K
  # z is N*J
  # R_bin is N*J*C
  # alphas is of length K
  # betas is of length C
  
  N <- nrow(pi)
  J <- dim(theta)[1]
  K <- dim(theta)[2]
  C <- dim(theta)[3]
  pi_sum <- matrix(0, N, K)
  theta_sum <- array(0, dim=c(J, K, C))
  theta_trace = pi_trace <- rep(NA, n_burnin+n_sample)
  
  for(iter in 1:(n_burnin+n_sample)) {
    if(iter %% 1000 == 0) print(iter)
    
    samples <- gibbs_sample_pol(theta, pi, z, R_bin, alphas, betas)
    theta <- samples$theta; pi <- samples$pi; z <- samples$z
    theta_trace[iter] <- theta[1, 1, 1]; pi_trace[iter] <- pi[1, 1]
    
    if(iter > n_burnin) {
      theta_sum <- theta_sum + theta
      pi_sum <- pi_sum + pi
    }
  }
  theta_est <- theta_sum / n_sample
  pi_est <- pi_sum / n_sample
  return(list(theta_est=theta_est, pi_est=pi_est, theta_trace=theta_trace, pi_trace=pi_trace))
}
