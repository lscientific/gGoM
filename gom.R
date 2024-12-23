library(RSpectra)

flatten <- function(R, C_list) { # the response of R[j] is 0, 1, ..., C_list[j]-1
  N <- nrow(R)
  J <- ncol(R)
  R_flattened <- matrix(NA, N, sum(C_list))
  for(i in 1:N) {
    end <- 0
    for(j in 1:J) {
      start <- end + 1
      end <- end + C_list[j]
      r <- rep(0, C_list[j])
      r[R[i, j]+1] <- 1
      R_flattened[i, start:end] <- r
    }
  }
  return(R_flattened)
}


spa <- function(mat) {
  N <- dim(mat)[1]
  K <- dim(mat)[2]
  indices = rep(0, K)
  vertices = matrix(rep(0, K*K), K, K)
  Y <- mat
  
  for (k in 1:K) {
    row_norms <- apply(Y, 1, function(x) sqrt(sum(x^2)))
    idx <- which.max(row_norms)
    indices[k] = idx
    vertices[k, ] = mat[idx, ]
    u = Y[idx, ] / row_norms[idx]
    Y = Y %*% (diag(K) - matrix(u, nc=1) %*% matrix(u, nr=1))
  }
  
  return(list(indices=indices, vertices=vertices))
}

pruning <- function(mat, r=10, q=0.4, e=0.2) {
  N <- dim(mat)[1]
  row_norms = apply(mat, 1, function(x) sqrt(sum(x^2)))
  S0 <- which(row_norms > quantile(row_norms, 1-q))
  
  x = c()
  for (s in S0){
    mat_s <- t(apply(mat, 1, function(x) x - mat[s,]))
    norm_s <- apply(mat_s, 1, function(x) sqrt(sum(x^2)))
    d = norm_s[sort(norm_s, index.return=T)$ix[2:(r+1)]]
    x <- c(x, mean(d))
  }
  S <- S0[which(x > quantile(x, 1-e))]
  
  return(S)
}

gomSVD <- function(U, V, d, J, r=10, q=0.4, e=0.2, eps=0, truncate=T) {
  t1 <- Sys.time()
  
  N <- dim(U)[1]
  K <- dim(U)[2]
  S <- pruning(U, r, q, e)
  X = U[-S,]
  indices_X <- spa(X)$indices
  indices_U <- ((1:N)[-S])[indices_X]
  vertices <- X[indices_X, ]
  P1 <- U %*% solve(vertices)
  P2 <- t(apply(P1, 1, function(x) ifelse(x<0, 0, x)))
  P_hat <- t(apply(P2, 1, function(x) x / sum(x)))
  
  R_hat = U %*% diag(d) %*% t(V)
  if (truncate) {
    R_hat[R_hat > 1-eps] <- 1 - eps
    R_hat[R_hat < eps] <- eps
  }
  
  T_hat <- t(solve(t(P_hat) %*% P_hat) %*% t(P_hat) %*% U %*% diag(d) %*% t(V))
  if (truncate) {
    T_hat[T_hat > 1-eps] <- 1 - eps
    T_hat[T_hat < eps] <- eps
  }
  
  t2 <- Sys.time()
  
  return(list(P_hat=P_hat, T_hat=T_hat, R_hat=R_hat, 
              S=S, indices_U=indices_U, vertices=vertices, t=t2-t1))
}

gomSVD_Pois <- function(U, V, d, J, r=10, q=0.4, e=0.2, eps=0, truncate=T) {
  t1 <- Sys.time()
  
  N <- dim(U)[1]
  K <- dim(U)[2]
  S <- pruning(U, r, q, e)
  X = U[-S,]
  indices_X <- spa(X)$indices
  indices_U <- ((1:N)[-S])[indices_X]
  vertices <- X[indices_X, ]
  P1 <- U %*% solve(vertices)
  P2 <- t(apply(P1, 1, function(x) ifelse(x<0, 0, x)))
  P_hat <- t(apply(P2, 1, function(x) x / sum(x)))
  
  R_hat = U %*% diag(d) %*% t(V)
  if (truncate) {
    R_hat[R_hat > 1-eps] <- 1 - eps
    R_hat[R_hat < eps] <- eps
  }
  
  T_hat <- t(solve(t(P_hat) %*% P_hat) %*% t(P_hat) %*% U %*% diag(d) %*% t(V))
  if (truncate) {
    T_hat[T_hat < eps] <- eps
  }
  t2 <- Sys.time()
  
  return(list(P_hat=P_hat, T_hat=T_hat, R_hat=R_hat, 
              S=S, indices_U=indices_U, vertices=vertices, t=t2-t1))
}

perm <- function(x, p, K) {
  x_perm <- rep(NA, length(x))
  for (i in 1:length(x)) {
    for (k in 1:K) {
      x_perm[x == k] <- p[k]
    }
  }
  return(x_perm)
}

perm_mat <- function(c){
  n <- length(c)
  mat <- matrix(0, n, n)
  for(i in 1:n){
    mat[i, c[i]] <- 1 
  }
  return(mat)
}

rescale <- function(T_mat, C_list) {
  J <- nrow(T_mat)
  end <- 0
  for (j in 1:length(C_list)) {
    start <- end + 1
    end <- end + C_list[j]
    mat <- T_mat[start:end, ]
    T_mat[start:end, ] <- apply(mat, 2, function(x) x / sum(x))
  }
  return(T_mat)
}


ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  rho^exponent
}

