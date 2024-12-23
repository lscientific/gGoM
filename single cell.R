source("./GoM.R")

#### load the data ####
# install_github('kkdey/singleCellRNASeqMouseDeng2014') 
library(singleCellRNASeqMouseDeng2014)
deng.counts <- exprs(Deng2014MouseESC)
deng.meta_data <- pData(Deng2014MouseESC) # cell type
deng.gene_names <- rownames(deng.counts)

print(max(deng.counts)) # maximum gene expression count
mean(deng.counts != 0) # sparsity

# data matrix R
R <- t(deng.counts)
N <- nrow(R)
J <- ncol(R)

cell_type <- as.character(deng.meta_data$cell_type)
embryo <- deng.meta_data$embryo_id
table(cell_type) # cell types

#### estimation ####
t1 <- Sys.time()
K <- 6
cat("N =", N, "J =", J, "K =", K, "\n")
# svd
svd_res <- svds(R, K)
U <- svd_res$u
V <- svd_res$v
# empirical incoherence numbers
sqrt(N) * max(apply(U, 1, function(x) sqrt(sum(x^2)))) / sqrt(6)
sqrt(J) * max(apply(V, 1, function(x) sqrt(sum(x^2)))) / sqrt(6)

d <- svd_res$d

# GoM
r <- 10
q <- 0.4
e <- 0.2
eps <- 0
res <- gomSVD(U, V, d, J, r=r, q=q, e=e, eps=eps, truncate = T)
t2 <- Sys.time()
print(t2-t1)

P_hat <- res$P_hat
T_hat <- res$T_hat

print(head(cbind(cell_type, P_hat)))

