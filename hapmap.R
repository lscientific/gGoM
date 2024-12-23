source('./gom.R')
library(rARPACK)
library(fossil)
library(ggplot2)
library(GGally)

path <- "/Users/lingchen/Downloads/" # path to downloaded hapmap data
# downloaded from https://figshare.com/s/9b4d5964af498d167e85
# raw data at https://www.broadinstitute.org/medical-and-population-genetics/hapmap-3
data_hapmap <- read.csv(paste0(path, "hapmap.traw"), sep="\t")
eth <- read.csv(paste0(path, "hapmap_eth.csv"))[, 2]
eth_num <- as.numeric(as.factor(eth))
print(dim(data_hapmap))

data_gene <- t(data_hapmap[, -(1:6)])
data_hapmap <- data_hapmap[, 1:6]

# consider 4 ethnicity groups
subset <- which(eth %in% c('CEU', 'YRI', 'MEX', 'ASW'))
length(subset)

# data matrix R
R <- data_gene[subset, ]
rm(data_gene)
print(dim(R)) # 467x274128

N <- nrow(R)
J <- ncol(R)
cat("N =", N, ", J =", J, "\n")

# SVD
t1 <- Sys.time()
K <- 3
svd_res <- svds(R/2, K)
U <- svd_res$u
V <- svd_res$v
d <- svd_res$d
rm(svd_res)
mu1_hat <- sqrt(N) * max(apply(U, 1, function(x) sqrt(sum(x^2)))) / sqrt(K)
mu2_hat <- sqrt(J) * max(apply(V, 1, function(x) sqrt(sum(x^2)))) / sqrt(K)
# empirical incoherence numbers
print(mu1_hat)
print(mu2_hat)

# GoM
r <- 10
q <- 0.4
e <- 0.2
eps <- 0
res <- gomSVD(U, V, d, J, r=r, q=q, e=e, eps=eps)
t2 <- Sys.time()
print(t2 - t1) # 8.77 secs


#### residual correlation ####
library(reshape2)
R_hat <- P_hat %*% t(T_hat)

cormat_dep <- cov((R/2)[,1:100] - R_hat[,1:100])
cor_dep <- melt(cormat_dep)
colnames(cor_dep)[3] <- 'covariance'
ggplot(data = cor_dep, aes(x=Var1, y=Var2, fill=covariance)) + geom_tile() +
  scale_fill_gradient(low = "white", high = "black") + xlab("") + ylab("") + 
  theme(panel.grid.major = element_blank())
# ggsave("./figures/HapMap_cov_heatmap.png", width=6, height=5)
## this figure corresponds to the left panel in Figure 8


#### plot ####
# ternary
library(ggtern)
data_tern_hapmap <- as.data.frame(P_hat)
data_tern_hapmap$Population <- as.factor(eth[subset])
p1 <- ggtern(data=data_tern_hapmap, aes(V1, V2, V3, color=Population)) + 
  geom_point(size=0.7, alpha=1) + 
  xlab("CEU") + ylab("YRI") + zlab("3rd ancestral") + 
  theme(axis.title=element_text(size=8.5), 
        legend.title= element_text(size=9), 
        legend.text=element_text(size=9), 
        legend.key.size=unit(0.6, 'cm'), 
        legend.position = "bottom", 
        tern.axis.title.L = element_text(hjust = 0.4, vjust = 0.5),
        tern.axis.title.R = element_text(hjust = 0.65, vjust = 0.5),
        legend.box.margin = margin(-35, 0, 0, 0), 
        plot.margin = margin(-15, -15, 0, -30)) +
  scale_color_manual(values=c('#c19883', '#FC9239', '#78bf9d', '#028AD0'))
p1

# ggsave("./figures/hapmap_ternary.png", width=5, height=4)
## This figure corresponds to the left panel in Figure 5
