source("./gom.R")

path <- "/Users/lingchen/Downloads/" # path to downloaded ANES 2022 data
# data is downloaded from https://electionstudies.org/data-center/
dat <- read.csv(paste0(path, "anes_pilot_2022_csv_20221214.csv"))[, -1]
#dat <- dat[which(!is.na(dat$weight)), ] # use QC-ed ones
pid <- dat$pid_x # party ID
idx_noparty <- which(is.na(pid)) # indices without party ID
pid3 <- dat$pid3

#### data pre-processing ####
# manually exclude items including: rate (1-100), text, 
# "how important is being ... to your identity?", and followup questions
fields_unchosen <- c("howreg", "whenreg", "regdiff", "howreg_os", 
                     "vharder_0", "vharder_1", "vharder_2", 
                     "vharder_3", "vharder_4", "vharder_5", "vharder_6", "vharder_7", 
                     "vharder_8", "vharder_9", "vharder_10", "vharder_11", "vharder_12", 
                     "turnout22", "turnout22ns", "turnout22w", "pipevote22a", "pipevote22b", "votehard",
                     "waittime", "triptime", "house22t", "house22p", "senate22t", "senate22p", "gov22t",
                     "gov22p", "turnout20", "turnout20b", "pipevote20", "vote20", 
                     "pidstr", "pidlean", "maga", "pid_x", "roefollow", "abpros1srt", "abpros2str",
                     "vharder_t", "pid2d", "pid2r", "abort_4pta_os", "abort_4ptb_os", 
                     "attnflag", "ftblack1", "ftblack2", "ftwhite1", "ftwhite2", "fthisp", 
                     "ftasian", "ftfbi", "ftscotus", "fttrump", "ftbiden", "ftdem",
                     "femwell", "nonfemwell", "femimp", "nonfemimp", 
                     "ftrep", "ftteach", "ftfem", "ftnfem", "ftjourn", "ftmen", "ftwomen", 
                     "fttrans", "jan6therm", "identabort_os", "impabortid", "imppartyid", 
                     "impfundid", "impevanid", "impamerid", "imphispid", "impwhite", "impblackid", 
                     "impnatid", "impasianid", "impnhpiid", "whpriv2", "blpriv2", "hipriv2", 
                     "aspriv2", "libcon")
R <- dat[, 11:253]
R_sub <- R[, !(colnames(R) %in% fields_unchosen)]
dim(R_sub)

# need to combine pid1r & pid1d, abort_4pta & abort_4ptb, abortimp_a & abortimp_b, transbath1 & transbath2
# transsport1 & transsport2, rr1a & rr1b, rr2a & rr2b, rr3a & rr3b, rr4a & rr4b, 

pre_processing <- function(x) {
  if ((x[1] == -1) & (x[2] == -1)) return (-1)
  if ((x[1] == -1) & (x[2] != -1)) return (x[2])
  if ((x[1] != -1) & (x[2] == -1)) return (x[1])
}

R_sub$pid <- apply(R_sub[, c("pid1r", "pid1d")], 1, pre_processing)
R_sub$abort_4pt <- apply(R_sub[, c("abort_4pta", "abort_4ptb")], 1, pre_processing)
R_sub$abortimp <- apply(R_sub[, c("abortimp_a", "abortimp_b")], 1, pre_processing)
R_sub$transbath <- apply(R_sub[, c("transbath1", "transbath2")], 1, pre_processing)
R_sub$transsport <- apply(R_sub[, c("transsport1", "transsport2")], 1, pre_processing)
R_sub$rr1 <- apply(R_sub[, c("rr1a", "rr1b")], 1, pre_processing)
R_sub$rr2 <- apply(R_sub[, c("rr2a", "rr2b")], 1, pre_processing)
R_sub$rr3 <- apply(R_sub[, c("rr3a", "rr3b")], 1, pre_processing)
R_sub$rr4 <- apply(R_sub[, c("rr4a", "rr4b")], 1, pre_processing)

R_sub <- R_sub[, !(colnames(R_sub) %in% c("pid1r", "pid1d", "abort_4pta", "abort_4ptb", 
                                          "abortimp_a", "abortimp_b", "transbath1", "transbath2", 
                                          "transsport1", "transsport2", "rr1a", "rr1b", 
                                          "rr2a", "rr2b", "rr3a", "rr3b", "rr4a", "rr4b"))]

idx_missing <- apply(R_sub, 1, function(x) {
  any(x == -7)
})
idx_missing <- unique(c((1:nrow(R_sub))[idx_missing], idx_noparty))
length(idx_missing)

R_sub <- R_sub[-idx_missing, ]
R_sub <- R_sub - 1 # to binary
table(apply(R_sub, 2, function(x) length(unique(x)))) # number of categories


#### estimation ####
K <- 3

# number of categories
C_list <- apply(R_sub, 2, function(x) {
  length(unique(x))
})

df_num <- data.frame(item=names(C_list), ans_num = C_list)

t1 <- Sys.time()
# flatten
R_flattened <- flatten(R_sub, C_list)
dim(R_flattened)

N <- nrow(R_flattened)
J <- ncol(R_flattened)
cat("N =", N, "J =", J, "K =", K, "\n")

svd_res <- svds(R_flattened, K)
U <- svd_res$u
V <- svd_res$v
d <- svd_res$d

# empirical incoherence numbers
sqrt(N) * max(apply(U, 1, function(x) sqrt(sum(x^2)))) / sqrt(K)
sqrt(J) * max(apply(V, 1, function(x) sqrt(sum(x^2)))) / sqrt(K)

# GoM estimation
r <- 10
q <- 0.4
e <- 0.2
eps <- 0
res <- gomSVD(U, V, d, J, r=r, q=q, e=e, eps=eps)
P_hat <- res$P_hat
T_hat <- res$T_hat
T_hat <- rescale(T_hat, C_list)
t2 <- Sys.time()
print(t2-t1)


#### residual correlation ####
R_hat <- P_hat %*% t(T_hat)

cormat_dep <- cov(R_flattened[,1:50] - R_hat[,1:50])
cor_dep <- melt(cormat_dep)
colnames(cor_dep)[3] <- 'covariance'
ggplot(data = cor_dep, aes(x=Var1, y=Var2, fill=covariance)) + geom_tile() +
  scale_fill_gradient(low = "white", high = "black") + xlab("") + ylab("") + 
  theme(panel.grid.major = element_blank())

# ggsave("./figures/ANES_cov_heatmap.png", width=6, height=5)
## This figure corresponds to the right panel in Figure 8

#### plot ####
library(ggtern)
library(latex2exp)
library(GGally)
library(ggpubr)

# ternary
data_tern <- as.data.frame(P_hat)
data_tern$pid3 <- as.factor(pid3[-idx_missing])

idx <- which(!(pid3[-idx_missing] %in% c(1, 2, 3))) # remove other pids
data_tern <- data_tern[-idx, ]
data_tern$Party <- sapply(data_tern$pid3, function(x) {
  if (x == 1) return ("Democrat")
  if (x == 2) return ("Republican")
  if (x == 3) return ("Independent")
})

p2 <- ggtern(data=data_tern, aes(V2, V1, V3)) + 
  geom_point(size=0.7, aes(color=Party), alpha=1) + 
  xlab("Conservative") + ylab("Indifferent") + zlab("Liberal") + 
  theme(axis.title=element_text(size=8.5), 
        legend.title= element_text(size=9), 
        legend.text=element_text(size=9), 
        legend.key.size=unit(0.6, 'cm'), 
        legend.position = "bottom", 
        tern.axis.title.L = element_text(hjust = 0.2, vjust = 0.5),
        tern.axis.title.R = element_text(hjust = 0.8, vjust = 0.5),
        legend.box.margin = margin(-35, 0, 0, 0), 
        plot.margin = margin(-15, -15, 0, -30)) +
  scale_color_manual(values=c('#377EC2', "#7FBFBB", '#e26b57'))
p2

# ggsave("./figures/tern_pid.png", width=5, height=4)
## this figure corresponds to the right panel in Figure 5


#### interpretation ####
idx_first <- c(1, cumsum(C_list)[-length(C_list)] + 1) # index of the first response
names(idx_first)[1:(length(C_list)-1)] <- names(idx_first[2:length(C_list)])
T_hat_sub <- T_hat[idx_first, ]

# select items with big difference among the extreme profiles
T_diff <- apply(T_hat_sub, 1, function(x) max(c(abs(x[1]-x[2]), abs(x[1]-x[3]), abs(x[3]-x[2]))))
selected <- sort(sort(T_diff, index.return=T, decreasing = T)$ix[1:30])
head(T_hat_sub[selected,])
colnames(R_sub)[selected]
# head(T_hat_sub[selected, ])
selected <- selected[-c(22, 23, 24, 26, 27, 28)]
selected
selected <- c(1, 2, 8, selected)
colnames(R_sub)[selected]

items_selected <- c("Do you follow what’s going on in government and public affairs? (Most of the time)", 
                    "Are you registered to vote, or not? (Yes, registered to vote at my current address)", 
                    "During the past 12 months, have you posted a message or comment online about a political issue or campaign? (Yes)", 
                    "Would you vote for Donald Trump, Joe Biden, someone else, or probably not vote? (Donald Trump)", 
                    "How important is illegal immigration in the country today? (Extremely important)", 
                    "How important is what's being taught in public schools in the country today? (Extremely important)", 
                    "Which political party would do a better job handling Covid-19? (Democrats)", 
                    "Which political party would do a better job handling illegal immigration? (Democrats)", 
                    "Which political party would do a better job handling jobs and employment? (Democrats)", 
                    "Which political party would do a better job handling the cost of living and rising prices? (Democrats)", 
                    "Which political party would do a better job handling climate change? (Democrats)", 
                    "Which political party would do a better job handling abortion? (Democrats)", 
                    "Which political party would do a better job handling gun policy? (Democrats)", 
                    "Which political party would do a better job handling taxes? (Democrats)", 
                    "Which political party would do a better job handling health care? (Democrats)", 
                    "Which political party would do a better job handling voting rights? (Democrats)", 
                    "Which political party would do a better job handling voter fraud? (Democrats)", 
                    "Do you favor, oppose, or neither favor nor oppose the Supreme Court’s decision to overturn Roe v. Wade? (Favor)",
                    "Does Roe v Wade being overturned make you feel angry? (Not at all)",
                    "Does Roe v Wade being overturned make you feel afraid? (Not at all)",
                    "Does Roe v Wade being overturned make you feel worried? (Not at all)",
                    "Does Roe v Wade being overturned make you feel outraged? (Not at all)",
                    "How important to you is protecting the right to own guns? (Extremely important)",
                    "In the 2020 presidential election, who do you believe to be the legitimate winner? (Biden)",
                    "How much responsibility do you think Donald Trump bears for the violence that occurred on 1/6/2021 at the US Capitol? (A great deal)",
                    "Regarding abortion, do you usually think of yourself as pro-choice, pro-life, or something else? (Pro-choice)", 
                    "In American society, do you think that being White comes with advantages, disadvantages, or doesn't it matter? (Advantages)")

T_hat_sub <- T_hat_sub[selected, ]
hm <- heatmap(T_hat_sub)
items_cluster <- hm$rowInd
items_coded <- factor(items_selected[items_cluster], levels=items_selected[items_cluster])

T_data <- data.frame(value=round(as.vector(T_hat_sub[items_cluster,]), 3), Item=rep(items_coded, 3), 
                     Profile=factor(c(rep("Indifferent", length(selected)), 
                                      rep("Conservative", length(selected)), 
                                      rep("Liberal", length(selected)))))

ggplot(T_data, aes(x = Profile, y = Item, fill = value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = value), color = "black", size = 2) +
  scale_fill_gradient(low = "white", high = "coral2") +
  scale_y_discrete(expand=c(0, 0)) + scale_x_discrete(expand=c(0, 0)) + ylab("") + xlab("") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 9), 
        axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8), 
        panel.border = element_blank(), panel.grid = element_blank(),
        panel.spacing = element_blank(), line = element_blank(), 
        panel.background = element_blank()) 

# ggsave("./figures/election_heatmap.pdf", width=9.5, height=5)
## this figure corresponds to Figure 7

