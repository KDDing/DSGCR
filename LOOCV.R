rm(list=ls())

setwd("...")
RDS <- read.table("***.txt")
SDA1 <- read.table("***.txt")
SDA2 <- read.table("***.txt")
SDA3 <- read.table("***.txt")
RDG <- read.table("***.txt")
RGG <- read.table("***.txt")
RGG <- as.matrix(RGG)
RDS <- as.matrix(RDS)
SDA1 <- as.matrix(SDA1)
SDA2 <- as.matrix(SDA2)
SDA3 <- as.matrix(SDA3)
SDA <- (SDA1+SDA2+SDA3)/3


setwd("...")
source("BMA.R")
source("Lapla.R")
source("CrossValidation.R")
source("Sim_Rank.R")

SGG1 <- Sim_rank(RGG,2)
SGG2 <- Sim_rank(RGG,3)
SGG3 <- Sim_rank(RGG,4)
SGG <- (1/3)*(SGG1+SGG2+SGG3)

S_ge_geSet <- Cal_ge_geSet(SGG, t(RDG))
SDT <- Cal_drug_tar(S_ge_geSet,t(RDG))

runs <- 1
drug_num <- dim(RDS)[1]
p_in_all <- P_positive(RDS)
K <- dim(p_in_all)[1]
num_te <- floor(dim(p_in_all)[1]/K)

TPR_ALL_N <- matrix(nrow=(drug_num*(drug_num-1)/2-(dim(p_in_all)[1]-num_te))+1,ncol=runs*floor(K))
FPR_ALL_N <- matrix(nrow=(drug_num*(drug_num-1)/2-(dim(p_in_all)[1]-num_te))+1,ncol=runs*floor(K))
PRE_ALL_N <- matrix(nrow=(drug_num*(drug_num-1)/2-(dim(p_in_all)[1]-num_te))+1,ncol=runs*floor(K))

lamdaA <- 0.1
lamdaS <- 0.3
lamdaT <- 0.01

    
for (r in 1:runs) {
        
  RDS_v <- RDS
  tep_pos_set <- sample(dim(p_in_all)[1],dim(p_in_all)[1])
  tep_pos_set_all <- rep(0,num_te)
  for (i in 1:floor(K)) {
    t_p <- 1
    RDS_v <- RDS
    for (j in ((i-1)*num_te+1):(i*num_te)) {
      RDS_v[p_in_all[tep_pos_set[j],1],p_in_all[tep_pos_set[j],2]] <- 0
      RDS_v[p_in_all[tep_pos_set[j],2],p_in_all[tep_pos_set[j],1]] <- 0
    }
          
    LDT <- Lapla_nor(SDT)
    LDS <- Lapla_nor(RDS_v)
    LDA <- Lapla_nor(SDA)
    LDC <- Lapla_nor(SDC)
    I <- diag(drug_num)
    F <- (RDS_v)%*%solve(I+lamdaA*(I-LDA)+lamdaS*(I-LDS)+lamdaT*(I-LDT))
    F <- (F+t(F))*0.5
          
    data_ROC_n <- Get_Test_Score(F,RDS_v,RDS)
    FTP <- Get_fpr_tpr_pre(data_ROC_n)
    TPR_ALL_N[,((r-1)*floor(K)+i)] <- FTP$tpr_n
    FPR_ALL_N[,((r-1)*floor(K)+i)] <- FTP$fpr_n
    PRE_ALL_N[,((r-1)*floor(K)+i)] <- FTP$pre_n
  }
}
tpr_p <- rowMeans(TPR_ALL_N)
fpr_p <- rowMeans(FPR_ALL_N)
pre_p <- rowMeans(PRE_ALL_N)
library("caTools")
AUC <- trapz(fpr_p,tpr_p)
 







