rm(list=ls())

setwd("...")
RDS <- read.table("***.txt")
SDA1 <- read.table("***.txt")
SDA2 <- read.table("***.txt")
SDA3 <- read.table("***.txt")
SDC <- read.table("***.txt")
RDG <- read.table("***.txt")
RGG <- read.table("***.txt")
RGG <- as.matrix(RGG)
RDS <- as.matrix(RDS)
SDA1 <- as.matrix(SDA1)
SDA2 <- as.matrix(SDA2)
SDA3 <- as.matrix(SDA3)
SDA <- (SDA1+SDA2+SDA3)/3
SDC <- as.matrix(SDC)

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

drug_num <- dim(RDS)[1]
lamdaA <- 0.1
lamdaS <- 0.3
lamdaT <- 0.01

RDS_v <- RDS

LDT <- Lapla_nor(SDT)
LDS <- Lapla_nor(RDS_v)
LDA <- Lapla_nor(SDA)
LDC <- Lapla_nor(SDC)
          
I <- diag(drug_num)
F <- solve(I+lamdaA*(I-LDA)+lamdaS*(I-LDS)+lamdaT*(I-LDT))%*%(RDS_v)
F <- F+t(F)          
data_ROC_n <- Get_Test_Score_List(F,RDS_v,RDS)
          
