Cal_ge_geSet <- function(SGG, RGD) {
  gene_num <- dim(RGD)[1]
  drug_num <- dim(RGD)[2]
  S_ge_geSet <- matrix(nrow=gene_num,ncol=drug_num)
  for (i in 1:gene_num) {
    for (j in 1:drug_num) {
      S_ge_geSet[i,j] <- max((SGG[i,] & RGD[,j]) *SGG[i,])
    }
  }
  return(S_ge_geSet)
}

Cal_drug_tar <- function(S_ge_geSet,RGD) {
  drug_num <- dim(S_ge_geSet)[2]
  gene_num <- dim(S_ge_geSet)[1]
  SDT <- diag(drug_num)
  for (i in 1:drug_num) {
    for (j in 1:drug_num) {
      if ((sum(RGD[,i]==1)+sum(RGD[,j]==1))!=0)
        SDT[i,j] <- (sum(S_ge_geSet[,i]*RGD[,j])+sum(S_ge_geSet[,j]*RGD[,i]))/(sum(RGD[,i]==1)+sum(RGD[,j]==1))
    }
  }
  return(SDT)
}