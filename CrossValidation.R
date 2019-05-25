P_positive <- function(DD_mat) {
  drug_num <- dim(DD_mat)[1]
  DD_mat_O <- DD_mat
  p_in_all <- matrix(nrow=(sum(DD_mat_O)-drug_num)/2,ncol=3)
  x <- 1
  p <- 1
  for (i in 2:drug_num) {
    for (j in 1:(i-1)) {
      if (DD_mat_O[i,j]==1) {
        p_in_all[p,1] <- i
        p_in_all[p,2] <- j
        p_in_all[p,3] <- x
        p <- p+1
      }
      x <- x+1
    }
  }
  return(p_in_all)
}

Get_Test_Score <- function(F, DD_mat, DD_mat_O) {
  drug_num <- dim(F)[1]
  train_n <- c()
  n_pos <- 1
  probability <- c()
  label <- c()
  list_pos <- 1
  for (m in 2:drug_num) {
    for (n in 1:(m-1)) {
      train_n[n_pos] <- F[m,n]
      n_pos <- n_pos + 1
      if (DD_mat[m,n]==0) {
        probability[list_pos] <- F[m,n]
        label[list_pos] <- DD_mat_O[m,n]
        list_pos <- list_pos + 1
      }
    }
  }
  data_ROC_n <- data.frame(prob=probability,obs=label)
  data_ROC_n <- data_ROC_n[order(-data_ROC_n$prob),]
  return (data_ROC_n)
}

Get_Test_Score_List <- function(F, DD_mat, DD_mat_O) {
  drug_num <- dim(F)[1]
  train_n <- c()
  n_pos <- 1
  probability <- c()
  label <- c()
  row_n <- c()
  col_n <- c()
  list_pos <- 1
  for (m in 2:drug_num) {
    for (n in 1:(m-1)) {
      train_n[n_pos] <- F[m,n]
      n_pos <- n_pos + 1
      if (DD_mat[m,n]==0) {
        probability[list_pos] <- F[m,n]
        label[list_pos] <- DD_mat_O[m,n]
        row_n[list_pos] <- m
        col_n[list_pos] <- n
        list_pos <- list_pos + 1
      }
    }
  }
  data_ROC_n <- data.frame(prob=probability,r=row_n,c=col_n,obs=label)
  data_ROC_n <- data_ROC_n[order(-data_ROC_n$prob),]
  return (data_ROC_n)
}

Get_fpr_tpr_pre <- function(data_ROC_n) {
  n <- dim(data_ROC_n)[1]
  tpr_n <- fpr_n <- pre_n <- rep(0,n+1)
  for ( j in 1:n) {
    threhold <- data_ROC_n$prob[j]
    tp_n <- sum(data_ROC_n$prob > threhold & data_ROC_n$obs == 1)
    fp_n <- sum(data_ROC_n$prob > threhold & data_ROC_n$obs == 0)
    tn_n <- sum(data_ROC_n$prob <= threhold & data_ROC_n$obs == 0)
    fn_n <- sum(data_ROC_n$prob <= threhold & data_ROC_n$obs == 1)
    tpr_n[j] <- tp_n/(tp_n+fn_n)
    fpr_n[j] <- fp_n/(tn_n+fp_n)
    pre_n[j] <- tp_n/(tp_n+fp_n)
  }
  pre_n[1] <- 0
  fpr_n[n+1] <- 1
  tpr_n[n+1] <- 1
  fpr_tpr_pre <- data.frame(fpr_n,tpr_n,pre_n)
}




