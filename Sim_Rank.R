#normalize matrix along the row vector
row_nor <- function(input) {
  rows <- dim(input)[1]
  cols <- dim(input)[2]
  output <- matrix(nrow = rows,ncol = cols)
  for (i in 1:rows) {
    s <- sum(input[i,])
    for (j in 1:cols) {
      if (s!=0) {
        output[i,j] <- input[i,j]/s
      } else {
        output[i,j] <- 0
      }
    }
  }
  return(output)
}

#normalize matrix along the column vector
col_nor <- function(input) {
  rows <- dim(input)[1]
  cols <- dim(input)[2]
  output <- matrix(nrow = rows, ncol = cols)
  for (i in 1:rows) {
    s <- sum(input[,i])
    for (j in 1:cols) {
      output[i,j] <- input[i,j]/s
    }
  }
  return(output)
}

#decomposition of atomic relation for odd-length meta-path
split_mat <- function(input) {
  rows <- dim(input)[1]
  cols <- dim(input)[2]
  mid <- sum(input)
  output1 <- matrix(nrow=rows,ncol=mid)
  output1[is.na(output1)] <- 0
  output2 <- matrix(nrow=mid,ncol=cols)
  output2[is.na(output2)] <- 0
  start <- 1
  for (i in 1:rows) {
    num <- sum(input[i,])
    if (num != 0) {
      for (j in start:(start+num-1)) {
        output1[i,j] <- 1
      }
      start <- start + num
    }
  }
  start <- 1
  for (i in 1:cols) {
    num <- sum(input[,i])
    if (num != 0) {
      for (j in start:(start+num-1)) {
        output2[j,i] <- 1
      }
      start <- start + num
    }
  }
  output <- list(output1,output2)
  return(output)
}

#normalization
final_nor <- function(input1, input2) {
  rows1 <- dim(input1)[1]
  cols1 <- dim(input1)[2]
  rows2 <- dim(input2)[1]
  cols2 <- dim(input2)[2]
  output <- matrix(nrow = rows1, ncol = cols2)
  for (i in 1:rows1) {
    for (j in 1:cols2) {
      f1 <- sum(input1[i,]*input2[,j])
      f2 <- sqrt(sum(input1[i,]^2) * sum(input2[,j]^2))
      if (f2!=0) {
        output[i,j] <- f1/f2
      } else {
        output[i,j] <- 0
      }
      
    }
  }
  return(output)
}

Sim_rank <- function(RGG, len) {
  gene_num <- dim(RGG)[1]
  comp_mat1 <- diag(gene_num)
  sim_mat <- diag(gene_num)
  if (len%%2==1) {
    so <- floor(len/2)
    temp <- row_nor(RGG)
    for (j in 1:so) {
      comp_mat1 <- row_nor(comp_mat1)
      comp_mat1 <- comp_mat1%*%temp
    }
    AB <- split_mat(RGG)
    A <- AB[[1]]
    B <- AB[[2]]
    comp_mat1 <- row_nor(comp_mat1)
    A <- row_nor(A)
    comp_mat1 <- comp_mat1%*%A
    sim_mat <- final_nor(comp_mat1,t(comp_mat1))
  }
  else if (len%%2==0) {
    so <- floor(len/2)
    temp <- row_nor(RGG)
    for (j in 1:so) {
      comp_mat1 <- row_nor(comp_mat1)
      comp_mat1 <- comp_mat1%*%temp
    }
    sim_mat <- final_nor(comp_mat1,t(comp_mat1))
  }
  return(sim_mat)
}








