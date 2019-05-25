Lapla_nor <- function(S) {
  rL <- rowSums(S)
  for (i in 1:length(rL)) {
    rL[i] <- rL[i]^(-1.0/2.0)
  }
  DS <- diag(rL)
  L <- DS%*%S%*%DS
  return(L)
}


Lapla_matrix <- function(S) {
  rL <- rowSums(S)
  D <- diag(rL)
  L <- D-S
  return(L)
}


