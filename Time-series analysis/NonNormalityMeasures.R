library(psych)
library(expm)

henrici_departure <- function(W){
  return(sqrt(norm(W,"F")^2 - sum(abs(eigen(W)$values)^2)))
}

normalised_henrici_departure <- function(W){
  return((sqrt(norm(W,"F")^2 - sum(abs(eigen(W)$values)^2)))/norm(W,"F"))
}

asymmetry_measure <- function(W){
  N = dim(W)[1]
  num=0
  den=0
  for (j in (2:N)){
    for (i in (1:(j-1))){
      num = num+W[i,j]-W[j,i]
      den = den+W[i,j]+W[j,i]
    }
  }
  return (abs(num)/den)
}

mackkay_index <- function(W){
  return(tr((t(W)%*%W)%^%2)-tr((t(W)%^%2)%*%(W%^%2)))
}

normalised_mackkay_index <- function(W){
  return((tr((t(W)%*%W)%^%2)-tr((t(W)%^%2)%*%(W%^%2)))/(tr((t(W)%*%W)%^%2)))
}