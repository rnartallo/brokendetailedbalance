library(pracma)
library(maotai)

EPR_OUP <- function(B,D){
  S = lyapunov(B,2*D)
  Q = 0.5*(B%*%S-D)
  Phi= -tr(inv(D)%*%B%*%Q)
  return(Phi)
}

EPR_NOUP <- function(W,sigma){
  N = dim(W)[1]
  I = eye(N)
  B = (I-W)
  S = lyapunov(B,2*sigma*I)
  Q = B%*%S-sigma*I
  Phi = -tr((1/sigma)*B%*%Q)
  return(Phi)
}