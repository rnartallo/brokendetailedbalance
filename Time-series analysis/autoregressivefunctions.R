library(nnls)
library(mAr)


nonnegativenetworkautoregressivemodel <- function(X,time_step){
  N = dim(X)[1]
  T = dim(X)[2]
  matrixX = t(X[,1:(T-1)])
  W = matrix(0,nrow=N,ncol=N)
  xi = matrix(0,nrow=N,ncol=(T-1))
  for (i in 1:N){
    a = nnls(matrixX,X[i,2:T]-(1-time_step)*X[i,1:(T-1)])
    W[,i] = (a$x)/time_step
    xi[i,] = a$residuals
  }
  C = cov(t(xi))
  sigma = mean(diag(C/(2*time_step)))
  return(list(W=W,resid=xi,C=C,sigma=sigma))
}

autoregressivemodel <- function(X,time_step){
  N = dim(X)[1]
  T = dim(X)[2]
  AR1 = mAr.est(t(X),1,nonnegative=TRUE)
  B = (1/time_step)*(eye(N)-AR1$AHat)
  C = cov(AR1$resid)
  D = C/(2*time_step)
  return(list(B=B,resid=AR1$resid,C=C,D=D))
}


nonnegativeautoregressivemodel <- function(X){
  N = dim(X)[1]
  T = dim(X)[2]
  matrixX = t(X[,1:(T-1)])
  A = matrix(0,nrow=N,ncol=N)
  xi = matrix(0,nrow=N,ncol=(T-1))
  for (i in 1:N){
    a = nnls(matrixX,X[i,2:T])
    A[i,] = a$x
    xi[i,] = a$residuals
  }
  return(list(A=A,resid=xi,noisecov=cov(t(xi))))
}