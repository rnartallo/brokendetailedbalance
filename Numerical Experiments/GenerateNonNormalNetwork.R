generate_non_normal_network <- function(N,mu,p0,recip=10){
  W = matrix(0,nrow=N,ncol = N)
  A = matrix(0,nrow=N,ncol = N)
  for (n in 2:N){
    deg = rowSums(A)
    rands = runif(n,min=0,max=1)
    for (m in 1:(n-1)){
      if (sum(deg)>0){
        p = deg[m]/sum(deg) + mu
      }
      if (sum(deg)==0){
        p = p0
      }
      if (rands[m]<p){
        weight = runif(1,min=0,max=1)
        W[n,m] = weight
        W[m,n] = weight/recip
        A[n,m] = 1
        A[m,n] = 1
      }
    }
  }
  return(list(W=W,A=A))
}

parameterise_network <- function(W,eps){
  return((1-eps)*0.5*(W+t(W))+eps*W)
}