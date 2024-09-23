generate_directed_erdos_renyi <- function(N,p){
  W = matrix(0, nrow = N,ncol=N)
  for (i in 1:N){
    for (j in 1:N){
      u = runif(1)
      if ((u<=p)&&(i!=j)){
        W[i,j]=1
      }
    }
  }
  return(W)
}