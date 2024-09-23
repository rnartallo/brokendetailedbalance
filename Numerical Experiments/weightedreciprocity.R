weighted_reciprocity <- function(W){
  N = dim(W)[1]
  W_symm = matrix(0,nrow=N,ncol=N)
  W_right = matrix(0,nrow=N,ncol=N)
  Wcap = 0
  Wcap_symm = 0
  for (i in 1:N){
    for (j in 1:N){
      W_symm[i,j] = min(W[i,j],W[j,i])
      if (i!=j){
        Wcap = Wcap+W[i,j]
        Wcap_symm = Wcap_symm+W_symm[i,j]
      }
    }
  }
  return(Wcap_symm/Wcap)
}