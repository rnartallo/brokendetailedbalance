Ising_model <- function(W,temp,time_steps,random){
  N = dim(W)[1]
  X = matrix(0,nrow=N,ncol = time_steps)
  #Random initial state
  if (random){
    X[,1]=sample(c(-1,1),size = N,replace = TRUE)
  }
  else{
    X[,1]=X[,1]+1
  }
  for (t in 2:time_steps){
    exponents = W%*%X[,t-1]
    rands = runif(N,0,1)
    for (alpha in 1:N){
      p = (1+exp(-(2/temp)*(exponents[alpha])))^(-1)
      if (rands[alpha]<p){
        X[alpha,t] = 1
      }
      else{
        X[alpha,t]= -1
      }
    }
  }
  return(X)
}