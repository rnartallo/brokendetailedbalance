transitionmatrix <- function(condprob,steady_states,W,states,T){
  no_states = dim(states)[1]
  P = matrix(0,nrow=no_states,ncol = no_states)
  for (i in 1:no_states){
    for (j in 1:no_states){
      P[i,j] = condprob(states[j,],states[i,],W,T)*steady_states[i]
    }
  }
  return(P)
}
