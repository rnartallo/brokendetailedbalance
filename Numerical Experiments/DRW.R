directed_random_walk <- function(W,time_steps,initial_state){
  N = dim(W)[1]
  P = W/rowSums(W)
  states = numeric(time_steps)
  states[1]=initial_state
  for (t in 2:time_steps){
    u = runif(1)
    cumulativeP = cumsum(P[states[t-1],])
    transition_found = FALSE
    check = N
    while (!transition_found){
      if ((u<=cumulativeP[check]) & (u>cumulativeP[check-1])){
        transition_found = TRUE
        move = check
      }
      else{
        check = check-1
        if (check==1){
          transition_found = TRUE
          move = 1
        }
      }
    }
    states[t] = move
  }
  return(states)
}

estimate_steady_state_drw <-function(states,N){
  state_count = numeric(N)
  total_count = length(states)
  for (s in 1:total_count){
    state_count[states[s]]=state_count[states[s]]+1
  }
  return(state_count/total_count)
}

joint_transition_prob <-function(pi,W){
  T = W/rowSums(W)
  N = dim(W)[1]
  P = matrix(0,nrow=N,ncol = N)
  for (i in 1:N){
    for (j in 1:N){
      P[i,j]=T[i,j]*pi[i]
    }
  }
  return(P)
}

EPR_drw <- function(P){
  N = dim(P)[1]
  M = 0.5*(P+t(P))
  S=0
  for (i in 1:N){
    for (j in 1:N){
      if (P[i,j]>0){
        S=S+0.5*(P[i,j]*log(P[i,j]/M[i,j]))
      }
      if (P[j,i]>0){
        S=S+0.5*(P[j,i]*log(P[j,i]/M[i,j]))
      }
    }
  }
  return(S)
}
