count_states <- function(spins,states){
  S = dim(states)[1]
  time_steps = dim(spins)[2]
  count = c(0,0,0,0)
  for (t in 1:time_steps){
    for (s in 1:S){
      if (setequal(spins[,t],states[s,])){
        count[s]=count[s]+1
      }
    }
  }
  return(count)
}
