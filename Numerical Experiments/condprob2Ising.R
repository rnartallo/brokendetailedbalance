conditional_transition_prob_2Ising<-function(A,B,W,T){
  #This calculates the probability of transitioning to A=(a1,a2) given that you are in B=(b1,b2)
  #The Ising model has network W and temperature
  a1=A[1]; a2=A[2]; b1=B[1]; b2=B[2]; w12=W[1,2]; w21=W[2,1];
  p1 = 0.5*(1+a1)*((1+exp((-2/T)*w12*b2))^(-1))+0.5*(1-a1)*(1-(1+exp((-2/T)*w12*b2))^(-1))
  p2 = 0.5*(1+a2)*((1+exp((-2/T)*w21*b1))^(-1))+0.5*(1-a2)*(1-(1+exp((-2/T)*w21*b1))^(-1))
  return(p1*p2)
}

conditional_transition_prob_NIsing<-function(A,B,W,T){
  #This calculates the probability of transitioning to A=(a1,a2,...,aN) given that you are in B=(b1,b2,...,bN)
  #The Ising model has network W and temperature
  N = length(A)
  WB = W%*%B
  p = numeric(N)
  for (i in 1:N){
    p[i] = 0.5*(1+A[i])*((1+exp((-2/T)*WB[i]))^(-1))+0.5*(1-A[i])*(1-(1+exp((-2/T)*WB[i]))^(-1))
  }
  return(prod(p))
}