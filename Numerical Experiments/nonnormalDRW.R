source('NonNormalityMeasures.R')
source('GenerateNonNormalNetwork.R')
library(pracma)

N =10;
#Generate network
Network = generate_non_normal_network(N,0.5,1,recip = 4)
W = Network$W
heatmap(W,Rowv=NA,Colv=NA)

P = W/rowSums(W)
e<-eigen(P)
ss <- e$vectors[,1]/sum(e$vectors[,1])

eps_range = linspace(0,1,100)
Phi=numeric(100)
for (i in 1:100){
  print(i)
  W_hat = parameterise_network(Network$W,eps_range[i])
  T = W_hat/rowSums(W_hat)
  e<-eigen(t(T))
  steady_state <- Re(e$vectors[,1]/sum(e$vectors[,1]))
  P = joint_transition_prob(steady_state,W_hat)
  Phi[i] = EPR_drw(P)
}

write.csv(Phi,'DRWPhi10Perron.csv')





# Investigating if we actually reach the steady state
N =100;
#Generate network
Network = generate_non_normal_network(N,0.5,1,recip = 4)
W = Network$W
T = W/rowSums(W)
e<-eigen(t(T))
steady_state_true <- Re(e$vectors[,1]/sum(e$vectors[,1]))
time_lengths = linspace(10,500000,100)
states = directed_random_walk(W,500000,1)
errors = numeric(100)
for (t in 1:length(time_lengths)){
  steady_state = estimate_steady_state_drw(states[1:time_lengths[t]],N)
  errors[t] = norm(steady_state-steady_state_true,type="2")
}
plot(time_lengths,errors, type = "l", lty = 1,log="xy")


