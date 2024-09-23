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






