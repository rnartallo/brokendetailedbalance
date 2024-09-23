source('IsingModel.R')
source('countstates2Ising.R')
source('transitionmatrix.R')
source('entropyrate.R')
source('GenerateNonNormalNetwork.R')
source('NonNormalityMeasures.R')
source('condprob2Ising.R')
library(pracma)

N =5;
#Generate network
Network = generate_non_normal_network(N,0.2,1)
W = Network$W
heatmap(W,Rowv=NA,Colv=NA)
states_enumerated = expand.grid(replicate(N, c(-1,1), simplify = FALSE))
SE = list()
for (j in 1:(2^N)){
  SE[[j]] = as.numeric(states_enumerated[j,])
}
SE = unlist(SE)
eps_range = linspace(0,1,100)
S =numeric(100)
H = numeric(100)
nn = numeric(100)
E = numeric(100)
T=1; time_steps = 1000000; burnin=10000;
for (i in 1:100){
  print(i)
  W_hat = parameterise_network(Network$W,eps_range[i])
  nn[i] = henrici_departure(W_hat)
  spins = Ising_model(W_hat,temp=T,time_steps = time_steps,random=TRUE)
  spins = spins[,(burnin+1):time_steps]
  clusters = kmeans(t(spins),2^N)
  states = clusters$centers
  steady_states = (clusters$size)/(time_steps-burnin)
  P = transitionmatrix(condprob = conditional_transition_prob_NIsing,steady_states,W_hat,states,T)
  H[i] = henrici_departure(P)
  E[i] = EntropyRate(P)
}
write.csv(H,'2HofP.csv')
write.csv(E,'2E.csv')
