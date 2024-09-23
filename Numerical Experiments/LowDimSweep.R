source('IsingModel.R')
source('countstates2Ising.R')
source('transitionmatrix.R')
source('entropyrate.R')
source('GenerateNonNormalNetwork.R')
source('NonNormalityMeasures.R')
source('condprob2Ising.R')
library(pracma)

sizes= c(2,3,4,5,6,7,8,9,10)
Sz = length(sizes)
eps_range = linspace(0,1,100)
S =matrix(0,nrow=100,ncol =Sz)
H = matrix(0,nrow=100,ncol =Sz)
nn = matrix(0,nrow=100,ncol =Sz)
E = matrix(0,nrow=100,ncol =Sz)
for (M in 1:Sz){
  N = sizes[M]
  #Generate network
  Network = generate_non_normal_network(N,0.5,1)
  W = Network$W
  heatmap(W,Rowv=NA,Colv=NA)
  T=1; time_steps = 1000000; burnin=10000;
  for (i in 1:100){
    print(paste(i,N))
    W_hat = parameterise_network(Network$W,eps_range[i])
    nn[i,M] = henrici_departure(W_hat)
    spins = Ising_model(W_hat,temp=T,time_steps = time_steps,random=TRUE)
    spins = spins[,(burnin+1):time_steps]
    clusters = kmeans(t(spins),2^N)
    states = clusters$centers
    steady_states = (clusters$size)/(time_steps-burnin)
    P = transitionmatrix(condprob = conditional_transition_prob_NIsing,steady_states,W_hat,states,T)
    H[i,M] = henrici_departure(P)
    E[i,M] = EntropyRate(P)
  }
}
write.csv(H,'LowSweepH_fixed.csv')
write.csv(E,'LowSweepE_fixed.csv')