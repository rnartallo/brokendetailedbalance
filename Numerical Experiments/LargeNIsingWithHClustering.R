source('GenerateNonNormalNetwork.R')
source('IsingModel.R')
source('entropy_clusters.R')
source('NonNormalityMeasures.R')
library(MASS)
library(pracma)

N =100;
#Generate network
Network = generate_non_normal_network(N,0.5,1)
W = Network$W
heatmap(W,Rowv=NA,Colv=NA)



eps_range = linspace(0,1,100)
S =numeric(100)
H = numeric(100)
T=0.8; time_steps = 1000000; burnin=10000;
for (i in 1:100){
  print(i)
  W_hat = parameterise_network(Network$W,eps_range[i])
  spins = Ising_model(W_hat,temp=T,time_steps = time_steps,random=TRUE)
  spins = spins[,(burnin+1):time_steps]
  print('Spins done')
  L = entropy_hclusters(spins,10)
  print('Entropy and clusters done')
  H[i] = henrici_departure(L$P)
  S[i] = L$S
}
write.csv(H,'100Hwithhclust.csv')
write.csv(S,'100Ewithhclust.csv')
