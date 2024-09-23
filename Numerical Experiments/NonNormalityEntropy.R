source('GenerateNonNormalNetwork.R')
source('IsingModel.R')
source('TEFunctionalNetwork.R')
source('entropy_clusters.R')
source('NonNormalityMeasures.R')
library(MASS)
library(pracma)


Network = generate_non_normal_network(100,0.5,1)
W = Network$W
heatmap(W,Rowv=NA,Colv=NA)
eps_range = linspace(0,1,100)
S =numeric(100)
H = numeric(100)
nn = numeric(100)
for (i in 1:100){
  print(i)
  W_hat = parameterise_network(Network$W,eps_range[i])
  nn[i] = henrici_departure(W_hat)
  spins = Ising_model(W_hat,temp=1,time_steps = 1000000,random=TRUE)
  spins = spins[,10000:1000000]
  L = entropy_clusters(spins,10)
  H[i] = henrici_departure(L$P)
  S[i] = L$S
}
write.csv(W,'W3.csv')
write.csv(H,'H3.csv')
write.csv(S,'S3.csv')