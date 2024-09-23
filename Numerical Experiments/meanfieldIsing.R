source('NonNormalityMeasures.R')
source('GenerateNonNormalNetwork.R')
library(maotai)
library(matlib)
library(pracma)
library(matrixcalc)

N =1000;
# Generate network
Network = generate_non_normal_network(N,0.2,1,N)
W = Network$W
heatmap(W,Rowv=NA,Colv=NA)

# Start with all units active
m_prev = numeric(N)+1
D = matrix(0,nrow = N,ncol=N)
iterations = 10; burnin=10000;
eps_range = linspace(0,1,100)
Phi=matrix(0,nrow=100,ncol=iterations)
for (i in 1:100){
  print(i)
  W_hat = parameterise_network(Network$W,eps_range[i])
  for (t in 1:iterations){
    v = W_hat%*%m_prev
    m_new = tanh(v)
    for (j in 1:N){
      for (l in 1:N){
        D[j,l] = W_hat[j,l]*(1-m_new[j]^2)*(1-m_prev[l]^2)
        Phi[i,t] = Phi[i,t] + (W_hat[j,l]-W_hat[l,j])*D[j,l]
      }
    }
    m_prev = m_new
  }
}
heatmap(Phi,Rowv=NA,Colv=NA)
plot(Phi[,10],xlab='Epsilon index',ylab='EPR')
title('N=100 Ising model with mean-field')
write.csv(Phi[,10],'N1000IsingMF_snn.csv')
