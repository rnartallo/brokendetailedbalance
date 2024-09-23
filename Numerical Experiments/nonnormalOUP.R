source('NonNormalityMeasures.R')
source('GenerateNonNormalNetwork.R')
library(maotai)
library(matlib)
library(pracma)
library(matrixcalc)

N =10;
#Generate network
Network = generate_non_normal_network(N,0.5,1)
W = Network$W
heatmap(W,Rowv=NA,Colv=NA)

Theta = 0.5
gamma = 0.005
sigma = 0.5
I = diag(N)
B = Theta*(I-gamma*W)
print(all(Re(eigen(B)$values)>0))
S = lyapunov(B,2*sigma*I)
heatmap(S,Rowv=NA,Colv=NA)




Q = B%*%S-sigma*I
S_inv = inv(S)
Phi = -tr((1/sigma)*B%*%Q)

eps_range = linspace(0,1,100)
Phi=numeric(100)
for (i in 1:100){
  print(i)
  W_hat = parameterise_network(Network$W,eps_range[i])
  B_hat = Theta*(I-gamma*W_hat)
  print(all(Re(eigen(B_hat)$values)>0))
  S = lyapunov(B_hat,2*sigma*I)
  Q = B_hat%*%S-sigma*I
  Phi[i] = -tr((1/sigma)*B_hat%*%Q)
}
write.csv(Phi,'OUPhi1000.csv')

