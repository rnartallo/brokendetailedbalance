source('NonNormalityMeasures.R')
source('GenerateNonNormalNetwork.R')
source('DRW.R')
library(pracma)
library(maotai)
library(matlib)
library(expm)

N =1000;
#Generate network
Network = generate_non_normal_network(N,0.5,1,recip = 4)
W = Network$W
heatmap(W,Rowv=NA,Colv=NA)

T = W/rowSums(W)
e<-eigen(T)
ss <- e$vectors[,1]/sum(e$vectors[,1])

bidirectional_transition_matrix <- function(N,T){
  B = matrix(0,nrow=N,ncol=N)
  for (i in 1:N){
    for (j in 1:N){
      if ((T[i,j]>0) & (T[j,i]>0)){
        B[i,j] = T[i,j]
      }
    }
  }
  return(B)
}

B = bidirectional_transition_matrix(N,T)

# varying non-normality
eps_range = linspace(0,1,100)
Phi=numeric(100)
for (k in 1:100){
  print(k)
  W_hat = parameterise_network(Network$W,eps_range[k])
  T = W_hat/rowSums(W_hat)
  B = bidirectional_transition_matrix(N,T)
  e<-eigen(t(T))
  ss <- Re(e$vectors[,1]/sum(e$vectors[,1]))
  sigma = 0
  for (i in 1:N){
    for (j in i:N){
      if (B[i,j]!=0){
        sigma = sigma + B[j,i]*ss[j]*log((B[j,i]*ss[j])/(B[i,j]*ss[i]))
      }
    }
  }
  Phi[k] = sigma
}

plot(eps_range,Phi,xlab='Epsilon',ylab='EPR')
write.csv(Phi,"CTRWPhi1000.csv")



