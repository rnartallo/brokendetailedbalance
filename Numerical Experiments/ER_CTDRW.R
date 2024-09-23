source('DirectedErdosRenyi.R')
source('DRW.R')
source('bidirectionaltransitionmatrix.R')

library(pracma)
library(maotai)
library(matlib)
library(expm)
library(igraph)
library(matrixStats)

N=100
p_range = linspace(0,1,101)
trials = 100

CTRW=matrix(0,nrow=92,ncol=trials)
DTRW=matrix(0,nrow=92,ncol=trials)


for (k in 10:101){
  print(k)
  trial = 1
  while (trial<=trials){
    sigma=0
    W_hat = generate_directed_erdos_renyi(N,p_range[k])
    c = components(graph_from_adjacency_matrix(W_hat,mode="directed"), mode="strong")
    no_cc = c$no
    if (no_cc==1){
      T = W_hat/rowSums(W_hat)
      
      # Continuous time
      B = bidirectional_transition_matrix(N,T)
      e<-eigen(t(T))
      ss <- Re(e$vectors[,1]/sum(e$vectors[,1]))
      for (i in 1:(N)){
        for (j in i:N){
          if (B[i,j]!=0){
            sigma = sigma + B[j,i]*ss[j]*log((B[j,i]*ss[j])/(B[i,j]*ss[i]))
          }
        }
      }
      CTRW[k-9,trial] = sigma
      
      # Discrete time
      #P = joint_transition_prob(ss,W_hat)
      #DTRW[k-9,trial] = EPR_drw(P)
      
      trial = trial +1
    }
  }
}
write.csv(CTRW,"CTRW_ER_100_100trials.csv")
#write.csv(DTRW,"DTRW_ER_100_100trials.csv")


plot(rowMedians(CTRW))
plot(rowMedians(DTRW))
