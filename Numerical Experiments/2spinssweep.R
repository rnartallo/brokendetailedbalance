source('condprob2Ising.R')
source('IsingModel.R')
source('countstates2Ising.R')
source('transitionmatrix.R')
source('entropyrate.R')
library(pracma)

#Enumerate states
states = matrix(0,nrow=4,ncol = 2)
states[1,1] = 1; states[1,2]=1
states[2,1] = 1; states[2,2]=-1
states[3,1] = -1; states[3,2]=1
states[4,1] = -1; states[4,2]=-1
#Build weight matrix
N=2
W = matrix(0,nrow=N,ncol = N)
T=1; time_steps = 1000000; burnin=10000;
weight_range = logseq(0.1,4,40)
Entropies = matrix(0,nrow=40,ncol = 40)
for (w12 in 1:40){
  print(w12)
  for (w21 in 1:40){
    W[1,2] = weight_range[w12]
    W[2,1] = weight_range[w21]
    spins = Ising_model(W,temp=T,time_steps = time_steps,random=TRUE)
    spins = spins[,(burnin+1):time_steps]
    steady_states = count_states(spins,states)/(time_steps-burnin)
    P = transitionmatrix(condprob = conditional_transition_prob_2Ising,steady_states,W,states,T)
    E = EntropyRate(P)
    Entropies[w12,w21] = E
    }
}
write.csv(Entropies,'2SpinSweepLog_fixed.csv')
