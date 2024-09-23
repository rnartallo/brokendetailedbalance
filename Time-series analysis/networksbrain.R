source("autoregressivefunctions.R")
source("NonNormalityMeasures.R")
source("EntropyProductionRate.R")
source("trophic_coherence.R")
source("weightedreciprocity.R")

mvts <- read.csv('Rest100PPT.csv',header = FALSE)

N=80
W = matrix(0,nrow=N,ncol=N)

for (ppt in 1:100){
  print(ppt)
  # Extract a single participant
  ppt_mvts = as.matrix(mvts[(80*ppt-79):(80*ppt),1:274])
  X=matrix(data = ppt_mvts,nrow=80,ncol = 274)
  
  nnNARM <- nonnegativenetworkautoregressivemodel(X,0.1)
  
  W = W + (nnNARM$W)/N
} 

write.matrix(W,'RestEffectiveNetwork.csv')

troph_details = trophic_coherence_dg(W)$h
write.matrix(troph_details,'RestEffectiveLevels.csv')




mvts <- read.csv('Social100PPT.csv',header = FALSE)

N=80
troph_participants = matrix(0,nrow=N,ncol=100)

for (ppt in 1:100){
  print(ppt)
  # Extract a single participant
  ppt_mvts = as.matrix(mvts[(80*ppt-79):(80*ppt),1:274])
  X=matrix(data = ppt_mvts,nrow=80,ncol = 274)
  
  nnNARM <- nonnegativenetworkautoregressivemodel(X,0.1)
  
  troph_participants[,ppt] = trophic_coherence_dg(nnNARM$W)$h
} 

write.csv(troph_participants,'SocialTrophicPPT.csv')

