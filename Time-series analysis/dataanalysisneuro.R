source("autoregressivefunctions.R")
source("NonNormalityMeasures.R")
source("EntropyProductionRate.R")
source("trophic_coherence.R")
source("weightedreciprocity.R")

mvts <- read.csv('Rest100PPT.csv',header = FALSE)
mvts <- read.csv('Motor100PPT.csv',nrows=8000)


trophicd <- numeric(100)
henrici <- numeric(100)
irrecip <- numeric(100)
nhenrici <- numeric(100)
Phi <-numeric(100)
N=80
for (ppt in 1:100){
  print(ppt)
  # Extract a single participant
  ppt_mvts = as.matrix(mvts[(80*ppt-79):(80*ppt),1:274])
  X=matrix(data = ppt_mvts,nrow=80,ncol = 274)
  
  nnNARM <- nonnegativenetworkautoregressivemodel(X,0.01)
  ARM <- autoregressivemodel(X,0.01)
  
  W = nnNARM$W
  sigma = nnNARM$sigma
  
  trophicd[ppt] = trophic_coherence(W)$oneminusF0
  henrici[ppt] = henrici_departure(W)
  nhenrici[ppt] = normalised_henrici_departure(W)
  irrecip[ppt] = 1- weighted_reciprocity(W)
  
  B = ARM$B
  D = ARM$D
  print(all(Re(eigen(B)$values)>0))

  Phi[ppt] = EPR_OUP(B,D)
} 

social = list(trophicd = trophicd,henrici=henrici,nhenrici=nhenrici,irrecip = irrecip,Phi = Phi)
motor = list(trophicd = trophicd,henrici=henrici,nhenrici=nhenrici,irrecip = irrecip,Phi = Phi)
rest = list(trophicd = trophicd,henrici=henrici,nhenrici=nhenrici,irrecip = irrecip,Phi = Phi)


