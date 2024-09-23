source("autoregressivefunctions.R")
source("NonNormalityMeasures.R")
source("EntropyProductionRate.R")
source("trophic_coherence.R")
source("weightedreciprocity.R")

stock <- subset(read.csv('NYSE/stock.csv',header = TRUE),select = -c(Date))
stock_names <- subset(read.csv('NYSE/list_stocknames.txt'),select = -c(X0))

X=t(data.matrix(frame = stock))



nnNARM <- nonnegativenetworkautoregressivemodel(X,0.1)
ARM <- autoregressivemodel(X,0.1)
Phi = EPR_OUP(ARM$B,ARM$D)
W = nnNARM$W
write.matrix(W,'financenetwork.csv',sep = ',')
N = dim(W)[1]
sigma = nnNARM$sigma


trophiccoh = trophic_coherence_dg(W)
trophiccoh$oneminusF0
write.matrix(trophiccoh$h,'financelevels.csv',sep = ',')


# trophic d 0.4047484
# henrici 1.533352
# n henrici 0.1602408
# irreciprocity 0.9753964



