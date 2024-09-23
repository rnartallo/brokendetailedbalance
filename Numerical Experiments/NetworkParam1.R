library(pracma)
library(psych)
source('GenerateNonNormalNetwork.R')
source('NonNormalityMeasures.R')

# Generates non-normal networks
Network = generate_non_normal_network(100,0.5,1)

# Measure the departure from non-normality
d_HW = henrici_departure(Network$W)
d_DeltaW = asymmetry_measure(Network$W)

# Parameterise network
d_Delta = numeric(1000)
d_H = numeric(1000)
eps_range = linspace(0,1,1000)
#Show that parameterising network in this fashion parametrises non-normality
for (i in 1:1000){
  W_hat = parameterise_network(Network$W,eps_range[i])
  d_Delta[i] = asymmetry_measure(W_hat)/d_DeltaW
  d_H[i] = henrici_departure(W_hat)/d_HW
}
write.csv(d_H,'d_H.csv')
write.csv(d_Delta,'d_Delta.csv')


