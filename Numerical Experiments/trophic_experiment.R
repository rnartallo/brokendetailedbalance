source('trophic_coherence.R')
source('DirectedErdosRenyi.R')
W = matrix(0,nrow=3,ncol=3)
W[1,2]=2
W[2,1]=1

N = dim(W)[1]
col_sum = colSums(W)
row_sum = colSums(t(W))
u = col_sum+row_sum;
v = col_sum-row_sum;
Lambda = diag(u)-W-t(W);

# Find connected components and modify the equation for one node from each cc
cc_info = sort_connected_components(W)
no_cc = cc_info$no_cc
cc = cc_info$cc

# Modify the linear system to guarantee a unique solution 
for (n in (1:no_cc)){
  found=FALSE
  check = 1
  while (!found){
    if (cc[check] == n){
      eq1 = numeric(N)
      eq1[check] = 1
      Lambda[check,] = eq1
      v[check] = 0
      found = TRUE
    }
    check = check+1
  }
}
h = solve(Lambda,v)
for (n in (1:no_cc)){
  mask = cc==n
  h[mask] = h[mask]-min(h[mask])
}
