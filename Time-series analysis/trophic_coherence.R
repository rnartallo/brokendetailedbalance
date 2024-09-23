library(matlib) 
library(pracma)
library(igraph)

trophic_coherence <- function(W){
  N = dim(W)[1]
  col_sum = colSums(W)
  row_sum = colSums(t(W))
  u = col_sum+row_sum;
  v = col_sum-row_sum;
  Lambda = diag(u)-W-t(W);
  #Place constant equation for first node
  eq1 = numeric(N)
  eq1[1] = 1
  Lambda[1,] = eq1
  v[1] = 0
  h = solve(Lambda,v)
  num =0
  denom =0
  for (i in 1:N){
    for (j in 1:N){
      num = num + W[i,j]*(h[j]-h[i]-1)^2
      denom = denom + W[i,j]
    }
  }
  return(list(oneminusF0 = 1-num/denom, h = h))
}

trophic_coherence_dg <- function(W){
  N = dim(W)[1]
  col_sum = colSums(W)
  row_sum = colSums(t(W))
  u = col_sum+row_sum;
  v = col_sum-row_sum;
  Lambda = diag(u)-W-t(W);
  
  # Find connected components and modify the equation for one node from each cc
  #cc_info = sort_connected_components(W)
  #no_cc = cc_info$no_cc
  #cc = cc_info$cc
  c = components(graph_from_adjacency_matrix(1*(W>0),mode="directed"), mode="weak")
  no_cc = c$no
  cc = c$membership
  
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
  num =0
  denom =0
  for (i in 1:N){
    for (j in 1:N){
      num = num + W[i,j]*(h[j]-h[i]-1)^2
      denom = denom + W[i,j]
    }
  }
  return(list(oneminusF0 = 1-num/denom, h = h))
}

sort_connected_components<-function(W){
  N = dim(W)[1]
  L = diag(rowSums(W))-W
  ns = nullspace(L)
  no_cc = dim(ns)[2]
  cc = kmeans(ns[,1],no_cc)$cluster
  return(list(cc=cc,no_cc = no_cc))
}