hkmeans <- function(spins,k){
  N = dim(spins)[1]
  T = dim(spins)[2]
  #Split original cluster into 2
  clusters = kmeans(t(spins),2)
  no_clusters = 2
  state_series = clusters$cluster
  #Find cluster with maximal internal distance
  distances = clusters$withinss
  max_diff = which.max(distances)
  cluster_to_split_indices = which(state_series==max_diff)
  ones_to_replace = which(state_series==1)
  #Relabel cluster to split as cluster 1
  state_series[ones_to_replace] = max_diff
  temp = distances[max_diff]
  distances[max_diff] = distances[1]
  distances[1] = temp
  state_series[cluster_to_split_indices] = 1
  new_spins = spins[,cluster_to_split_indices]
  while (no_clusters<k){
    clusters = kmeans(t(new_spins),2)
    no_clusters = no_clusters + 1
    partial_state_series = clusters$cluster
    twos_to_replace = which(state_series==2)
    state_series[twos_to_replace] = no_clusters
    state_series[cluster_to_split_indices] = partial_state_series
    distances = append(distances,distances[2])
    distances[1] = max(clusters$withinss)
    distances[2] = min(clusters$withinss)
    max_diff = which.max(distances)
    cluster_to_split_indices = which(state_series==max_diff)
    ones_to_replace = which(state_series==1)
    #Relabel cluster to split as cluster 1
    state_series[ones_to_replace] = max_diff
    temp = distances[max_diff]
    distances[max_diff] = distances[1]
    distances[1] = temp
    state_series[cluster_to_split_indices] = 1
    new_spins = spins[,cluster_to_split_indices]
  }
  return(state_series)
}
