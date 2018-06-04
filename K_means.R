K_means=function(MAT, K, measure_type, R, SEED){
  #' @param MAT data matrix
  #' @param K number of clusters
  #' @param measure_type choice between "Correlation" or "Euclidean" distance in clustering
  #' @param R number of restarts of the k-means algorithm
  #' @param SEED Random number generator SEED
  #' @return list containing the genes grouped in K clusters and K centroids
  #' 
  #' As an EXTERNAL termination condition we impose that Restarts are reached,
  #' as an INTERNAL condition that none of the genes changes cluster in the entire FOR cycle
  
  source('measure.R')
  set.seed(SEED)
  N <- nrow(MAT)
  M <- ncol(MAT) # samples
  
  for(r in (1:R)) {
    
    # avoid empty clusters
    clustering <- sample(K,N,replace=TRUE)
    while(length(unique(clustering))!=K){
      clustering <- sample(K,N,replace=TRUE)
    } 
    
    centroids <- matrix(0,K,M) 
    for(j in (1:K)){
      ind.clust <- which(clustering==j)
      if(length(ind.clust)==1){centroids[j,]=MAT[ind.clust,]} # only 1 element in the cluster
      else{
      centroids[j,]=apply(MAT[ind.clust,],2,mean)}
      # centroids contain on row j coordinates of j-th centroid
    }
    
    min.cost.fun <- NA
    condizione <- TRUE
    while(condizione){
      
      counter <- 0  # to exit
      for(i in (1:N)){

        dist.min <- measure(MAT[i,],centroids[1,],measure_type)
        for(h in (1:K)){
          # calculate the distance of the i-th gene from the h-th cluster
          tmp <- measure(MAT[i,],centroids[h,],measure_type)
          if(tmp<=dist.min){
            dist.min <- tmp
            ind.dist.min <- h
          }
        }
        if(clustering[i]!=ind.dist.min){
          counter <- counter+1 # in this iteration +1 shift
          old <- clustering[i]
          clustering[i] <- ind.dist.min # assigning to the new cluster
          
          ind.clust <- which(clustering==ind.dist.min) # update the arriving cluster
          if(length(ind.clust)==1){centroids[ind.dist.min,]=MAT[i,]}
          else{
          centroids[ind.dist.min,]=apply(MAT[ind.clust,],2,mean)}
          
          ind.clust <- which(clustering==old) # update the starting cluster
          if(length(ind.clust)==1){centroids[old,]=MAT[which(clustering==old),]}
          else if(length(ind.clust)>1){centroids[old,]=apply(MAT[ind.clust,],2,mean)}
          
        
        }
        
      }
      
      # I've updated cluster passing the N genes
      if(counter<=0.01*N){condizione <- FALSE} # less than 1% of genes move; quit while
      
    }
    # exiting while, final config of clusters, calculate cost function
    cost.fun <- 0
    for(q in (1:N)){
      cost.fun <- cost.fun + measure(MAT[q,], centroids[clustering[q],], measure_type)
    }
    if(is.na(min.cost.fun)){
      min.cost.fun <- cost.fun
      best.clustering <- clustering
      best.centroids <- centroids
    }
    else{
      if(cost.fun<=min.cost.fun){
        min.cost.fun <- cost.fun
        best.clustering <- clustering
        best.centroids <- centroids
      }
    }
    
  }
  
  output <- list(clustering=best.clustering,centroids=best.centroids)
  return(output)
}