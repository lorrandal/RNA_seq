TSdiff=function(matr, dataframe, phi, pFP, test, correction, P0)
  #' @param matr normalized and filtered counts matrix
  #' @param dataframe dataframe with information on the names of the samples
  #' @param phi a vector receiving the time-specific dispersion parameters
  #' @param pFP the highest tolerated pFP probability of committing a Type I error
  #' @param test the selection method (string): "t_test", "Wald"
  #' @param correction the correction method (string): 'NO_correction', 'Bonferroni', 'FDR'
  #' @param P0 the percentage P0 of genes not differentially expressed
  #' @return returns a list of 2 items: Element 1: a vector with EntrezIDs of genes
  #'  differentially expressed in at least one temporal instant, Element 2: a 12x4 matrix that
  #'  shows, for each line, the TP, FP, TN and FN expected at time t_j (j=1,...,12) 
  #'  in correspondence to the threshold of required significance.
{
  source('NBdiff.test.R')
  expressed <- c()
  matr.out <- matrix(,12,4)
  G <- nrow(matr)
  G0 <- P0*G
  for (i in (1:12))
  {
    # create X
    ind <- which(dataframe$Time==i & dataframe$Condition=="CTRL")
    X <- matr[,ind]
    
    # create Y
    ind <- which(dataframe$Time==i & dataframe$Condition=="INS")
    Y <- matr[,ind]
    
    # perform t_test/wald
    p.values <- NBdiff.test(X, Y, phi[i+1], test)
    
    if (correction=="NO_correction")
    {
      alpha <- pFP
    }
    if (correction=="Bonferroni")
    {
      alpha <- pFP/G
    }
    if (correction=="FDR")
    {
      grid_support <- sort(unique(p.values))
      delta <- grid_support[2] - grid_support[1]
      for (hh in (3:length(grid_support)))
      {
        tmp <- grid_support[hh]-grid_support[hh-1]
        if (tmp < delta)
        {delta <- tmp}
      }
      grid_support <- grid_support+delta/2 # update grid_support
      
      # parameters for FDR calculation
      FDR <- 0
      j <- 1
      while(FDR<=pFP & j<=length(grid_support))
      {
        alpha <- grid_support[j]
        num.sel <- length(which(p.values < alpha))
        FDR <- (G0 * alpha) / num.sel
        j <- j + 1
      }
    }
    ind.expressed <- which(p.values<alpha)
    expressed[ind.expressed]=TRUE 
    
    # create matrix TP FP TN FN
    selected <- length(ind.expressed) # at time istant i
    FP <- G0 * alpha
    TN <- G0 * (1 - alpha)
    TP <- selected - FP 
    FN <- G - selected - TN
    values <- c(TP, FP, TN, FN)
    matr.out[i,] <- values
  }
  # create ENTREZ_ID vector
  ind.expressed <- which(expressed==TRUE)
  entrez.expressed <- rownames(matr)[ind.expressed]
  
  # rename rows and columns of dataframe
  colnames(matr.out) <- c("TP", "FP", "TN", "FN")
  row_names <- c()
  for (ii in (1:12)) 
  {row_names[ii] <- paste("Time", ii,sep="")
  }
  rownames(matr.out) <- row_names
  
  # create the list to return
  list.output <- list(entrez.expressed, matr.out)
  return(list.output)
}
