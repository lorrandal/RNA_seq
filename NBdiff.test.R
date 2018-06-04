NBdiff.test=function(X, Y, phi, type.test)
  #' @param X matrix of counts with the samples of the "CTRL" condition at time t_j
  #' @param Y matrix of counts with the samples of the "INS" condition at time t_j
  #' @param phi dispersion parameter
  #' @param type.test the selection method (string): "t_test", "Wald"
  #' @return a vector with the p-values associated with EntrezID of the genes, ordered by rownames(X)=rownames(Y).
{
 ind <- nrow(X)
 p.values <- c()
 if(type.test=="t_test") ## t_test WELCH
 {
  for (i in (1:ind))
  {
    a <- X[i,]
    b <- Y[i,]
    
    if (var(a)==0 | var(b)==0 | var(a)==var(b)) {
      p.values[i] <- NA } 
    else      
      p.values[i] <- t.test(a,b)$p.value
  }   
   
 }
 if (type.test=="Wald") ## WALD
 {for (i in (1:ind))
 { a <- X[i,]
 b <- Y[i,]
 nA <- length(a)
 nB <- length(b)
 sumA <- sum(a)
 sumB <- sum(b)

 if (var(a)==0 | var(b)==0 | var(a)==var(b)) 
 {p.values[i] <- NA}
 else
 { y <- (sumA-(nA/nB)*sumB)/sqrt(sumA+phi*sumA^2/nA+(nA^2/nB^2)*sumB+phi*sumB^2/nB)
 p.values[i] <- 2*pnorm(-abs(y))}
 }}  
 
 if(type.test!="t_test" & type.test!="Wald"){return(NULL)}
 
 return(p.values)  
}

