TMMnorm=function(matr, index, param)
  #' @param matr matrix to be normalised
  #' @param index index of the reference sample chosen to standardise
  #' @param param trim parameter to use for the average
  #' @return The same matrix with the normalized data using scaling normalization with TMM
  #' 
  #'  Print a .svg file: each plot has 2 adjacent boxes showing the MvA plots (in log2-scale)
  #'  of the sample before and after normalization 
{  
  col_index <- matr[, index] 
  matr_norm <- matr
  SF <- 0;

# Calculate the normalized matrix using trimmed mean
  for (i in (2:length(matr[1,])))
       {
       SF[i-1] <- mean(log2(col_index)-log2(matr[, i]), trim=param, na.rm = TRUE)
       matr_norm[,i] <- matr_norm[,i]*(2^SF[i-1])
  }
  
  
## Save some specific values for MvA plot purposes 
# TODO: collapse the code belove in a single loop   

 ## 2
  m2.pre <- (log2(col_index)-log2(matcount[,2]))
  a2.pre <- (log2(col_index)+log2(matcount[,2]))/2
  
  m2.norm <- (log2(col_index)-log2(matr_norm[,2]))
  a2.norm <- (log2(col_index)+log2(matr_norm[,2]))/2
  TMM.norm2 <- mean(m2.norm, trim=param, na.rm=TRUE)
 ## 3
  m3.pre <- (log2(col_index)-log2(matcount[,3]))
  a3.pre <- (log2(col_index)+log2(matcount[,3]))/2
  
  m3.norm <- (log2(col_index)-log2(matr_norm[,3]))
  a3.norm <- (log2(col_index)+log2(matr_norm[,3]))/2
  TMM.norm3 <- mean(m3.norm,trim=param,na.rm=TRUE)
 ## 4
  m4.pre <- (log2(col_index)-log2(matcount[,4]))
  a4.pre <- (log2(col_index)+log2(matcount[,4]))/2
  
  m4.norm <- (log2(col_index)-log2(matr_norm[,4]))
  a4.norm <- (log2(col_index)+log2(matr_norm[,4]))/2
  TMM.norm4 <- mean(m4.norm,trim=param,na.rm=TRUE)
 ## 5
  m5.pre <- (log2(col_index)-log2(matcount[,5]))
  a5.pre <- (log2(col_index)+log2(matcount[,5]))/2
  
  m5.norm <- (log2(col_index)-log2(matr_norm[,5]))
  a5.norm <- (log2(col_index)+log2(matr_norm[,5]))/2
  TMM.norm5 <- mean(m5.norm,trim=param,na.rm=TRUE)
 ## 6  
  m6.pre <- (log2(col_index)-log2(matcount[,6]))
  a6.pre <- (log2(col_index)+log2(matcount[,6]))/2
  
  m6.norm <- (log2(col_index)-log2(matr_norm[,6]))
  a6.norm <- (log2(col_index)+log2(matr_norm[,6]))/2
  TMM.norm6 <- mean(m6.norm,trim=param,na.rm=TRUE)

  labels_pdf <- colnames(matr_norm)[2:6]
  
  ## print svg
  
  svg("plots_TMMnorm.svg")
  par(mfrow=c(3,2))
  
  plot(a2.pre, m2.pre, main=paste0(labels_pdf[1],'_Pre-Norm'), xlab = 'A', ylab = 'M')
  lines(rep(0,18), col=2)
  lines(rep(SF[1], 18), col=3)
  legend("topright", "TMM", col=3, lty=1, bty="n")
  plot(a2.norm, m2.norm, main=paste0(labels_pdf[1],'_Norm'), xlab = 'A', ylab = 'M')
  lines(rep(TMM.norm2, 18), col=2)
  legend("topright", "TMM", col=2, lty=1, bty="n")
  
  plot(a3.pre,m3.pre, main=paste0(labels_pdf[2],'_Pre-Norm'), xlab = 'A', ylab = 'M')
  lines(rep(0,18),col=2)
  lines(rep(SF[2],18),col=3)
  legend("topright","TMM",col=3,lty=1,bty="n")
  plot(a3.norm,m3.norm, main=paste0(labels_pdf[2],'_Norm'), xlab = 'A', ylab = 'M')
  lines(rep(TMM.norm3,18),col=2)
  legend("topright","TMM",col=2,lty=1,bty="n")
  
  plot(a4.pre,m4.pre, main=paste0(labels_pdf[3],'_Pre-Norm'), xlab = 'A', ylab = 'M')
  lines(rep(0,18),col=2)
  lines(rep(SF[3],18),col=3)
  legend("topright","TMM",col=3,lty=1,bty="n")
  plot(a4.norm,m4.norm, main=paste0(labels_pdf[3],'_Norm'), xlab = 'A', ylab = 'M')
  lines(rep(TMM.norm4,18),col=2)
  legend("topright","TMM",col=2,lty=1,bty="n")
  
  # plot(a5.pre,m5.pre, main=paste0(labels_pdf[4],'_Pre-Norm'), xlab = 'A', ylab = 'M')
  # lines(rep(0,18),col=2)
  # lines(rep(SF[4],18),col=3)
  # legend("topright","TMM",col=3,lty=1,bty="n")
  # plot(a5.norm,m5.norm, main=paste0(labels_pdf[4],'_Norm'), xlab = 'A', ylab = 'M')
  # lines(rep(TMM.norm5,18),col=2)
  # legend("topright","TMM",col=2,lty=1,bty="n")
  # 
  # plot(a6.pre,m6.pre, main=paste0(labels_pdf[5],'_Pre-Norm'), xlab = 'A', ylab = 'M')
  # lines(rep(0,18),col=2)
  # lines(rep(SF[5],18),col=3)
  # legend("topright","TMM",col=3,lty=1,bty="n")
  # plot(a6.norm,m6.norm, main=paste0(labels_pdf[5],'_Norm'), xlab = 'A', ylab = 'M')
  # lines(rep(TMM.norm6,18),col=2)
  # legend("topright","TMM",col=2,lty=1,bty="n")
  
  dev.off()
  
  
return(matr_norm)
} 