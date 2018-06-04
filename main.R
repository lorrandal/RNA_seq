######################################################################
## Main Script
######################################################################

rm(list=ls())
load("matcount_1095199.RData")
load("Disp_1095199.RData")

source('TMMnorm.R')
matr.norm <- TMMnorm(matcount, 1, 0.1)

# Eliminate genes that always have a count number <=1 
# on all 78 normalized samples

ind.rm <- c()
k <- 1
for (i in (1:nrow(matr.norm)))
{
  j <- 1
  while (j<=78)
  {
    if(matr.norm[i, j]<=1)
    {
      j <- j+1
      if(j==79)
      { ind.rm[k] <- i
      k <- k+1}
    }
    else j=79
      
  }
}

matr.filt <- matr.norm[-ind.rm,]

# create a 3-column data frame associating the column
# names of the normalized matrix to the times
# and conditions to be tested

samples <- colnames(matr.filt)
ind.C <- charToRaw("C")
condition <- c()
time_istant <- c()

for (i in (1:length(samples)))
{
  tmp <- charToRaw(samples[i])
  if (tmp[1]==ind.C)
    {
    condition[i] <- "CTRL"
    if (length(tmp)==13) 
        {time_istant[i] <- as.integer(rawToChar(tmp[7]))}
    else {time_istant[i] <- as.integer(rawToChar(tmp[7:8]))}
    }
  else
  {
    condition[i] <- "INS"
    if (length(tmp)==12)
    {time_istant[i] <- as.integer(rawToChar(tmp[6]))}
    else {time_istant[i] <- as.integer(rawToChar(tmp[6:7]))}
  }
  
}
time_istant <- time_istant-1

d.frame <- data.frame(Sample=samples, Time=time_istant, Condition=condition)


# load and execute the function to select the genes differentially expressed
# over time using the parameters that allow you to obtain the best solution

source('TSdiff.R')
list_TSdiff <- TSdiff(matr.filt, d.frame, disp, 0.05/12, 'Wald', 'Bonferroni', 0.99) 
#### if input != t_test or Wald, list_TSdiff[[1]] is empty; if selected = zero (t_test+Bonferroni)
#### quit or kmeans enters in loop

if(length(list_TSdiff[[1]])<=1){
  print(paste0("Insufficient number of selected genes"))
}

# load files "GOannotations.txt" and "geneinfo.txt"
GOannotations <- read.delim('GOannotations.txt', header=TRUE, sep="\t")
geneinfo <- read.delim('geneinfo.txt', header=TRUE ,sep="\t")

# create a data frame of 2 columns associating the name of the genes 
# to a flag indicating if the gene is (=1) or not (=0) selected
row.entrez <- rownames(matr.filt)
row.flag <- row.entrez%in%list_TSdiff[[1]]
row.flag[which(row.flag==TRUE)] <- 1
row.flag[which(row.flag==FALSE)] <- 0
flag.frame <- data.frame(EntrezID=row.entrez, Selcted=row.flag)

# load and execute the function for enrichment analysis
source('EnrichGO.R')
Enrich.frame <- EnrichGO(flag.frame, GOannotations, geneinfo)

# extract differentially expressed genes from normalised data,
# average replicates, log-transform (with offsets) and calculate treatment-control
ind.diff <- which(flag.frame$Selcted==1)
matr.diff <- matr.filt[ind.diff,]

MAT <- matrix(0,length(ind.diff),13)
rownames(MAT) <- rownames(matr.diff)
for (i in (1:13)){ # using d.frame, I already have the disjointed names of the variables 
  ind.CTRL <- which(d.frame$Time==(i-1) & d.frame$Condition=='CTRL')
  ind.INS <- which(d.frame$Time==(i-1) & d.frame$Condition=='INS')
  # calculate average of replicates
  if(is.null(nrow(matr.diff))){
    mean_CTRL <- mean(matr.diff[ind.CTRL])
    mean_INS <- mean(matr.diff[ind.INS])
    } else{
  mean_CTRL <- apply(matr.diff[,ind.CTRL],1,mean)
  mean_INS <- apply(matr.diff[,ind.INS],1,mean)
  MAT[,i] <- log(mean_INS+1)-log(mean_CTRL+1)}
}

# load and run the k-means clustering function, 
# selecting an adequate number of clusters and algorithm
# restarts and an adequate distance measurement
source('K_means.R')
K <- 8
R <- 300
SEED <- 11

### 2 methods: "Euclidean", "Correlation"
list_Kmeans <- K_means(MAT, K, 'Euclidean', R, SEED)

# plot all the clusters in a single figure, arranged in two lines,
# with the time profiles of the individual elements in grey and the centroid in red

# create layout, plots in 2 rows
if(K%%2==0) {num.colonne=K/2} else
{num.colonne=(K%/%2+1)}
  
layout(matrix(c(1:(2*num.colonne)),2,num.colonne,byrow=TRUE))
for (i in (1:K)) {
  plot(list_Kmeans$centroids[i,],type="l",col="red",main=paste("Cluster",i),xlab="time_istant",ylab="Profile")
  num.el <- length(which(list_Kmeans$clustering==i))
  if (num.el!=0) {
    elementi <- which(list_Kmeans$clustering==i)
    for (j in (1:num.el)) 
      lines(MAT[elementi[j],],col="grey")
  }
}

