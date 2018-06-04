EnrichGO=function(flag.frame, GOannotations, geneinfo)
  #' @param flag.frame the data frame on the names of the selected/unselected genes
  #' @param GOannotations the data frame from "GOannotations.txt"
  #' @param geneinfo the data frame from "geneinfo.txt"
  #' @return returns in output a data frame, ordered according to increasing p-values, which reports for each row:
  #   GOID: GO term identifier
  #   Term: Name of the GO term
  #   p-value: nominal p-value of the Fisher test.
  #   FDR: The corresponding FDR, using P0=1
  #   Annotated: number of genes annotated in the GO term
  #   Selected: Number of selected genes that are annotated in the GO term
  #   EntrezID: EntrezID identifiers of the genes selected in the GO term (separator = ";")
  #   Gene Name: Names of genes selected in the GO term (separator = ";")
  #   Gene Symbol: symbols of the genes selected in the GO term(separator = ';')
  #
  # The function performs enrichment analysis using Fisher's test, considering an FDR
  # correction for multiple tests, on all GO terms reported in the file "GOannotations.txt".

{
  # initializing vectors for output
  GOID <- as.character(GOannotations[,1])
  L <- length(GOID) #per inizializzare gli altri vettori
  term <- as.character(GOannotations[,2])
  p.values <- rep(NA,L)
  FDR <- rep(NA,L)
  annotated <- rep(NA,L)
  selected <- rep(NA,L)
  EntrezID <- rep(NA,L)
  gene.name <- rep(NA,L)
  gene.symbol <- rep(NA,L)
  
  N <- nrow(flag.frame) # total number of genes (filtered) in the samples
  ID.tot.sel <- flag.frame[which(flag.frame[,2]==1),1] #ID selected genes
  tot.sel <- length(ID.tot.sel) # total selected genes
  tot.no.sel <- N-tot.sel # total unselected genes
  
  
  for (i in (1:nrow(GOannotations)))
  {
   ID.annotated <- strsplit(as.character(GOannotations[i,3]),";") # split 3 cols 
   ID.annotated <- unique(ID.annotated[[1]]) # eliminate duplicated annotated
   n.annotated <- length(which(flag.frame[,1] %in% ID.annotated)) # Hanlde ID not present in the samples
   annotated[i] <- n.annotated # a+c fisher
   
   ID.sel.annotated <- ID.annotated[which(ID.annotated%in%ID.tot.sel)] # selected among annotated
   EntrezID[i] <- paste(ID.sel.annotated,collapse = ";") # for output
   a <- length(ID.sel.annotated) # a fisher test, num sel annotated
   selected[i] <- a 
   
   c <- n.annotated - a #fisher
   b <- tot.sel - a #fisher
   d <- tot.no.sel-c #fisher
   
   # fisher test
   mat.fish <- matrix(c(a,c,b,d),2,2)
   res <- fisher.test(mat.fish, alternative = "greater")
   pval <- res$p.value
   p.values[i] <- pval
   
   #gene name & symbol
   ind.name <- which(geneinfo[,1] %in% ID.sel.annotated)
   gene.name[i] <- paste(as.character(geneinfo[ind.name,2]), collapse=";")
   gene.symbol[i] <- paste(as.character(geneinfo[ind.name,3]), collapse=";")
   
   #FDR
   if(a==0){FDR[i]=0} else{
   FDR[i]=(pval*L)/a}
   
  }
  
  
  output <- data.frame(GOID=GOID, Term=term, p_value=p.values, FDR=FDR, Annotated=annotated,
                    Selected=selected, EntrezID=EntrezID, Gene_Name=gene.name, Gene_Symbol=gene.symbol)
  output <- output[order(output[,3]),]
  
  return(output)
}