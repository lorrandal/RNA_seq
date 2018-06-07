# RNA-seq

In this project, I elaborate an analysis pipeline for RNA-seq data through the creation of appropriate functions called from a single ```main``` script . 

The aim is therefore not to use functions from the most famous NGS data analysis packages but to create custom ones in order to work at a low level on the primitives.

Each step of the pipeline will be represented by a different function.

### **Pipeline scheme:**
 - Data normalization
 - Selection of differentially expressed genes
 - Biological interpretation with functional annotation
 - Clustering Expression Profiles (k-means)
## DATA

The data reproduce an experimental design in which the dynamics of **beta cells** gene expression is monitored under two conditions:

- activation of insulin secretion ("**INS**" samples)
- inhibition of insulin secretion ("**CTRL**" samples)

The data are stored in the file ```matcount_1095199.RData```, matrix 10000x78:

**N** = 10000 genes monitored

**M** = 13 time instants

3 **biological replicates** for each temporal instant

where:

```rownames(matcount) = EntrezID``` of the monitored genes
(http://www.ncbi.nlm.nih.gov/gene)

```colnames(matcount) = CONDITION.TIME.REPLICATES```

CONDITION = {"CTRL","INS"}

TIME = {"T0","T1",...,"T12"}

REPLICATES = {"REPL1","REPL2","REPL3"}

It is possible to use this viewer for a better visualization of the notebook:
https://nbviewer.jupyter.org/github/lorrandal/RNA_seq/blob/master/RNA_seq.ipynb

