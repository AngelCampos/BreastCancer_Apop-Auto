################################################################################
# Correlation PDS vs gene
## Author: Miguel Angel Garcia Campos github: https://github.com/AngelCampos
################################################################################

# Load data
load("fPDS.Rdata")
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/METABRIC/0_ExpressionData/expData_METABRIC.Rdata")
gs <- as.matrix(read.delim("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/0_Genesets/GOBP_Genesets_Apoptosis_Autophagy.txt", header = F, row.names = 1)[,-1])

# Prepare objects

genes <- intersect(rownames(exp.data), na.omit(unique(as.vector(gs))))
exp.data <- exp.data[genes,]

# Calculate gene-PDS correlations
PDScor <- list()
for(i in rownames(fPDS)){
  gset <- na.omit(unique(as.vector(gs[i,]))); gset <- gset[gset != ""]
  C <- NULL
  for (j in gset) {
    if(j %in% rownames(exp.data)){
      c <- abs(cor(fPDS[i,], exp.data[j,]))
      C <- c(C,c)
    }
  }
  names(C) <- gset[gset %in% rownames(exp.data)]
  C <- sort(C, decreasing = T)
  PDScor[[i]] <- C
}
save(PDScor, file="PDScorrelation.Rdata")
rm(i,j,c,C) # Trash

# Median  relevance in deregulation analysis ####
M <- NULL
for(i in genes){
  m <- median(unlist(PDScor)[grep(pattern = paste("^", i,"$", sep = ""), substring(names(unlist(PDScor)),12))])
  M[i] <- m
}
M <- sort(M, decreasing = T)
save(M, file = "GeneRelevanceInPDS.Rdata")

# Finding a particular gene correlation values ####
gene <- "ATG9B"
unlist(PDScor)[grep(pattern = paste("^",gene,"$", sep = ""), substring(names(unlist(PDScor)),12))]

