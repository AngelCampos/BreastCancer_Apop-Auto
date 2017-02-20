################################################################################
# 
## Author: Miguel Angel Garcia Campos github: https://github.com/AngelCampos
################################################################################
# Load correlations
wd <- getwd(); setwd(".."); setwd("1_PDS correlation/")
consensus_5 <- read.delim(file = "ConsensusPC_0.05.txt", header = F)
consensus_1 <- read.delim(file = "ConsensusPC_0.01.txt", header = F)

# Load Genesets
setwd(".."); setwd(".."); setwd("TCGA/0_Genesets/")
member <- as.matrix(read.delim("GOids_Apoptosis_Autophagy_Membership2.txt", F, row.names = 1))
genesets <- as.matrix(read.delim("GOBP_Genesets_Apoptosis_Autophagy.txt", F, row.names = 1))
setwd(wd); rm(wd)

# Which genes are shared by pathways?
genes_5 <- list()
pathSize <- list()
for(i in 1:nrow(consensus_5)){
  x <- na.omit(genesets[rownames(member)[member[,1] %in% consensus_5[i,1]],-1])
  x <- x[x != ""]
  y <- na.omit(genesets[rownames(member)[member[,1] %in% consensus_5[i,2]],-1])
  y <- y[y != ""]
  genes_5[[i]] <- intersect(x,y)
  pathSize[[i]] <- c(length(x), length(y))
  rm(x, y)
}
save(genes_5, file = "Genes_overlap_0.05.Rdata")
genes_1 <- list()
for(i in 1:nrow(consensus_1)){
  x <- na.omit(genesets[rownames(member)[member[,1] %in% consensus_1[i,1]],-1])
  x <- x[x != ""]
  y <- na.omit(genesets[rownames(member)[member[,1] %in% consensus_1[i,2]],-1])
  y <- y[y != ""]
  genes_1[[i]] <- intersect(x,y)
  rm(x, y)
}
save(genes_1, file = "Genes_overlap_0.01.Rdata")
freqGenes_5 <- sort(table(unlist(genes_5)), decreasing = T)

# 123 apoptosis
# 28 autophagy
apopIds <- rownames(member)[member[,2] == "Apoptosis"]
apopIds <- intersect(rownames(genesets), apopIds)
apopGenes <- as.vector(genesets[apopIds, -1])
apopGenes <- sort(table(apopGenes[apopGenes != ""]), decreasing = T)

autoIds <- rownames(member)[member[,2] == "Autophagy"]
autoIds <- intersect(rownames(genesets), autoIds)
autoGenes <- as.vector(genesets[autoIds, -1])
autoGenes <- sort(table(autoGenes[autoGenes != ""]), decreasing = T)