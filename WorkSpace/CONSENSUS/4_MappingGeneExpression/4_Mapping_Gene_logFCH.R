################################################################################
# 
## Author: Miguel Angel Garcia Campos github: https://github.com/AngelCampos
################################################################################
# Load/install CRAN packages
ipak <- function(pkg){ # https://gist.github.com/stevenworthington/3178163
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("pathview", "devtools", "gplots")
ipak(packages)
# Installing GithHub Package Annotables (uncomment next line)
# devtools::install_github("stephenturner/annotables")
library(annotables) #GitHub Package

# Loading gene data
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/5_Differential_Expression/HvsLumB.Rdata")
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/5_Differential_Expression/HvsBasal.Rdata")
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/5_Differential_Expression/HvsHer2.Rdata")
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/5_Differential_Expression/HvsLumA.Rdata")
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/CONSENSUS/2_Genes_in_PDSs/Genes_overlap_0.05.Rdata")

# Get logFCH for selected genes
genes <- unique(unlist(genes_5))
genes <- sort(table(unlist(genes_5)), decreasing = T)
genes <- names(genes[genes >=2])
LumA <- HvsLumA[genes,"logFC"]
LumB <- HvsLumB[genes,"logFC"]
Her2 <- HvsHer2[genes,"logFC"]
Basal <- HvsBasal[genes,"logFC"]
LFCtable <- t(rbind(LumA, LumB, Her2, Basal))
rownames(LFCtable) <- genes

# Heatmap LFC
png(filename = "logFC_CrossTalkGenes.png")
heatmap.2(LFCtable, col=greenred(10), trace = "none",density.info = "none",
          keysize = 1.2)
dev.off()

# Annotating Gene symbols -> ENTREZ IDS
hsaAnnot <- grch38
ENTREZ <- NULL
for(i in 1:length(genes)){
  x <- as.character(hsaAnnot[which(genes[i] == hsaAnnot[,"symbol"])[1],"entrez"])
  ENTREZ[i] <- x
}
rownames(LFCtable) <- ENTREZ

# Mapping in KEGG diagrams
# Apoptosis KEGGid 04210
# Autophagy KEGGid 04140
for(i in 1:ncol(LFCtable)){
  pv.out <- pathview(gene.data = LFCtable[, i], pathway.id = "04210", 
                     species = "hsa", out.suffix = colnames(LFCtable)[i])
  pv.out <- pathview(gene.data = LFCtable[, i], pathway.id = "04140", 
                     species = "hsa", out.suffix = colnames(LFCtable)[i])
}