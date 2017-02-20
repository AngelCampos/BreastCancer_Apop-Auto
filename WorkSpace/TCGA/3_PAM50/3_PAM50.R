################################################################################
# Input variables for the subtype prediction script (PAM50)
# Original author: JS Parker et al., J Clin Onc 27(8):1160-7 (2009)
# Original code: https://genome.unc.edu/pubsup/breastGEO/PAM50.zip
################################################################################
if (!require("HGNChelper")) {
  install.packages("HGNChelper", dependencies = TRUE)
  library(HGNChelper)
}
if (!require("genefu")){
  source("https://bioconductor.org/biocLite.R")
  biocLite("genefu", ask = F)
  library(genefu)
}

# Load data ####################################################################
load("C:/Users/WildFang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/0_ExpressionData/expMatrix_TCGA_onlyFemales.Rdata")
load("C:/Users/WildFang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/0_ExpressionData/normals_TCGA_onlyFemales.Rdata")
normals <- as.logical(onlyFemNormals)
names(normals) <- colnames(onlyFemNormals)
names(normals) <- gsub(replacement = "-", pattern = "\\.", names(normals))
exp.mat <- onlyFemData
colnames(exp.mat) <- names(normals)
rm(onlyFemData, onlyFemNormals)

# GENEFU PAM50 method ##########################################################
# Annotation
checkedGenes <- checkGeneSymbols(rownames(exp.mat), unmapped.as.na = F)[,3]
annot <- data.frame(HUGO.gene.symbol = checkedGenes, 
                    NCBI.gene.symbol = checkedGenes,
                    EntrezGene.ID = checkedGenes,
                    Gene.Symbol = checkedGenes,
                    probe = checkedGenes)
rownames(annot) <- rownames(exp.mat)

# PAM50 subtyping using PAM50
exp.mat <- exp.mat[intersect(rownames(pam50$centroids), rownames(exp.mat)),]
data <- t(exp.mat)
PAM50pred <- molecular.subtyping(sbt.model = "pam50", data = t(exp.mat),
                             annot = annot, do.mapping= F)

genefuPAM50 <- as.character(PAM50pred$subtype)
names(genefuPAM50) <- names(PAM50pred$subtype)
save(genefuPAM50, file = "genefu_PAM50.Rdata") # Save PAM50 results (genefu)

# Tumoral samples wrongly classified
table(genefuPAM50[!normals])
# Normals samples classification
table(genefuPAM50[normals])

# Save session
save.image(file = "CheckPoint_3.Rdata") # Check-Point!
load("CheckPoint_3.Rdata")