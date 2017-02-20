################################################################################
# PAM50 subtyping with genefu and HGNChelper
################################################################################
library(HGNChelper)
library(genefu)

# Prepare matrix for PAM50 #####################################################
wd <- getwd(); setwd(".."); setwd("0_ExpressionData/")
load("expData_METABRIC.Rdata")
load("normals_METABRIC.Rdata")
setwd(wd)
colnames(exp.data) <- gsub(replacement = "-", pattern = "\\.", colnames(exp.data))
colnames(exp.data) <- gsub(replacement = "-", pattern = " \\.", colnames(normals))
exp.mat <- exp.data[,colnames(normals)]
rm(exp.data, wd)

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
exp.mat <- exp.mat[which(checkedGenes %in% rownames(pam50$centroids), arr.ind = T),]
data <- t(exp.mat)
PAM50pred <- molecular.subtyping(sbt.model = "pam50", data = t(exp.mat),
                                 annot = annot, do.mapping= F)
genefuPAM50 <- as.character(PAM50pred$subtype)
names(genefuPAM50) <- names(PAM50pred$subtype)
save(genefuPAM50, file = "genefu_PAM50.Rdata") # Save results

# Compare OFFICIAL PAM50 vs genefu algorithm results
wd <- getwd()
setwd("..")
load("0_ExpressionData/official_PAM50.Rdata")
setwd(wd)
rm(wd)

# Intersecting samples
genefuPAM50_f <- genefuPAM50[intersect(names(officialPAM50), names(genefuPAM50))]

# Proportion of identical calls
(sum(officialPAM50 == genefuPAM50_f) + sum(genefuPAM50[normals] == "Normal")) / length(genefuPAM50)

# Tumoral samples wrongly classified
sum(genefuPAM50_f == "Normal") # 25 "Normal-like" tumors
# Proportion of tumoral subtypes
table(genefuPAM50_f)/ length(genefuPAM50_f)

# Save session
save.image(file = "CheckPoint_3_3.Rdata") # Check-Point!
load("CheckPoint_3_3.Rdata")