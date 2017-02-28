################################################################################
# Data preparation for ARACNE TCGA
## Author: Miguel Angel Garcia Campos github: https://github.com/AngelCampos
################################################################################
# Load data ####
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/0_ExpressionData/expMatrix_TCGA_onlyFemales.Rdata")
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/0_ExpressionData/normals_TCGA_onlyFemales.Rdata")
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/5_Differential_Expression/Subgroup_filtered.Rdata")
gs <- as.matrix(read.delim("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/0_Genesets/GOBP_Genesets_Apoptosis_Autophagy.txt", header = F))

# Selecting data from Apoptosis and Autophagy genes ####
genes <- unique(as.vector(gs[,-2:-1]))
genes <- sort(intersect(rownames(exp.matrix_TCGA), genes))
exp.matrix_TCGA <- exp.matrix_TCGA[genes,names(subg)]
normals_TCGA <- normals_TCGA[names(subg)]
sampleIds <- colnames(exp.matrix_TCGA)
exp.matrix_TCGA <- rbind(sampleIds, exp.matrix_TCGA)

# Subgroup specific expression matrices ####
# Subgroups Labels: 1=LumA, 2=Healthy, 3=Basal, 4= LumB, 5= Her2
subTable <- table(subg)
names(subTable) <- c("LumA","Healthy","Basal","LumB","Her2")
lumA <- exp.matrix_TCGA[,subg == 1]
healthy <- exp.matrix_TCGA[,subg == 2]
basal <- exp.matrix_TCGA[,subg == 3]
lumB <- exp.matrix_TCGA[,subg == 4]
her2 <- exp.matrix_TCGA[,subg == 5]

# Writing INPUT files for ARACNE ####
write(genes, file = "probesets.txt")
write.table(healthy, file = "healthy_TCGA.txt", sep = "\t", quote=F,col.names=F)
write.table(lumA, file = "lumA_TCGA.txt", sep = "\t", quote = F, col.names = F)
write.table(lumB, file = "lumB_TCGA.txt", sep = "\t", quote = F, col.names = F)
write.table(her2, file = "her2_TCGA.txt", sep = "\t", quote = F, col.names = F)
write.table(basal, file = "basal_TCGA.txt", sep = "\t", quote = F, col.names =F)
write.table(subTable, file = "SubtypesN.txt", sep = "\t", col.names = F, row.names = F, 
            quote = F)
