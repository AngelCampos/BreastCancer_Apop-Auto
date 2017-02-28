################################################################################
# Data preparation for ARACNE METABRIC
## Author: Miguel Angel Garcia Campos github: https://github.com/AngelCampos
################################################################################
# Load data ####
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/METABRIC/0_ExpressionData/expData_METABRIC.Rdata")
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/METABRIC/0_ExpressionData/normals_METABRIC.Rdata")
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/METABRIC/5_Differential_Expression/Subgroup_filtered.Rdata")
gs <- as.matrix(read.delim("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/METABRIC/0_Genesets/GOBP_Genesets_Apoptosis_Autophagy.txt", header = F))
colnames(normals) <- colnames(exp.data)

# Selecting data from Apoptosis and Autophagy genes ####
genes <- unique(as.vector(gs[,-2:-1]))
genes <- sort(intersect(rownames(exp.data), genes))
exp.data <- exp.data[genes,names(subg)]
normals <- normals[,names(subg)]
sampleIds <- colnames(exp.data)
exp.data <- rbind(sampleIds, exp.data)

# Subgroup specific expression matrices ####
# Subgroups Labels: 1=Healthy, 2=LumA, 3=LumB, 4= Her2, 5= Basal
subTable <- table(subg)
names(subTable) <- c("Healthy","LumA","LumB","Her2","Basal")
healthy <- exp.data[,subg == 1]
lumA <- exp.data[,subg == 2]
lumB <- exp.data[,subg == 3]
her2 <- exp.data[,subg == 4]
basal <- exp.data[,subg == 5]

# Writing INPUT files for ARACNE ####
write(genes, file = "probesets.txt")
write.table(healthy, file = "healthy_METABRIC.txt", sep = "\t", quote=F,col.names=F)
write.table(lumA, file = "lumA_METABRIC.txt", sep = "\t", quote = F, col.names = F)
write.table(lumB, file = "lumB_METABRIC.txt", sep = "\t", quote = F, col.names = F)
write.table(her2, file = "her2_METABRIC.txt", sep = "\t", quote = F, col.names = F)
write.table(basal, file = "basal_METABRIC.txt", sep = "\t", quote = F, col.names =F)
write.table(subTable, file = "SubtypesN.txt", sep = "\t", col.names = F, row.names = F, 
            quote = F)