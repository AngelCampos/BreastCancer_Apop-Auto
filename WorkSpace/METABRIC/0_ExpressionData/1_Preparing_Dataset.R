################################################################################
# Preparing TCGA dataset 
## Author: Miguel Angel Garcia Campos github: https://github.com/AngelCampos
################################################################################
# Loading dataset and clinical info ############################################
n.expdata <- read.delim("normals_ExpressionMatrix.txt", header = T,row.names = 1)
sClinical <- as.matrix(read.delim("table_S2_revised.txt", header = T))
eData.All <- as.matrix(read.delim("Exp.matrix_collapsed.txt", header = T, 
                                  sep = "\t", row.names = 1))

# Logical vector for normal samples (1 = Normal, 0 = Tumor)
normals <- matrix(rep(1, length(colnames(n.expdata))), nrow = 1) 
colnames(normals) <- colnames(n.expdata)
file <- "discovery_ExpressionMatrix.txt" # Count tumoral samples (discovery)
con <- file(description=file, open="r"); line = readLines(con, n = 1)
ids <- unlist(strsplit(line, "\t")); 
tumors <- matrix(rep(0, length(ids)), nrow = 1)
colnames(tumors) <- ids; normals <- cbind(normals, tumors)

# Clinical data
officialPAM50 <- sClinical[,"Pam50Subtype"]
names(officialPAM50) <- sClinical[,"METABRIC_ID"]
save(officialPAM50, file = "official_PAM50.Rdata")

# Expression data 
# Write working files #########################################################
setwd(wd)
write.table(normals, "normals_METABRIC.txt", sep ="\t", quote = F, col.names = NA)
write.table(pam50, "pam50_METABRIC_official.txt", sep = "\t", quote = F, col.names = NA)
save(normals, file = "normals_METABRIC.Rdata")
save(pam50, file = "pam50_METABRIC_official.Rdata")
save(eData.All, file = "expData_METABRIC.Rdata")
