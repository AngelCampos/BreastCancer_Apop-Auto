################################################################################
# Preparing TCGA dataset 
## Author: Miguel Angel Garcia Campos github: https://github.com/AngelCampos
################################################################################

# Loading dataset and clinical info ############################################
sClinical <- read.delim("data_bcr_clinical_data_sample.txt", header = T,skip= 5)
pClinical <- read.delim("data_bcr_clinical_data_patient.txt", header= T,skip= 5)
data <- as.matrix(read.delim("expMatrix_TCGA_cBioPortal.txt", row.names = 1, 
                             header = T))
normals <- read.delim("normals_TCGA_cBioPortal.txt", header = T)
rownames(normals) <- "NORMAL.LOGICAL"

# Removing male data ###########################################################
#Male sample ids
malePatients <- as.character(pClinical$PATIENT_ID[pClinical$GENDER =="MALE"])
maleSamples <- as.character(sClinical$SAMPLE_ID[sClinical$PATIENT_ID %in% malePatients])
maleSamples <- gsub(pattern = "-", ".", x = maleSamples)

#Female sample ids
femalePatients <- as.character(pClinical$PATIENT_ID[pClinical$GENDER=="FEMALE"])
femaleSamples <- as.character(sClinical$SAMPLE_ID[sClinical$PATIENT_ID %in% femalePatients])
femaleSamples <- gsub(pattern = "-", ".", femaleSamples)

# Without male data and normals
onlyFemData <- data[,!colnames(data) %in% maleSamples]
onlyFemNormals <- (normals[,!colnames(normals) %in% maleSamples])
normals_TCGA <- as.logical(onlyFemNormals); names(normals_TCGA) <- colnames(onlyFemNormals)
exp.matrix_TCGA <- log2(onlyFemData+1)

# Save results #################################################################
save(exp.matrix_TCGA, file = "expMatrix_TCGA_onlyFemales.Rdata")
save(normals_TCGA, file = "normals_TCGA_onlyFemales.Rdata")