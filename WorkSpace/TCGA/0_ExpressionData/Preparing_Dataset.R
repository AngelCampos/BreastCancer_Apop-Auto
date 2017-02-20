################################################################################
# Preparing TCGA dataset 
## Author: Miguel Angel Garcia Campos github: https://github.com/AngelCampos
################################################################################

# Loading dataset and clinical info ############################################
sClinical <- read.delim("data_bcr_clinical_data_sample.txt", header = T,skip= 5)
pClinical <- read.delim("data_bcr_clinical_data_patient.txt", header= T,skip= 5)
data <- read.delim("expMatrix_TCGA_cBioPortal.txt", row.names = 1, header = T)
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
onlyFemNormals <- normals[,!colnames(normals) %in% maleSamples]

#Write results #################################################################
write.table(onlyFemData, "expMatrix_TCGA_onlyFemales.txt",quote = F, sep = "\t",
            col.names = NA )
write.table(onlyFemNormals, "normals_TCGA_onlyFemales.txt", quote = F, 
	        sep = "\t", col.names = NA)
save(onlyFemData, file = "expMatrix_TCGA_onlyFemales.Rdata")
save(onlyFemNormals, file = "normals_TCGA_onlyFemales.Rdata")