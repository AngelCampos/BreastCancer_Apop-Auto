################################################################################
# 
## Author: Miguel Angel Garcia Campos github: https://github.com/AngelCampos
################################################################################

# Load data
PCmetabric_5 <- read.delim("Pathways_Correlation_METABRIC_0.05.txt", header = F)
PCmetabric_1 <- read.delim("Pathways_Correlation_METABRIC_0.01.txt", header = F)
PCtcga_5 <- read.delim("Pathways_Correlation_TCGA_0.05.txt", header = F)
PCtcga_1 <- read.delim("Pathways_Correlation_TCGA_0.01.txt", header = F)

# Compare coincidences TCGA-METABRIC #################################
union_1 <- rbind(PCmetabric_1, PCtcga_1)
length(which(duplicated(union_1[,-3])))
# coincidences = 1
union_5 <- rbind(PCmetabric_5, PCtcga_5)
length(which(duplicated(union_5[,-3]))) 
# coincidences = 9

# Hypergeometric test to assertain posibility of obtaining coinciding results 
# by chance
# phyper() arguments:
# q: number of succeses in sample: Number of coincidences = 1 and 9
# m: number of succeses in population: Selected correlations = 12 and 59
# n: number of failures in population: Possible correlations = 1168-12 & 1168-59
# k: sample size: Selected correlations = 12 and 59

# For 1% of best correlations
a <- phyper(1, 12, 1168-12, 12, lower.tail = F)
# For 5% of best correlations
b <- phyper(9, 59, 1168-59, 59, lower.tail = F)

# Concensus results
c1 <- rbind(union_1[which(duplicated(union_1[,-3], fromLast = T)),],
      union_1[which(duplicated(union_1[,-3], fromLast = F)),])
c5 <- rbind(union_5[which(duplicated(union_5[,-3], fromLast = T)),],
            union_5[which(duplicated(union_5[,-3], fromLast = F)),])
meanc1 <- NULL
for(i in seq(1, nrow(c1), 2)){
  meanc1 <- c(meanc1, mean(c(c1[i,3], c1[i+1, 3])))
}
meanc5 <- NULL
for(i in seq(1, nrow(c5), 2)){
  meanc5 <- c(meanc5, mean(c(c5[i,3], c5[i+1, 3])))
}

concensus_1 <- cbind(union_1[which(duplicated(union_1[,-3])),-3], meanc1)
concensus_5 <- cbind(union_5[which(duplicated(union_5[,-3])),-3], meanc5)

# Write results
write.table(concensus_1, file = "ConsensusPC_0.01.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(concensus_5, file = "ConsensusPC_0.05.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write(a, file = "pValue_0.01.txt")
write(b, file = "pValue_0.05.txt")
