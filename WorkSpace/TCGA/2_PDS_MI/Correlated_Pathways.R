################################################################################
# Correlation between pathways in Mutual Information 
## Author: Miguel Angel Garcia Campos github: https://github.com/AngelCampos
################################################################################
# Load packages ################################################################
library(gplots)
library(RColorBrewer)

# Load data ####################################################################
wd <- getwd()
setwd("..")
setwd("1_Deregulation_Analysis/")
load("fPDS.Rdata")
load("ZscorePDS_TCGA.Rdata")
setwd("..")
setwd("0_ExpressionData/")
load("normals_TCGA_onlyFemales.Rdata")
setwd("..")
setwd("0_Genesets/")
member <- as.matrix(read.delim("GOids_Apoptosis_Autophagy_Membership.txt", 
                               header = F, row.names = 1))
setwd(wd)
rm(wd)
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/4_Subgroups/Subgroups.Rdata")
subg <- mytree; rm(mytree)

# Membership ###################################################################
M <- member[rownames(fPDS),2]
auto <- names(M)[M=="Autophagy"]
apop <- names(M)[M=="Apoptosis"]

# Spearman correlation #########################################################
SCmatrix <- matrix(0, length(auto), length(apop))
rownames(SCmatrix) <- rownames(fPDS[auto,])
colnames(SCmatrix) <- rownames(fPDS[apop,])
for(i in auto){
  x1 <- fPDS[i,]
  for(j in apop){
    x2 <- fPDS[j,]
    c <- abs(cor(x1, x2, method = "spearman"))
    SCmatrix[i,j] <- c
  }
}
rm(x1,i)
png(filename = "PDS_Correlations.png")
heatmap.2(SCmatrix, trace = "none", labCol = "", labRow = "", 
          ylab = "Autophagic Pathways", xlab = "Apoptotic Pathways",
          density.info = "none", keysize = 1.2)
dev.off()

# Best correlations ############################################################
BEST99 <- which(SCmatrix >= quantile(SCmatrix, 0.99), arr.ind = T)
BEST95 <- which(SCmatrix >= quantile(SCmatrix, 0.95), arr.ind = T)
WORST <- which(SCmatrix <= quantile(SCmatrix, 0.10), arr.ind = T)

# Selected correlations plot
png(filename = "Selected Correlations.png")
plot(sort(SCmatrix), main= "Selected correlations (0.99 Quantile)")
abline(h = quantile(SCmatrix, 0.99), col = "red")
abline(h = quantile(SCmatrix, 0.95), col = "blue")
dev.off()
#Pathways best correlations ####################################################
#BEST 1%
PW <- NULL
for(i in 1:nrow(BEST99)){
  x <- c(rownames(SCmatrix)[BEST99[i,1]], colnames(SCmatrix)[BEST99[i,2]],
         SCmatrix[BEST99[i,1],BEST99[i,2]])
  PW <- rbind(PW, x); rm(x)
}
PW <- PW[order(PW[,3], decreasing = T),]
# With names 
PWnames <- NULL
for (i in 1:nrow(PW)){
  x <- c(member[PW[i,1],1], member[PW[i,2],1], PW[i,3])
  PWnames <- rbind(PWnames, x); rm(x)
}
PWnames <- gsub(pattern = "_",replacement = " ", x = PWnames)
save(PWnames, file = "PathsCorrelation_0.01.Rdata")

#BEST 5%
PW <- NULL
for(i in 1:nrow(BEST95)){
  x <- c(rownames(SCmatrix)[BEST95[i,1]], colnames(SCmatrix)[BEST95[i,2]],
         SCmatrix[BEST95[i,1],BEST95[i,2]])
  PW <- rbind(PW, x); rm(x)
}
PW <- PW[order(PW[,3], decreasing = T),]
# With names 
PWnames <- NULL
for (i in 1:nrow(PW)){
  x <- c(member[PW[i,1],1], member[PW[i,2],1], PW[i,3])
  PWnames <- rbind(PWnames, x); rm(x)
}
PWnames <- gsub(pattern = "_",replacement = " ", x = PWnames)
save(PWnames, file = "PathsCorrelation_0.05.Rdata")

# Plotting PDS of best and worst correlations ##################################
# Subtype colors
scolors <- as.factor(subg)
scolors <- brewer.pal(5, "Set1")[subg]

# 4 best correlations
png(filename = "BestCorrelations.png", width = 750, height = 750)
par(mfrow=c(2,2))
par(pch= 19)
plot(fPDS[PW[1,1],], fPDS[PW[1,2],], col = scolors)
plot(fPDS[PW[2,1],], fPDS[PW[2,2],], col = scolors)
plot(fPDS[PW[3,1],], fPDS[PW[3,2],], col = scolors)
plot(fPDS[PW[4,1],], fPDS[PW[4,2],], col = scolors)
dev.off()


# Worst correlations
#Worst
badPW <- NULL
for(i in 1:nrow(WORST)){
  x <- c(rownames(SCmatrix)[WORST[i,1]], colnames(SCmatrix)[WORST[i,2]],
         SCmatrix[WORST[i,1],WORST[i,2]])
  badPW <- rbind(badPW, x); rm(x)
}
badPW <- badPW[order(badPW[,3], decreasing = F),]

png(filename = "WorstCorrelations.png", width = 750, height = 750)
par(mfrow=c(2,2))
par(pch= 19)
plot(fPDS[badPW[1,1],], fPDS[badPW[1,2],], col = scolors)
plot(fPDS[badPW[2,1],], fPDS[badPW[2,2],], col = scolors)
plot(fPDS[badPW[3,1],], fPDS[badPW[3,2],], col = scolors)
plot(fPDS[badPW[4,1],], fPDS[badPW[4,2],], col = scolors)
dev.off()

# Plotting correlated PDSs
sPDS <- fPDS[unique(as.character(PW[,1:2])),] # PDS matrix of correlated PDSs

png(filename = "AllPDS.png", width = 1200, height = 300)
par(pch=16)
plot(sPDS[1,], rep(1, ncol(sPDS)), col = scolors, 
     xlim = c(0,1), ylim = c(0,nrow(sPDS)), main= "Correlated PDSs Samples",
     xlab = "PDS", ylab = "Pathways")
for(i in 2:nrow(sPDS)){
  points(sPDS[i,], rep(i, ncol(sPDS)), col = scolors, xlim = c(0,1), ylim = c(0,ncol(sPDS)))
}
dev.off()