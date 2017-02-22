################################################################################
# Differential expression with LIMMA
## Author: Miguel Angel Garcia Campos github: https://github.com/AngelCampos
################################################################################
# Load/Install Packages
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("limma")) {
  install.packages("limma", dependencies = TRUE)
  library(limma)
}
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}

# Load files
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/4_Subgroups/Subgroups.Rdata")
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/0_ExpressionData/normals_TCGA_onlyFemales.Rdata")
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/0_ExpressionData/expMatrix_TCGA_onlyFemales.Rdata")
subg <- mytree; rm(mytree)

# Removing miss-classified samples
# normals_TCGA classified as tumors
missmatchN <- subg[as.logical(normals_TCGA)] != "2"
missmatchT <- subg[!as.logical(normals_TCGA)] == 2
selected <- c(names(missmatchN)[!missmatchN], names(missmatchT)[!missmatchT])
exp.matrix_TCGA <- exp.matrix_TCGA[,selected]
normals_TCGA <- normals_TCGA[selected]
subg <- subg[selected]

# Differential Expression Analysis (LIMMA) & Volcano plots #####################
# Design matrix
design = matrix(rep(0,ncol(exp.matrix_TCGA)*length(unique(subg))), nrow = ncol(exp.matrix_TCGA)) #Design matrix
colnames(design) = c("LumA", "Healthy", "Basal", "LumB"
                     , "Her2")
for (i in 1:length(unique(subg))){
  design[subg == i,i] <- 1 # Genera la matriz de diseño a partir de la lista de subgrupos
}

# Matriz de contrastes
cont.matrix = makeContrasts(Healthy - LumA, Healthy - LumB, 
                            Healthy - Basal, Healthy - Her2,
                            levels = design)
fit = lmFit(exp.matrix_TCGA, design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

#Get the results for every contrast
for (i in 1:4){
  name <- paste("tt", i, sep = "")
  assign(name, topTable(fit2, coef = i, adjust = "fdr", number = nrow(exp.matrix_TCGA)))
}

contrasts = c("Healthy - LumA", "Healthy - LumB", 
              "Healthy - Basal", "Healthy - Her2")

coef = c("Healthy - LumA", "Healthy - LumB", 
         "Healthy - Basal", "Healthy - Her2")

# Volcano plots PDF report
pdf(file="volcano_plots.pdf", height=7, width=10)

# Parametros gráficos
M = 1; B = 5; y0 = -10; y1 = 400; x0 = -5; x1 = 5

for (i in 1:4){
  volcanoplot(fit2, col = "blue", ylim = c(y0,y1), xlim = c(x0,x1),
              coef = i, main = contrasts[i], cex.lab = 1.3)
  par(new = T)
  abline(v = -M, col = "darkgray")
  par(new = T)
  abline(v = M, col="darkgray")
  par(new = T)
  abline(h = B, col = "black")
  par(new = T)
  ind1 = abs(fit2$coef[,i])>M
  ind2 = fit2$lods[,i] >B
  ind3 = (fit2$coef[,i]>M & fit2$lods[,i]>B)
  ind4 = (fit2$coef[,i]< -M & fit2$lods[,i]>B)
  x = as.matrix(fit2$coef[ind1,i])
  y = as.matrix(fit2$lods[ind1,i])
  plot(x, y, col="magenta", ylim = c(y0,y1), xlim = c(x0,x1),
       main = "", pch = "*", xlab = "", ylab = "", cex.lab = 1.2)
  x = as.matrix(fit2$coef[ind2,i])
  y = as.matrix(fit2$lods[ind2,i])
  par(new=T)
  plot(x, y, col = "orange", ylim = c(y0,y1), xlim = c(x0,x1),
       main = "", pch = 19, xlab = "", ylab = "", cex.lab = 1.2)
  x = as.matrix(fit2$coef[ind3,i])
  y = as.matrix(fit2$lods[ind3,i])
  par(new = T)
  plot(x, y, col = "red",  ylim =c (y0,y1), xlim = c(x0,x1),
       main = "", pch = 19, xlab = "", ylab = "",cex.lab = 1.2)
  x = as.matrix(fit2$coef[ind4,i])
  y = as.matrix(fit2$lods[ind4,i])
  par(new = T)
  plot(x, y, col = "green3",  ylim = c(y0,y1), xlim = c(x0,x1),
       main = "", pch = 19, xlab = "", ylab = "", cex.lab = 1.2)
}
dev.off()

## Write results
HvsLumA <- tt1; HvsLumB <- tt2; HvsBasal <- tt3; HvsHer2 <- tt4; 
save(HvsLumA, file = "HvsLumA.Rdata"); save(HvsLumB, file = "HvsLumB.Rdata")
save(HvsHer2, file = "HvsHer2.Rdata"); save(HvsBasal, file = "HvsBasal.Rdata")