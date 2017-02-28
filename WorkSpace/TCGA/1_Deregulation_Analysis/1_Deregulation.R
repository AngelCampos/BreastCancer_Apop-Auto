################################################################################
# Deregulation Analysis Pathifier - TCGA
## Author: Miguel Angel Garcia Campos github: https://github.com/AngelCampos
################################################################################
# Installing and/or loading required packages ##################################
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("pathifier")) {
  source()
  install.packages("pathifier", dependencies = TRUE)
  library(pathifier)
}

# Loading data #################################################################
wd <- getwd(); setwd(".."); setwd(".."); setwd("0_ExpressionData/TCGA/")
exp.mat <- as.matrix(read.delim("expMatrix_TCGA_onlyFemales.txt",row.names = 1))
exp.mat <- log2(exp.mat+1)
normals <- as.matrix(read.delim("normals_TCGA_onlyFemales.txt", header = T,
                                row.names = 1))
setwd(".."); setwd(".."); setwd("0_Genesets/")
gs <- as.matrix(read.delim("GOBP_Genesets_Apoptosis_Autophagy.txt", header = F,
                sep = "\t"))
rownames(gs) <- gs[,1]; setwd(wd)
sampleIds <- colnames(exp.mat)

# Prepare DATA for Pathifier ###################################################
# Calculate MIN_STD
N.exp.mat <- exp.mat[,as.logical(normals)]
rsd <- apply(N.exp.mat, 1, sd)
min_std <- quantile(rsd, 0.25)

# Filter low value genes. At least 10% of samples with values over min_std
min_exp <- quantile(as.vector(exp.mat), 0.1) # min_std = Percentile 10 of data
over <- apply(exp.mat, 1, function(x) x > min_exp) # Values higher than min_exp
G.over <- apply(over, 2, mean)
G.over <- names(G.over)[G.over > 0.1]
exp.mat <- exp.mat[G.over,]

# Filter relevant genes (Apoptosis Autophagy)
genes <- unique(as.vector(gs[,c(-1,-2)]))
genes <- intersect(genes, rownames(exp.mat))
exp.mat <- exp.mat[genes,]

# Set every value bellow min_exp to min_exp
exp.mat[exp.mat < min_exp] <- min_exp
genes <- rownames(exp.mat) # Checking genes

# Run Pathifier ################################################################
geneSets <- list()
for (i in 1:nrow(gs)){
  a <- as.vector(gs[i,3:ncol(gs)])
  a <- head(unique(na.omit(a)), -1)
  a <- matrix(a, ncol = 1)
  geneSets[[length(geneSets)+1]] <- a
  rm(a,i)
}

# Generate a list that contains the names of the genesets used
pathwaynames <- as.list(gs[,1])

# Generate a list that contains the previos two lists: genesets and their names
PATHWAYS <- list()
PATHWAYS$geneSets <- geneSets
PATHWAYS$pathwaynames <- pathwaynames

# Extract information from binary phenotypes. 1 = Normal, 0 = Tumor
normals <- as.logical(normals)
data <- exp.mat
allgenes <- rownames(data)

# Generate a list that contains previous data: gene expression, normal status,
# and name of genes
DATASET <- list(); DATASET$allgenes <- allgenes
DATASET$normals <- normals
DATASET$data <- data

# Run Pathifier
PDS <- quantify_pathways_deregulation(data = DATASET$data,
                                      allgenes = DATASET$allgenes,
                                      syms = PATHWAYS$geneSets,
                                      pathwaynames = PATHWAYS$pathwaynames, 
                                      normals = DATASET$normals,
                                      maximize_stability = F,
                                      attempts = 2,
                                      logfile="logfile.txt",
                                      min_std = min_std,
                                      min_exp = min_exp)


save.image("Check-Point.Rdata") # Check Point! ###########################
save(PDS, file = "TCGA_PDS.Rdata") # Save PDS object

# Post Pathifier filter ########################################################
# Get scores
PDSmatrix <- t(mapply(FUN = c, PDS$scores))
colnames(PDSmatrix) <- sampleIds

# Filtering pathways with lower median scores on normal samples (Best 75%)
rowMed <- apply(PDSmatrix[,normals], 1, median) # Median by row NORMALS
rowMed2 <- apply(PDSmatrix[,!normals], 1, median) # Med by row TUMORS

t <- 0.75 # Treshold
fPDS <- PDSmatrix[rowMed <= quantile(rowMed, t),]
out <- PDSmatrix[rowMed > quantile(rowMed, t),]
save(fPDS, file = "fPDS.Rdata")

ZscorePDS <- t(scale(t(fPDS))) # Zscore Transform
colnames(ZscorePDS) <- sampleIds
save(ZscorePDS, file = "ZscorePDS_TCGA.Rdata")

rowMed3 <- apply(fPDS[,normals], 1, median)
rowMed4 <- apply(fPDS[,!normals], 1, median)

save.image() # CheckPoint! ##################

# Plotting #############
###############################################################################
## Creating Custom Palette
###############################################################################

# creates a own color palette passing from blue, green yellow to dark red
my_palette <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(n = 1000))

###############################################################################
## Clustering Methods
###############################################################################

row.distance = dist(fPDS, method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")

col.distance = dist(t(fPDS), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")

###############################################################################
## Assign Column labels (Optional)
###############################################################################

colLabels <- as.character(normals)
colLabels[colLabels == "TRUE"] <- "#377EB8"
colLabels[colLabels == "FALSE"] <- "#E41A1C"

###############################################################################
## Plotting the Heatmap!! (where all colorful things happen...)
###############################################################################

png("PDS_Heatmap_TCGA.png", # Name of png file
    width = 6 * 500,      # Easier scaling 6*500 = 3000 pixels
    height = 6 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 10)        # font size

heatmap.2(fPDS,
          main = "PDS-Heatmap-TCGA",   # heat map title
          density.info = "none",  # turns off density plot inside color legend
          trace = "none",         # turns off trace lines inside the heat map
          margins = c(10,21),     # widens margins around plot
          col = my_palette,       # use on color palette defined earlier
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          Colv = as.dendrogram(col.cluster),  # apply selected clustering method
          labRow = "",
          labCol = "",
          keysize = 0.8,          # size of color key
          #Additional Options
          ## Color labeled columns
          ColSideColors = colLabels
)
## Legend for ColumnSide color labeling
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Normals", "Tumors"), # category labels
       col = c("dodgerblue", "firebrick1"),  # color key
       lty= 1,          # line style
       lwd = 5, unit    # line width
)
dev.off()               # close the PNG device

###############################################################################
## Zscore Heatmap
# Clustering Methods
row.distance = dist(ZscorePDS, method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")

col.distance = dist(t(ZscorePDS), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")

png("PDS-Zscore_Heatmap_TCGA.png", # Name of png file
    width = 6 * 500,      # Easier scaling 6*500 = 3000 pixels
    height = 6 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 10)        # font size

heatmap.2(ZscorePDS,
          main = "PDS-Zscore-Heatmap-TCGA",   # heat map title
          density.info = "none",  # turns off density plot inside color legend
          trace = "none",         # turns off trace lines inside the heat map
          margins = c(10,21),     # widens margins around plot
          col = my_palette,       # use on color palette defined earlier
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          Colv = as.dendrogram(col.cluster), # apply selected clustering method
          labRow = "",
          labCol = "",
          keysize = 0.8,          # size of color key
          #Additional Options
          ## Color labeled columns
          ColSideColors = colLabels
)
## Legend for ColumnSide color labeling
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Normals", "Tumors"), # category labels
       col = c("dodgerblue", "firebrick1"),  # color key
       lty= 1,          # line style
       lwd = 5, unit    # line width
)
dev.off()               # close the PNG device