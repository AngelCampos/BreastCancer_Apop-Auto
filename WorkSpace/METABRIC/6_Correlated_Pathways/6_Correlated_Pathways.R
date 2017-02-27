################################################################################
# PDS heatmaps of correlated pathways 
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

#Load data
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/2_PDS_MI/PathsCorrelation_0.05.Rdata")
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/4_Subgroups/Subgroups.Rdata")
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/1_Deregulation_Analysis/fPDS.Rdata")
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/0_ExpressionData/normals_TCGA_onlyFemales.Rdata")
load("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/0_Genesets/GOids_Apoptosis_Autophagy_Membership2.txt")
M <- as.matrix(read.delim("C:/Users/Fang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/TCGA/0_Genesets/GOids_Apoptosis_Autophagy_Membership2.txt", header = F))

# All samples heatmap ONLY selected pathways
rownames(M) <- M[,2]
selected <- M[unique(c(PWnames[,c(1,2)])),1]
fPDS <- fPDS[selected, ]

### PLOT ####
# creates a own color palette passing from blue, green yellow to dark red
my_palette <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(n = 1000))

## Clustering Method
row.distance = dist(fPDS, method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")

col.distance = dist(t(fPDS), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")

## Assign Column labels (Optional)
colLabels <- as.character(normals)
colLabels[colLabels == "TRUE"] <- "#377EB8"
colLabels[colLabels == "FALSE"] <- "#E41A1C"

## Plotting the Heatmap!! (where all colorful things happen...)
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




# Subgroups heatmap


