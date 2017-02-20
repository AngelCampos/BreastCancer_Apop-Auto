################################################################################
#
## Author: Miguel Angel Garcia Campos github: https://github.com/AngelCampos
################################################################################
# Load packages
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("stats")) {
  install.packages("stats", dependencies = TRUE)
  library(stats)
}
if (!require("gridExtra")) {
  install.packages("gridExtra", dependencies = TRUE)
  library(gridExtra)
}

# Load data
load("C:/Users/WildFang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/METABRIC/1_Deregulation_Analysis/fPDS.Rdata")
load("C:/Users/WildFang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/METABRIC/0_ExpressionData/normals_METABRIC.Rdata")
load("C:/Users/WildFang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/METABRIC/0_ExpressionData/official_PAM50.Rdata")
load("C:/Users/WildFang/Google Drive/Science/BreastCancer_AA Project/WorkSpace/METABRIC/3_PAM50/genefu_PAM50.Rdata")

# Hierarchichal Clustering #####################################################
# Clustering Method
row.distance = dist(fPDS, method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")
col.distance = dist(t(fPDS), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")

# Subgroups (5)
cut <- 5
mytree <- cutree(col.cluster, cut)
save(mytree, file = "Subgroups.Rdata" )
# Checking if all normals are in the same group
sum(mytree[normals] == mytree[normals[1]]) / length(mytree[normals]) # TRUE!

# Heatmap ######################################################################
# Color parameters for Heatmap
# Creating Custom Palette
my_palette <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(n = 1000))

# Column labels
ccol <- brewer.pal(name = "Set1", cut)
colLabels <- ccol[mytree]

## Plotting the Heatmap!! (where all colorful things happen...) ################
png("Subgroups_METABRIC.png", # Name of png file
    width = 6 * 500,      # Easier scaling 6*500 = 3000 pixels
    height = 6 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 10)        # font size

heatmap.2(fPDS,
          main = "PDS-Subgroups-METABRIC",   # heat map title
          density.info = "none",  # turns off density plot inside color legend
          trace = "none",         # turns off trace lines inside the heat map
          margins = c(10,21),     # widens margins around plot
          col = my_palette,       # use on color palette defined earlier
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          Colv = as.dendrogram(col.cluster),  # apply selected clustering method
          labRow = "",
          labCol = "",
          keysize = 0.8,          # size of color key
          ColSideColors = colLabels
)

## Legend for ColumnSide color labeling
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Healthy", "SG1", "SG2", "SG3", "SG4"), # category labels
       col = ccol,  # color key
       lty= 1,          # line style
       lwd = 5, unit    # line width
)
dev.off()               # close the PNG device

# Hypergeometric Test of PAM50 Subtype match ###################################
names(mytree) <- gsub(pattern = "\\.", replacement = "-", x = names(mytree))
healthy <- rep("Normal", sum(normals))
names(healthy) <- names(mytree)[as.logical(normals)]
officialPAM50 <- c(healthy, officialPAM50)
officialPAM50 <- officialPAM50[names(mytree)]
subg <- mytree

# Obtener numero de muestras por subtipo en cada subgrupo
y0 <- NULL
for (i in c("Basal", "LumA", "LumB", "Her2", "Normal")){
  g <- as.logical(officialPAM50 == i) & as.logical(subg == 1)
  y <- sum(g)
  y0 <- cbind(y0, y) 
}
y1 <- NULL
for (i in c("Basal", "LumA", "LumB", "Her2", "Normal")){
  g <- as.logical(officialPAM50 == i) & as.logical(subg == 2)
  y <- sum(g)
  y1 <- cbind(y1, y) 
}
y2 <- NULL
for (i in c("Basal", "LumA", "LumB", "Her2", "Normal")){
  g <- as.logical(officialPAM50 == i) & as.logical(subg == 3)
  y <- sum(g)
  y2 <- cbind(y2, y) 
}
y3 <- NULL
for (i in c("Basal", "LumA", "LumB", "Her2", "Normal")){
  g <- as.logical(officialPAM50 == i) & as.logical(subg == 4)
  y <- sum(g)
  y3 <- cbind(y3, y) 
}
y4 <- NULL
for (i in c("Basal", "LumA", "LumB", "Her2", "Normal")){
  g <- as.logical(officialPAM50 == i) & as.logical(subg == 5)
  y <- sum(g)
  y4 <- cbind(y4, y) 
}

sub_table <- rbind(y0, y1, y2, y3, y4)
rownames(sub_table) <- c("SG0","SG1","SG2","SG3","SG4")
colnames(sub_table) <- c("Basal","LumA", "LumB", "Her2", "Normal")

# Hypergeometric Test 
pop <- sum(sub_table) # Total population
# Basal
pv1 <- NULL
for (i in 1:5){
  m <- sum(sub_table[,1])
  n <- pop - m
  k <- sum(sub_table[i,])
  x <- sub_table[i,1]
  pv <- phyper(x, m, n, k, lower.tail = F)
  pv1 <- signif(c(pv1, pv), 5)
}
# LumA
pv2 <- NULL
for (i in 1:5){
  m <- sum(sub_table[,2])
  n <- pop - m
  k <- sum(sub_table[i,])
  x <- sub_table[i,2]
  pv <- phyper(x, m, n, k, lower.tail = F)
  pv2 <- signif(c(pv2, pv), 5)
}
# LumB
pv3 <- NULL
for (i in 1:5){
  m <- sum(sub_table[,3])
  n <- pop - m
  k <- sum(sub_table[i,])
  x <- sub_table[i,3]
  pv <- phyper(x, m, n, k, lower.tail = F)
  pv3 <- signif(c(pv3, pv), 5)
}
# Her2
pv4 <- NULL
for (i in 1:5){
  m <- sum(sub_table[,4])
  n <- pop - m
  k <- sum(sub_table[i,])
  x <- sub_table[i,4]
  pv <- phyper(x, m, n, k, lower.tail = F)
  pv4 <- signif(c(pv4, pv), 5)
}
# Normal
pv5 <- NULL
for (i in 1:5){
  m <- sum(sub_table[,5])
  n <- pop - m
  k <- sum(sub_table[i,])
  x <- sub_table[i,5]
  pv <- phyper(x, m, n, k, lower.tail = F)
  pv5 <- signif(c(pv5, pv), 5)
}


# Multiple comparison adjustment
pv1 <- p.adjust(p = pv1, method = "bonferroni")
pv2 <- p.adjust(p = pv2, method = "bonferroni")
pv3 <- p.adjust(p = pv3, method = "bonferroni")
pv4 <- p.adjust(p = pv4, method = "bonferroni")
pv5 <- p.adjust(p = pv5, method = "bonferroni")

pv_table <- cbind(pv1, pv2, pv3, pv4, pv5)
rownames(pv_table) <- c("SG0","SG1","SG2","SG3","SG4")
colnames(pv_table) <- c("Basal", "LumA", "LumB", "Her2", "Normal")

# Plotting results in table
png(filename = "p-values.png", width = 1800, 
    height = 600, res = 300, pointsize = 1)
grid.table(pv_table)
dev.off()

# CheckPoint! ##################################################################
save.image("CheckPoint_4.Rdata")
load("CheckPoint_4.Rdata")