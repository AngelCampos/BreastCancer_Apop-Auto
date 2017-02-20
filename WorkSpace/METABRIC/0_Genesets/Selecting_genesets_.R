################################################################################
# Select Apoptosis and Autophagy Processes & Prepare 4 PATHIFIER
## Author: Miguel Angel Garcia Campos - github: https://github.com/AngelCampos
################################################################################

genesets <- as.matrix(read.delim("PATHIFIER_genesets_apop_autop.gmt", header = F, sep = "\t", row.names = 1))
GODB <- as.matrix(read.delim("GO_biologicalprocess_adjusted.txt", header = F, sep = "\t", row.names = 1))

# Selecting GO terms to extract
ids <- substring(rownames(genesets),1,10)
names <- substring(rownames(genesets), 12)
N <- NULL; for(i in 1:length(names)){
  x <- substring(names[i],0, nchar(names[i])-16); N[i] <- x}
AAids <- cbind(ids,N)
rm(N, x, i, names)

# Some ids may be deprecated, we have to filter them out
sum(ids %in% rownames(GODB))/ length(ids) # Percentage of present ids in GODB
present <- which(ids %in% rownames(GODB)) # Filtered ids vector
gs <- GODB[ids[present], -1]

# Removing genesets containing less than 3 genes
S <- NULL; for(i in 1:nrow(gs)){s <- length(unique(gs[i,-1]))-1; S[i] <- s}
gs <- gs[S >= 3,]

# Write tables
write.table(AAids, file = "GOids_Apoptosis_Autophagy.txt", quote = F, row.names = F, col.names = F, sep= "\t")
write.table(x= gs, file= "GOBP_Genesets_Apoptosis_Autophagy.txt", quote = F, row.names = T, col.names = F, sep = "\t")
