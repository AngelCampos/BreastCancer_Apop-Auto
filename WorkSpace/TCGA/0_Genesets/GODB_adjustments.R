################################################################################
# Adjustments to Gene Ontology DB downloaded from Bader et al.,
# http://baderlab.org/GeneSets
# Merico, Daniele, et al. "Enrichment map: a network-based method for gene-set 
# enrichment visualization and interpretation." PloS one 5.11 (2010)
# Script author: Miguel Angel Garcia-Campos - https://github.com/AngelCampos
################################################################################

# Load Gene Ontology DB, remove genesets with no GO identifier ####
# "con" can be changed to get the latest GO database at baderlab.org/EM_Genesets 
con <- "Human_GOBP_AllPathways_no_GO_iea_January_01_2017_symbol.gmt"
no_col <- max(count.fields(con, sep="\t"))-1
GODB <- read.delim(con, col.names = 1:no_col, header = F, fill = T)

# Functions ######
# Extract GO ids in format "GO:XXXXXXXX"
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# Filter BioProcesses ####
BP <- grep(GODB[,1], pattern = "GOBP")
PDB <- GODB[BP,]

# Modify labels to remove white spaces ####
labels <- as.character(PDB[,1])
x <- gsub(" ", "_", labels)
y <- gsub("%", "_", x)
PDB[,1] <- y

# Extract labels and assign to rownames
etiquetas <- substrRight(labels, 10)
rownames(PDB) <- etiquetas

# Write adjusted PDB to TXT file
write.table(x = PDB, file = "GO_biologicalprocess_adjusted.txt", sep = "\t",
            row.names = TRUE, col.names = FALSE )