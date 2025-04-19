# Packages needed -----------------------------------------------------
library(tidyverse)
library(SNPRelate)

# Create diff file --------------------------------------------------------
# Generate a matrix with the average dissimilaties between individuals, 
# estimated from SNPRelate

# Open the GDS file
genofile <- snpgdsOpen("data/Tangara_2023_83-II/Tangara_2023_83-II.gds")


# Create coords file ------------------------------------------------
# Open the GDS file
genofile <- snpgdsOpen("data/Tangara_2023_83-II/Tangara_83-II.gds")

# Generate dissimilarity matrix
Diss_obj <- snpgdsDiss(genofile, autosome.only = F)
Diss_matrix <- Diss_obj$diss
write.table(Diss_matrix, 
            file="data/Tangara_2023_83-II/05_EEMS/Tangara_2023_83-II_SNPRelateDiss.txt", 
            sep = "\t",
            row.names=FALSE, 
            col.names=FALSE)

# Otra opcion que uso Cassia en su paper:
# Generate similarity matrix 
IBS_obj <- snpgdsIBS(genofile, autosome.only = F)
IBS_matrix <- IBS_obj$ibs

# Substract by 1 to obtain a dissimilarity index
DisIBS_matrix <- 1 - IBS_matrix
write.table(DisIBS_matrix, 
            file="data/Tangara_2023_83-II/05_EEMS/Tangara_2023_83-II_SNPRelateDisIBS.txt", 
            sep = "\t",
            row.names=FALSE, 
            col.names=FALSE)

# create outer file -------------------------------------------
remotes::install_github("JeffWeinell/misc.wrappers")

# La cree directamente en QGIS! Para evitar bardo

