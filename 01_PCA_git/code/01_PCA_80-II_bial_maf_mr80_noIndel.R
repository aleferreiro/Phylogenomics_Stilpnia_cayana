# Packages needed -----------------------------------------------------------

library(tidyverse)
library(SNPRelate)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(patchwork)


# 1. VCF conversion and filter --------------------------------------

## Load VCF file

Tangara_vcf_path <- "00_IPYRAD/data/assembly_80/Stilpnia_cayana_80_bial_maf_mr80_noIndel.recode.vcf"

## Convert to GDS file 

# This software uses the format GDS to perform analysis, then we reformat
# In this new format only bial_maf_mr80_noIndellelic snps were retained.
snpgdsVCF2GDS(Tangara_vcf_path, "01_PCA/data/80/Stilpnia_cayana_80_bial_maf_mr80_noIndel.gds", method="copy.num.of.ref")

# Summary
snpgdsSummary("01_PCA/data/80/Stilpnia_cayana_80_bial_maf_mr80_noIndel.gds") 

# 2. Load data ---------------------------------------------------------------

## Load gds data -----------------------------------------------------------
# Open the GDS file generated in the previous step
genofile <- snpgdsOpen("01_PCA/data/80/Stilpnia_cayana_80_bial_maf_mr80_noIndel.gds")

## Load population data to plot --------------------------------------------
Sampling_sites <- read.csv("00_Sampling_sites/Sampling_sites83.csv",
                           sep = ";")

# 3. PCA ------------------------------------------------------------
# Run PCA
pca_mr01 <- snpgdsPCA(genofile, num.thread=1, autosome.only = F,
                      missing.rate = NaN) # remuevo SNPs sin datos en al menos 10% de las muestras

# variance proportion (%)
pc.percent_mr01 <- pca_mr01$varprop*100
head(round(pc.percent_mr01, 2))

write.csv(as.data.frame(pc.percent_mr01), "01_PCA/result/80/bial_maf_mr80_noIndel/Scay_variance_proportion.csv")

# make a data.frame
tab_mr01 <- data.frame(sample.id = pca_mr01$sample.id,
                       EV1 = pca_mr01$eigenvect[,1],    # the first eigenvector
                       EV2 = pca_mr01$eigenvect[,2],    # the second eigenvector
                       EV3 = pca_mr01$eigenvect[,3],    # the third eigenvector
                       EV4 = pca_mr01$eigenvect[,4],    # the fourth eigenvector
                       stringsAsFactors = FALSE)

# Guardo la matriz con los valores de los PC
write.csv(tab_mr01, "01_PCA/result/80/bial_maf_mr80_noIndel/Scay_EVvalues.csv")

# Generate a data frame to use hereafter with the sample id, the pop and the subsp 
# of the 80 samples used
data_pop <- Sampling_sites |> 
  filter(Sample %in% tab_mr01$sample.id) |> # me quedo con las muestras usadas en el PCA
  dplyr::select(Sample, Populations2, Subespecie, Longitude, Latitude) |> 
  rename(Populations = Populations2 )

write.csv(data_pop,
          file = "00_Sampling_sites/Populations_Stilpnia80.csv")
# Columnas con poblaciones
tab_pop_mr10 <- bind_cols(tab_mr01, data_pop)

# By Populations
# Order populations north to south
tab_pop_mr10$Populations <- factor(tab_pop_mr10$Populations, levels = c("Guyana",
                                                                        "Para",
                                                                        "Amapa",
                                                                        "Marajo",
                                                                        "Caatinga",
                                                                        "East",
                                                                        "West"))

PCAplot_mr10 <- ggplot(tab_pop_mr10, aes(EV1, EV2, colour = Populations)) +
  geom_point(size = 3) +
  labs(x = "PC1 (13.56%)",
       y = "PC2 (2.01%)", 
       colour = "Populations",
       title = "80 samples, biallelic, maf (0.05), mr80 (20,278 SNPs)") +
  scale_color_manual(
    values = c("#d4b9da", # Guyana
               "#df65b0", # Para
               "#e7298a", # Amapa
               "#a50f15", # Marajo
               "#54278f", # Caatinga
               "#dab600", # East
               "#6baed6")) + # West 
  # geom_label_repel(aes(label = sample.id),
  #                  box.padding   = 0.35, 
  #                  point.padding = 0.5,
  #                  segment.color = 'grey50',
  #                  max.overlaps = 10,
  #                  show.legend = F) +
  theme_bw()
PCAplot_mr10
ggsave("01_PCA/plot/80/bial_maf_mr80_noIndel/01_PCAplot.pdf", 
       PCAplot_mr10,
       width = 7,
       height = 4.17)
ggsave("01_PCA/plot/80/bial_maf_mr80_noIndel/01_PCAplot.png", 
       PCAplot_mr10,
       width = 7,
       height = 4.17)

PCAplot_mr10 <- ggplot(tab_pop_mr10, aes(EV3, EV4, colour = Populations)) +
  geom_point(size = 3) +
  labs(x = "PC3 (1.82%)",
       y = "PC4 (1.77%)", 
       colour = "Populations",
       title = "80 samples, biallelic, maf (0.05), mr80 (20,278 SNPs)") +
  scale_color_manual(
    values = c("#d4b9da", # Guyana
               "#df65b0", # Para
               "#e7298a", # Amapa
               "#a50f15", # Marajo
               "#54278f", # Caatinga
               "#dab600", # East
               "#6baed6")) + # West 
  # geom_label_repel(aes(label = sample.id),
  #                  box.padding   = 0.35, 
  #                  point.padding = 0.5,
  #                  segment.color = 'grey50',
  #                  max.overlaps = 10,
  #                  show.legend = F) +
  theme_bw()
PCAplot_mr10
ggsave("01_PCA/plot/80/bial_maf_mr80_noIndel/01_PCAplot_PC3y4.pdf", 
       PCAplot_mr10,
       width = 7,
       height = 4.17)
ggsave("01_PCA/plot/80/bial_maf_mr80_noIndel/01_PCAplot_PC3y4.png", 
       PCAplot_mr10,
       width = 7,
       height = 4.17)