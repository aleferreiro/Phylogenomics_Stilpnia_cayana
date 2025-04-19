
# Packages needed ---------------------------------------------------------

library(tidyverse)
library(vcfR)
library(adegenet)
library(ape)
library(ggplot2)


# 1. Load the VCF file ----------------------------------------------------

vcf <- read.vcfR("00_IPYRAD/data/assembly_80/Stilpnia_cayana_80_bial_maf_mr80_noIndel.recode.vcf")  # Replace with actual VCF file path


# 2. Convert VCF to a genlight object -------------------------------------

genlight_obj <- vcfR2genlight(vcf)



# 3. Convert VCF to a genlight object -------------------------------------

gen_dist <- dist(as.matrix(genlight_obj))



# 4. # Perform Principal Coordinates Analysis (PCoA) ----------------------

pcoa_result <- pcoa(gen_dist)

# PCs contributions
var_explained <- pcoa_result$values$Relative_eig * 100
var_explained

# Extract scores for first two axes
pcoa_df <- data.frame(
  Sample = rownames(pcoa_result$vectors),
  PC1 = pcoa_result$vectors[,1],
  PC2 = pcoa_result$vectors[,2],
  Populations = Sampling_sites
)

# 5. Plot -----------------------------------------------------------------
# Load population data
Sampling_sites <- read.csv("00_Sampling_sites/Sampling_sites83.csv", sep = ";") |> 
  filter(Sample %in% pcoa_df$Sample) 
  
# Columnas con poblaciones
pcoa_df_pop <- pcoa_df |> 
  mutate(Populations = Sampling_sites$Populations2)

# Order populations
pcoa_df_pop$Populations <- factor(pcoa_df_pop$Populations, levels = c("Guyana",
                                                                      "Para",
                                                                      "Amapa",
                                                                      "Marajo",
                                                                      "Caatinga",
                                                                      "East",
                                                                      "West"))
# Plot with ggplot
PCoAplot <- ggplot(pcoa_df_pop, aes(PC1, PC2, colour = Populations)) +
  geom_point(size = 3) +
  labs(x = "PC1 (18.31%)",
       y = "PC2 (2.01%)", 
       colour = "Populations",
       title = "PCoA") +
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
PCoAplot
ggsave("01_PCoA/plot/01_Scay_83_80samp_bial_maf_mr80_noIndel_PCoAplot.pdf", 
       PCoAplot,
       width = 7,
       height = 4.17)
ggsave("01_PCoA/plot/01_Scay_83_80samp_bial_maf_mr80_noIndel_PCoAplot.png", 
       PCoAplot,
       width = 7,
       height = 4.17)
