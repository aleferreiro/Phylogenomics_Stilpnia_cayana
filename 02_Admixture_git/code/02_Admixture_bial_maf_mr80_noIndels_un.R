# Packages needed ---------------------------------------------------------

library(LEA)
library(tidyverse)
library(ggplot2)
library(forcats)
library(ggthemes)
library(patchwork)
library(ggplot2)

# 1. Prepare genomic dataset -------------------------------------------------

# Covert vcf with unlinked SNPs to geno format
vcf2geno("00_IPYRAD/data/assembly_80/Stilpnia_cayana_80_bial_maf_mr80_noIndel_un.vcf", 
         "02_Admixture/data/80/bial_maf_mr80_noIndels_un/Scay_80_bial_maf_mr80_noIndels_un.geno")

# 2. Load data ------------------------------------------------------------

## Load geno file ---------------------------------------------------------
genofile <- read.geno("02_Admixture/data/80/bial_maf_mr80_noIndels_un/Scay_80_bial_maf_mr80_noIndels_un.geno")


## Load population data ----------------------------------------------------
# Data frames with populations of the 83 samples
data_pop <- read.csv("00_Sampling_sites/Sampling_sites83.csv") |> 
  filter(Sample != "UFG4028",
         Sample != "UFG5668",
         Sample != "UFG5753")

# 1. Admixture models ----------------------------------------------------------
# main options
# K = number of ancestral populations
# entropy = TRUE computes the cross-entropy criterion,
# CPU = 4 is the number of CPU used (hidden input)
project_Scay = NULL
project_Scay = snmf("02_Admixture/data/80/bial_maf_mr80_noIndels_un/Scay_80_bial_maf_mr80_noIndels_un.geno",
                             K = 1:10,
                             entropy = TRUE,
                             repetitions = 10,
                             project = "new")

# An snmfProject can be load in a different session.
project_Scay = load.snmfProject("02_Admixture/data/80/bial_maf_mr80_noIndels_un/Scay_80_bial_maf_mr80_noIndels_un.snmfProject") 

# 2. Plot to select K -----------------------------------------------------
# plot cross-entropy criterion for all runs in the snmf project
pdf(file="02_Admixture/plot/80/bial_maf_mr80_noIndels_un/Stilpnia_80_bial_un_k_select.pdf")
plot(project_Scay, col = "blue", pch = 19, cex = 1.2)
dev.off()

png(file="02_Admixture/plot/80/bial_maf_mr80_noIndels_un/Stilpnia_80_bial_un_k_select.png")
plot(project_Scay, col = "blue", pch = 19, cex = 1.2)
dev.off()


# 3. STRUCTURE-like plot --------------------------------------------------

## a. K2 plot -----------------------------------------------------------------

# select the best run for K = 2 clusters
best2 = which.min(cross.entropy(project_Scay, K = 2))

# Create matrix of ancestry proportions:
qmatrix_k2 = LEA::Q(project_Scay, K = 2, run = best2)

# Create a data frame with info for each sample
k2_df <- as.data.frame(qmatrix_k2) |> 
  mutate(sample_id = data_pop$Sample,
    pop = data_pop$Populations,
    pop2 = data_pop$Populations2,
    subsp = data_pop$Subespecie) |> 
  pivot_longer(1:2)


k2plot <-
  ggplot(k2_df, aes(sample_id, value, fill = name)) +
  geom_col(color = NA, lwd = 0.1) +
  facet_grid(~factor(pop2, 
                     levels = c("Guyana", "Para", 
                                "Amapa", "Marajo", 
                                "Caatinga", 
                                "East", "West")),
             switch = "x", 
             scales = "free", 
             space = "free") +
  theme_minimal() + labs(x = "", title = "K=2", y = "") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    # axis.text = element_blank(), sin nombre de las muestras
    axis.text.x = element_text(size = 6,angle = 90, vjust = 0.5, hjust=1),
    panel.grid = element_blank(),
    strip.text.x = element_blank() # remove facets names
  ) +
  scale_fill_grey(guide = FALSE)
k2plot

# Save individual plot
ggsave("02_Admixture/plot/80/bial_maf_mr80_noIndels_un/Stilpnia_80_bial_maf_mr80_noIndels_un_K2withSamples.pdf", 
       k2plot, width = 12.58, height = 4.33)


## b. K3 plot -----------------------------------------------------------------

# select the best run for K = 3 clusters
best3 = which.min(cross.entropy(project_Scay, K = 3))

# Create matrix of ancestry proportions:
qmatrix_k3 = LEA::Q(project_Scay, K = 3, run = best3)

k3_df <- as.data.frame(qmatrix_k3) |> 
  mutate(sample_id = data_pop$Sample,
         pop = data_pop$Populations,
         pop2 = data_pop$Populations2,
         subsp = data_pop$Subespecie) |> 
  pivot_longer(1:3)


k3plot <-
  ggplot(k3_df, aes(sample_id, value, fill = name)) +
  geom_col(color = NA, lwd = 0.1) +
  facet_grid(~factor(pop2, 
                     levels = c("Guyana", "Para", 
                                "Amapa", "Marajo", 
                                "Caatinga", 
                                "East", "West")), 
             switch = "x", 
             scales = "free", 
             space = "free") +
  theme_minimal() + labs(x = "", title = "K=3", y = "") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    strip.text.x = element_blank() # remove facets names
  ) +
  scale_fill_grey(start = 0.2, end = 0.8, guide = FALSE)
k3plot

# Save individual plot
# ggsave( "plots/Tangara_2023_83-II/02a_randSNP_Admixture_K3plot.pdf", k3plot, width = 12.58, height = 4.33)

## c. K4 plot -----------------------------------------------------------------

# select the best run for K = 3 clusters
best4 = which.min(cross.entropy(project_Scay, K = 4))

# Create matrix of ancestry proportions:
qmatrix_k4 = LEA::Q(project_Scay, K = 4, run = best4)

k4_df <- as.data.frame(qmatrix_k4) |> 
  mutate(sample_id = data_pop$Sample,
         pop = data_pop$Populations,
         pop2 = data_pop$Populations2,
         subsp = data_pop$Subespecie) |> 
  pivot_longer(1:4)


k4plot <-
  ggplot(k4_df, aes(sample_id, value, fill = name)) +
  geom_col(color = NA, lwd = 0.1) +
  facet_grid(~factor(pop2, 
                     levels = c("Guyana", "Para", 
                                "Amapa", "Marajo", 
                                "Caatinga", 
                                "East", "West")), 
             switch = "x", 
             scales = "free", 
             space = "free") +
  theme_minimal() + labs(x = "", title = "K=4", y = "") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    strip.text.x = element_blank() # remove facets names
  ) +
  scale_fill_grey(start = 0.2, end = 0.8, guide = FALSE)

k4plot

# Save plot individuall
#ggsave( "plots/Tangara_2023_83-II/02a_randSNP_Admixture_K4plot.pdf", k3plot, width = 12.58, height = 4.33)

## d. Pop axis plot ------------------------------------------------------

# x-axis plot to with colors of each population
kplot_axis_pop <- 
  ggplot(k2_df, aes(sample_id, y = 1, fill = pop2)) +
  geom_tile() +
  facet_grid(~factor(pop2, 
                     levels = c("Guyana", "Para", 
                                "Amapa", "Marajo", 
                                "Caatinga", 
                                "East", "West" )), 
             switch = "x", scales = "free", space = "free") +
  theme_void() +
  scale_fill_manual(values = c("#e7298a", # Amapa
                               "#54278f", # Caatinga
                               "#dab600", # East
                               "#d4b9da", # Guyana
                               "#a50f15", # Marajo
                               "#df65b0", # Para
                               "#6baed6")) + # West 
  theme(legend.position = 'none')
kplot_axis_pop

# ggsave( "plots/Tangara_2023_83-II/02_Admixture_popAxis.pdf", k2plot_axis_pop, width = 12.58, height = 1)

## e. Multi K plots --------------------------------------------------------

plot_admix <- k2plot + k3plot + k4plot + 
  kplot_axis_pop +
  plot_layout(ncol = 1)
plot_admix

ggsave("02_Admixture/plot/80/bial_maf_mr80_noIndels_un/Stilpnia_83_bial_maf_mr80_noIndels_un_Admixture_plot.pdf", 
       plot_admix, width = 10, height = 5.73)
ggsave("02_Admixture/plot/80/bial_maf_mr80_noIndels_un/Stilpnia_80_bial_maf_mr80_noIndels_un_Admixture_plot.png", 
       plot_admix, width = 10, height = 5.73)


