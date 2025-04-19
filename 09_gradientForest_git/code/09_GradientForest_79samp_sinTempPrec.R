# Intro -------------------------------------------------------------------

# This script shows the code used to perform the Gradient Forest analysis 

# This analysis allows to evaluate and model the effects of different 
# predictors over the effect of the turnover of genomic diversity


# Packages needed ---------------------------------------------------------

library(SNPRelate)
library(sf)
library(tidyverse)
library(tmap)
library(tmaptools)
library(rnaturalearth)
library(rnaturalearthdata)
library(terra)
library(flexsdm)
library(geodata)
library(vegan)
library(geosphere)
library(gradientForest)

### Loading packages
library(raster)
library(rgdal)
#install.packages("gradientForest", repos="http://R-Forge.R-project.org")
library(gradientForest)
library(vegan)
library(caret)
library(sp)
library(tseries)

# Base maps for visual inspection -----------------------------------------
# Cargo limites políticos de Sudamerica

# Países
Sudam_raw <- ne_countries(scale = 10, 
                          continent = "south america", 
                          returnclass = "sf")

French_guiana <- ne_countries(scale = 10,
                              type = "map_units",
                              geounit = "French Guiana",
                              returnclass = "sf")
Sudam <- st_union(Sudam_raw, French_guiana)
# Mapa base
bb <- c(-83, # longitud mínima 
        -30, # latitud mínima
        -33, # longitud máxima
        12) # latitud máxima


mapa_base <- tm_shape(Sudam, bbox = bb) + tm_polygons(col = "lightgray") +
  tmap_options(check.and.fix = TRUE)


# Study sites -------------------------------------------------------------

Localities_df <- read.csv("00_Sampling_sites/Sampling_sites83.csv", sep = ";") |> 
  filter(!Sample == "UFG4028") |> 
  filter(!Sample == "UFG5668") |> 
  filter(!Sample == "UFG5753") |> 
  filter(!Sample == "UFG5639")
         
Localities_sf <- st_read("00_Sampling_sites/Sampling_sites83.gpkg")|> 
  filter(!Sample == "UFG4028") |> 
  filter(!Sample == "UFG5668") |> 
  filter(!Sample == "UFG5753") |> 
  filter(!Sample == "UFG5639")

# 1. Response variable ---------------------------------------------------------
## 1.1. vcftools preproccess (Linux)--------------------------------------------

# 1- cleaning the unliked vcf to retain biallelic SNPs and only sites without mising data

# Due that with 83 samples the retained SNP without missing data was too low, 
# we remove the three samples with high missing data (UFG4028, UFG5668, UFG5753, UFG5639)
vcftools --vcf 00_IPYRAD/Unlinked_data/Tangara_2023_83-II_randSNP.vcf --remove-indv UFG4028 --remove-indv UFG5668 --remove-indv UFG5753 --remove-indv UFG5639 --out 09_gradientForest/data/79samp_sinTempPrec/Stilpnia_79_unlink --recode
# Then apply filters
vcftools --vcf 09_gradientForest/data/79samp_sinTempPrec/Stilpnia_79_unlink.recode.vcf --min-alleles 2 --max-alleles 2 --max-missing 1 --out 09_gradientForest/data/79samp_sinTempPrec/Stilpnia_79_unlink_bial_mr0 --recode
# After filtering we kept 705 out of a possible 12791 Sites

# 2- Converting the SNP data to the format accepted by Gradient Forest
vcftools --vcf  09_gradientForest/data/79samp_sinTempPrec/Stilpnia_79_unlink_bial_mr0.recode.vcf --012 --out 09_gradientForest/data/79samp_sinTempPrec/Stilpnia_79_unlink_bial_mr0.gradforest.snp

# 3- Editting the name of the columns and exporting in final format (linux console)
cut -f2- Stilpnia_79_unlink_bial_mr0.gradforest.snp.012 | sed 's/-1/NA/g' >Stilpnia_79_unlink_bial_mr0.gradforest.snp.temp
cut -f1 <Stilpnia_79_unlink_bial_mr0.gradforest.snp.012.pos | tr -d '\t' | tr '\n' '\t' | sed 's/[[:space:]]*$//' > header
paste <(echo "ID" | cat - Stilpnia_79_unlink_bial_mr0.gradforest.snp.012.indv) <(echo "" | cat header - Stilpnia_79_unlink_bial_mr0.gradforest.snp.temp) > Stilpnia_79_unlink_bial_mr0.gradforest.snp.forR
rm header Stilpnia_79_unlink_bial_mr0.gradforest.snp.temp


## 1.2. Load genetic data ------------------------------------------------------
genotype_df <- read.table("09_gradientForest/data/79samp_sinTempPrec/Stilpnia_79_unlink_bial_mr0.gradforest.snp.forR",
                          header = T, row.names = 1)

# 2. Predictor variables ----------------------------------------------------

## 2.1. Climate ---------------------------------------------------------------

# consisted of uncorrelated (Pearson's r < .75) bioclimatic variables 
# obtained from WorldClim (Bio)
paths_CurrentBios <- list.files(path = "C:/Capas_GIS/ClimVars/WorldClim_1_4/bio_2-5m_bil/",
                                pattern = ".bil$",full.names = TRUE)
CurrentBios <- rast(paths_CurrentBios) |> 
  subset(c(1:9,12:17)) |> # Quito biovars 8, 9, 18 y 19
  crop(Sudam) |> # Corto por area de estudio
  mask(Sudam)

VIF <- correct_colinvar(CurrentBios, 
                        method = c('vif', th = '10'), 
                        proj = NULL)
CurrentBios_filter <- VIF$env_layer

Clima_df <- extract(CurrentBios_filter, Localities_sf) 

## 2.2. Tree cover -------------------------------------------------------
# Percentage of tree cover obtained from Geospatial Information Authority 
# of Japan, Chiba University and collaborating organization—version 1; 
# https://globalmaps.github.io/ptc. html)
Tree_cover <- rast("C:/Capas_GIS/Deforestacion/Percent_tree_coverage/gm_ve_v1.tif")
Tree_cover_df <- extract(Tree_cover, Localities_sf) |> 
  dplyr::select(gm_ve_v1) |> 
  rename(TC = gm_ve_v1)
## 2.3. Vegetation abundance ---------------------------------------------

# We used the maximum green vegetation fraction (MGVF; Broxton et al., 2014),
# as proxy of vegetation abundance.

MGVF <- rast("C:/Capas_GIS/Vegetation_index/MGVF_Broxton_2014/Average.tif/Average.tif")

MGVF_df <- extract(MGVF, Localities_sf) |> 
  dplyr::select(Average) |> 
  rename(MGVF = Average)



## 2.4. Elevation ---------------------------------------------------------------
Elevation <- rast("C:/Capas_GIS/Elevation/WorldClim_v21/wc2.1_30s_elev.tif")

Elevation_df <- extract(Elevation, Localities_sf) |> 
  dplyr::select(wc2.1_30s_elev) |> 
  rename(Elevation = wc2.1_30s_elev)

## 2.5. Spatial structure ------------------------------------------------

# PCNMs or MEMs (principal coordinates of neighbor matrices or Moran's eigenvector maps)
# Used to account for the influence of spatial processes and unmeasured 
# environmental variation
# Matrix with coordinates
Coords <- Localities_df |>
  dplyr::select(Longitude, Latitude) |> 
  as.matrix()
Matrix_distgeo <- distm(Coords) 

# Generate the PCNM matrix using vegan::pcnm
# you could stop here if you want all of them
pcnm_matrix <- pcnm(Matrix_distgeo)
pcnm_keep <- round(length(which(pcnm_matrix$value > 0))/2)
#keep half of positive ones as suggested by some authors
MEM_df <- scores(pcnm_matrix)[,1:pcnm_keep] |> 
  as.data.frame()


# 3. Gradient Forest analysis ---------------------------------------------
Stilpnia_env <- cbind(Clima_df, Tree_cover_df, 
                       MGVF_df, Elevation_df, MEM_df) |> 
  dplyr::select(!1)
Stilpnia_snp <- genotype_df

#a maximum number of splits can be defined following the developers suggestion
Stilpnia_maxLevel <- log2(0.368*nrow(Stilpnia_env)/2) 

Stilpnia_gf <- gradientForest(cbind(Stilpnia_env, Stilpnia_snp),
                     predictor.vars=colnames(Stilpnia_env),
                     response.vars=colnames(Stilpnia_snp), 
                     ntree=1000,
                     maxLevel=Stilpnia_maxLevel, 
                     trace = T, 
                     corr.threshold=0.50) 
# there will be warnings about having less than five values for 
# response variables, which is because we have only three: 0, 1, or 2. 
# You can ignore them.

# Save the gradient forest results
Stilpnia_gf
Stilpnia_gf$dens
write.table(Stilpnia_gf$result, "09_gradientForest/result/79samp_sinTempPrec/Stilpnia_gf_SNPs_R2.txt", sep="\t", quote=F)  
write.table(Stilpnia_gf$overall.imp, "09_gradientForest/result/79samp_sinTempPrec/Stilpnia_gf_vars_overall_imp.txt", sep="\t", quote=F)  
write.table(Stilpnia_gf$overall.imp2, "09_gradientForest/result/79samp_sinTempPrec/Stilpnia_gf_vars_overall_imp2.txt", sep="\t", quote=F)  
write.table(Stilpnia_gf$imp.rsq, "09_gradientForest/result/79samp_sinTempPrec/Stilpnia_gf_imp.rsq.txt", sep="\t", quote=F)  
write.table(Stilpnia_gf$res, "09_gradientForest/result/79samp_sinTempPrec/Stilpnia_gf_results.txt", sep="\t", quote=F)  
write.table(Stilpnia_gf$res.u, "09_gradientForest/result/79samp_sinTempPrec/Stilpnia_gf_res.u.txt", sep="\t", quote=F)  

### Plot main results -------------------------------------------------------------------------
# the mean accuracy importance and the mean importance weighted by species R2
png("09_gradientForest/plot/79samp_sinTempPrec/Stilpnia_splits_mean_importance.png")
plot(Stilpnia_gf, plot.type = "O")
dev.off()

pdf("09_gradientForest/plot/79samp_sinTempPrec/Stilpnia_splits_mean_importance.pdf")
plot(Stilpnia_gf, plot.type = "O")
dev.off()

# The splits density plot shows binned split importance and location on each gradient (spikes), kernel density of splits (black lines), of observations (red lines) and of splits standardised by observations density (blue lines). Each distribution integrates to predictor importance
png("09_gradientForest/plot/79samp_sinTempPrec/Stilpnia_splits_density_plot.png")
Stilpnia_most_important <- names(importance(Stilpnia_gf)) #[1:10]
plot(Stilpnia_gf, plot.type = "S", imp.vars = Stilpnia_most_important,
     leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6, cex.lab = 0.7,
     line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))
dev.off()

# The species cumulative plot, which for each species shows cumulative importance distributions of splits improvement scaled by R2 weighted importance, and standardised by density of observations.
png("09_gradientForest/plot/79samp_sinTempPrec/Stilpnia_species_cumulative_plot.png")
plot(Stilpnia_gf, plot.type = "C", imp.vars = Stilpnia_most_important,
     show.overall = F, legend = T, leg.posn = "topleft",
     leg.nspecies = 5, cex.lab = 0.7, cex.legend = 0.4,
     cex.axis = 0.6, line.ylab = 0.9, par.args = list(mgp = c(1.5,
                                                              0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, 0.3, 0, 0)))
dev.off()

# The predictor cumulative plot (plot.type="C", show.species=F), which for each predictor shows cumulative importance distributions of splits improvement scaled by R2 weighted importance, and standardised by density of observations, averaged over all species.
png("09_gradientForest/plot/79samp_sinTempPrec/Stilpnia_Predictor_cumulative_plot.png")
plot(Stilpnia_gf, plot.type = "C", imp.vars = Stilpnia_most_important,
     show.species = F, common.scale = T, cex.axis = 0.6,
     cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5,
                                                             0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0,0.3, 0, 0)))
dev.off()

# The R2 measure of the fit of the random forest model for each species, ordered in various ways.
png("09_gradientForest/plot/79samp_sinTempPrec/Stilpnia_overall_performance.png")
plot(Stilpnia_gf, plot.type = "P", show.names = T, horizontal = F,
     cex.axis = 1, cex.labels = 0.7, line = 2.5)
dev.off()

# Number of SNPs with R2 positive
number_SNPs_with_positive_R2 <- Stilpnia_gf$species.pos.rsq

# Mean SNP R2
mean_SNP_r2 <- mean(Stilpnia_gf$result)

# 4.Gradient Forest predictions -----------------------------------------------

## 4.1. Spatial prediction -------------------------------------------------

# We must resample the raster of predictors to get all the same resol and extent
qtm(CurrentBios_filter)
Tree_cover_rcl <- resample(Tree_cover, CurrentBios_filter)
MGVF_rcl <- resample(MGVF, CurrentBios_filter)
Elevation_rcl <- resample(Elevation, CurrentBios_filter)

# Load a the polygon of the study area usen for ENM 
M_Stilpnia <- st_read("07_ENM/ENMTML/M_Sptilpnia.gpkg")
qtm(M_Stilpnia)

# Create an stack with predictors layers and crop them by area M
predictors <- c(CurrentBios_filter, Tree_cover_rcl, MGVF_rcl, Elevation_rcl) |> 
  crop(M_Stilpnia)

# Extract predictors values for each cell
Stilpnia_clim_land <- extract(predictors, 1:ncell(predictors)) |> 
  mutate(ID = 1:ncell(predictors), .before = bio1) |> 
  rename(TC = gm_ve_v1) |> 
  rename(MGVF = Average) |> 
  rename(Elevation = wc2.1_30s_elev) |> 
  na.omit() # Remove cells with NA in any of the predictors

# Predict the allelic turnover across the study area 
# using our fitted model (gf) and the climate values 
# on the landscape (clim.land)
Stilpnia_pred <- predict(Stilpnia_gf, Stilpnia_clim_land[,-1])  #note the removal of the cell ID column with [,-1]

# These predictions then need to be converted to a color scale for mapping. 
# One way is to use principal components analysis (PCA) on the predictions
Stilpnia_PCs <- prcomp(Stilpnia_pred, center=T, scale.=F)

# And use the first 3 axes to define red, green, blue color scales, 
# respectively. 
r <- Stilpnia_PCs$x[, 1] # PC1 as red
g <- Stilpnia_PCs$x[, 2] # PC2 as green
b <- Stilpnia_PCs$x[, 3] # PC3 as blue
r <- (r - min(r))/(max(r) - min(r)) * 255 # Re-scale to plot in the red scale
g <- (g - min(g))/(max(g) - min(g)) * 255 # Re-scale to plot in the green scale
b <- (b - min(b))/(max(b) - min(b)) * 255 # Re-scale to plot in the blue scale

# Create a SpatRast with the same resol and extent of the predictors cropped
study_area <- (predictors$wc2.1_30s_elev)
study_area[]<-as.numeric(study_area[]>0)
qtm(study_area)

# From this SpatRast create a raster for each PC/color
rastR <- rastG <- rastB <- study_area

# Asigno
rastR[Stilpnia_clim_land$ID] <- r
qtm(rastR)
rastG[Stilpnia_clim_land$ID] <- g
qtm(rastG)
rastB[Stilpnia_clim_land$ID] <- b
qtm(rastB)
Stilpnia_rgb.rast <- c(rastR, rastG, rastB)

## 4.2. Map plot --------------------------------------------------------
png("09_gradientForest/plot/79samp_sinTempPrec/Stilpnia_GF_Map.png")
plotRGB(Stilpnia_rgb.rast, bgalpha=0)
points(Localities_df$Longitude, Localities_df$Latitude, 
       pch = 20, col = "yellow")
dev.off()

pdf("09_gradientForest/plot/79samp_sinTempPrec/Stilpnia_GF_Map.pdf")
plotRGB(Stilpnia_rgb.rast, bgalpha=0)
points(Localities_df$Longitude, Localities_df$Latitude, 
       pch = 20, col = "yellow")
dev.off()
