# Paquetes ----------------------------------------------------------------

library(tidyverse)
library(sf)
library(terra)
library(geosphere)      # Calculate geographic distance
# library(gdistance)      # Calculate geographic distance 2 and resistance distances
# library(leastcostpath)  # Calculate resistance distances
library(SNPRelate)      # Genetic distance
library(olsrr)          # Residuals evaluation
library(fmsb)           # VIF function
library(yhat)           # Commonality analysis
library(ecodist)
library(descr)          # addon for Commonality analysis
source("08_MMRR_and_CA/code/CA_logistic.R") # script for calculate comm coeff for log reg
source("08_MMRR_and_CA/code/CA_semiStdCoef.R")

# Plotting
library(genio)
library(ggplot2) # plots figures
library(gt) # plot tables
library(tmap) # plot maps
library(tmaptools) # plot maps


# Functions needed --------------------------------------------------------

# Z-transformation (standardization  by  subtracting  the  mean  and  dividing  by  the  standard  deviation)
matrix_standarization <- function(x){ 
  x=data.matrix(x)
  x=x[upper.tri(x, diag = FALSE)]
  x=(x-mean(x))/sqrt(var(x)) 
}

# Logaritmic conversion
lg <- function(x){ exp(x)/(1 + exp(x)) }


# 1. Geographic distance -------------------------------------------------
# Se utilizaron distancias geográficas geodésicas entre individuos (DGEO), estimadas
# con la función distm() del paquete R “geosphere”, para probar aislamiento
# por distancia.

# Cargo las coordenadas de cada individuo
Localities83 <- read.csv("00_Sampling_sites/Populations_Stilpnia83.csv")

Coords <- read.csv("00_Sampling_sites/Populations_Stilpnia83.csv") |>
  dplyr::select(Longitude, Latitude) |> 
  as.matrix()

Coords_sf <- Coords |> 
  as.data.frame() |> 
  st_as_sf(coords = c("Longitude","Latitude"),
           crs = 4326)

Matrix_distgeo <- distm(Coords) 
write.table(Matrix_distgeo, 
            file="08_MMRR_and_CA/data/Matrix_distgeo.txt", 
            row.names=FALSE, col.names=FALSE)
Matrix_distgeo <- read_matrix("08_MMRR_and_CA/data/Matrix_distgeo.txt")

# Standarization
DistGeog <- matrix_standarization(Matrix_distgeo)

# 2. Present resistance -------------------------------------------------
# Resitance was calculated as 1 - conductance. We used 
# suitability surfaces (ENM) as conductance surfaces. Then to calculate resistance,
# we substracted 1 to the suitability layer.

# First generate input files to process them and obtain the matrix of 
# resistance between localities in the program Circuitscape v5.0.0 that run 
# in the CCAD cluster. 

# First we generate a raster of the occurence data.
# Load current suitability layer

## 2.1.Input data for Circuitscape ---------------------------------------------
### 2.1.1. Resistance raster ----------------------------------------------------
ENM_Present <- rast("07_ENM/ENMTML/Result9/Projection/Present/Ensemble/W_MEAN/Stilpnia_cayana.tif") 
# Set negative values to 0
ENM_Present[ENM_Present < 0] <- 0
qtm(ENM_Present)

# Convert conductance to resistance
Present_Resist <- 1-ENM_Present 
qtm(Present_Resist) # plot to check
# Save to use it in CS
writeRaster(Present_Resist,
            "08_MMRR_and_CA/data/Circuitscape/Present_Resist.asc",
            overwrite = T)

### 2.1.2. Focal nodes raster ----------------------------------------------------

# Generate a sample id for each record
Coords_sf$id <- Localities83$X 
Coords_spatvect <- vect(Coords_sf)

# Create a raster of the presence records, where the value of the cell is 
# the sample id
Present_Coords <- terra::rasterize(Coords_spatvect, 
                                   Present_Resist, 
                                   field = "id")
qtm(Coords_rast)
# IMPORTANTE GENERARLO CON EL RASTER QUE SE USARA EN EL CS

# Save it to use with CircuitScape
writeRaster(Present_Coords,
            "08_MMRR_and_CA/data/Circuitscape/Present_Coords.asc",
            overwrite = T)


### 2.1.3. ini file ---------------------------------------------------------------

# Make an .ini file
Present_CS_ini <- c("[circuitscape options]",
                    "data_type = raster",
                    "scenario = pairwise",
                    "point_file = Present_Coords.asc",
                    "habitat_file = Present_Resist.asc",
                    "habitat_map_is_resistances = true",
                    "output_file = Present_CS.out")
# Write it to your working directory
writeLines(Present_CS_ini,"08_MMRR_and_CA/data/Circuitscape/Present_CS.ini")

# Then, I copied the 3 iput files generated in the directory of serafin, 
# where I run Circuitscape (/home/aferreiro/Circuitscape/Stilpnia/Present/)
# From there I initialize Julia, load Circuitscape, and compute the pairwise 
# distance matrix

## 2.2. Post-proccessing of CS outputs ------------------------------------------

# Import the effective resistance distances obtained from CS
Present_CS_distances <- read.csv("08_MMRR_and_CA/data/Circuitscape/Present_CS_resistances.out", header=TRUE, sep=" ", dec=".")
names(Present_CS_distances) <- str_remove(string = names(Present_CS_distances),
                                       pattern = "X") 
names(Present_CS_distances) <- str_remove(string = names(Present_CS_distances),
                                       pattern = "\\.0")
colnames(Present_CS_distances)[1] <- "id"
Present_CS_distances

# To see which point was associated with each sample we use
Cell_id_present <- extract(Present_Coords, Coords_spatvect)

# The CS will only generate the result for localities. Then I can use the vector created,
# by extracting cells id for each locality to generate a vector to duplicate rows and columns
Sit <- Cell_id_present$last |> as.character()

# Then, I Duplicate rows and columns basen on this Sit vector
Matrix_ResistPresent <- Present_CS_distances[match(Sit,Present_CS_distances$id), Sit] |> 
  as.matrix()

# Export matrix
write.table(Matrix_ResistPresent, file="08_MMRR_and_CA/data/Matrix_ResistPresent_CS.txt", row.names=FALSE, col.names=FALSE)
# Read matrix
Matrix_ResistPresent <- read_matrix("08_MMRR_and_CA/data/Matrix_ResistPresent_CS.txt")

# Standarization
ResistPresent <- matrix_standarization(Matrix_ResistPresent)


# 3. Historic -----------------------------------------------------------

## 3.1. Input data for Circuitscape ---------------------------------------------

### 3.1.1. Resistance raster ----------------------------------------------------
# Load historic suitability surfaces
ENM_MH <- rast("07_ENM/ENMTML/Result9/Projection/MH/Ensemble/W_MEAN/Stilpnia_cayana.tif") 
ENM_LGM <- rast("07_ENM/ENMTML/Result9/Projection/LGM/Ensemble/W_MEAN/Stilpnia_cayana.tif") 
ENM_LIG <- rast("07_ENM/ENMTML/Result9/Projection/LIG/Ensemble/W_MEAN/Stilpnia_cayana.tif") 

# Set negative values to 0
ENM_MH[ENM_MH < 0] <- 0
ENM_LGM[ENM_LGM < 0] <- 0
ENM_LIG[ENM_LIG < 0] <- 0
qtm(ENM_MH)
qtm(ENM_LGM)
qtm(ENM_LIG)

# Edit LIG and MH to associate a 0 values to the LGM coastline
# Load LGM coastline
LGM_coast <- st_read("C:/Users/ale_f/OneDrive/Documentos/Capas_GIS/Capas_Sudam/lgm_LINEA_DE_COSTA/LGMcoastline_RayAdams_2001.gpkg")
LGM_coast <- st_transform(LGM_coast,
                          crs = st_crs(ENM_LGM))

# Convert NA to -1 values
rcl_1 <- matrix(c(NA, -1),
                ncol = 2,)
rcl_2 <- matrix(c(-1, 0),
                ncol = 2,)
ENM_MH_2 <- ENM_MH |> 
  classify(rcl_1) |> # Change NA values for other
  crop(LGM_coast) |> # Crop using LGM coast line
  mask(LGM_coast) |> 
  classify(rcl_2)
qtm(ENM_MH_2)

ENM_LIG_2 <- ENM_LIG |> 
  classify(rcl_1) |> # Change NA values for other
  crop(LGM_coast) |> # Crop using LGM coast line
  mask(LGM_coast) |> 
  classify(rcl_2)
qtm(ENM_LIG_2)

# Also crop LGM raster to match extent
ENM_LGM_2 <- ENM_LGM |> 
  crop(ENM_MH_2) |> # Crop using LGM coast line
  mask(ENM_MH_2)
qtm(ENM_LGM_2)

# Average historic suitability
ENM_Historic <- terra::mean(ENM_MH_2, ENM_LGM_2, ENM_LIG_2)
qtm(ENM_Historic)

# Convert conductance to resistance
Historic_Resist <- 1-ENM_Historic 
qtm(Historic_Resist) # plot to check
# Save to use it in CS
writeRaster(Historic_Resist,
            "08_MMRR_and_CA/data/Circuitscape/Historic_Resist.asc")

### 3.1.2. Focal nodes raster ----------------------------------------------------

# Create the coords file
Historic_Coords <- terra::rasterize(Coords_spatvect, 
                                   Historic_Resist, 
                                   field = "id")

# Save it to use with CircuitScape
writeRaster(Historic_Coords,
            "08_MMRR_and_CA/data/Circuitscape/Historic_Coords.asc",
            overwrite = T)

### 3.1.3. ini file ----------------------------------------------------

# Make an .ini file
Historic_CS_ini <- c("[circuitscape options]",
                     "data_type = raster",
                     "scenario = pairwise",
                     "point_file = Historic_Coords.asc",
                     "habitat_file = Historic_Resist.asc",
                     "habitat_map_is_resistances = true",
                     "output_file = Historic_CS.out")
# SI LO CORRES EN EL CLUSTER DEBES ADAPTAR LOS PATHS DEL .ini!!!

# Write it to your working directory
writeLines(Historic_CS_ini,"08_MMRR_and_CA/data/Circuitscape/Historic_CS.ini")

# Then, I copied the 3 input files generated in the directory of serafin, 
# where I run Circuitscape (/home/aferreiro/Circuitscape/Stilpnia/Present/)
# From there I initialize Julia, load Circuitscape, and compute the pairwise 
# distance matrix


## 3.2. Post-proccessing of CS output -------------------------------------------


# Import the effective resistance distances obtained from CS
Historic_CSdistances <- read.csv("08_MMRR_and_CA/data/Circuitscape/Historic_CS_resistances.out", 
                                 header=TRUE, sep=" ", dec=".")
names(Historic_CSdistances) <- str_remove(string = names(Historic_CSdistances),
                                       pattern = "X") 
names(Historic_CSdistances) <- str_remove(string = names(Historic_CSdistances),
                                       pattern = "\\.0")
colnames(Historic_CSdistances)[1] <- "id"

# To see which point was associated with each sample we use
Cell_id_historic <- extract(Historic_Coords, Coords_spatvect)

# The CS will only generate the result for localities. Then I can use the vector created,
# by extracting cells id for each locality to generate a vector to duplicate rows and columns
Sit_historic <- Cell_id_historic$last |> as.character()


# Then, I Duplicate rows and columns basen on this Sit vector
Matrix_ResistHistoric <- Historic_CSdistances[match(Sit_historic, Historic_CSdistances$id), Sit_historic] |> 
  as.matrix()

# Export matrix
write.table(Matrix_ResistHistoric, file="08_MMRR_and_CA/data/Matrix_ResistHistoric_CS.txt", row.names=FALSE, col.names=FALSE)
# Read matrix
Matrix_ResistHistoric <- read_matrix("08_MMRR_and_CA/data/Matrix_ResistHistoric_CS.txt")

# Standarization
ResistHistoric <- matrix_standarization(Matrix_ResistHistoric)
# 4. Amazon river barrier -----------------------------------------------------

# Load pop info, and replace in order to only get 3 groups (North, South and Marajo)
Localities_df <- read.csv("00_Sampling_sites/Populations_Stilpnia83.csv", sep = ",") |> 
  mutate(pop = replace(Populations, 
                       Populations == "East" |
                         Populations == "West" |
                         Populations == "Caatinga" |
                         Populations == "Marajo", 
                       "South"))


# Generate a matrix that tells if  two sites are in the same group
Matrix_Barrier <- outer(Localities_df$pop, Localities_df$pop, `!=`) 
+Matrix_Barrier # convert to numerical matrix instead of boolean
diag(Matrix_Barrier) <- 0 # convert main diagonal values to 0.
write.table(Matrix_Barrier, 
            file = "08_MMRR_and_CA/data/Matrix_Barrier.txt", 
            row.names=FALSE, col.names=FALSE)
Matrix_Barrier <- read_matrix("08_MMRR_and_CA/data/Matrix_Barrier.txt")

# Standarization
Barrier <- matrix_standarization(Matrix_Barrier)

# 5. Genetic distance -----------------------------------------------------

# Calculated as 1-Identity by descent coefficient, in SNPRelate.
# Open the GDS file
genofile <- snpgdsOpen("08_MMRR_and_CA/data/Tangara_2023_83-II.gds")

# Generate similarity matrix 
IBS_obj <- snpgdsIBS(genofile, 
                     sample.id = Localities_df$Sample, # Saco muestras de West y Marajo
                     autosome.only = F)
IBS_matrix <- IBS_obj$ibs

# Substract by 1 to obtain a dissimilarity index
DisIBS_matrix <- 1 - IBS_matrix
write.table(DisIBS_matrix, 
            file="08_MMRR_and_CA/data/Matrix_gendist_IBS.txt", 
            sep = "\t",
            row.names=FALSE, 
            col.names=FALSE)

# Cargo matriz con distancias geneticas obtenida de SNPRelate, 
# ver script de input data de EEMS, para ver como se obtiene
Matrix_gendist_IbD <- read_matrix("04_EEMS/data/Tangara_2023_83-II_randSNP_DisIBS.txt")

DistGenetica_IBD <- matrix_standarization(Matrix_gendist_IbD)





# 6. Data inspection -------------------------------------------------
## 6.1. Multicolinearity --------------------------------------------------

# Generate a Pearson correlation matrix
Corr_mat <- data.frame(DistGenetica_IBD, DistGeog, ResistPresent, ResistHistoric, Barrier)
CorrMat_df <- cor(Corr_mat, method = "pearson") |> 
  as.data.frame()

# Estimate VIF
vif_df <- vif(lm(DistGenetica_IBD~DistGeog+ResistPresent+ResistHistoric+Barrier,data=Corr_mat))
vif_df <- c(NA, vif_df)

# Multicolinearity table
Multicol_table <- CorrMat_df |> 
  mutate(VIF = vif_df, .before = 1)

# Plot correlation between variables
CorrMat_gtTable <- Multicol_table |>
  gt(rownames_to_stub = T) |> 
  tab_header(
    title = md("Multicolinearity matrix")
  ) |> 
  tab_spanner(
    label = "Pearson correlation",
    columns = c(3:7)
  ) |> 
  data_color(
    columns = VIF,
    rows = everything(),
    method = "bin",
    palette = c("white","red"),
    bins = c(0, 10, 10E9),
    na_color = "white"
  ) |> 
  data_color(
    columns = c(DistGenetica_IBD:Barrier),
    rows = everything(),
    method = "bin",
    palette = c("white","red", "white"),
    bins = c(0 ,0.7, 0.999, 1)
  ) |> 
  cols_align(
    align = "center",
    columns = everything()
  ) |> 
  fmt_number(
    columns = everything(),
    decimals = 3
  ) |> 
  sub_missing(
    missing_text = "-"
  )

CorrMat_gtTable


# Save
gtsave(data = CorrMat_gtTable,
       filename = "08_MMRR_and_CA/plot/LRDM/Multicol_Mat.png")


## 6.2. Plot to see distribution of Y -------------------------------------------

DistGenetica_IBD_df <- DistGenetica_IBD |> 
  as.data.frame()


histogram_DistGen <- ggplot(DistGenetica_IBD_df, aes(x = DistGenetica_IBD)) +
  geom_histogram(bins = 80) +
  xlab(label = "Genetic distance") +
  ylab(label = "Frequency" ) +
  theme_classic()

ggsave("08_MMRR_and_CA/plot/00_histogram_DistGen.png",
       histogram_DistGen)

## 6.3. Residuals evaluation -----------------------------------------------------

# Test to check normality of residuals
lm_Stilpnia <- lm(DistGenetica_IBD~DistGeog+ResistPresent+ResistHistoric+Barrier,data=Corr_mat)

# Graph for detecting violation of normality assumption.

QQplot_Stilpnia <- ols_plot_resid_qq(lm_Stilpnia)

ggsave("08_MMRR_and_CA/plot/00_QQplot.png",
       QQplot_Stilpnia)

# Residual Normality Test
# Test for detecting violation of normality assumption.
norm_test_Stilpnia <- ols_test_normality(lm_Stilpnia) 

# Correlation between observed residuals and expected residuals under normality.
ols_test_correlation(lm_Stilpnia)
# Residual vs Fitted Values Plot
# It is a scatter plot of residuals on the y axis and fitted values on the x axis
# to detect non-linearity, unequal error variances, and outliers.
# Characteristics of a well behaved residual vs fitted plot:
# - The residuals spread randomly around the 0 line indicating that the relationship is linear.
# - The residuals form an approximate horizontal band around the 0 line indicating homogeneity of error variance.
# - No one residual is visibly away from the random pattern of the residuals indicating that there are no outliers.

ols_plot_resid_fit(lm_Stilpnia)
ols_plot_resid_fit_spread(lm_Stilpnia)

# Residual Histogram
# Histogram of residuals for detecting violation of normality assumption.
hist_resid_Stilpnia <- ols_plot_resid_hist(lm_Stilpnia)
ggsave("08_MMRR_and_CA/plot/00_histogram_residuals.png",
       plot = hist_resid_Stilpnia)

# # Collinearity diagnostics
# # VIF = Variance inflation factors measure the inflation in the variances of the parameter estimates due to collinearities that exist among the predictors
# # tolerance = Percent of variance in the predictor that cannot be accounted for by other predictors.
# ols_coll_diag(lm_Stilpnia)

# # Diagnostics Panel: Panel of plots for regression diagnostics
# pdf("Asthenes_Diagnostics_neutral.pdf")
# ols_plot_diagnostics(lm_Stilpnia)
# dev.off()

# 7. LRDM and CA  --------------------------------------------------------

# As we see in the section before, genetic distances showed a multimodal 
# distribution, thus violating the assumptions of normality of residuals 
# in OLS regression

# In this case it should be used nonparametric regressions. But, as 
# variance-partitioning procedures such as Commonality Analysis cannot 
# be applied # in the case of nonparametric regressions, we transformed 
# all genetic distances into binary variables and performed logistic 
# regression on distance matrices (LRDM).

# Using zero as a threshold, z-transformed genetic distances were thus recoded 
# into binary data with 0 for pairs of individuals with negative z-scores and 
# 1 for pairs of individuals with positive z-scores
DistGenetica_IBD_bin <- +(DistGenetica_IBD > "0") 

## 7.1. CA -----------------------------------------------------------------

temp2 <- data.frame(DistGenetica_IBD_bin, 
                    DistGeog, ResistPresent, ResistHistoric, Barrier)

# CA for logistic regression using function from Roberts & Nimon (2012)
CA_Stilpnia <- cc4log(temp2, # dataset containing the dependent and independent variables
                      dv = "DistGenetica_IBD_bin", # dependent variable named in the dataset
                      iv = list("DistGeog","ResistPresent","ResistHistoric","Barrier"), # list of independent variables named in the dataset
                      PR2 = "N") # 
CA_Stilpnia

## 7.2. LRDM ---------------------------------------------------------------

### 7.2.1. Basic model -------------------------------------------------------
glmt_Stilpnia <- glm(DistGenetica_IBD_bin ~ DistGeog+ResistPresent+ResistHistoric+Barrier, 
                     family=binomial("logit"), data=temp2)
glmt_Stilpnia_summ <- summary(glmt_Stilpnia)

# To evaluate overall model fit, we used the Nagelkerke’s Index as a pseudo-R2,
# with a range (from 0 to 1) identical to the range of OLS multiple R2
pseudoR_Stilpnia <- NagelkerkeR2(glmt_Stilpnia)
pseudoR_Stilpnia$R2

### 7.2.2. Plot -------------------------------------------------------
### Plot presented in Prunier 2015 
minTOT <- min(min(DistGeog),min(ResistPresent),min(ResistHistoric),min(Barrier))
maxTOT <- max(max(DistGeog),max(ResistPresent),max(ResistHistoric),max(Barrier))
Zscore=seq(minTOT,maxTOT,0.1)
p1 <- rep(mean(DistGeog),length(Zscore));
p2 <- rep(mean(ResistPresent),length(Zscore));
p3 <- rep(mean(ResistHistoric),length(Zscore));
p4 <- rep(mean(Barrier),length(Zscore));

# Plot and save simultaneusly
png("08_MMRR_and_CA/plot/LRDM/prob_success.png")
for (i in c(1,2,3,4)){
  if(i==1){
    Newdata = data.frame(DistGeog = Zscore, ResistPresent = p2,
                         ResistHistoric = p3, Barrier = p4)
    prGLM <- predict(glmt_Stilpnia,newdata=Newdata,type="link",se=TRUE)
    plot(Zscore,lg(prGLM$fit),type="l",ylim=c(0,1),col="black",lty=1,lwd=2,
         xlab="Z-scores", ylab="Probability of success")
    text(-2, 0, "DistGeo", col="black")
    par(new=T)
    }
  else if(i==2){
    Newdata = data.frame(DistGeog = p1, ResistPresent = Zscore,
                         ResistHistoric = p3, Barrier = p4)
    prGLM <- predict(glmt_Stilpnia,newdata=Newdata,type="link",se=TRUE)
    plot(Zscore,lg(prGLM$fit),type="l",ylim=c(0,1),col="darkgreen",lty=1,lwd=2,
         xlab="Z-scores", ylab="Probability of success")
    text(-2, 1, "ResPresent", col="darkgreen")
    par(new=T)
    }
  else if(i==3){
    Newdata = data.frame(DistGeog = p1, ResistPresent = p2,
                         ResistHistoric = Zscore, Barrier = p4)
    prGLM <- predict(glmt_Stilpnia,newdata=Newdata,type="link",se=TRUE)
    plot(Zscore,lg(prGLM$fit),type="l",ylim=c(0,1),col="red",lty=1,lwd=2,
         xlab="Z-scores", ylab="Probability of success")
    text(-2, 0.7, "ResHistoric", col="red")
    par(new=T)}
  else if(i==4){
    Newdata = data.frame(DistGeog = p1, ResistPresent = p2,
                         ResistHistoric = p3, Barrier = Zscore)
    prGLM <- predict(glmt_Stilpnia,newdata=Newdata,type="link",se=TRUE)
    plot(Zscore,lg(prGLM$fit),type="l",ylim=c(0,1),col="darkorchid4",lty=1,lwd=2,
         xlab="Z-scores", ylab="Probability of success")
    text(-2, 0.4, "Barrier", col="darkorchid4")
  }}	
dev.off()

### 7.2.3. Calculate beta and odd ratios ---------------------------------------
# Semi-Standardized beta weights and Odd-ratios
stdBeta=rbind(0,0,0,0)
stdBeta2=rbind(0,0,0,0)
for (i in c(2:5)){
  Newdata = data.frame(DistGeog, ResistPresent, ResistHistoric, Barrier)
  prGLM <- predict(glmt_Stilpnia,newdata=Newdata,type="response",se=TRUE)
  mprob=mean(prGLM$fit)
  mprob2=0.5
  stdBeta[i-1]=(1/(1+exp(-(log(mprob/(1-mprob))+0.5*glmt_Stilpnia$coefficients[i]))))-(1/(1+exp(-(log(mprob/(1-mprob))-0.5*glmt_Stilpnia$coefficients[i]))))
  stdBeta2[i-1]=(1/(1+exp(-(log(mprob2/(1-mprob2))+0.5*glmt_Stilpnia$coefficients[i]))))-(1/(1+exp(-(log(mprob2/(1-mprob2))-0.5*glmt_Stilpnia$coefficients[i]))))
}
betas <- c(stdBeta)
odd_ratios <- c(exp(stdBeta))

### 7.2.4. Level of significance -----------------------------------------------

source("08_MMRR_and_CA/code/CAonDM_Rstudio.R")

mydata1 <- list(Matrix_gendist_IbD,     # a list object containing all original 
                Matrix_distgeo,         # pairwise squared matrices. The dependent
                Matrix_ResistPresent,   # matrix (Y) is to be the first attribute
                Matrix_ResistHistoric,
                Matrix_Barrier)
names(mydata1) <- c("DistGen", "DistGeo", "ResistPresent", "ResistHistoric", "Barrier")

calrdm_Stilpnia <- CAlrdm(mydata =  mydata1,
                          regrnperm = 10000, # number of permutations for tests of significance
                          bootn = 1000, # number of bootstrap iterations for assessing parameters' robustness to the removal of a random subset of objects
                          bootprop = (1/83)) # proportion of sampled objects (populations or individuals) to be randomly removed at each bootstrap iteration

# Level of significance
calrdm_Stilpnia$modelfit

# 8. CA Boostrap replicates ----------------------------------------------------

## 8.1. Bootstrap to calculate commonalities by Variable -----------------------

# Calculation of bootstrap estimates to obtain CI (Peterman et al. 2014)
# Number of boostrap replicates
nperm <- 1000
# Number of predictors
n_predictors <- 4
# Create empty matrix that will contain bootstrapped values for parameters
Unique_bs <- matrix(data = 0, nrow = nperm, ncol = n_predictors)
Common_bs <- matrix(data = 0, nrow = nperm, ncol = n_predictors)
Total_bs <- matrix(data = 0, nrow = nperm, ncol = n_predictors)
# Create an empty data frame to use it as final data frame with with estimates 
# and CI for each pred
Stilpnia_resBootstrap1 <- data.frame(n = rep(0,n_predictors),
                                     Unique = rep(0,n_predictors),
                                     Unique_CI_low = rep(0,n_predictors), 
                                     Unique_CI_upper = rep(0,n_predictors),
                                     Common = rep(0,n_predictors),
                                     Common_CI_low = rep(0,n_predictors), 
                                     Common_CI_upper = rep(0,n_predictors),
                                     Total = rep(0,n_predictors),
                                     Total_CI_low = rep(0,n_predictors), 
                                     Total_CI_upper = rep(0,n_predictors))
# Number of genetic samples
n <- ncol(Matrix_gendist_IbD)
# Number of samples to use in each bootstrap replicate
sn <- n-1 

# Loop to create 1000 bootstrapped CA 
for (i in 1:nperm){ # on each bootstrap:
  rarray = sort(sample(n,sn,replace=F)) # crate a vector of the samples to include in the bs sampling
  mmDG = matrix_standarization(Matrix_gendist_IbD[rarray,rarray]) # subsampling genetic data
  mmGEO = matrix_standarization(Matrix_distgeo[rarray,rarray]) # subsampling predictor
  mmCUR = matrix_standarization(Matrix_ResistPresent[rarray,rarray]) # subsampling predictor
  mmPAST = matrix_standarization(Matrix_ResistHistoric[rarray,rarray]) # subsampling predictor
  mmRIV = matrix_standarization(Matrix_Barrier[rarray,rarray]) # subsampling predictor
  
  mmDG <- +(mmDG>"0") # binarization of genetic distances
  prout=data.frame(mmDG,mmGEO, mmCUR, mmPAST, mmRIV)
  COMM_Stilpnia=cc4log(prout,"mmDG",list("mmGEO", "mmCUR", "mmPAST", "mmRIV"), "N") 
  
  # Fill the row correspondent in the boot file with the results obtained
  Unique_bs[i,] = COMM_Stilpnia$CCTotalbyVar[c(1:n_predictors),1]
  Common_bs[i,] = COMM_Stilpnia$CCTotalbyVar[c(1:n_predictors),2]
  Total_bs[i,] = COMM_Stilpnia$CCTotalbyVar[c(1:n_predictors),3]
}
# After this loop i will obtain a boot file with 1000 estimations of commonalities
# and I can use it to estimate CI
for (i in 1:n_predictors){
  q_Unique = quantile(Unique_bs[,i], c(.025,.975)) # Calc CI limits
  q_Common = quantile(Common_bs[,i], c(.025,.975))
  q_Total = quantile(Total_bs[,i], c(.025,.975))
  Stilpnia_resBootstrap1[i,1]=i # Write CI limits in the df created previously
  Stilpnia_resBootstrap1[i,3]=q_Unique[1]
  Stilpnia_resBootstrap1[i,4]=q_Unique[2]
  Stilpnia_resBootstrap1[i,6]=q_Common[1]
  Stilpnia_resBootstrap1[i,7]=q_Common[2]
  Stilpnia_resBootstrap1[i,9]=q_Total[1]
  Stilpnia_resBootstrap1[i,10]=q_Total[2]
}

# Add to the data point estimates of CA parameters  
Stilpnia_resBootstrap1[,2] = CA_Stilpnia$CCTotalbyVar[c(1:n_predictors),1]
Stilpnia_resBootstrap1[,5] = CA_Stilpnia$CCTotalbyVar[c(1:n_predictors),2]
Stilpnia_resBootstrap1[,8] = CA_Stilpnia$CCTotalbyVar[c(1:n_predictors),3]

# Set the rownames for each predictor
rownames(Stilpnia_resBootstrap1) = rownames(CA_Stilpnia$CCTotalbyVar)[1:n_predictors]
Stilpnia_resBootstrap1

# Add bootstrap results to the final table, AND PUT ALL VALUES IN THE SAME COLUMN
Results_table <- Stilpnia_resBootstrap1 |> 
  dplyr::select(!1) |> # remove n column
  mutate(Beta = betas, .before = 1) |> 
  mutate(Odd_ratio = odd_ratios, .before = 2) |> 
  mutate(P = glmt_Stilpnia_summ$coefficients[2:5,4], .before = 3) |> 
  mutate(across(where(is.numeric), round, 3)) |> # redondeo a dos decimales
  unite(col = "Unique_CI", Unique_CI_low:Unique_CI_upper, sep = ";") |> # Unite CI limits columns
  unite(col = "Common_CI", Common_CI_low:Common_CI_upper, sep = ";") |> 
  unite(col = "Total_CI", Total_CI_low:Total_CI_upper, sep = ";") |>
  mutate(Unique_CI = paste0(Unique_CI, ")"),
         Common_CI = paste0(Common_CI, ")"),
         Total_CI = paste0(Total_CI, ")")) |> 
  unite(col = "Unique", Unique:Unique_CI, sep = " (") |> 
  unite(col = "Common", Common:Common_CI, sep = " (") |> 
  unite(col = "Total", Total:Total_CI, sep = " (")
  
  
Results_table

# 9. Fancy table -------------------------------------------------------------

# Generate a fancy table
CA_gt_table <- Results_table |> 
  gt(rownames_to_stub = T) |> 
  tab_header(
    title = md("LRDM and CA of _Stilpnia cayana_"),
    subtitle = md("Pseudo-R^2^ = **0.5765** ; p-value: **0.000**")
  ) |>
  opt_align_table_header(align = "center") |>
  cols_label(
    Beta = md("β"),
    Odd_ratio = md("Ψ"),
    P = md("p-value"),
    Unique = md("_U_"),
    Common = md("_C_"),
    Total = md("_T_")
  ) |>
  fmt_number(
    columns = everything(),
    decimals = 3
    ) |>
  # tab_stubhead(
  #   label = md("Predictors")
  #   ) |> 
  # data_color(
  #   columns = Beta_weights,
  #   palette = "RdBu",
  #   domain = c(-1,1)
  #   ) |>
  cols_width(
    Beta ~ px(100),
    Odd_ratio ~ px(100),
    P ~ px(100),
    Unique ~ px(100),
    Common ~ px(100),
    Total ~ px(100)
  ) |> 
  cols_align(
    align = "center",
    columns = everything()
  )
CA_gt_table

# Save
gtsave(data = CA_gt_table,
       filename = "08_MMRR_and_CA/plot/LRDM/CA_table.png")

