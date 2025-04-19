# Paquetes ----------------------------------------------------------------

library(sf)
library(ENMTML)
library(tidyverse)
library(raster)
library(terra)
library(tmap)
library(tmaptools)

# Mapa base para chequeo de capas -----------------------------------------
# Cargo limites políticos de Sudamerica
library(rnaturalearth)
library(rnaturalearthdata)
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




# 1. Directories and data creation ---------------------------------------

# First will be created a folder with a working directory
dir.create("ENMTML")


## 1.1. Ocurrence file -----------------------------------------------------

# Leer csv con occs de tesis de Samira
Occs_Tesis <- read.csv(file = "data/Occs_Tangara_TesisSamira.csv",
                       sep = ";")
Occs_raw <- Occs_Tesis |> 
  dplyr::select(Lat, Long) |> 
  dplyr::mutate(sp = rep("Stilpnia_cayana", 263),
         .before = Lat)

Occs_raw_sf <- st_as_sf(Occs_raw,
                        coords = c("Long","Lat"),
                        crs = 4326)
# Will be saved some ENMTML data sets to ENMTML_example folder
utils::write.table(Occs_raw, "ENMTML/occ_stilpnia.txt", sep = '\t', row.names = FALSE)



## 1.2. Predictors ---------------------------------------------------------

paths_CurrentBios <- list.files(path = "C:/Users/ale_f/OneDrive/Documentos/Capas_GIS/ClimVars/PaleoClim/1_Current/2_5min",
                                pattern = ".tif$",full.names = TRUE)
CurrentBios <- rast(paths_CurrentBios) |> 
  subset(c(1:9,12:17)) |> 
  crop(Sudam) |> # Corto por area de estudio
  mask(Sudam)
qtm(CurrentBios)
# Guardo capas con variables predictoras
dir.create("ENMTML/PredVars")
lapply(names(CurrentBios), 
       function(x){
         writeRaster(CurrentBios[[x]], 
                     paste0("ENMTML/PredVars/",x,".tif"),
                     overwrite=TRUE,
                     NAflag = -9999)})

## 1.3. Projection variables ----------------------------------------------------

### 1.3.1. Present ----------------------------------------------------------



# Guardo variables para proyectar
dir.create("ENMTML/ProyVars")
dir.create("ENMTML/ProyVars/Present")
lapply(names(CurrentBios), 
       function(x){
         writeRaster(CurrentBios[[x]], 
                     paste0("ENMTML/ProyVars/Present/",x,".tif"),
                     overwrite=TRUE, NAflag = -9999)})

### 1.3.2. MH -----------------------------------------------------
paths_MHBios <- list.files(path = "C:/Users/ale_f/OneDrive/Documentos/Capas_GIS/ClimVars/PaleoClim/2_MH",
                           pattern = ".tif$",full.names = TRUE)
# Las cargo en R, selecciono las bios a usar y corto por el area de estudio
MHBios <- rast(paths_MHBios) |>
  subset(c(1:9,12:17)) |> # selecciono bios a usar
  crop(Sudam) |> # Corto por area de estudio
  mask(Sudam)
# Creo carpeta con variables de proyección para kuenm
dir.create("ENMTML/ProyVars/MH")
lapply(names(MHBios), 
       function(x){
         writeRaster(MHBios[[x]], 
                     paste0("ENMTML/ProyVars/MH/",x,".tif"),
                     overwrite=TRUE, NAflag = -9999)})


### 1.3.3. LGM --------------------------------------------------------------
paths_LGMBios <- list.files(path = "C:/Users/ale_f/OneDrive/Documentos/Capas_GIS/ClimVars/PaleoClim/3_LGM/2_5min",
                            pattern = ".tif$",full.names = TRUE)
# Las cargo en R, selecciono las bios a usar y corto por el area de estudio
LGMBios <- rast(paths_LGMBios) |>
  subset(c(1:9,12:17)) |> # selecciono bios a usar
  crop(Sudam) |> # Corto por area de estudio
  mask(Sudam)
# Creo carpeta con variables de proyección 
dir.create("ENMTML/ProyVars/LGM")
lapply(names(LGMBios), 
       function(x){
         writeRaster(LGMBios[[x]], 
                     paste0("ENMTML/ProyVars/LGM/",x,".tif"),
                     overwrite=TRUE, NAflag = -9999)})


### 1.3.4. LIG --------------------------------------------------------------
paths_LIGBios <- list.files(path = "C:/Users/ale_f/OneDrive/Documentos/Capas_GIS/ClimVars/PaleoClim/4_LIG",
                            pattern = ".tif$",full.names = TRUE)
# Las cargo en R, selecciono las bios a usar y corto por el area de estudio
LIGBios <- rast(paths_LIGBios) |>
  subset(c(1:9,12:17)) |> # selecciono bios a usar
  crop(Sudam) |> # Corto por area de estudio
  mask(Sudam)
# Creo carpeta con variables de proyección para kuenm
dir.create("ENMTML/ProyVars/LIG")
lapply(names(LIGBios), 
       function(x){
         writeRaster(LIGBios[[x]], 
                     paste0("ENMTML/ProyVars/LIG/",x,".tif"),
                     overwrite=TRUE, NAflag = -9999)})

## 1.4. Area M -------------------------------------------------------------
M_Stilpnia <-
  flexsdm::calib_area(
    data = Occs_raw,
    x = 'Long',
    y = 'Lat',
    method =  c('bmcp', width = 300000),
    crs = crs(CurrentBios)
  ) # create a calibration area with 300 km buffer around MCP
M_Stilpnia_rast <- rasterize(M_Stilpnia, CurrentBios)
qtm(M_Stilpnia_rast)
writeRaster(M_Stilpnia_rast, 
            "ENMTML/M_Stilpnia.tif")

M_Stilpnia_sf <- st_as_sf(M_Stilpnia)

mapa2 <- mapa_base +
  tm_shape(M_Stilpnia_sf) + tm_borders(col = "green", lwd = 2)
mapa2
# Gurado como shape
st_write(M_Stilpnia_sf, "ENMTML/M_Sptilpnia.gpkg")
M_Stilpnia_sf <- st_read("07_ENM/ENMTML/M_Sptilpnia.gpkg")

# 2. Construction ENM with ENMTML ----------------------------------------

# ENMTML provides a variety of tools to build different models
# depending on the modeling objectives.
# Here will be provided a single modeling procedure.
# For more example and exploration of models
# see <https://github.com/andrefaa/ENMTML>

# Will be fitted models for five virtual species with
# current and future conditions. Please read ENMTML arguments.

# The next object contains the directory and file path data and folders that will be used
d_occ <- "ENMTML/occ_stilpnia.txt" # file path with species occurrences
d_env <- "ENMTML/PredVars/"        # directory path with current environmental conditions (raster in tiff format)
d_pro <- "ENMTML/ProyVars/"        # directory path with folders with future environmental conditions (raster in tiff format)
d_arm <- "ENMTML/M_Stilpnia.tif"   # file path with shapefile used to constrain models
d_res <- "ENMTML/Result7"          # Directory where results will be written

ENMTML(
  pred_dir = d_env,
  proj_dir = d_pro,
  result_dir = d_res,
  occ_file = d_occ,
  sp = 'sp',
  x = 'Long',
  y = 'Lat',
  min_occ = 10,
  thin_occ = c(method='USER-DEFINED', distance='25'),
  eval_occ = NULL,
  colin_var = c(method='VIF'),
  imp_var = FALSE,
  sp_accessible_area = c(method='MASK', filepath=d_arm),
  pseudoabs_method = c(method='GEO_CONST', width='150'),
  pres_abs_ratio = 1,
  part=c(method= 'BLOCK'),
  save_part = T,
  save_final = T,
  algorithm = c("BRT","MXD","SVM","RDF","GAU"),
  thr = c(type='MAX_TSS'),
  msdm = NULL,
  ensemble = c(method='W_MEAN', metric='TSS'),
  extrapolation = T,
  cores = 2
)

# ENMTML function will create a folder named Result a directory
# prior to the directory specified in the pred_dir argument


# 3. Check partition and extent mask ---------------------------------------------------

ExtentMask <- rast("07_ENM/ENMTML/Result7/Extent_Masks/Stilpnia_cayana.tif")
qtm(ExtentMask)

# Load raster and occs partitioned

SpatPart <- rast("07_ENM/ENMTML/Result7/BLOCK/Stilpnia_cayana.tif")
OccsPart <- read_tsv("07_ENM/ENMTML/Result7/BLOCK/OccBlocks.txt") |> 
  st_as_sf(coords = c("x","y"),
           crs = 4326)



# Plot to see partitio in space
Partition_map <- mapa_base +
  tm_shape(SpatPart) + tm_raster(palette = get_brewer_pal("Accent", n = 4, plot = F)) +
  tm_shape(OccsPart) + tm_symbols(size = 0.2, 
                                  col = "Partition",
                                  shape = "PresAbse",
                                  palette = get_brewer_pal("Accent", n = 4, plot = F)) +
  tm_shape(Sudam, bbox = bb) + tm_borders(lwd = 2) +
  tm_shape(M_Stilpnia_sf, bbox = bb) + tm_borders(lwd = 2, col = "red") +
  tm_layout(legend.show = F)

Partition_map

tmap_save(Partition_map, "ENMTML/Result6/Plots/Partition_map.pdf")


# 4. Check evaluation metrics ----------------------------------------------

Eval_metrics <- read.table("ENMTML/Result7/Evaluation_Table.txt")

# 5. Plot results ---------------------------------------------------------

## 5.1. Present ------------------------------------------------------------
Present_rast <- rast("ENMTML/Result7/Projection/Present/Ensemble/W_MEAN/MAX_TSS/Stilpnia_cayana.tif")

Present_map <-  mapa_base +
  tm_shape(Present_rast) + tm_raster(style = "cont",
                                       palette = viridisLite::viridis(20, direction = -1),
                                       title= "Suitability",
                                       legend.format =list(text.separator="-"),
                                       legend.reverse = T) +
  tm_shape(Occs_raw_sf) + tm_symbols(size = 0.1, col = "red", alpha = 0.8) +
  tm_shape(Sudam, bbox = bb) + tm_borders(lwd = 2) +
  tm_layout(title = "Present",
            title.position = c("right", "top"),
            legend.position = c("right", "bottom"))
Present_map

tmap_save(Present_map, "ENMTML/Result7/Plots/SDM_Present_GAU.pdf")
  
## 5.2. MH -----------------------------------------------------------------
MH_rast <- rast("ENMTML/Result7/Projection/MH/Ensemble/W_MEAN/Stilpnia_cayana.tif")

MH_map <-  tm_shape(Sudam, bbox = bb) + tm_polygons(col = "lightgrey") +
  tm_shape(MH_rast) + tm_raster(style = "cont",
                                    palette = viridisLite::viridis(20, direction = -1),
                                    title= "Suitability",
                                    legend.format =list(text.separator="-"),
                                    legend.reverse = T) +
  tm_shape(Sudam, bbox = bb) + tm_borders(lwd = 2) +
  tm_layout(title = "Mid-Holocene",
            title.position = c("right", "top"),
            legend.show = F)
MH_map

tmap_save(ens_MH_map, "ENMTML/Result2/Plots/SDM_MH_ens.pdf")


## 5.3. LGM ----------------------------------------------------------------

LGM_rast <- rast("ENMTML/Result7/Projection/LGM/Ensemble/W_MEAN/MAX_TSS/Stilpnia_cayana.tif")

LGM_map <-  tm_shape(Sudam, bbox = bb) + tm_polygons(col = "lightgrey") +
  tm_shape(LGM_rast) + tm_raster(style = "cont",
                               palette = viridisLite::viridis(20, direction = -1),
                               title= "Suitability",
                               legend.format =list(text.separator="-"),
                               legend.reverse = T) +
  tm_shape(Sudam, bbox = bb) + tm_borders(lwd = 2) +
  tm_layout(title = "Last Glacial Maximum",
            title.position = c("right", "top"),
            legend.show = F)
LGM_map

tmap_save(ens_LGM_map, "ENMTML/Result2/Plots/SDM_LGM_ens.pdf")


## 5.4. LIG -----------------------------------------------------------------

LIG_rast <- rast("ENMTML/Result7/Projection/LIG/Ensemble/W_MEAN/Stilpnia_cayana.tif")

LIG_map <-  tm_shape(Sudam, bbox = bb) + tm_polygons(col = "lightgrey") +
  tm_shape(LIG_rast) + tm_raster(style = "cont",
                                palette = viridisLite::viridis(20, direction = -1),
                                title= "Suitability",
                                legend.format =list(text.separator="-"),
                                legend.reverse = T) +
  tm_shape(Sudam, bbox = bb) + tm_borders(lwd = 2) +
  tm_layout(title = "Last Interglacial",
            title.position = c("right", "top"),
            legend.show = F)
LIG_map

tmap_save(ens_LIG_map, "ENMTML/Result2/Plots/SDM_LIG_ens.pdf")

# 5.5. Map arrange --------------------------------------------------------

PaleoENM_map <- tmap_arrange(LIG_map, LGM_map, MH_map, Present_map,
                             ncol = 2,
                             nrow = 2,
                             asp = NULL)
PaleoENM_map

tmap_save(PaleoENM_map, "ENMTML/Result7/Plots/PaleoENM_ENMTML.pdf")


# 6. Check extrapolation (MOP) --------------------------------------------

# 0 values indicate higher extrapolation


## 6.1. Present ------------------------------------------------------------

MOP_Present <- rast("ENMTML/Result/Projection/Present/Extrapolation/Stilpnia_cayana_MOP.tif")

MOP_Present_map <-  tm_shape(Sudam, bbox = bb) + tm_polygons(col = "lightgrey") +
  tm_shape(MOP) + tm_raster(style = "cont",
                                palette = get_brewer_pal("RdBu", n = 5, plot = F),
                                title= "MOP",
                                legend.format =list(text.separator="-"),
                                legend.reverse = T) +
  tm_shape(Occs_raw_sf) + tm_symbols(size = 0.1, col = "red", alpha = 0.8) +
  mapa_base +
  tm_layout(title = "Present",
            title.position = c("left", "top"),
            legend.position = c("right", "bottom"))
MOP_Present_map

tmap_save(MOP_Present_map, "ENMTML/Result/ENMTML_MOP_Present.pdf")

## 6.2. MH ------------------------------------------------------------

MOP_MH <- rast("ENMTML/Result/Projection/MH/Extrapolation/Stilpnia_cayana_MOP.tif")

MOP_MH_map <-  tm_shape(Sudam, bbox = bb) + tm_polygons(col = "lightgrey") +
  tm_shape(MOP) + tm_raster(style = "cont",
                            palette = get_brewer_pal("RdBu", n = 5, plot = F),
                            title= "MOP",
                            legend.format =list(text.separator="-"),
                            legend.reverse = T) +
  tm_shape(Occs_raw_sf) + tm_symbols(size = 0.1, col = "red", alpha = 0.8) +
  mapa_base +
  tm_layout(title = "MH",
            title.position = c("left", "top"),
            legend.position = c("right", "bottom"))
MOP_MH_map

tmap_save(MOP_MH_map, "ENMTML/Result/ENMTML_MOP_MH.pdf")

## 6.3. LGM ------------------------------------------------------------

MOP_LGM <- rast("ENMTML/Result/Projection/LGM/Extrapolation/Stilpnia_cayana_MOP.tif")

MOP_LGM_map <-  tm_shape(Sudam, bbox = bb) + tm_polygons(col = "lightgrey") +
  tm_shape(MOP) + tm_raster(style = "cont",
                            palette = get_brewer_pal("RdBu", n = 5, plot = F),
                            title= "MOP",
                            legend.format =list(text.separator="-"),
                            legend.reverse = T) +
  tm_shape(Occs_raw_sf) + tm_symbols(size = 0.1, col = "red", alpha = 0.8) +
  mapa_base +
  tm_layout(title = "LGM",
            title.position = c("left", "top"),
            legend.position = c("right", "bottom"))
MOP_LGM_map

tmap_save(MOP_LGM_map, "ENMTML/Result/ENMTML_MOP_LGM.pdf")

## 6.4. LIG ------------------------------------------------------------

MOP_LIG <- rast("ENMTML/Result/Projection/LIG/Extrapolation/Stilpnia_cayana_MOP.tif")

MOP_LIG_map <-  tm_shape(Sudam, bbox = bb) + tm_polygons(col = "lightgrey") +
  tm_shape(MOP) + tm_raster(style = "cont",
                            palette = get_brewer_pal("RdBu", n = 5, plot = F),
                            title= "MOP",
                            legend.format =list(text.separator="-"),
                            legend.reverse = T) +
  tm_shape(Occs_raw_sf) + tm_symbols(size = 0.1, col = "red", alpha = 0.8) +
  mapa_base +
  tm_layout(title = "LIG",
            title.position = c("left", "top"),
            legend.position = c("right", "bottom"))
MOP_LIG_map

tmap_save(MOP_LIG_map, "ENMTML/Result/ENMTML_MOP_LIG.pdf")

