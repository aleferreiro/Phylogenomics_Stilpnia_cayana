# Script to generate input data for ENM of Stilpnia cayana with 
# flexsdm, and ENMTML

# Packages ----------------------------------------------------------------
library(sf)
library(tidyverse)
library(terra)
library(tmap)
library(flexsdm) 

# Basemap for layers check -----------------------------------------
# Load South American country limits
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


# Basemap
mapa_base <- tm_shape(Sudam, bbox = c(-83, # minimum longitude 
                                      -30, # minimum latitude
                                      -33, # maximum longitude
                                      12   # maximum latitude
                                      )) + tm_borders(lwd = 2) +
  tmap_options(check.and.fix = TRUE)
mapa_base

# 1. Occurrence data ------------------------------------------------------

## 1.1. Duarte (2022) database -------------------------------------------------- 
# Load csv with occs of Duarte database
Occs_Tesis <- read.csv(file = "data/Occs_Tangara_TesisSamira.csv",
                       sep = ";")
Occs_Duarte <- Occs_Tesis |> 
  select(Lat, Long) |> 
  mutate(sp = rep("Stilpnia_cayana", 263),
         .before = Lat)

Occs_Duarte_sf <- st_as_sf(Occs_raw,
                        coords = c("Long","Lat"),
                        crs = 4326)

## 1.2. GBIF occurrences ---------------------------------------------------

# We used the rgbif package to download data fron GBIF database
library(rgbif)
# To reproducibility
set.seed(1234)
# Download presence records of Stilpnia cayana using RGBIF
name_backbone("Stilpnia cayana")
gbif_download <- occ_download(pred("taxonKey", 10895981),
                              pred("hasGeospatialIssue", FALSE),
                              pred("hasCoordinate", TRUE),
                              pred("occurrenceStatus","PRESENT"), 
                              format = "SIMPLE_CSV")
occ_download_wait(gbif_download) # check if download is finished
gbif_Occs_raw <- occ_download_get(gbif_download) %>%
  occ_download_import(path = "07_ENM/data/")
# Get doi for citation of the dataset
gbif_citation("0117483-240626123714530")

# Select column of interest
gbif_Occs <- gbif_Occs_raw |> select(species, decimalLongitude, decimalLatitude, 
                                     countryCode, stateProvince, locality, 
                                     year, gbifID) |> # selecciono columnas relevantes
  rename(longitude = decimalLongitude,  latitude = decimalLatitude) |> 
  distinct(longitude, latitude, .keep_all = TRUE) |> # elimino duplicados
  filter(latitude < 15,
         year > 1970)     # remove outliers records


# Convert df to sf
gbif_Occs_sf <- st_as_sf(gbif_Occs,
                         coords = c("longitude", "latitude"),
                         crs = st_crs(4326))

# Map presence records
mapa_base + 
  tm_shape(gbif_Occs_sf) + tm_symbols(size = 0.3, col = "blue")


st_write(gbif_Occs_sf,
         "07_ENM/data/Occs_gbif.gpkg")

write.csv()

# Reduce to species, lat, long data
Occs_gbif <- gbif_Occs |> 
  select(1:3)

## 1.3. GBIF + Duarte dataset -------------------------------------------
Occs_Combined <- rbind(Occs_Duarte, Occs_gbif) 

## 1.3. Accesible area -------------------------------------------------------------
# Como area M o de estudio voy a usar un poligono generado con un buffer de 
# 100km alrededor del MCP que rodea los registros de presencia

# Convierto a crs proyectado para trabajar en metros
Occs_t_sf_proj <- st_transform(Occs_t_sf,
                               crs = 31972)

# Crear poligono mínimo convexo (en inglés, MCP)
MCP_proj <- st_convex_hull(st_union(Occs_t_sf_proj)) 

# Crear poligono con buffer de 300km del PMC que representa el Area M
area_m_proj <- st_buffer(MCP_proj, 
                         dist = 300000) 

area_m <- st_transform(area_m_proj, crs = 4326)

# Recorto usando el shape del LGM
# Linea de costa del LGM
LGMcoast <- st_read("C:/Users/ale_f/OneDrive/Documentos/Capas_GIS/Capas_Sudam/lgm_LINEA_DE_COSTA/LGMcoastline_RayAdams_2001.gpkg")

M_Stilpnia <- st_intersection(LGMcoast, area_m)
M_Stilpnia <- st_sf(M_Stilpnia)

# Ploteo el mapa
mapa_area_m <- mapa_occs_t +
  tm_shape(LGMcoast) + tm_borders(lty = "dotted", lwd = 2) +
  tm_shape(M_Stilpnia) + tm_borders(col = "green", lwd = 2)
mapa_area_m

#Guardo el vector con el area M
st_write(M_Stilpnia,
         "data/M_Spiltnia.gpkg")

M_Stilpnia <- st_read("data/M_Spiltnia.gpkg")

# Area M obtenida usando flexsdm
M_Stilpnia_flex_raw <-  calib_area(
  data = Occs_t,
  x = 'Longitude',
  y = 'Latitude',
  method =  c('bmcp', width = 300000),
  crs = crs(CurrentBios)
) # create a calibration area with 300 km buffer around MCP

# Recorto usando el shape del LGM para quitar areas marinas
M_Stilpnia_flex_sf_uncutted <- st_as_sf(M_Stilpnia_flex_raw)  
M_Stilpnia_flex_sf <- st_intersection(LGMcoast, M_Stilpnia_flex_sf_uncutted)
M_Stilpnia_flex_sf <- st_sf(M_Stilpnia_flex_sf)

mapa_area_m_flex <- mapa_occs_t +
  tm_shape(LGMcoast) + tm_borders(lty = "dotted", 
                                  lwd = 2) +
  tm_shape(M_Stilpnia_flex_sf) + tm_borders(col = "green", 
                                            lwd = 2)
mapa_area_m_flex


## 1.1. Filtrado de ocurrencias -------------------------------------
### 1.1.1. Con spThin -----------------------------------------------
# Quito sitios de presencia que estan más proximos a los 10km
library(spThin)
set.seed(1234)
Occs_thin <- spThin::thin(Occs_raw,
                          lat.col = "Lat",
                          long.col = "Long",
                          spec.col = "sp",
                          thin.par = 10, reps = 1,
                          locs.thinned.list.return = T, 
                          write.files = F)
Occs_thinned <- Occs_thin[[1]]
specie <- rep("Stilpnia_cayana", n = 263)
Occs_t <- cbind(specie,Occs_thinned)
Occs_t_sf <- st_as_sf(Occs_t,
                      coords = c("Longitude","Latitude"),
                      crs = 4326)

# Ploteo para ver donde estan las presencias que quedan luego del filtrado
mapa_occs_t <- mapa_base +
  tm_shape(Occs_t_sf) + tm_symbols(size = 0.2, col = "black")

write.csv(Occs_t,
          "kuenm/Occs_t.csv")

### 1.1.2. Con flexsdm -------------------------------------------------
Occs_filtered <- Occs_raw |> 
  occfilt_env(
    data = .,
    x = "x",
    y = "y",
    id = "id",
    nbins = 8,
    env_layer = somevar
  ) %>%
  left_join(abies_p, by = c("id", "x", "y"))

plot(layer1, col="gray80", legend=FALSE, axes=FALSE)
plot(crop(ca, layer1), add=TRUE)
points(abies_p[,c("x", "y")], col = "#00000480")
points(abies_pf[,c("x", "y")], col = "#5DC86180")


# 2. Variables para modelar el nicho -------------------------------------------------

# Utilizamos variables de la base de datos de PaleoClim (CHELSA 1979-2013; Karger et al, 2017) 
### 2.1. Carga de variables del presente -------------------------------------------------------- 
# Estas variables serán utilizadas con el objetivo de generar los modelos de nicho
# Genero vector con direccion de cada capa bioclimatica
paths_CurrentBios <- list.files(path = "C:/Users/ale_f/OneDrive/Documentos/Capas_GIS/ClimVars/PaleoClim/1_Current/2_5min",
                                pattern = ".tif$",full.names = TRUE)
# Las cargo en R y las corto por el area de estudio
CurrentBios_all <- rast(paths_CurrentBios) |>
  crop(M_Stilpnia) |> # Corto por area de estudio
  mask(M_Stilpnia)

## 2.2. Seleccion de variables ---------------------------------------------

# Selecciono variables correlacionadas mediante Variance Inflation Factor

VIF <- correct_colinvar(CurrentBios_all, 
                                 method = c('vif', th = '10'), 
                                 proj = NULL)
write.csv(VIF$vif_table,
          "data/VIF_values.csv")

# Variables seleccionadas
CurrentBios <- VIF$env_layer

dir.create("kuenm")
dir.create("kuenm/M_variables")
dir.create("kuenm/M_variables/Set_1")
# Guardo capas de calibracion
lapply(names(CurrentBios), 
       function(x){
         writeRaster(CurrentBios[[x]], 
                     paste0("kuenm/M_variables/Set_1/",x,".asc"),
                     overwrite=TRUE,
                     NAflag = -9999)})
# Guardo capas de proyeccion
dir.create("kuenm/G_variables")
dir.create("kuenm/G_variables/Set_1")
dir.create("kuenm/G_variables/Set_1/Present")
lapply(names(CurrentBios), 
       function(x){
         writeRaster(CurrentBios[[x]], 
                     paste0("kuenm/G_variables/Set_1/Present/",x,".asc"),
                     overwrite=TRUE, NAflag = -9999)})

## 2.3. Partición de los datos ---------------------------------------------

# Para poder evaluar los modelos necesito datos independientes que me permitan
# evaluar que tan bien funcionan. Como no tengo fuentes independientes voy
# a particionar los datos de modo que calibre los modelos con una parte y con
# la otra los evalúo. Usaré la función part_sblock() delpaquete flexsdm para
# particionar en bloques espaciales los datos.

# Previamente debo agregar una columna que identifique las presencias
Occs_pres <- Occs_t |> 
  mutate(pr_ab = rep("1", 236))

# Corro la función para seleccionar la mejor partición
set.seed(10)

occ_part <- part_sblock(
    data = Occs_pres,
    env_layer = CurrentBios,
    pr_ab = "pr_ab",
    x = "Longitude",
    y = "Latitude",
    n_part = 4,
    min_res_mult = 3,
    max_res_mult = 200,
    num_grids = 30,
    min_occ = 50,
    prop = 1
  )
Occs_p <- occ_part$part

# Transform best block partition to a raster layer with same resolution and extent than
# predictor variables
block_layer <- get_block(env_layer = CurrentBios, best_grid = occ_part$grid)

# Plot best spatial blocking partitions
cl <- c("#64146D", "#9E2962", "#F47C15", "#FCFFA4")
plot(block_layer, col=cl, legend=FALSE, axes=FALSE)
points(Occs_p[,c("x", "y")])

# Number of presences per block
Occs_p %>%
  dplyr::group_by(.part) %>%
  dplyr::count()


## 2.4. Muestreo de pseudo-ausencias y puntos de background ----------------

# Spatial blocks where species occurs
# Sample background points throughout study area with random method, 
# allocating 10X the number of presences a background
set.seed(10)
bg <- lapply(1:4, function(x) {
  sample_background(
    data = Occs_p,
    x = "x",
    y = "y",
    n = sum(Occs_p$.part == x) * 10,
    method = "random",
    rlayer = block_layer,
    maskval = x,
    calibarea = M_Stilpnia
  )
}) %>%
  bind_rows()
bg <- sdm_extract(data = bg, x = "x", y = "y", env_layer = block_layer)

# Sample a number of pseudo-absences equal to the presence in each partition
set.seed(10)
psa <- lapply(1:4, function(x) {
  sample_pseudoabs(
    data = Occs_p,
    x = "x",
    y = "y",
    n = sum(Occs_p$.part == x),
    method = "random",
    rlayer = block_layer,
    maskval = x,
    calibarea = M_Stilpnia
  )
}) %>%
  bind_rows()
psa <- sdm_extract(data = psa, x = "x", y = "y", env_layer = block_layer)

cl <- c("#280B50", "#9E2962", "#F47C15", "#FCFFA4")
plot(block_layer, col="gray80", legend=FALSE, axes=FALSE)
points(bg[,c("x", "y")], col=cl[bg$.part], cex=0.8) # Background points
points(psa[,c("x", "y")], bg=cl[psa$.part], cex=0.8, pch=21) # Pseudo-absences


class(Occs_p$pr_ab) = "numeric"
class(psa$pr_ab)
# Bind a presences and pseudo-absences
Occs_pa <- bind_rows(Occs_p, psa)
Occs_pa # Presence-Pseudo-absence database
bg # Background points

# Extract environmental values for presence and pseudo-absence/background
Occs_pa <- Occs_pa %>%
  sdm_extract(
    data = .,
    x = "x",
    y = "y",
    env_layer = CurrentBios,
    filter_na = TRUE
  )
bg <- bg %>%
  sdm_extract(
    data = .,
    x = "x",
    y = "y",
    env_layer = CurrentBios,
    filter_na = TRUE
  )

# 3. Projections layers ------------------------------------------------

## 3.1. Mid Holocene (MH) ----------------------------------------------
# Pleistocene: mid-Holocene, Northgrippian (8.326-4.2 ka), v1.0 
# Genero vector con direccion de cada capa bioclimatica
paths_MHBios <- list.files(path = "C:/Users/ale_f/OneDrive/Documentos/Capas_GIS/ClimVars/PaleoClim/2_MH",
                           pattern = ".tif$",full.names = TRUE)
# Las cargo en R, selecciono las bios a usar y corto por el area de estudio
MHBios <- rast(paths_MHBios) |>
  subset(names(CurrentBios)) |> # selecciono bios a usar
  crop(M_Stilpnia) |> # Corto por area de estudio
  mask(M_Stilpnia)

# Creo carpeta con variables de proyección para kuenm
dir.create("kuenm/G_variables/Set_1/MH")
lapply(names(MHBios), 
       function(x){
         writeRaster(MHBios[[x]], 
                     paste0("kuenm/G_variables/Set_1/MH/",x,".asc"),
                     overwrite=TRUE, NAflag = -9999)})
## 3.2. Last Glacial Maximum (LGM) --------------------------------------------
# (ca. 21 ka)

# Genero vector con direccion de cada capa bioclimatica
paths_LGMBios <- list.files(path = "C:/Users/ale_f/OneDrive/Documentos/Capas_GIS/ClimVars/PaleoClim/3_LGM/2_5min",
                            pattern = ".tif$",full.names = TRUE)
# Las cargo en R, selecciono las bios a usar y corto por el area de estudio
LGMBios <- rast(paths_LGMBios) |>
  subset(names(CurrentBios)) |> # selecciono bios a usar
  crop(M_Stilpnia) |> # Corto por area de estudio
  mask(M_Stilpnia)
  
# Creo carpeta con variables de proyección para kuenm
dir.create("kuenm/G_variables/Set_1/LGM")
lapply(names(LGMBios), 
       function(x){
         writeRaster(LGMBios[[x]], 
                     paste0("kuenm/G_variables/Set_1/LGM/",x,".asc"),
                     overwrite=TRUE,                      NAflag = -9999)})

## 3.3. Last Interglacial (LIG) ---------------------------------------------
# Genero vector con direccion de cada capa bioclimatica
paths_LIGBios <- list.files(path = "C:/Users/ale_f/OneDrive/Documentos/Capas_GIS/ClimVars/PaleoClim/4_LIG",
                            pattern = ".tif$",full.names = TRUE)
# Las cargo en R, selecciono las bios a usar y corto por el area de estudio
LIGBios <- rast(paths_LIGBios) |>
  subset(names(CurrentBios)) |> # selecciono bios a usar
  crop(M_Stilpnia) |> # Corto por area de estudio
  mask(M_Stilpnia)

# Creo carpeta con variables de proyección para kuenm
dir.create("kuenm/G_variables/Set_1/LIG")
lapply(names(LIGBios), 
       function(x){
         writeRaster(LIGBios[[x]], 
                     paste0("kuenm/G_variables/Set_1/LIG/",x,".asc"),
                     overwrite=TRUE,                      NAflag = -9999)})
