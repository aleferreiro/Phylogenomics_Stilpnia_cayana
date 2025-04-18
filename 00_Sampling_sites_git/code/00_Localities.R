# Packages needed -----------------------------------------------------
library(tidyverse)
library(sp)
library(measurements)
library(sf)
library(tmap)
library(tmaptools)
library(rnaturalearth)
library(rnaturalearthdata)

# 0. Load raw sampling data ---------------------------------------------------

Coords_raw <- read.csv("00_Sampling_sites/Samples_tangara_NGS.csv", sep = ";") 
round(Coords_raw$lat3) # redondeo para sacar los decimales de los segundos sino no anda la funcion
round(Coords_raw$long3) # redondeo para sacar los decimales de los segundos sino no anda la funcion

# Cleaning the data frame
Coords_clean <- Coords_raw |>
  filter(Sample != "UFG3804") |> # Quito el outgroup
#  filter(Sample != "UFG4028", 
#         Sample != "UFG5660",
#         Sample != "UFG5667",  # !BAJA CALIDAD DE READS!!
#         Sample != "UFG5668",
#         Sample != "UFG5692",
#         Sample != "UFG5753") |> # Quito muestras de baja calidad de reads
  unite(Latitude_GMS, lat1, lat2, lat3, sep = " ") |> 
  unite(Longitude_GMS, long1, long2, long3, sep = " ") # Renombro columna de coordenadas en GMS


# 1. Convert GMS in DEC degrees --------------------------------------------
Coords <- Coords_clean |> 
  mutate(Latitude = conv_unit(Latitude_GMS, "deg_min_sec", "dec_deg"),
         Longitude = conv_unit(Longitude_GMS, "deg_min_sec", "dec_deg")) 
# Coords for Marajo samples are 0 but from the southern side of the Ecuador. 
# I multiply them by -1 to turn negative
Coords[20:21,11] <- -0.860833333333333
Coords[22,11] <- -0.852222222222222

# Create a column that differentiates samples that were at south or 
# north of Amazon river
Coords_NS <- mutate(Coords, popNS = if_else(Latitude >= 0,"North","South"))

# Then convert to sf object 
Coords_sf <- st_as_sf(Coords_NS, coords = c("Longitude", "Latitude"))
head(Coords_sf)

## Add columns with Bioregions (NOT USED) --------------

# Para crear una columna que me diga en que ecorregion esta cada registro es
# un poco más complejo porque tengo que operar con objetos espaciales

# Primero cargo un vector con las ecoregiones de Sudamerica de 2017
Ecoregions <- st_read("C:/Capas_GIS/Bioregiones/Ecoregions2017/Ecoregions2017_sudam.gpkg")
# y las regiones biogeograficas segun Morrone, 2022
Bioregions <- st_read("C:/Capas_GIS/Bioregiones/BioRegionesNeotrop_Morrone2022/Geographic_vector_all/NeotropicMap_Geo.shp")
st_crs(Coords_sf) <- st_crs(Ecoregions)
sf_use_s2(FALSE) # Desactivo algo si no tengo error

# Usando la funcion st_join me agregara al Coords_sf una columna con los nombres
# de la ecoregion en la que se encuentran
st_crs(Coords_sf) <- st_crs(Ecoregions)
Coords_sf_ecoreg <- st_join(Coords_sf, Ecoregions["ECO_NAME"]) |> 
  rename(popEcoreg = ECO_NAME) # Le cambio el nombre a la columna

# Usando la funcion st_join me agregara al Coords_sf una columna con los nombres
# de la ecoregion en la que se encuentran
st_crs(Coords_sf_ecoreg) <- st_crs(Ecoregions)
Coords_sf_ecoreg2 <- st_join(Coords_sf_ecoreg, Ecoregions["BIOME_NAME"]) |> 
  rename(popBiome = BIOME_NAME) # Le cambio el nombre a la columna

# sumo tb otra columna con el nombre la provincia biogeografica (Morrone, 2022)
st_crs(Coords_sf_ecoreg) <- st_crs(Bioregions)
Coords_sf_ecoYbioreg <- st_join(Coords_sf_ecoreg2, Bioregions["Provincias"]) |> 
  rename(popBioreg = Provincias) # Le cambio el nombre a la columna


# 2. Add population for each sample -------------------------------------

## 2.1. Five Population column -------------------------------------------------------

# Load polygons with each pop generated in QGIS
pops <- st_read("data/pops_Tangara_cayana.gpkg") 

# Add pop column
st_crs(Coords_sf) <- st_crs(pops)
Coords_sf_pop <- st_join(Coords_sf, pops["pop"])

# Modify the pop column in order to get the 5 populations finally used in 
# G-PhoCS and fastsimcoal
pop5 <- Coords_sf_pop$pop
pop5[pop5 == "NorthEast Cerrado"] <- "East"
pop5[pop5 == "SouthEast Cerrado"] <- "East"
pop5[pop5 == "West Cerrado"] <- "West"
pop5[pop5 == "Amapa"] <- "North"
pop5[pop5 == "Guyana"] <- "North"
pop5[pop5 == "Para"] <- "North"
pop5[is.na(pop5)] <- "West"
pop5

Coords_sf_pop5 <- Coords_sf_pop |> 
  mutate(Populations = pop5, .before = geometry) 

## 2.2. Population column with north pop divided -------------------------------

Coords_sf_pop2 <- Coords_sf_pop5 |> 
  mutate(Populations2 = replace(pop,
                              pop == "NorthEast Cerrado" | pop == "SouthEast Cerrado",
                              "East"), .before = geometry) |> 
  mutate(Populations2 = replace(Populations2,
                                Populations2 == "West Cerrado",
                                "West")) 
 #  mutate(Populations2 = replace(Populations2,
 #                                Populations2 == "Guyana",
 #                                "North (Guyana)")) |>
 #  mutate(Populations2 = replace(Populations2,
 #                               Populations2 == "Para",
 #                               "North (Para)")) |>
 #  mutate(Populations2 = replace(Populations2,
 #                                Populations2 == "Amapa",
 #                                "North (Amapa)"))
  

## 2.3. Population column for 6 populations ---------------------------------

# Used in Model 11 of fastsimcoal2.

Coords_sf_pop3 <- Coords_sf_pop2 |> 
  mutate(Populations3 = replace(pop,
                                pop == "NorthEast Cerrado" | pop == "SouthEast Cerrado",
                                "East"), .before = geometry) |> 
  mutate(Populations3 = replace(Populations3,
                                Populations3 == "West Cerrado",
                                "West")) |> 
  mutate(Populations3 = replace(Populations3,
                                Populations3 == "Para" | Populations3 == "Amapa",
                                "ParaAmapa"))
# Save as spatial vector
st_write(Coords_sf_pop3, "00_Sampling_sites/Sampling_sites.gpkg")

# Convert to data frame
Coords_df_pop3 <- Coords_sf_pop3 |> 
  sfheaders::sf_to_df(fill = T) |> 
  mutate(sfg_id = NULL,
         point_id = NULL) |> 
  rename(Latitude = y,
         Longitude = x)
# Save data frame
write.csv(Coords_df_pop3, "00_Sampling_sites/Sampling_sites.csv")


# 3. Genomics samples analyzed ------------------------------------------

#  sampling sites of only those samples used in the analysis
Localities_83 <- read.csv("00_Sampling_sites/Populations_Stilpnia.csv")
Coords_sf_pop83 <- Coords_sf_pop3 |> 
  filter(Sample %in% Localities_83$Sample) # filtro por muestras usadas en analisis
  
# Guardo el archivo tanto en formato espacial
st_write(Coords_sf_pop83, "00_Sampling_sites/Sampling_sites83.gpkg")

# Lo vuelvo a data frame
Coords_df_pop83 <- Coords_sf_pop83 |> 
  sfheaders::sf_to_df(fill = T) |> 
  mutate(sfg_id = NULL,
         point_id = NULL) |> 
  rename(Latitude = y,
         Longitude = x)

write.csv(Coords_df_pop83, "00_Sampling_sites/Sampling_sites83.csv")



# 4. Map with sampling sites (Opcional) --------------------------------------

## 4.1. Mapa base ---------------------------------------------------------------------------
# Cargo limites políticos de Sudamerica

# Países
Sudam_raw <- ne_countries(scale = 10, 
                          continent = "south america", 
                          returnclass = "sf")
# Add french guiana
French_guiana <- ne_countries(scale = 10,
                              type = "map_units",
                              geounit = "French Guiana",
                              returnclass = "sf")
Sudam <- st_union(Sudam_raw, French_guiana)
qtm(Sudam)

# Muestreo neotropical savannas del shape de Ecoregiones
Savannas <- Ecoregions |> 
  subset(BIOME_NAME == "Tropical & Subtropical Grasslands, Savannas & Shrublands") |> 
  subset(!(ECO_NAME == "Dry Chaco" | 
             ECO_NAME == "Humid Chaco" |
             ECO_NAME == "Uruguayan savanna"))

# Cargo cuenca del Amazonas
Amazona <- st_read("C:/Capas_GIS/Capas_Sudam/Hidrografia_Sudam/Hidrografia_Sudamérica.shp") |> 
  filter(NOMBRE == "Amazonas") |> 
  slice(1)
Amazona[1,1] <- "Amazon river"

qtm(Amazona)

# Load geographic range of Stilpnia cayana from BirdLife
Scay_range1 <- st_read("C:/Capas_GIS/Species_range/Stilpnia_cayana_BirdLife2024/SppDataRequest.shp")
Scay_range2 <- st_read("C:/Capas_GIS/Species_range/Stilpnia_cayana_BelenBukowski/Tangara_cayana_22722884.shp")

qtm(Scay_range2)

# Mapa base
mapa_base <- tm_shape(Sudam, bbox = c(-83, # longitud mínima 
                                      -28, # latitud mínima
                                      -33, # longitud máxima
                                      12   # latitud máxima
                                      )) + tm_polygons(col = "#f0f0f0", lwd = 1) +
  tm_shape(Scay_range2) + tm_polygons(col = "#969696") +
  tm_shape(Sudam) + tm_borders() +
  tm_shape(Amazona) + 
    tm_lines(col = "#2c7fb8", lwd = 4) + 
    tm_text(text = "NOMBRE", xmod = 2.3, ymod = 0.7) +
  tmap_options(check.and.fix = TRUE)
mapa_base

## 2.2. Inset map ----------------------------------------------------------
# convert bbox to a poligon
bbox = c(-83, # longitud mínima 
         -28, # latitud mínima
         -33, # longitud máxima
         12)   # latitud máxima

bb <- tmaptools::bb_poly(bbox)

inset_map <- tm_shape(Sudam, bbox = c(-85, # longitud mínima 
                                      -55, # latitud mínima
                                      -30, # longitud máxima
                                      14)   # latitud máxima
) + tm_polygons(col = "#f0f0f0", lwd = 1) +
  tm_shape(bb) + tm_borders(lwd = 2) +
  tmap_options(check.and.fix = TRUE)
inset_map

tmap_save(tm = inset_map,
          filename = "00_Sampling_sites/plots/Inset_map.pdf")

## 2.3. Map by population ---------------------------------------------------------------

# Order populations north to south
Coords_sf_pop83$Populations2 <- factor(Coords_sf_pop83$Populations2, 
                                     levels = c("Guyana",
                                                "Para",
                                                "Amapa",
                                                "Marajo",
                                                "Caatinga", 
                                                "East", 
                                                "West"))
Coords_sf_pop83_map <- Coords_sf_pop83 |>
  select(!Populations) |> 
  rename(Populations = Populations2)

# Map
Mapa_Sampling_sites <- mapa_base + 
  tm_shape(Savannas) + tm_polygons(col = "#addd8e",
                                   lwd = 0,
                                   #border.col = "",
                                   legend.show = F) +
  tm_shape(Coords_sf_pop83_map) +  tm_symbols(col = "Populations",
                                          size = 0.7,
                                          shape = 21,
                                          palette = c("#d4b9da", # Guyana
                                                      "#ff7f2a", # Para
                                                      "#e7298a", # Amapa
                                                      "#a50f15", # Marajo
                                                      "#54278f", # Caatinga
                                                      "#dab600", # East
                                                      "#6baed6"),# West
                                          border.col = "black") +
  tm_layout(title = "",
            title.position = c("left", "top"),
            legend.position = c(0.8,0.65),
            frame = F)
Mapa_Sampling_sites

tmap_save(tm = Mapa_Sampling_sites,
          filename = "00_Sampling_sites/plots/Map_sampling_sites.pdf")
tmap_save(tm = Mapa_Sampling_sites,
          filename = "00_Sampling_sites/plots/Map_sampling_sites.png")


## 2.4. Smaller map by population (Fig 2) --------------------------------------

small_map_sampSite <- tm_shape(Sudam, bbox = c(-83, # longitud mínima 
                                               -25, # latitud mínima
                                               -33, # longitud máxima
                                               12   # latitud máxima
                                               )) + tm_polygons(col = "#f0f0f0", lwd = 1) +
  # tm_shape(Savannas) + tm_polygons(col = "#addd8e",
  #                                  lwd = 0,
  #                                  #border.col = "",
  #                                  legend.show = F) +
  tm_shape(Coords_sf_pop83_map) +  tm_symbols(col = "Populations",
                                              size = 3,
                                              shape = 21,
                                              palette = c("#d4b9da", # Guyana
                                                          "#df65b0", # Para
                                                          "#e7298a", # Amapa
                                                          "#a50f15", # Marajo
                                                          "#54278f", # Caatinga
                                                          "#dab600", # East
                                                          "#6baed6"),# West
                                              border.col = "black") +
  tm_shape(Amazona) + 
  tm_lines(col = "#2c7fb8", lwd = 4) + 
  tm_text(text = "NOMBRE", xmod = 2.3, ymod = 0.7) +
  tm_layout(title = "",
            title.position = c("left", "top"),
            legend.position = c(0.8,0.65),
            frame = F) 
small_map_sampSite

tmap_save(tm = small_map_sampSite,
          filename = "00_Sampling_sites/plots/SmallMap_sampling_sites.pdf")



## Map by supsp ------------------------------------------------------------
# 
# # Order populations north to south
# Coords_sf_pop5$Subespecie <- factor(Coords_sf_pop5$Subespecie, 
#                                     levels = c("cayana", 
#                                                "huberi", 
#                                                "flava",
#                                                "sincipitalis", 
#                                                "margaritae",
#                                                "chloroptera"))
# 
# # Map
# Mapa_Sampling_sites_subsp <- mapa_base + 
#   tm_shape(Savannas) + tm_polygons(col = "#addd8e",
#                                    lwd = 0,
#                                    legend.show = F) +
#   tm_shape(Coords_sf_pop5) +  tm_symbols(col = "Subespecie", 
#                                          size = 0.7, 
#                                          shape = 21,
#                                          palette = c("#df65b0", # cayana
#                                                      "#a50f15", # huberi
#                                                      "#54278f", # flava
#                                                      "#dab600", # sincipitalis
#                                                      "#6baed6", # margaritae
#                                                      "#007e19"),# chloroptera
#                                          border.col = "black") +
#   tm_layout(title = "",
#             title.position = c("left", "top"),
#             legend.position = c(0.87,0.8),
#   )
# Mapa_Sampling_sites_subsp
# 
# tmap_save(tm = Mapa_Sampling_sites_subsp,
#           filename = "00_Sampling_sites/plots/Map_sampling_sites_subsp.pdf")
# tmap_save(tm = Mapa_Sampling_sites_subsp,
#           filename = "00_Sampling_sites/plots/Map_sampling_sites_subsp.png")
# 

## 2.6. Smaller map with amazon river and sampling sites --------------------------------------

small_map_sampSite <- tm_shape(Sudam, bbox = c(-83, # longitud mínima 
                                               -25, # latitud mínima
                                               -33, # longitud máxima
                                               12   # latitud máxima
)) + tm_polygons(col = "#f0f0f0", lwd = 1) +
  # tm_shape(Savannas) + tm_polygons(col = "#addd8e",
  #                                  lwd = 0,
  #                                  #border.col = "",
  #                                  legend.show = F) +
  tm_shape(Coords_sf_pop83_map) +  tm_symbols(col = "black") +
  tm_shape(Amazona) + 
  tm_lines(col = "#2c7fb8", lwd = 4) + 
  tm_text(text = "NOMBRE", size = 2, xmod = 2.3, ymod = 0.7) +
  tm_layout(title = "",
            title.position = c("left", "top"),
            legend.position = c(0.8,0.65),
            frame = F) 
small_map_sampSite

tmap_save(tm = small_map_sampSite,
          filename = "00_Sampling_sites/plots/SmallMap_sampling_sites_RAO.pdf")


## 2.7. Savanna map for oral presentation ----------------------------------

savanna_map <- tm_shape(Sudam, bbox = c(-83,   # longitud mínima 
                                        -28,   # latitud mínima
                                        -33,   # longitud máxima
                                        12 # latitud máxima
  )) + tm_polygons(col = "#f0f0f0", lwd = 1) +
  tm_shape(Savannas) + tm_polygons(col = "#addd8e",
                                   lwd = 0,
                                   #border.col = "",
                                   legend.show = F) +
  # tm_shape(Coords_sf_pop83_map) +  tm_symbols(col = "black",
  #                                             size = 0.5,
  #                                             shape = 21) +
  tm_shape(Sudam) + tm_borders() +
  tm_shape(Amazona) + 
  tm_lines(col = "#2c7fb8", lwd = 4) + 
  tm_text(text = "NOMBRE", xmod = 2.3, ymod = 0.7) +
  tmap_options(check.and.fix = TRUE)
savanna_map

tmap_save(tm = savanna_map,
          filename = "00_Sampling_sites/plots/Savannas_RAO.pdf")


## 2.8. Stilpnia cayana map -----------------------------------------------------

Scayana_map <- tm_shape(Sudam, bbox = c(-83,   # longitud mínima 
                                        -28,   # latitud mínima
                                        -33,   # longitud máxima
                                        12) # latitud máxima
                        ) + 
    tm_polygons(col = "#f0f0f0", lwd = 1) +
  tm_shape(Savannas) + 
    tm_polygons(col = "#addd8e",
                                   lwd = 0,
                                   #border.col = "",
                                   legend.show = F) +
  tm_shape(Scay_range2) + 
    tm_polygons(col = "#969696", alpha = 0.4) +
  # tm_shape(Coords_sf_pop83_map) +  tm_symbols(col = "Populations",
  #                                             size = 0.7,
  #                                             shape = 21,
  #                                             palette = c("#d4b9da", # Guyana
  #                                                         "#ff7f2a", # Para
  #                                                         "#e7298a", # Amapa
  #                                                         "#a50f15", # Marajo
  #                                                         "#54278f", # Caatinga
  #                                                         "#dab600", # East
  #                                                         "#6baed6"),# West
  #                                             border.col = "black") +
  tm_shape(Amazona) + 
    tm_lines(col = "#2c7fb8", lwd = 4) + 
  #  tm_text(text = "NOMBRE", xmod = 2.3, ymod = 0.7) +
  tm_layout(title = "",
            title.position = c("left", "top"),
            legend.position = c(0.8,0.65),
            frame = F)
Scayana_map

tmap_save(tm = Scayana_map,
          filename = "00_Sampling_sites/plots/Scayana_map_RAO.pdf")


# 2.9. Sampling sites -----------------------------------------------------

Sampling_sites_map <- tm_shape(Sudam, bbox = c(-83,   # longitud mínima 
                                        -28,   # latitud mínima
                                        -33,   # longitud máxima
                                        12) # latitud máxima
  ) + 
    tm_polygons(col = "#f0f0f0", lwd = 1) +
  tm_shape(Savannas) + 
    tm_polygons(col = "#addd8e",
              lwd = 0,
              #border.col = "",
              legend.show = F) +
  tm_shape(Scay_range2) + 
    tm_polygons(col = "#969696", alpha = 0.4) +
  # tm_shape(Coords_sf_pop83_map) +  tm_symbols(col = "Populations",
  #                                             size = 0.7,
  #                                             shape = 21,
  #                                             palette = c("#d4b9da", # Guyana
  #                                                         "#ff7f2a", # Para
  #                                                         "#e7298a", # Amapa
  #                                                         "#a50f15", # Marajo
  #                                                         "#54278f", # Caatinga
  #                                                         "#dab600", # East
  #                                                         "#6baed6"),# West
  #                                             border.col = "black") +
  tm_shape(Amazona) + 
    tm_lines(col = "#2c7fb8", lwd = 4) + 
  #  tm_text(text = "NOMBRE", xmod = 2.3, ymod = 0.7) +
  tm_shape(Coords_sf_pop83_map) + 
    tm_symbols(col = "black", size = 0.5) +
  tm_layout(title = "",
            title.position = c("left", "top"),
            legend.position = c(0.8,0.65),
            frame = F)
Sampling_sites_map

tmap_save(tm = Sampling_sites_map,
          filename = "00_Sampling_sites/plots/Sampling_sites_map_RAO.pdf")

