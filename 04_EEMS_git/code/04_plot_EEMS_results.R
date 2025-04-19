## Check that the current directory contains the rEEMSplots source directory (from GitHub)
if (dir.exists("rEEMSplots")) {
 install.packages("rEEMSplots", repos = NULL, type = "source")
} else {
  stop("Move to the directory that contains the rEEMSplots source to install the package.")
}

library("rEEMSplots")	
library("rworldmap")
library("rworldxtra")


projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=merc +datum=WGS84"
eems.plots(mcmcpath = 'data/Tangara_2023_83-II/05_EEMS/Tangara_2023_83-II_V1',
           plotpath = 'data/Tangara_2023_83-II/05_EEMS/Tangara_2023_83-II_V1',
           add.grid = FALSE, add.outline = FALSE, add.demes = TRUE,
           longlat = TRUE,
           projection.in = projection_none,
           projection.out = projection_mercator,
           add.map = TRUE,
           col.map = "black",
           lwd.map = 2,
           col.demes = "red")
