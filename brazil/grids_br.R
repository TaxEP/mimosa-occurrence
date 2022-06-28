#======================================================================================================#

library(raster)
library(rgdal)

# Set wd (Yago)
setwd("B:/yagob/GoogleDrive/Academia/Parallel_projects/Mimosa_occurrence")

# Set wd (Monique)
setwd("G:/.shortcut-targets-by-id/19Bt9xRgbQsy9ySW31FgR7E5aEscE0jKG/Mimosa_occurrence")

#======================================================================================================#

# Defining grid cells' size
grids_size <- c(3, 3)

# Reading the shapefile of the Brazilian terrestrial territory
br <- readOGR("shapefiles/BR/BR_UF_2020.shp")

# Projecting
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # EPSG:4326 - WGS 84
br <- spTransform(br, crswgs84)

# Creating spatial grid cells on the Brazilian terrestrial territory 
grid <- raster(extent(br), resolution = grids_size, crs = proj4string(br))
gridPolygon <- rasterToPolygons(grid)
gridPolygon$id <- 1:nrow(gridPolygon)
intersectGridClipped <- raster::intersect(gridPolygon, br)
intersectGrid <- gridPolygon[gridPolygon$id %in% intersectGridClipped$id, ]

# Saving spatial grids in a shapefile
writeOGR(intersectGrid, ".", "shapefiles/grids_br2/grids_br2", driver="ESRI Shapefile")
