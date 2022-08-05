#======================================================================================================#

library(raster)
library(rgdal)

# Set wd (Yago)
setwd("B:/yagob/GoogleDrive/Academia/Parallel_projects/Mimosa_occurrence")

# Set wd (Monique)
setwd("G:/.shortcut-targets-by-id/19Bt9xRgbQsy9ySW31FgR7E5aEscE0jKG/Mimosa_occurrence")

#======================================================================================================#

#==================#
# HONEY COMB GRIDS # 
#==================#

# https://strimas.com/post/hexagonal-grids/#:~:text=Creating%20grids-,Hexagonal%20grids,grid%20of%20polygons%20with%20HexPoints2SpatialPolygons%20

# Defining grid cells' size
size <- 1

# Reading the shapefile of the Brazilian terrestrial territory
br <- readOGR("shapefiles/BR/BR_UF_2020.shp")

# Projecting
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # EPSG:4326 - WGS 84
br <- spTransform(br, crswgs84)

# There is no function in R that will directly generate a hexagonal 
# grid of polygons covering a given region; however, it can be accomplished 
# by first generating a hexagonal grid of points with spsample, then converting 
# this point grid to a grid of polygons with HexPoints2SpatialPolygons.
hex_points <- spsample(br, type = "hexagonal", cellsize = size)
hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = size)
#plot(br, col = "grey90", axes = T)
#plot(hex_grid, border = "black", add = T)

## Assigning grids' IDs

# Extract polygon IDs
pid <- sapply(slot(hex_grid, "polygons"), function(x) slot(x, "ID"))

# Create dataframe row names
p.df <- data.frame( id = 1:length(hex_grid), row.names = pid)   

# Coersing and testing class
p <- SpatialPolygonsDataFrame(hex_grid, p.df)
class(p)

# Saving spatial grids in a shapefile
writeOGR(p, ".", "shapefiles/grids_br/hc_grids", driver="ESRI Shapefile")

#======================================================================================================#







