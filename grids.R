#======================================================================================================#

library(raster)
library(rgdal)

# Defining grid cells' size
grids_size <- c(0.1, 0.1)

# Reading shapefiles
pncv <- readOGR("data/shapefiles/PNCV/2conferido_Zoneamento_PNCV_10_03_2020.shp") 
pnsc <- readOGR("data/shapefiles/PNSC/parna_serra_do_cipo-polygon.shp") 

# Checking shapefiles
plot(pncv)
plot(pnsc)

# Projecting
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") #EPSG:4326 - WGS 84
pncv <- spTransform(pncv, crswgs84)
pnsc <- spTransform(pnsc, crswgs84)

# Grid PNCV
grid_pncv <- raster(extent(pncv), resolution = grids_size, crs = proj4string(pncv))
gridPolygon_pncv <- rasterToPolygons(grid_pncv)
gridPolygon_pncv$id <- 1:nrow(gridPolygon_pncv)
# checking 
plot(pncv); plot(gridPolygon_pncv, add = T)

# Grid PNSC
grid_pnsc <- raster(extent(pnsc), resolution = grids_size, crs = proj4string(pnsc))
gridPolygon_pnsc <- rasterToPolygons(grid_pnsc)
gridPolygon_pnsc$id <- 1:nrow(gridPolygon_pnsc)
# checking 
plot(pnsc); plot(gridPolygon_pnsc, add = T)

# Creating spatial grids representing the campos rupestres areas
grid <- raster(extent(pnsc), resolution = grids_size, crs = proj4string(pnsc))
gridPolygon <- rasterToPolygons(grid)
gridPolygon$id <- 1:nrow(gridPolygon)
intersectGridClipped <- raster::intersect(gridPolygon, pnsc)
intersectGrid <- gridPolygon[gridPolygon$id %in% intersectGridClipped$id, ]

# Saving spatial grids shapefile
writeOGR(gridPolygon_pncv, ".", "data/shapefiles/grids/PNCV", driver = "ESRI Shapefile")
writeOGR(gridPolygon_pnsc, ".", "data/shapefiles/grids/PNSC", driver = "ESRI Shapefile")
