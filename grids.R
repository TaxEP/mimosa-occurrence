#======================================================================================================#

library(raster)
library(rgdal)

#Defining grid cells' size
grids_size <- c(0.6, 0.6)

#Reading shapefiles: campos rupestres and Brazilian terrestrial territory
# cr <- readOGR("shapefiles/shapefilecamporupestre/cr.shp") #Silveira (2016)
# br <- readOGR("shapefiles/br_unidades_da_federacao/BRUFE250GC_SIR.shp") #IBGE: https://downloads.ibge.gov.br/downloads_geociencias.htm

#checking shapefiles
br
go_mg
pncv
pnsc

#Projecting
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") #EPSG:4326 - WGS 84
proj4string(pnsc) <- crswgs84
br <- spTransform(br, crswgs84)

#Creating spatial grids representing the campos rupestres areas
grid <- raster(extent(pnsc), resolution = grids_size, crs = proj4string(pnsc))
gridPolygon <- rasterToPolygons(grid)
gridPolygon$id <- 1:nrow(gridPolygon)
intersectGridClipped <- raster::intersect(gridPolygon, pnsc)
intersectGrid <- gridPolygon[gridPolygon$id %in% intersectGridClipped$id, ]

#Saving spatial grids in a shapefile
#writeOGR(intersectGrid, ".", "shapefiles/grids_cr/grids_cr", driver="ESRI Shapefile")

# Getting coordinates for each grid
coords <- coordinates(intersectGrid)
dist_grids <- distm(coords)
rownames(dist_grids) <- rownames(intersectGrid@data)
colnames(dist_grids) <- rownames(intersectGrid@data)
write.csv(file = "results/dist_grids.csv", dist_grids)
