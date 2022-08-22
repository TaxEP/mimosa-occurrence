# Function to get equivalent color with different transparency
makeTransparent = function(..., alpha=0.5) {
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  return(newColor)
}

#=====================================================================================================#

library(tidyverse)
library(raster)
library(sp)
library(rgdal)
library(latticeExtra)
library(viridis)

setwd("B:/yagob/GoogleDrive/Academia/Parallel_projects/Mimosa_occurrence")

# Reading csv
pachy <- read.csv("datasets/mimosa-clean.csv")

# Standardizing names
pachy$gen_sp <- paste(pachy$genus, "_", pachy$species, sep = "")

# Reading and separating pachycarpae names
list_pachy <- scan("lista_pachycarpae.txt", what = "character")
list_pachy <- str_split(list_pachy, pattern = ",")[[1]]

# Filtering
pachy <- pachy[pachy$gen_sp %in% list_pachy, ]

# Reading shapefile: biomes
biomes <- readOGR("shapefiles/biomes/bioma_1milhao_uf2015_250mil_IBGE_albers_v4_revisao_pampa_lagoas.shp")

crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # EPSG:4326 - WGS 84


coords <- pachy
coordinates(coords) <- ~ longitude + latitude
proj4string(coords) <- crswgs84

plot(biomes)
points(coords, col = makeTransparent("blue", alpha = 0.5), lwd = 3, cex = 0.5)

#Loading an elevation raster 
alt <- raster('HYP_HR_SR/HYP_HR_SR.tif')
e <- as(extent(biomes), 'SpatialPolygons')
crs(e) <- crswgs84
r <- crop(alt, e)
rr <- mask(r, biomes)

plot(rr)

alt_plot  <- spplot(rr, col.regions = grey(1:100/100, alpha = 0.7), 
                    maxpixels = 2e10, colorkey = FALSE)

alt_plot + layer(panel.points(longitude, latitude, col="green", pch=19), data=coords@coords)


#=====================================================================================================#

# Reading grid cells
grids_br <- readOGR("shapefiles/grids_br/grids_br.shp") 

# Matrix
# Creating a presence/absence matrix
pachy_matrix <- matrix(data = NA, nrow = length(unique(pachy$id_grid)), 
                        ncol = length(unique(pachy$gen_sp)))
pachy_matrix <- as.data.frame(pachy_matrix)
colnames(pachy_matrix) <- unique(pachy$gen_sp)
rownames(pachy_matrix) <- unique(pachy$id_grid)
for(i in 1:nrow(pachy_matrix)){
  for(j in 1:ncol(pachy_matrix)){
    if(colnames(pachy_matrix)[j] %in% pachy$gen_sp[pachy$id_grid == rownames(pachy_matrix)[i]]){
      pachy_matrix[i, j] <- 1
    } else {
      pachy_matrix[i, j] <- 0  
    }
  }
}

pachy_df <- data.frame("id" = rownames(pachy_matrix), "sr" = rowSums(pachy_matrix))

# Assigning results for each grid cell by merging results and a spatial polygon
sr_poly <- merge(grids_br, pachy_df, by.x = "id")

summary(sr_poly@data)

# Species richness
sr_plot <- spplot(sr_poly,
                  zcol = "sr", 
                  xlim = c(-73.99045 , -28.99045), ylim = c(-33.72816 , 5.271841),
                  par.settings = list(fontsize = list(text = 21)),
                  at = seq(0, 31, length.out = 16),
                  colorkey = TRUE, 
                  sp.layout = list(list(biomes, fill = "gray90", lwd = 2)), 
                  col.regions = makeTransparent(rev(magma(16)), alpha = 0.85), 
                  scales = list(draw = FALSE))
