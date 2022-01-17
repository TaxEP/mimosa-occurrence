#=====================================================================================================#

library(ape)
library(tidyverse)

#=====================================================================================================#

#=======#
# INPUT #
#=======#

# Reading CR dataset
mimosa_clean <- read.csv("data/dataset/mimosa_coordClean.csv", na.strings = c("", NA))

# Reading grids
grids_PNCV <- readOGR("data/shapefiles/grids/PNCV.shp") 
grids_PNSC <- readOGR("data/shapefiles/grids/PNSC.shp")

# Projection
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Intersecting dataset and grids

# PNCV
coords_pncv <- mimosa_clean
coords_pncv_2 <- coords_pncv 
coordinates(coords_pncv_2) <- ~ longitude + latitude
proj4string(coords_pncv_2) <- crswgs84
coords_pncv_2 <- over(coords_pncv_2, grids_PNCV)
coords_pncv_2$id_2 <- 1:nrow(coords_pncv_2) 
coords_pncv$id_2 <- 1:nrow(coords_pncv_2)
coords_pncv_2 <- coords_pncv_2 %>% filter(!is.na(id))
coords_pncv <- coords_pncv %>% filter(id_2 %in% coords_pncv_2$id_2)
coords_pncv$id_grid <- coords_pncv_2$id
coords_pncv <- coords_pncv[ , -which(colnames(coords_pncv) == "id_2")]
# checking
# coordinates(coords_pncv) <- ~longitude + latitude
# proj4string(coords_pncv) <- crswgs84
# plot(grids_PNCV); points(coords_pncv)

# PNSC
coords_pnsc <- mimosa_clean
coords_pnsc_2 <- coords_pnsc 
coordinates(coords_pnsc_2) <- ~ longitude + latitude
proj4string(coords_pnsc_2) <- crswgs84
coords_pnsc_2 <- over(coords_pnsc_2, grids_PNSC)
coords_pnsc_2$id_2 <- 1:nrow(coords_pnsc_2) 
coords_pnsc$id_2 <- 1:nrow(coords_pnsc_2)
coords_pnsc_2 <- coords_pnsc_2 %>% filter(!is.na(id))
coords_pnsc <- coords_pnsc %>% filter(id_2 %in% coords_pnsc_2$id_2)
coords_pnsc$id_grid <- coords_pnsc_2$id
coords_pnsc <- coords_pnsc[ , -which(colnames(coords_pnsc) == "id_2")]
# checking
# coordinates(coords_pnsc) <- ~longitude + latitude
# proj4string(coords_pnsc) <- crswgs84
# plot(grids_PNSC); points(coords_PNSC)

# Writing csv
write.csv(coords_pncv, "data/dataset/coords_PNCV.csv", row.names = F) 
write.csv(coords_pnsc, "data/dataset/coords_PNSC.csv", row.names = F) 
