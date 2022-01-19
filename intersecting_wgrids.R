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
coords_pncv <- mimosa_clean # loading the dataset 
coords_pncv_2 <- coords_pncv # Loading a copy that will undergo modifications below
coordinates(coords_pncv_2) <- ~ longitude + latitude # setting spatial coordinates as a SpatialPointsDataFrame
proj4string(coords_pncv_2) <- crswgs84 # projecting
coords_pncv_2 <- over(coords_pncv_2, grids_PNCV) # intersecting the dataset over the grids. this will generate a two-column data frame: NA values represent records that do not fall in the grids' area
coords_pncv_2$id_2 <- 1:nrow(coords_pncv_2) # generating a column that will assign a ID value for each record
coords_pncv$id_2 <- 1:nrow(coords_pncv_2) # the same is done for the original dataset
coords_pncv_2 <- coords_pncv_2 %>% filter(!is.na(id)) # removing NA values
coords_pncv <- coords_pncv %>% filter(id_2 %in% coords_pncv_2$id_2) # removing from the original dataset records that do not correspond to the records removed above
coords_pncv$id_grid <- coords_pncv_2$id # assigning a number that identify in which grid the record occurs
coords_pncv <- coords_pncv[ , -which(colnames(coords_pncv) == "id_2")] # removing the id_2 column, as it will not be required any further
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
