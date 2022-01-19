#Function to get equivalent color with different transparency
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

library(picante)
library(viridis)
library(tidyverse)
library(rgdal)
library(ape)

#=====================================================================================================#

#=======#
# INPUT #
#=======#

# Reading matrix
matrix_PNCV <- read.csv(file = "data/dataset/mimosa_matrix_PNCV.csv", row.names = 1)

# Reading tree
mimosa_tree <- read.nexus("data/trees/pruned_tree-mimosa_PNCV.nex")

# Reading grids
grids_PNCV <- readOGR("data/shapefiles/grids/PNCV.shp") 

# Reading PNCV shapefile
PNCV <- readOGR("data/shapefiles/PNCV/2conferido_Zoneamento_PNCV_10_03_2020.shp") 

# Projecting
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # EPSG:4326 - WGS 84
PNCV <- spTransform(PNCV, crswgs84)
grids_PNCV <- spTransform(grids_PNCV, crswgs84)
#=====================================================================================================#

#====#
# PD #
#====#

# Setting random number generator and running phylogenetic diversity analysis
set.seed(7)
pd_stats <- ses.pd(matrix_PNCV, mimosa_tree, include.root = TRUE, null.model = "taxa.label")

# Creating a column specifying grid ids
pd_stats$id <- rownames(pd_stats)

# Running a linear regression and incorporating residuals
pd_stats$residuals <- lm(pd.obs ~ ntaxa, pd_stats)$res

# Rounding values in order to standardize plot size (except for p-values, id and ntaxa)
pd_stats[ , 
          -which(colnames(pd_stats) %in% c("ntaxa",
                                           "pd.obs.p", 
                                           "id"))] <- round(pd_stats[ , 
                                                                      -which(colnames(pd_stats) %in% c("ntaxa",
                                                                                                       "pd.obs.p", 
                                                                                                       "id"))], 
                                                            2)

#Saving results
#write.csv(file = "results/Mimosa/pd_stats.csv", pd_stats, row.names = FALSE)

#Loading results
#pd_stats <- read.csv("results/Mimosa/pd_stats.csv")

#Assigning values for each grid by merging results and spatial grids
pd_poly <- merge(grids_PNCV, pd_stats, by.x = "id")

#=====================================================================================================#

#=========#
# FIGURES #
#=========#

# Label limits based on the range of each attribute. This was retrieved using summary(pd_poly@data)
summary(pd_poly@data)

# Getting grids extension
extent(grids_PNCV)

#PD
pd_plot <- spplot(pd_poly,
                  zcol = "pd.obs",
                  xlim = c(-47.90742, -47.00742 ), ylim = c(-14.21845, -13.61845),
                  par.settings=list(fontsize = list(text = 21)),
                  at = seq(0, 128.87, length.out = 16),
                  colorkey = TRUE, 
                  sp.layout = list(list(PNCV, fill = "gray")), 
                  col.regions = makeTransparent(rev(magma(16)), alpha = 0.7), 
                  scales = list(draw = FALSE))

#SR
sr_plot <- spplot(pd_poly,
                  zcol = "ntaxa", 
                  xlim = c(-47.90742, -47.00742 ), ylim = c(-14.21845, -13.61845),
                  par.settings=list(fontsize = list(text = 21)),
                  at = seq(0, 21, length.out = 16),
                  colorkey = TRUE, 
                  sp.layout = list(list(PNCV, fill = "gray")), 
                  col.regions = makeTransparent(rev(magma(16)), alpha = 0.7), 
                  scales = list(draw = FALSE))

# Residuals
pdres_plot <- spplot(pd_poly,
                     zcol = "residuals",
                     xlim = c(-47.90742, -47.00742 ), ylim = c(-14.21845, -13.61845),
                     par.settings=list(fontsize = list(text = 21)),
                     at = seq(-17.01001, 18.91001, length.out = 16),
                     colorkey = TRUE, 
                     sp.layout = list(list(PNCV, fill = "gray")), 
                     col.regions = makeTransparent(c(viridis(16)[1:7], "khaki1", rev(heat.colors(16)[1:8])), alpha = 0.7), 
                     scales = list(draw = FALSE))
