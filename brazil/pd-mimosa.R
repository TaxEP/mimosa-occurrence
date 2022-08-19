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

# Show colors
library(scales)
show_col(c("darkgreen","darkblue","darkred","darkorange"))
show_col(c("#231151FF",
           "#B63679FF",
           "#4AC16DFF",
           "#CFE11CFF"))

#=====================================================================================================#

library(picante)
library(viridis)
library(tidyverse)
library(rgdal)
library(ape)
library(raster)

# Set wd (Yago)
setwd("B:/yagob/GoogleDrive/Academia/Parallel_projects/Mimosa_occurrence")

# Set wd (Monique)
setwd("G:/.shortcut-targets-by-id/19Bt9xRgbQsy9ySW31FgR7E5aEscE0jKG/Mimosa_occurrence")

#=====================================================================================================#

#=======#
# INPUT #
#=======#

# Reading presence-absence matrix
mimosa_matrix <- read.csv(file = "datasets/mimosa_matrix.csv", row.names = 1)

# Reading tree
mimosa_tree <- read.nexus("trees/mimosa-pruned_tree.nex")

# Reading grid cells
grids_br <- readOGR("shapefiles/grids_br/grids_br.shp") 

# Reading shapefile: Brazilian terrestrial territory
br <- readOGR("shapefiles/BR/BR_UF_2020.shp")

# Reading shapefile: biomes
biomes <- readOGR("shapefiles/biomes/bioma_1milhao_uf2015_250mil_IBGE_albers_v4_revisao_pampa_lagoas.shp")

#=====================================================================================================#

#====#
# PD #
#====#

# Removing from the matrix species that do not occur in the tree 
mimosa_matrix <- mimosa_matrix[ , colnames(mimosa_matrix) %in% mimosa_tree$tip.label]

# Setting seed and running a phylogenetic diversity analysis
set.seed(7)
pd_stats <- ses.pd(mimosa_matrix, mimosa_tree, include.root = TRUE, null.model = "taxa.label")

# Creating a field that specifies grid cell id
pd_stats$id <- rownames(pd_stats)

# Running a linear regression and identifying residual values
pd_stats$residuals <- lm(pd.obs ~ ntaxa, pd_stats)$res

# Saving results as *.csv
write.csv(pd_stats, "datasets/mimosa_pd.csv", row.names = TRUE)

# Loading results
pd_stats <- read.csv("datasets/mimosa_pd.csv", row.names = 1)

# Assigning results for each grid cell by merging results and a spatial polygon
pd_poly <- merge(grids_br, pd_stats, by.x = "id")

# Defining cells with statistical significant PD
pd_poly$significance <- NA
for(i in 1:nrow(pd_poly@data)){
  if(is.na(pd_poly$pd.obs.p[i])){
    pd_poly$significance[i] <- NA
  } else if(pd_poly$pd.obs.p[i] >= 0.99 ){
    pd_poly$significance[i] <-  ">= 0.99"
  } else if(pd_poly$pd.obs.p[i] >= 0.975 & pd_poly$pd.obs.p[i] < 0.99){
    pd_poly$significance[i] <-  ">= 0.975"
  } else if(pd_poly$pd.obs.p[i] <= 0.025 & pd_poly$pd.obs.p[i] > 0.01){
    pd_poly$significance[i] <-  "<= 0.025"
  } else if(pd_poly$pd.obs.p[i] <= 0.01){
    pd_poly$significance[i] <-  "<= 0.01"
  } else {
    pd_poly$significance[i] <- "Not significant"
  }
}
pd_poly$significance <- factor(pd_poly$significance, levels = c("<= 0.01", 
                                                                "<= 0.025",
                                                                "Not significant",
                                                                ">= 0.975",
                                                                ">= 0.99"))

#=====================================================================================================#

#=========#
# FIGURES #
#=========#

# Scaling PD and species richness
pd_poly$pd.obs <- pd_poly$pd.obs/max(pd_poly$pd.obs, na.rm = T)
pd_poly$ntaxa<- pd_poly$ntaxa/max(pd_poly$ntaxa, na.rm = T)

# We can define label limits using the range of each attribute
summary(pd_poly@data)

# We must know the extent of the pd_poly object in order to define plot limits
extent(pd_poly)

# Phylogenetic diversty
pd_plot <- spplot(pd_poly,
                  zcol = "pd.obs",
                  xlim = c(-73.99045 , -28.99045), ylim = c(-33.72816 , 5.271841),
                  par.settings = list(fontsize = list(text = 21)),
                  at = seq(0, 0.336, length.out = 16),
                  colorkey = TRUE, 
                  sp.layout = list(biomes, fill = "gray90", lwd = 2), 
                  col.regions = makeTransparent(rev(magma(16)), alpha = 0.85), 
                  scales = list(draw = FALSE))

# Phylogenetic diversty (scaled)
pd_plot.scaled <- spplot(pd_poly,
                  zcol = "pd.obs",
                  xlim = c(-73.99045 , -28.99045), ylim = c(-33.72816 , 5.271841),
                  par.settings = list(fontsize = list(text = 21)),
                  at = seq(0, 1.000000001, length.out = 16),
                  colorkey = TRUE, 
                  sp.layout = list(biomes, fill = "gray90", lwd = 2), 
                  col.regions = makeTransparent(rev(magma(16)), alpha = 0.85), 
                  scales = list(draw = FALSE))

# Species richness
sr_plot <- spplot(pd_poly,
                  zcol = "ntaxa", 
                  xlim = c(-73.99045 , -28.99045), ylim = c(-33.72816 , 5.271841),
                  par.settings = list(fontsize = list(text = 21)),
                  at = seq(0, 55, length.out = 16),
                  colorkey = TRUE, 
                  sp.layout = list(list(biomes, fill = "gray90", lwd = 2)), 
                  col.regions = makeTransparent(rev(magma(16)), alpha = 0.85), 
                  scales = list(draw = FALSE))

# Species richness (scaled)
sr_plot.scaled <- spplot(pd_poly,
                  zcol = "ntaxa", 
                  xlim = c(-73.99045 , -28.99045), ylim = c(-33.72816 , 5.271841),
                  par.settings = list(fontsize = list(text = 21)),
                  at = seq(0, 1.000000001, length.out = 16),
                  colorkey = TRUE, 
                  sp.layout = list(list(biomes, fill = "gray90", lwd = 2)), 
                  col.regions = makeTransparent(rev(magma(16)), alpha = 0.85), 
                  scales = list(draw = FALSE))

# Correlation: phylogenetic diversity and species richness
corPdsr <- as.character(formatC(cor(pd_stats$ntaxa, 
                                    pd_stats$pd.obs, use = "na.or.complete"))) 
corPdsr_plot <- ggplot(data = pd_stats, mapping = aes(jitter(ntaxa), pd.obs))+
  geom_jitter()+
  labs(title = "Mimosa", subtitle = paste("r =", corPdsr))+
  xlab("Richness")+
  ylab("Phylogenetic diversity")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(face = "italic"))


# Standardized effect size
cols_ses <- c(viridis(16)[1:12],
              "khaki1",
              rev(heat.colors(12)[1:6]))
pdses_plot <- spplot(pd_poly,
                     zcol = "pd.obs.z",
                     xlim = c(-73.99045 , -28.99045), ylim = c(-33.72816 , 5.271841),
                     par.settings = list(fontsize = list(text = 21)),
                     at = seq(-5.21, 3.36, length.out = 16),
                     colorkey = TRUE, 
                     sp.layout = list(list(biomes, fill = "gray90")), 
                     col.regions = cols_ses, 
                     scales = list(draw = FALSE))

#Removing outliers: SES
sort(pd_poly$pd.obs.z)
cols_sesout <- c(viridis(16)[1:10],
                 "khaki1",
                 rev(heat.colors(16)[1:8]))
pdsesout_plot <- spplot(pd_poly,
                     zcol = "pd.obs.z",
                     xlim = c(-73.99045 , -28.99045), ylim = c(-33.72816 , 5.271841),
                     par.settings = list(fontsize = list(text = 21)),
                     at = seq(-3.85, 3.36, length.out = 19),
                     colorkey = TRUE, 
                     sp.layout = list(list(biomes, fill = "gray90")), 
                     col.regions = cols_sesout, 
                     scales = list(draw = FALSE))

# Residual values
cols_res <- c(viridis(20)[1:14], 
              "khaki1",
              rev(heat.colors(10)[1:4]))
pdres_plot <- spplot(pd_poly,
                     zcol = "residuals",
                     xlim = c(-73.99045 , -28.99045), ylim = c(-33.72816 , 5.271841),
                     par.settings = list(fontsize = list(text = 21)),
                     at = seq(-0.223, 0.077, length.out = 16),
                     colorkey = TRUE, 
                     sp.layout = list(list(biomes, fill = "gray90")), 
                     col.regions = cols_res, 
                     scales = list(draw = FALSE))

# Removing outliers: Residuals
sort(pd_poly$residuals)
cols_resout <- c(viridis(16)[1:10], 
              "khaki1",
              rev(heat.colors(10)[1:6]))
pdresout_plot <- spplot(pd_poly,
                     zcol = "residuals",
                     xlim = c(-73.99045 , -28.99045), ylim = c(-33.72816 , 5.271841),
                     par.settings = list(fontsize = list(text = 21)),
                     at = seq(-0.123, 0.077, length.out = 16),
                     colorkey = TRUE, 
                     sp.layout = list(list(biomes, fill = "gray90")), 
                     col.regions = cols_resout, 
                     scales = list(draw = FALSE))

#PD statistical significance
pdP_plot <- spplot(pd_poly,  
                   zcol = "significance", 
                   xlim = c(-73.99045 , -28.99045), ylim = c(-33.72816 , 5.271841),
                   colorkey = TRUE, 
                   sp.layout = list(list(biomes, fill = "gray90")), 
                   col.regions = c("firebrick4",
                                   "firebrick3",
                                    makeTransparent("khaki1", alpha = 0.5),
                                   "purple3",
                                   "purple4"), 
                   scales = list(draw = FALSE))
