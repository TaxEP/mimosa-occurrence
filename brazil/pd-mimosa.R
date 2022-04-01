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

# Reading presence-absence matrix
mimosa_matrix <- read.csv(file = "dataset/mimosa_matrix.csv", row.names = 1)

# Reading tree
mimosa_tree <- read.nexus("trees/mimosa-pruned_tree.nex")

# Reading grid cells
grids_br <- readOGR("shapefiles/grids_br/grids_br.shp") 

# Reading shapefile: Brazilian terrestrial territory
br <- readOGR("shapefiles/BR/BR_UF_2020.shp")

#=====================================================================================================#

#====#
# PD #
#====#

# Setting seed and running a phylogenetic diversity analysis
set.seed(7)
pd_stats <- ses.pd(mimosa_matrix, mimosa_tree, include.root = TRUE, null.model = "taxa.label")

# Creating a field that specifies grid cell id
pd_stats$id <- rownames(pd_stats)

# Running a linear regression and identifying residual values
pd_stats$residuals <- lm(pd.obs ~ ntaxa, pd_stats)$res

# Assigning results for each grid cell by merging results and a spatial polygon
pd_poly <- merge(grids_br, pd_stats, by.x = "id")

#=====================================================================================================#

#=========#
# FIGURES #
#=========#

# We can define label limits using the range of each attribute
summary(pd_poly@data)

# Phylogenetic diversty
pd_plot <- spplot(pd_poly,
                  zcol = "pd.obs",
                  xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5),
                  par.settings = list(fontsize = list(text = 21)),
                  at = seq(0, 0.581, length.out = 16),
                  colorkey = TRUE, 
                  sp.layout = list(list(br, fill = "gray")), 
                  col.regions = rev(magma(16)), 
                  scales = list(draw = FALSE))

# Species richness
sr_plot <- spplot(pd_poly,
                  zcol = "ntaxa", 
                  xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5),
                  par.settings = list(fontsize = list(text = 21)),
                  at = seq(0, 58, length.out = 16),
                  colorkey = TRUE, 
                  sp.layout = list(list(br, fill = "gray")), 
                  col.regions = rev(magma(16)), 
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
pdses_plot <- spplot(pd_poly,
                     zcol = "pd.obs.z",
                     xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5),
                     par.settings = list(fontsize = list(text = 21)),
                     at = seq(-3.661, 3.67, length.out = 16),
                     colorkey = TRUE, 
                     sp.layout = list(list(br, fill = "gray")), 
                     col.regions = viridis(16), 
                     scales = list(draw = FALSE))

# Residual values
pdres_plot <- spplot(pd_poly,
                     zcol = "residuals",
                     xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5),
                     par.settings = list(fontsize = list(text = 21)),
                     at = seq(-0.091, 0.091, length.out = 16),
                     colorkey = TRUE, 
                     sp.layout = list(list(br, fill = "gray")), 
                     col.regions = c(viridis(16)[1:7], "khaki1", rev(heat.colors(16)[1:8])), 
                     scales = list(draw = FALSE))