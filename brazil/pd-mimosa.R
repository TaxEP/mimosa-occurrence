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

#=====================================================================================================#

#====#
# PD #
#====#

# Removing species that do not occur in the tree from the matrix
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

# We can define label limits using the range of each attribute
summary(pd_poly@data)

extent(pd_poly)

# Phylogenetic diversty
pd_plot <- spplot(pd_poly,
                  zcol = "pd.obs",
                  xlim = c(-73.99045 , -28.99045), ylim = c(-33.72816 , 5.271841),
                  par.settings = list(fontsize = list(text = 21)),
                  at = seq(0, 0.336, length.out = 16),
                  colorkey = TRUE, 
                  sp.layout = list(list(br, fill = "gray")), 
                  col.regions = rev(magma(16)), 
                  scales = list(draw = FALSE))

# Species richness
sr_plot <- spplot(pd_poly,
                  zcol = "ntaxa", 
                  xlim = c(-73.99045 , -28.99045), ylim = c(-33.72816 , 5.271841),
                  par.settings = list(fontsize = list(text = 21)),
                  at = seq(0, 55, length.out = 16),
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
                     xlim = c(-73.99045 , -28.99045), ylim = c(-33.72816 , 5.271841),
                     par.settings = list(fontsize = list(text = 21)),
                     at = seq(-0.223, 0.223, length.out = 16),
                     colorkey = TRUE, 
                     sp.layout = list(list(br, fill = "gray")), 
                     col.regions = c(viridis(16)[1:7], "khaki1", rev(heat.colors(16)[1:8])), 
                     scales = list(draw = FALSE))

cols3 <- c(viridis(16)[1], 
           makeTransparent(viridis(16)[2], alpha = .9),
           makeTransparent(viridis(16)[3], alpha = .8),
           makeTransparent(viridis(16)[4], alpha = .7),
           makeTransparent(viridis(16)[5], alpha = .6),
           makeTransparent(viridis(16)[6], alpha = .6),
           makeTransparent(viridis(16)[7], alpha = .6),
           makeTransparent(viridis(16)[8], alpha = .6),
           makeTransparent(viridis(16)[9], alpha = .6),
           makeTransparent(viridis(16)[10], alpha = .5),
           makeTransparent(viridis(16)[11], alpha = .5),
           makeTransparent(viridis(16)[12], alpha = .5),
           makeTransparent(viridis(16)[13], alpha = .5),
           makeTransparent(viridis(16)[14], alpha = .5),
           makeTransparent(viridis(16)[15], alpha = .5),
           makeTransparent("khaki1"),
           makeTransparent(rev(heat.colors(16)[5]), alpha = .6),
           makeTransparent(rev(heat.colors(16)[4]), alpha = .7),
           makeTransparent(rev(heat.colors(16)[3]), alpha = .8),
           makeTransparent(rev(heat.colors(16)[2]), alpha = .9),
           rev(heat.colors(16)[1]))

pdres_plot2 <- spplot(pd_poly,
                      zcol = "residuals",
                      xlim = c(-73.99045 , -28.99045), ylim = c(-33.72816 , 5.271841),
                      at = seq(-0.223, 0.077, length.out = 22),
                      par.settings=list(fontsize = list(text = 22)),
                      colorkey = TRUE, sp.layout = list(list(br, 
                                                             fill = "gray90")), 
                      col.regions = cols3, scales = list(draw = FALSE))

#PD statistical significance
pdP_plot <- spplot(pd_poly,  
                   zcol = "significance", 
                   xlim = c(-73.99045 , -28.99045), ylim = c(-33.72816 , 5.271841),
                   colorkey = TRUE, 
                   sp.layout = list(list(br, fill = "gray")), 
                   col.regions = c("firebrick4",
                                   "firebrick3",
                                    makeTransparent("khaki1", alpha = 0.5),
                                   "purple3",
                                   "purple4"), 
                   scales = list(draw = FALSE))
