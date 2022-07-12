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

library(rgdal)
library(recluster)
library(dendextend)
library(ape)
library(picante)

# Set wd (Yago)
setwd("B:/yagob/GoogleDrive/Academia/Parallel_projects/Mimosa_occurrence")

# Set wd (Monique)
setwd("G:/.shortcut-targets-by-id/19Bt9xRgbQsy9ySW31FgR7E5aEscE0jKG/Mimosa_occurrence")

#======================================================================================================#

#=======#
# INPUT #
#=======#

#======#
# Tree #
#======#

# Reading tree
mimosa_tree <- read.nexus("trees/mimosa-pruned_tree2.nex")

#==========#
# Matrix #
#==========#

# Loading matrix
mimosa_matrix <- read.csv(file = "datasets/mimosa_matrix2.csv", row.names = 1)

# Removing taxa that only have one recorded presence
mimosa_matrix <- mimosa_matrix[ , which(colSums(mimosa_matrix) > 1)]

# Removing from the matrix species that do not occur in the tree 
mimosa_matrix <- mimosa_matrix[ , colnames(mimosa_matrix) %in% mimosa_tree$tip.label]

# Removing empty sites 
mimosa_matrix <- mimosa_matrix[which(rowSums(mimosa_matrix) > 0), ]

#============#
# shapefiles #
#============#

# Loading grids and, the Brazilian terrestrial territory and the biomes
grids_br <- readOGR("shapefiles/grids_br/grids_br.shp") 
br <- readOGR("shapefiles/BR/BR_UF_2020.shp")
biomes <- readOGR("shapefiles/biomes/bioma_1milhao_uf2015_250mil_IBGE_albers_v4_revisao_pampa_lagoas.shp")

# Projecting br and biomes
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
br <- spTransform(br, crswgs84)
biomes <- spTransform(biomes, crswgs84)

#=======#
# Grids #
#=======#

# Data frame with all grid ids
grids_df <- as.data.frame(grids_br@data)

#======================================================================================================#

#=========#
# UNIFRAC #
#=========#

# UniFrac matrix
unifrac_matrix <- as.matrix(unifrac(mimosa_matrix, mimosa_tree))

#===============#
# Running UPGMA #
#===============#

upgma <- recluster.cons(mimosa_matrix, mimosa_tree, dist = "unifrac",
                        tr = 1000, p = 0.5, method = "average")
upgma_cons <- upgma$cons
upgma_cons <- di2multi(upgma_cons) #identifying polytomies
hc <- as.hclust(upgma$cons) #dendrogram

# Fusion levels (useful to define the number of clusters)
plot(
  hc$height,
  nrow(mimosa_matrix):2,
  type = "S",
  main = "Fusion levels - Chord - UPGMA",
  ylab = "k (number of clusters)",
  xlab = "h (node height)",
  col = "grey"
)
text(hc$height,
     nrow(mimosa_matrix):2,
     nrow(mimosa_matrix):2,
     col = "red",
     cex = 0.8)

#====================#
# Dendogram and plot #
#====================#

# Number of clusters
ncluster <- 8
dend <- as.dendrogram(hc) 

# Defining colors
colors <- makeTransparent(c("#41D91E", 
            "#FB9A99",
            "#003200",
            "#008805",
            "#E300F7",
            "#B15928",
            "#FF0005",
            "#FF7F00",
            "#6A3D9A",
            "#1F78B4",
            "#FFFF99",
            "#0000A3",
            "#FAD900",
            "#FFB559"), alpha = 0.7)

# Coloring clusters in the dendogram
dend <- color_branches(dend, k = ncluster, col = colors[1:ncluster])

# Indentifying the clusters of each grid and assigning them the respective color
groups_id <- data.frame("id" = labels(dend), "color" = get_leaves_branches_col(dend), 
                        "cluster_membership" = NA)
rb <- colors[1:ncluster]
names(rb) <- as.character(1:length(rb))
for(i in 1:nrow(groups_id)){
  groups_id$cluster_membership[i] <- names(rb)[unname(rb) == groups_id$color[i]]
}

# Merging everything into a polygon dataframe
unifrac_poly <- sp::merge(grids_br, groups_id, by.x = "id")
unifrac_poly$cluster_membership <- factor(unifrac_poly$cluster_membership, 
                                          levels = unique(groups_id$cluster_membership))

#=====================================================================================================#

#=========#
# FIGURES #
#=========#

# Plot
unifrac_plot <- spplot(unifrac_poly, 
                       zcol = "cluster_membership", 
                       xlim = c(-73.99045 , -28.99045), ylim = c(-33.72816 , 5.271841),  
                       colorkey = TRUE, 
                       sp.layout = list(list(biomes, fill = "gray90")), 
                       col.regions = colors[1:length(levels(unifrac_poly$cluster_membership))], 
                       scales = list(draw = FALSE))

# Coloring the dendrogram without transparency
dend_colors <- c("#41D91E", 
                "#FB9A99",
                "#003200",
                "#008805",
                "#E300F7",
                "#B15928",
                "#FF0005",
                "#FF7F00",
                "#6A3D9A",
                "#1F78B4",
                "#FFFF99",
                "#0000A3",
                "#FAD900",
                "#FFB559")

# Coloring clusters in the dendrogram
dend <- as.dendrogram(hc) 
dend <- color_branches(dend, k = ncluster, col = dend_colors[1:ncluster])

# Dendogram
labels(dend) <- NULL
dend <- assign_values_to_branches_edgePar(dend = dend, value = 2, edgePar = "lwd")
plot_horiz.dendrogram(dend, axes = F)
