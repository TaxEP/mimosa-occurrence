#=====================================================================================================#

library(rgdal)
library(recluster)
library(vegan)
library(dendextend)

# Set wd (Yago)
setwd("B:/yagob/GoogleDrive/Academia/Parallel_projects/Mimosa_occurrence")

# Set wd (Monique)
setwd("G:/.shortcut-targets-by-id/19Bt9xRgbQsy9ySW31FgR7E5aEscE0jKG/Mimosa_occurrence")

#=====================================================================================================#

#=======#
# INPUT #
#=======#

#======#
# TREE #
#======#

# Reading tree
mimosa_tree <- read.nexus("trees/mimosa-pruned_tree.nex")

#==========#
# Matrix #
#==========#

# Loading matrix
mimosa_matrix <- read.csv(file = "datasets/mimosa_matrix.csv", row.names = 1)

# Removing taxa that only have one recorded presence
mimosa_matrix <- mimosa_matrix[ , which(colSums(mimosa_matrix) > 1)]

# Removing from the matrix species that do not occur in the tree 
mimosa_matrix <- mimosa_matrix[ , colnames(mimosa_matrix) %in% mimosa_tree$tip.label]

# Removing empty sites 
mimosa_matrix <- mimosa_matrix[which(rowSums(mimosa_matrix) > 0), ]

#============#
# shapefiles #
#============#

#Loading grids and, the Brazilian terrestrial territory and the biomes
grids_br <- readOGR("shapefiles/grids_br/grids_br.shp") 
br <- readOGR("shapefiles/BR/BR_UF_2020.shp")
biomes <- readOGR("shapefiles/biomes/bioma_1milhao_uf2015_250mil_IBGE_albers_v4_revisao_pampa_lagoas.shp")

#Projecting br and biomes
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
# JACCARD #
#=========#

# Jaccard distance matrix
jaccard_matrix <- as.matrix(vegdist(mimosa_matrix, method = "jaccard", diag = TRUE))

#===============#
# Running UPGMA #
#===============#

upgma <- recluster.cons(mimosa_matrix, dist = "jaccard",
                        tr = 1000, p = 0.5, method = "average")
upgma_cons <- upgma$cons
upgma_cons <- di2multi(upgma_cons) # identifying polytomies
hc <- as.hclust(upgma$cons) # dendrogram

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
     cex = 0.5)

#====================#
# Dendrogram and plot #
#====================#

#Number of clusters
ncluster <- 10
dend <- as.dendrogram(hc) 

#Defining colors
colors <- makeTransparent(c("#41D91E", 
            "#FB9A99",
            "#003200",
            "#008805",
            "#E300F7",
            "#FFB559",
            "#FF0005",
            "#FF7F00",
            "#6A3D9A",
            "#1F78B4",
            "#FFFF99",
            "#0000A3",
            "#FAD900",
            "#B15928"), alpha = 0.7)

#Coloring clusters in the dendogram
dend <- color_branches(dend, k = ncluster, col = colors[1:ncluster])

#Indentifying the clusters of each grid and assigning them the respective color
groups_id <- data.frame("id" = labels(dend), "color" = get_leaves_branches_col(dend), 
                        "cluster_membership" = NA)
rb <- colors[1:ncluster]
names(rb) <- as.character(1:length(rb))
for(i in 1:nrow(groups_id)){
  groups_id$cluster_membership[i] <- names(rb)[unname(rb) == groups_id$color[i]]
}

#Merging everything into a polygon dataframe
jaccard_poly <- sp::merge(grids_br, groups_id, by.x = "id")
jaccard_poly$cluster_membership <- factor(jaccard_poly$cluster_membership, 
                                          levels = unique(groups_id$cluster_membership))

#=====================================================================================================#

#=========#
# FIGURES #
#=========#

# Plot
jaccard_plot <- spplot(jaccard_poly, 
                       zcol = "cluster_membership", 
                       xlim = c(-73.99045 , -28.99045), ylim = c(-33.72816 , 5.271841), 
                       colorkey = TRUE, 
                       sp.layout = list(list(biomes, fill = "gray90")), 
                       col.regions = colors[1:length(levels(jaccard_poly$cluster_membership))], 
                       scales = list(draw = FALSE))

# Coloring the dendrogram without transparency
dend_colors <- c("#41D91E", 
                 "#FB9A99",
                 "#003200",
                 "#008805",
                 "#E300F7",
                 "#FFB559",
                 "#FF0005",
                 "#FF7F00",
                 "#6A3D9A",
                 "#1F78B4",
                 "#FFFF99",
                 "#0000A3",
                 "#FAD900",
                 "#B15928")
dend <- as.dendrogram(hc) 
dend <- color_branches(dend, k = ncluster, col = dend_colors[1:ncluster])

# Dendogram
labels(dend) <- NULL
dend <- assign_values_to_branches_edgePar(dend = dend, value = 1, edgePar = "lwd")
plot_horiz.dendrogram(dend, axes = F)

