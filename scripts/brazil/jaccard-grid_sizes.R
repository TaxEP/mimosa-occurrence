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

#======================================================================================================#

library(raster)
library(rgdal)
library(ape)
library(recluster)
library(vegan)
library(dendextend)
library(tidyverse)


# Set wd (Yago)
setwd("B:/yagob/GoogleDrive/Academia/Parallel_projects/Mimosa_occurrence")

# Reading the shapefile of the Brazilian terrestrial territory
br <- readOGR("shapefiles/BR/BR_UF_2020.shp")

# Reading biomes
biomes <- readOGR("shapefiles/biomes/bioma_1milhao_uf2015_250mil_IBGE_albers_v4_revisao_pampa_lagoas.shp")

# Projecting
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # EPSG:4326 - WGS 84
br <- spTransform(br, crswgs84)
biomes <- spTransform(biomes, crswgs84)

#=====================================================================================================#

#============#
# Sizes' loop #
#============#

sizes <- seq(2, 3, 0.1)
for(k in 1:length(sizes)){
  
  # Reading the dataset
  mimosa_clean <- read.csv(file = "datasets/mimosa-pre_grid.csv", na.strings = c("", NA),
                           encoding = "UTF-8")
  
  # Reading tree
  mimosa_tree <- read.tree("trees/mimosa_clean-VASCONCELOS2020.txt")
  
  # Defining grid cells' size
  grids_size <- c(sizes[k], sizes[k])
  
  # Creating spatial grid cells on the Brazilian terrestrial territory 
  grid <- raster(extent(br), resolution = grids_size, crs = proj4string(br))
  gridPolygon <- rasterToPolygons(grid)
  gridPolygon$id <- 1:nrow(gridPolygon)
  intersectGridClipped <- raster::intersect(gridPolygon, br)
  grids_br <- gridPolygon[gridPolygon$id %in% intersectGridClipped$id, ]
  
  # Intersecting coordinates with the grid cells
  coords <- mimosa_clean
  coordinates(coords) <- ~ longitude + latitude
  proj4string(coords) <- crswgs84
  coords <- over(coords, grids_br)
  mimosa_clean$id_grid <- coords$id
  mimosa_clean <- mimosa_clean %>% filter(!is.na(id_grid))
  coordinates(mimosa_clean) <- ~ longitude + latitude
  proj4string(mimosa_clean) <- crswgs84
  
  #=====================================================================================================#
  
  #=================#
  # TREE AND MATRIX #
  #=================#
  
  #===============#
  # STANDARDIZING #
  #===============#
  
  # We need to standardize both the tree and the dataset labels so they are consistent
  
  #
  # TREE
  #
  
  # Extracting tip labels
  tip_labels <- data.frame("tip_labels" =  mimosa_tree$tip.label, 
                           "std_labels" = NA)
  
  # Standardizing tip labels
  tip_labels$std_labels <- trimws(tip_labels$tip_labels) # removing white spaces
  
  # Replacing names in the tree
  mimosa_tree$tip.label <- tip_labels$std_labels
  
  #
  # DATASET
  #
  
  # Extracting tip labels
  data_labels <- data.frame("labels" =  mimosa_clean$gen_sp, 
                            "std_labels" = NA)
  
  # Standardizing tip labels
  data_labels$std_labels <- gsub("var.", "", data_labels$labels) # removing 'var.'
  data_labels$std_labels <- gsub("subsp.", "", data_labels$std_labels) # removing 'subsp.'
  data_labels$std_labels <- gsub("  ", " ", data_labels$std_labels) # replacing double spaces by single spaces
  data_labels$std_labels <- gsub(" ", "_", data_labels$std_labels) # replacing spaces by underline
  data_labels$std_labels <- trimws(data_labels$std_labels) # removing white spaces
  
  # Replacing names in the gen_sp field
  mimosa_clean$gen_sp <- data_labels$std_labels 
  
  #==============#
  # PRUNING TREE #
  #==============#
  
  mimosa_pruned.tree <- drop.tip(mimosa_tree, mimosa_tree$tip.label[!mimosa_tree$tip.label %in% mimosa_clean$gen_sp])
  
  #==========================#
  # RESCALING BRANCH LENGTHS # So results are comparable in the future
  #==========================#
  
  edge.l <- c()
  for(i in 1:length(mimosa_pruned.tree$edge.length)){
    edge.l[i] <- mimosa_pruned.tree$edge.length[i]/sum(mimosa_pruned.tree$edge.length)
  }
  
  mimosa_pruned.tree$edge.length <- edge.l
  
  mimosa_tree <- mimosa_pruned.tree
  
  #========#
  # MATRIX #
  #========#
  
  # Creating a presence/absence matrix
  mimosa_matrix <- matrix(data = NA, nrow = length(unique(mimosa_clean$id_grid)), 
                          ncol = length(unique(mimosa_clean$gen_sp)))
  mimosa_matrix <- as.data.frame(mimosa_matrix)
  colnames(mimosa_matrix) <- unique(mimosa_clean$gen_sp)
  rownames(mimosa_matrix) <- unique(mimosa_clean$id_grid)
  for(i in 1:nrow(mimosa_matrix)){
    for(j in 1:ncol(mimosa_matrix)){
      if(colnames(mimosa_matrix)[j] %in% mimosa_clean$gen_sp[mimosa_clean$id_grid == rownames(mimosa_matrix)[i]]){
        mimosa_matrix[i, j] <- 1
      } else {
        mimosa_matrix[i, j] <- 0  
      }
    }
  }
  
  #=========#
  # JACCARD #
  #=========#
  
  # Removing taxa that only have one recorded presence
  mimosa_matrix <- mimosa_matrix[ , which(colSums(mimosa_matrix) > 1)]
  
  # Removing from the matrix species that do not occur in the tree 
  mimosa_matrix <- mimosa_matrix[ , colnames(mimosa_matrix) %in% mimosa_tree$tip.label]
  
  # Removing empty sites 
  mimosa_matrix <- mimosa_matrix[which(rowSums(mimosa_matrix) > 0), ]
  
  #===============#
  # Running UPGMA #
  #===============#
  
  upgma <- recluster.cons(mimosa_matrix, dist = "jaccard",
                          tr = 1000, p = 0.5, method = "average")
  upgma_cons <- upgma$cons
  upgma_cons <- di2multi(upgma_cons) # identifying polytomies
  hc <- as.hclust(upgma$cons) # dendrogram
  
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
  
  png(paste("figs/jaccard-grid_sizes/plot", as.character(k), ".png", sep = ""),
      height = 4, width = 4, units = 'in', res=300)
  plot(jaccard_plot) 
  dev.off()
  png(paste("figs/jaccard-grid_sizes/dend", as.character(k), ".png", sep = ""),
      height = 4, width = 4, units = 'in', res=300)  
  labels(dend) <- NULL
  dend <- assign_values_to_branches_edgePar(dend = dend, value = 1, edgePar = "lwd")
  plot_horiz.dendrogram(dend, axes = F) 
  dev.off()
}

