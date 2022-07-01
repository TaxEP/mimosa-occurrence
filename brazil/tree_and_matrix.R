#=====================================================================================================#

library(ape)
library(tidyverse)

# Set wd (Yago)
setwd("B:/yagob/GoogleDrive/Academia/Parallel_projects/Mimosa_occurrence")

# Set wd (Monique)
setwd("G:/.shortcut-targets-by-id/19Bt9xRgbQsy9ySW31FgR7E5aEscE0jKG/Mimosa_occurrence")

#=====================================================================================================#

#=======#
# INPUT #
#=======#

# Reading occurrence dataset
mimosa_br <- read.csv("datasets/mimosa-clean3.csv", na.strings = c("", NA))

# Reading tree
mimosa_tree <- read.tree("trees/mimosa_clean-VASCONCELOS2020.txt")

#=====================================================================================================#

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
data_labels <- data.frame("labels" =  mimosa_br$gen_sp, 
                          "std_labels" = NA)

# Standardizing tip labels
data_labels$std_labels <- gsub("var.", "", data_labels$labels) # removing 'var.'
data_labels$std_labels <- gsub("subsp.", "", data_labels$std_labels) # removing 'subsp.'
data_labels$std_labels <- gsub("  ", " ", data_labels$std_labels) # replacing double spaces by single spaces
data_labels$std_labels <- gsub(" ", "_", data_labels$std_labels) # replacing spaces by underline
data_labels$std_labels <- trimws(data_labels$std_labels) # removing white spaces

# Replacing names in the gen_sp field
mimosa_br$gen_sp <- data_labels$std_labels 

#=====================================================================================================#

#==============#
# PRUNING TREE #
#==============#

mimosa_pruned.tree <- drop.tip(mimosa_tree, mimosa_tree$tip.label[!mimosa_tree$tip.label %in% mimosa_br$gen_sp])

#==========================#
# RESCALING BRANCH LENGTHS # So results are comparable in the future
#==========================#

edge.l <- c()
for(i in 1:length(mimosa_pruned.tree$edge.length)){
  edge.l[i] <- mimosa_pruned.tree$edge.length[i]/sum(mimosa_pruned.tree$edge.length)
}

mimosa_pruned.tree$edge.length <- edge.l

# Saving as *.nex
write.nexus(mimosa_pruned.tree, file = "trees/mimosa-pruned_tree3.nex")

#=====================================================================================================#

#========#
# MATRIX #
#========#

# Creating a presence/absence matrix
mimosa_matrix <- matrix(data = NA, nrow = length(unique(mimosa_br$id_grid)), 
                        ncol = length(unique(mimosa_br$gen_sp)))
mimosa_matrix <- as.data.frame(mimosa_matrix)
colnames(mimosa_matrix) <- unique(mimosa_br$gen_sp)
rownames(mimosa_matrix) <- unique(mimosa_br$id_grid)
for(i in 1:nrow(mimosa_matrix)){
  for(j in 1:ncol(mimosa_matrix)){
    if(colnames(mimosa_matrix)[j] %in% mimosa_br$gen_sp[mimosa_br$id_grid == rownames(mimosa_matrix)[i]]){
      mimosa_matrix[i, j] <- 1
    } else {
      mimosa_matrix[i, j] <- 0  
    }
  }
}

# Saving matrix as *.csv
write.csv(mimosa_matrix, "datasets/mimosa_matrix3.csv", row.names = TRUE)
