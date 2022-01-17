library(ape)
library(tidyverse)

#=====================================================================================================#

#=======#
# INPUT #
#=======#

#Reading PNCV dataset
coords_PNCV <- read.csv("data/dataset/coords_PNCV.csv", na.strings = c("", NA))

#Reading tree
mimosa_tree <- read.tree("data/trees/mimosa_clean-VASCONCELOS2020.txt")

#=====================================================================================================#

#=====================#
# STANDARDIZING NAMES #
#=====================#

# Replacing space by underline in the dataset
coords_PNCV$gen_sp <- gsub(" ", "_", x = coords_PNCV$gen_sp)

# Extracting tip labels
tip_labels <- data.frame("tip_labels" =  mimosa_tree$tip.label, 
                         "std_labels" = NA)

# Standardizing tip labels (removing white spaces)
tip_labels$std_labels <- trimws(tip_labels$tip_labels)

# Replacing names in the tree
mimosa_tree$tip.label <- tip_labels$std_labels

#=================#
# COLLAPSING VARS #
#=================#

# Supressing infraspecific epithet
tip_labels$vars_sup <- NA
for(i in 1:nrow(tip_labels)){
  tip_labels[i, "vars_sup"] <- paste(strsplit(tip_labels[i, "std_labels"], "_")[[1]][c(1, 2)], 
                                     collapse = "_")
}
mimosa_tree$tip.label <- tip_labels$vars_sup

# Removing duplicated labels
duplicated_logical <- duplicated(mimosa_tree$tip.label)
mimosa_tree$tip.label <- as.character(1:length(mimosa_tree$tip.label))
for(i in 1:length(mimosa_tree$tip.label)){
  if(duplicated_logical[i] == TRUE)
    mimosa_tree <- drop.tip(mimosa_tree, mimosa_tree$tip.label[mimosa_tree$tip.label == as.character(i)])
}
mimosa_tree$tip.label <- tip_labels$vars_sup[-which(duplicated_logical)]

#=====================================================================================================#

#==============#
# PRUNING TREE #
#==============#

mimosa_pruned.tree <- drop.tip(mimosa_tree, mimosa_tree$tip.label[!mimosa_tree$tip.label %in% coords_PNCV$gen_sp])

#=======#
# NEXUS #
#=======#

write.nexus(mimosa_pruned.tree, file = "data/trees/pruned_tree-mimosa_PNCV.nex")

#=====================================================================================================#

#========#
# MATRIX #
#========#

# Composing a presence/absence matrix
mimosa_matrix_PNCV <- matrix(data = NA, nrow = length(unique(coords_PNCV$id_grid)), 
                        ncol = length(unique(coords_PNCV$gen_sp)))
mimosa_matrix_PNCV <- as.data.frame(mimosa_matrix_PNCV)
colnames(mimosa_matrix_PNCV) <- unique(coords_PNCV$gen_sp)
rownames(mimosa_matrix_PNCV) <- unique(coords_PNCV$id_grid)
for(i in 1:nrow(mimosa_matrix_PNCV)){
  for(j in 1:ncol(mimosa_matrix_PNCV)){
    if(colnames(mimosa_matrix_PNCV)[j] %in% coords_PNCV$gen_sp[coords_PNCV$id_grid == rownames(mimosa_matrix_PNCV)[i]]){
      mimosa_matrix_PNCV[i, j] <- 1
    } else {
      mimosa_matrix_PNCV[i, j] <- 0  
    }
  }
}

#Saving matrix
write.csv(mimosa_matrix_PNCV, "data/dataset/mimosa_matrix_PNCV.csv", row.names = TRUE)
