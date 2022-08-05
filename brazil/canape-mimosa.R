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
library(PDcalc) #https://rdrr.io/github/davidnipperess/PDcalc/man/phyloendemism.html

# Set wd (Yago)
setwd("B:/yagob/GoogleDrive/Academia/Parallel_projects/Mimosa_occurrence")

# Set wd (Monique)
setwd("G:/.shortcut-targets-by-id/19Bt9xRgbQsy9ySW31FgR7E5aEscE0jKG/Mimosa_occurrence")

#=====================================================================================================#

#=======#
# INPUT #
#=======#

# Reading matrix
mimosa_matrix <- read.csv(file = "datasets/mimosa_matrix-hc.csv", row.names = 1)

# Reading tree
mimosa_tree <- read.nexus("trees/mimosa-pruned_tree_hc.nex")

# Reading grids
grids_br <- readOGR("shapefiles/grids_br/hc_grids.shp") 

# Reading shapefile: Brazilian terrestrial territory
br <- readOGR("shapefiles/BR/BR_UF_2020.shp")

# Reading shapefile: biomes
biomes <- readOGR("shapefiles/biomes/bioma_1milhao_uf2015_250mil_IBGE_albers_v4_revisao_pampa_lagoas.shp")

#=====================================================================================================#

#====#
# PE #
#====#

# Removing from the matrix species that do not occur in the tree 
mimosa_matrix <- mimosa_matrix[ , colnames(mimosa_matrix) %in% mimosa_tree$tip.label]

# Calculating PE
pe <- data.frame("id" = rownames(mimosa_matrix), 
                 "pe" = phyloendemism(mimosa_matrix, mimosa_tree, weighted = T)) 

# PE randomization (independent swap)
rand.matrix.pe <- function(t, mat){
  phyloendemism(randomizeMatrix(mat, null.model = "independentswap"), t)[ , 1]
}
set.seed(7)
null.output <- replicate(999, rand.matrix.pe(mimosa_tree, mimosa_matrix))
rownames(null.output) <- rownames(mimosa_matrix)
ses.all <- (pe$pe - apply(null.output, MARGIN = 1, mean)) / apply(null.output,
                                                                  MARGIN = 1, sd)
p.val.all <- apply(cbind(pe$pe, null.output), MARGIN = 1, rank)[1, ] / 1000

# Concatenating results
pe$p.pe <- p.val.all
pe$ses.pe <- ses.all

# Significance
pe_poly <- merge(grids_br, pe, by.x = "id")
pe_poly$significance <- NA
for(i in 1:nrow(pe_poly@data)){
  if(is.na(pe_poly$p.pe[i])){
    pe_poly$significance[i] <- NA
  } else if(pe_poly$p.pe[i] >= 0.99 ){
    pe_poly$significance[i] <-  ">= 0.99"
  } else if(pe_poly$p.pe[i] >= 0.975 & pe_poly$p.pe[i] < 0.99){
    pe_poly$significance[i] <-  ">= 0.975"
  } else if(pe_poly$p.pe[i] <= 0.025 & pe_poly$p.pe[i] > 0.01){
    pe_poly$significance[i] <-  "<= 0.025"
  } else if(pe_poly$p.pe[i] <= 0.01){
    pe_poly$significance[i] <-  "<= 0.01"
  } else {
    pe_poly$significance[i] <- "Not significant"
  }
}
pe_poly$significance <- factor(pe_poly$significance, levels = c("<= 0.01", 
                                                                "<= 0.025",
                                                                "Not significant",
                                                                ">= 0.975",
                                                                ">= 0.99"))

# Plotting 

# PE
pe_plot <- spplot(pe_poly, zcol = "pe", colorkey = TRUE, 
                  sp.layout = list(list(biomes, fill = "gray90")), 
                  xlim = c(-73.99, -28.99), ylim = c(-33.72, 5.27),
                  col.regions =  makeTransparent(rev(magma(16)), alpha = 0.85),
                  scales = list(draw = FALSE))

# p-values
peP_plot <- spplot(pe_poly,  
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

#=====================================================================================================#

#=========================#
# PE: tree for comparison #
#=========================#

# Tree for comparison
comp_tree <- mimosa_tree
comp_tree$edge.length <- rep(1/length(comp_tree$edge.length), 
                             length(comp_tree$edge.length)) #equal branch lengths

# Calculating PE
pe_comp <- data.frame("id" = rownames(mimosa_matrix), 
                      "pe_comp" = phyloendemism(mimosa_matrix, comp_tree, weighted = T)) 

# PE randomization (independent swap)
rand.matrix.pe <- function(t, mat){
  phyloendemism(randomizeMatrix(mat, null.model = "independentswap"), t)[ , 1]
}
set.seed(7)
null.output <- replicate(999, rand.matrix.pe(comp_tree, mimosa_matrix))
rownames(null.output) <- rownames(mimosa_matrix)
ses.all <- (pe_comp$pe_comp - apply(null.output, MARGIN=1, mean)) / apply(null.output,
                                                                          MARGIN = 1, sd)
p.val.all <- apply(cbind(pe_comp$pe_comp, null.output), MARGIN = 1, rank)[1, ] / 1000

# Concatenating results
pe_comp$p.pe_comp <- p.val.all
pe_comp$ses.pe_comp <- ses.all

#=====================================================================================================#

#=====#
# RPE #
#=====#

# Calculating RPE
rpe <-merge(pe, pe_comp)
rpe$rpe <- rpe$pe/rpe$pe_comp

# RPE randomization (independent swap)
rand.matrix.rpe <- function(t, mat){
  comp_t <- t
  comp_t$edge.length <- rep(1/length(comp_t$edge.length), 
                            length(comp_t$edge.length))
  rand.mat <- randomizeMatrix(mat, null.model = "independentswap")
  phyloe <- phyloendemism(rand.mat, t)[ , 1]
  phyloe_comp <- phyloendemism(rand.mat, comp_t)[ , 1]
  phyloe/phyloe_comp  
}
set.seed(7)
null.output <- replicate(999, rand.matrix.rpe(mimosa_tree, mimosa_matrix))
rownames(null.output) <- rownames(mimosa_matrix)
ses.all <- (rpe$rpe - apply(null.output, MARGIN=1, mean)) / apply(null.output,
                                                                  MARGIN = 1, sd)
p.val.all <- apply(cbind(rpe$rpe, null.output), MARGIN = 1, rank)[1, ] / 1000

# Concatenating results
rpe$p.rpe <- p.val.all
rpe$ses.rpe <- ses.all

# Significance
rpe_poly <- merge(grids_br, rpe, by.x = "id")
rpe_poly$significance <- NA
for(i in 1:nrow(rpe_poly@data)){
  if(is.na(rpe_poly$p.rpe[i])){
    rpe_poly$significance[i] <- NA
  } else if(rpe_poly$p.rpe[i] >= 0.99 ){
    rpe_poly$significance[i] <-  ">= 0.99"
  } else if(rpe_poly$p.rpe[i] >= 0.975 & rpe_poly$p.rpe[i] < 0.99){
    rpe_poly$significance[i] <-  ">= 0.975"
  } else if(rpe_poly$p.rpe[i] <= 0.025 & rpe_poly$p.rpe[i] > 0.01){
    rpe_poly$significance[i] <-  "<= 0.025"
  } else if(rpe_poly$p.rpe[i] <= 0.01){
    rpe_poly$significance[i] <-  "<= 0.01"
  } else {
    rpe_poly$significance[i] <- "Not significant"
  }
}
rpe_poly$significance <- factor(rpe_poly$significance, levels = c("<= 0.01", 
                                                                "<= 0.025",
                                                                "Not significant",
                                                                ">= 0.975",
                                                                ">= 0.99"))

# Plotting

# RPE
rpe_plot <- spplot(rpe_poly, zcol = "rpe", colorkey = TRUE, 
                  sp.layout = list(list(biomes, fill = "gray90")), 
                  xlim = c(-73.99, -28.99), ylim = c(-33.72, 5.27),
                  col.regions =  makeTransparent(rev(magma(16)), alpha = 0.85),
                  scales = list(draw = FALSE))

rpeP_plot <- spplot(rpe_poly,  
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

#=====================================================================================================#

#========#
# CANAPE #
#========#

canape <- rpe
canape$step_one <- c()
for(i in 1:nrow(canape)){
  if(canape$p.pe[i] >= 0.95 | canape$p.pe_comp[i] >= 0.95){
    canape$step_one[i] <- TRUE
  } else{
    canape$step_one[i] <- FALSE
  }
}

canape$results <- NULL
for(i in 1:nrow(canape)){
  if(canape$step_one[i] == FALSE){
    canape$results[i] <- "Not significant"
  } else if(canape$p.rpe[i] >= 0.975 & canape$step_one[i] == TRUE){
    canape$results[i] <- "Paleo-endemism"
  } else if(canape$p.rpe[i] <= 0.025 & canape$step_one[i] == TRUE){
    canape$results[i] <- "Neo-endemism"
  } else if(canape$p.pe[i] >= 0.95 & canape$p.pe_comp[i] >= 0.95 & canape$step_one[i] == TRUE){
    canape$results[i] <- "Mixed endemism"
  }
  if(canape$p.pe[i] >= 0.99 & canape$p.pe_comp[i] >= 0.99 & canape$step_one[i] == TRUE){
    canape$results[i] <- "Super endemism"
  }
}  

# Plotting
canape_poly <- merge(grids_br, canape[ , which(colnames(canape) %in% c("id", "results"))], by.x = "id")
canape_poly$results <- factor(canape_poly$results, levels = c("Neo-endemism",
                                                              "Paleo-endemism",
                                                              "Not significant",
                                                              "Mixed endemism",
                                                              "Super endemism"))

#=====================================================================================================#

#=========#
# FIGURES #
#=========#

# We must know the extent of the pd_poly object in order to define plot limits
extent(canape_poly)

# CANAPE
canape_plot <- spplot(canape_poly, zcol = "results", colorkey = TRUE, 
                      sp.layout = list(list(biomes, fill = "gray90")), 
                      xlim = c(-73.99, -28.99), ylim = c(-33.72, 5.27),
                      col.regions =  c("#231151FF",
                                       "#B63679FF",
                                       makeTransparent("khaki1", alpha = 0.5),
                                       "#4AC16DFF",
                                       "#CFE11CFF"), scales = list(draw = FALSE))

#=====================================================================================================#

# Running CANAPE with another code (update R version if needed)

# install.packages("remotes")
remotes::install_github("joelnitta/canaper")
library(canaper)

# Preparing matrix and tree for the analysis
phylocom <- list(comm = mimosa_matrix, phy = mimosa_tree)

# Running a series of analyses that are necessary for CANAPE
set.seed(071421)
rand_test_results <- cpr_rand_test(phylocom$comm, phylocom$phy, null_model = "swap", n_iterations = 1000)

# Running CANAPE
canape <- cpr_classify_endem(rand_test_results)

# Assigning IDs
canape$id <- row.names(canape)

# Merging the results with the spatial polygons
canape_poly <- merge(grids_br, canape[ , which(colnames(canape) %in% c("id", "endem_type"))], by.x = "id")

# Ordering levels
canape_poly$endem_type <- factor(canape_poly$endem_type, levels = c("neo",
                                                              "paleo",
                                                              "not significant",
                                                              "mixed",
                                                              "super"))

# Plotting
canape_plot <- spplot(canape_poly, zcol = "endem_type", colorkey = TRUE, 
                      sp.layout = list(list(biomes, fill = "gray90")), 
                      xlim = c(-73.99, -28.99), ylim = c(-33.72, 5.27),
                      col.regions =  c("#231151FF",
                                       "#B63679FF",
                                       makeTransparent("khaki1", alpha = 0.5),
                                       "#4AC16DFF",
                                       "#CFE11CFF"), scales = list(draw = FALSE))
