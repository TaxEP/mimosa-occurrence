#====================================================================================================#

library(tidyverse) 
library(data.table) 
library(CoordinateCleaner) 
library(flora) 
library(rgdal) 
library(raster)

#Loading functions
source('scripts/functions.R') 

#====================================================================================================#

#==================#
# READING DATASETS #
#==================#

#======#
# GBIF #
#======#

#Reading GBIF data
gbif <- fread(file = "datasets/Mimosa/0012888-190415153152247_gbif_mimosa/occurrence.txt",
              na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")