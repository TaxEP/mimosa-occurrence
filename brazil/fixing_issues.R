#====================================================================================================#

library(tidyverse) 
library(data.table) 
library(CoordinateCleaner) 
library(flora) 
library(rgdal) 
library(raster)

# Set wd (Yago)
setwd("B:/yagob/GoogleDrive/Academia/Parallel_projects/Mimosa_occurrence")

# Set wd (Monique)
setwd("G:/.shortcut-targets-by-id/19Bt9xRgbQsy9ySW31FgR7E5aEscE0jKG/Mimosa_occurrence")

# Loading functions
source('functions.R')

#====================================================================================================#

#==================#
# READING DATASETS #
#==================#

#======#
# GBIF #
#======#

# Reading GBIF data
gbif <- fread(file = "datasets/gbif.txt",
              na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")

# Reducing data dimensionality by selecting only necessary columns
gbif <- gbif %>% dplyr::select(institutionCode,
                               collectionCode,
                               catalogNumber,
                               genus,
                               specificEpithet,
                               infraspecificEpithet,
                               basisOfRecord,
                               identifiedBy,
                               recordNumber,
                               recordedBy,
                               year,
                               stateProvince,
                               county,
                               municipality,
                               locality,
                               decimalLongitude,
                               decimalLatitude)

# Renaming attributes in order to match speciesLink column names
gbif <- gbif %>% rename("species" = specificEpithet,
                        "institutioncode" = institutionCode,
                        "collectioncode" = collectionCode,
                        "catalognumber" = catalogNumber,
                        "basisofrecord" = basisOfRecord,
                        "identifiedby" = identifiedBy,
                        "collector" = recordedBy,
                        "yearcollected" = year,  
                        "collectornumber" = recordNumber,
                        "stateprovince" = stateProvince,
                        "longitude" = decimalLongitude,
                        "latitude" = decimalLatitude,
                        "subspecies" = infraspecificEpithet)

# Giving an unique ID number for each record
gbif <- cbind(id = 1:nrow(gbif), gbif)

# Converting gbif$yearcollected (integer) into string
gbif$yearcollected <- as.character(gbif$yearcollected)

#=============#
# speciesLink #
#=============#

# Reading splink (Warning message because of double quotes. Not a problem in this context)
splink <- fread(file = "datasets/splink.txt", 
                na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")

# Selecting important attributes
splink <- splink %>% dplyr::select(institutioncode,
                                   collectioncode,
                                   catalognumber,
                                   genus,
                                   species,
                                   subspecies, 
                                   basisofrecord,
                                   identifiedby,
                                   collector,
                                   collectornumber,
                                   yearcollected,
                                   stateprovince,
                                   county,
                                   locality,
                                   longitude,
                                   latitude)

# Coercing coordinates into numeric values (NAs are introduced by coercion in observations
# with string coordinate values)
splink$longitude <- as.numeric(as.character(splink$longitude))
splink$latitude <- as.numeric(as.character(splink$latitude))

# Giving an unique ID number for each record
splink <- cbind(id = (nrow(gbif) + 1):(nrow(gbif) + nrow(splink)), splink)

#======================================================================================#

#========================#
# CONCATENATING DATASETS #
#========================#

# Merging gbif and splink (), and adding a column to define the original dataset
# for each observation
merge.with.source <- function(x, y, name.x = "X", name.y = "Y") {
  x.df <- cbind(x, datsrc.x = name.x)
  y.df <- cbind(y, datsrc.y = name.y)
  merged.df <- merge(x = x.df,
                     y = y.df,
                     all = TRUE)
  merged.df[is.na(merged.df$datsrc.x), "datsrc.x"] <- ""
  merged.df[is.na(merged.df$datsrc.y), "datsrc.y"] <- ""
  merged.df$datsrc <- paste(merged.df$datsrc.x, merged.df$datsrc.y, sep = "")
  merged.df$datsrc.x <- rm()
  merged.df$datsrc.y <- rm()
  return(merged.df)
}

mimosa <- merge.with.source(x = gbif,
                            y = splink,
                            name.x = "gbif",
                            name.y = "splink")

rm(merge.with.source, gbif, splink)

#======================================================================================#

#=========================================#
# INVESTIGATING THE "ADMINISTRATOR ISSUE" #
#=========================================#

# Where do records identified by "Administrador" come from?
administrador <- mimosa %>% filter(identifiedby %in% c("Administrador")) %>% dplyr::select(institutioncode, collectioncode,
                                                                                           catalognumber, 
                                                                                           species, identifiedby)

# (1) DO RECORDS IDENTIFIED BY ADMINISTRADOR OCCUR IN SPECIESLINK AS WELL?

# Reading splink (Warning message because of double quotes. Not a problem in this context)
splink <- fread(file = "datasets/splink.txt", 
                na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")

# Selecting important attributes
splink <- splink %>% dplyr::select(institutioncode,
                                   collectioncode,
                                   catalognumber,
                                   genus,
                                   species,
                                   subspecies, 
                                   basisofrecord,
                                   identifiedby,
                                   collector,
                                   collectornumber,
                                   yearcollected,
                                   stateprovince,
                                   county,
                                   locality,
                                   longitude,
                                   latitude)

adm_splink <-  splink %>% filter(catalognumber %in% administrador$catalognumber)
