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

# Loading customized functions
source('functions.R')

#====================================================================================================#

#==================#
# READING DATASETS #
#==================#

#======#
# GBIF #
#======#

# Reading GBIF data (94,743)
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

# Converting the yearcollected field (integer) into string
# in order to match splink data
gbif$yearcollected <- as.character(gbif$yearcollected)

#=============#
# speciesLink #
#=============#

# Reading splink (Warning message because of double quotes. Not a problem here) (60,817)
splink <- fread(file = "datasets/splink.txt", 
                na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")

# Reducing data dimensionality by selecting only necessary columns
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

# Coercing coordinates into numeric values (NAs are introduced by coercion whenever coordinates
# are presented as strings)
splink$longitude <- as.numeric(as.character(splink$longitude))
splink$latitude <- as.numeric(as.character(splink$latitude))

# Giving an unique ID number for each record
splink <- cbind(id = (nrow(gbif) + 1):(nrow(gbif) + nrow(splink)), splink)

#======================================================================================#

#========================#
# CONCATENATING DATASETS #
#========================#

# Merging gbif and splink (155,560), and adding a column to define the original dataset
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

#================#
# PRE-REFINEMENT #
#================#

# Removing duplicated records (based on the institution code, the collection code and
# the catalog number) (observations with NA for any of those attributes
# are not considered) (119,549)
a <- mimosa
a.prime <- a[!is.na(a$institutioncode) &
               !is.na(a$collectioncode) &
               !is.na(a$catalognumber), ]
a.na <- a[is.na(a$institutioncode) |
            is.na(a$collectioncode) |
            is.na(a$catalognumber), ]
a <- unique(a.prime, by = c("institutioncode", 
                            "collectioncode",
                            "catalognumber"))
a <- rbind(a, a.na)
mimosa <- a
rm(a, a.prime, a.na)

# Indicating the source of the municipality attribute
mimosa <- mimosa %>% rename("municipality_gbif" = municipality)

# Replacing 0 by NA in the coordinates
# Even though the Equator line crosses the Brazilian territory, plain zero coordinates 
# are often unreliable (see Zizka et al., 2019)
mimosa$latitude[mimosa$latitude == 0] <- NA
mimosa$longitude[mimosa$longitude == 0] <- NA

# Removing records without determiner's name (81,988)
mimosa <- mimosa %>% filter(!is.na(identifiedby))

# Removing records without identification at species level (77,951)
mimosa <- mimosa %>% filter(!is.na(species))

# Removing records that are not based on preserved specimens and removing the 
# basisofrecord field (77,678)
plyr::count(mimosa$basisofrecord) # important to check the strings correspondent 
# to a preserved specimen basis of record
mimosa <- mimosa %>% filter(basisofrecord %in% c("Exsic", "PRESERVED_SPECIMEN",
                                                 "s", "S", "PreservedSpecimen"))
mimosa <- mimosa %>% dplyr::select(-basisofrecord)

# Removing records without valid coordinates (38,523)
mimosa <- mimosa[!is.na(mimosa$latitude), ]
mimosa <- mimosa[!is.na(mimosa$longitude), ]

#======================================================================================#

#===========================#
# CLEANING SCIENTIFIC NAMES #
#===========================#

# Generating a column for scientific names (without authors and including the infraspecific epithet)
mimosa$gen_sp <- paste(mimosa$genus,
                      mimosa$species,
                      mimosa$subspecies,
                      sep = " ")

# Extracting scientific names into a vector to work with
taxa <- plyr::count(mimosa$gen_sp)
taxa <- as.character(taxa$x)

# Standardizing the formatting of all names
for(i in 1:length(taxa)){
  taxa[i] <- gsub("NA", "", taxa[i]) # removing the character NA
  taxa[i] <- trimws(taxa[i]) # removing white space
  taxa[i] <- gsub(pattern = "  ", x = taxa[i], replacement = " ") # replacing double spaces by single spaces
}

# Suggesting names using the package 'flora' (and retrieving a few additional information that may be useful)
taxa_suggested <- get.taxa(taxa, vegetation.type = TRUE, 
                          habitat = TRUE, domain = TRUE, life.form = TRUE)

# Reorganizing columns in order to facilitate the manual checking procedure (see below)
taxa_suggested <- taxa_suggested[ , c(1, 2, 4, 5, 6, 8, 11, 12, 13, 14, 10, 9, 7, 3)]

# Adding a column dedicated to observations
taxa_suggested$obs <- NA

# Loading the checklist from the cleaning procedure that includes filtering records based on specialists' ids
taxa_corrected <- read.csv("lists/taxa_corrected.csv", stringsAsFactors = F, 
                           na.strings = c("NA",""), encoding = "UTF-8")

# Crossing current names with the checklist loaded above
for(i in 1:nrow(taxa_suggested)){
  if(taxa_suggested$original.search[i] %in% taxa_corrected$original.search){
    taxa_suggested$accepted.name[i] <- unique(taxa_corrected$accepted.name[taxa_corrected$original.search
                                                                    == taxa_suggested$original.search[i]])
  } else{
    taxa_suggested$accepted.name[i] <- NA
  }
}

# Writing *.csv for the manual checking procedure
#write.csv(taxa_suggested, file = "lists/taxa_suggested2.csv", row.names = F, fileEncoding = "UTF-8")

# Loading the checklist after the manual checking procedure
taxa_corrected <- read.csv("lists/taxa_corrected2.csv", stringsAsFactors = F, 
                           na.strings = c("NA",""), encoding = "UTF-8")

# Establishing a dataset with information on genus, species and varieties (or subspecies)
# The order of observations is the same as for the taxa_corrected dataset
taxa_gensp <- tibble(gen = NA, sp = NA, infra = NA, .rows = nrow(taxa_corrected))
for(i in 1:nrow(taxa_corrected)){
  str <- strsplit(taxa_corrected$accepted.name[i], split = " ")[[1]]
  if(length(str) == 2){
    taxa_gensp$gen[i] <- str[1]
    taxa_gensp$sp[i] <- str[2]
  } else if(length(str) == 4){
    taxa_gensp$gen[i] <- str[1]
    taxa_gensp$sp[i] <- str[2]
    taxa_gensp$infra[i] <- paste(str[3], str[4], sep = " ")
  }
}

# Adding a field that contains the original name for each corrected name (i.e., the names that should be replaced)
taxa_gensp$replace <- taxa_corrected$original.search

# Removing invalid taxa (27,146)
mimosa$gen_sp <- gsub("NA", "", mimosa$gen_sp) # removing NA strings from records of taxa without information regarding the infraspecific epithet
mimosa$gen_sp <- trimws(mimosa$gen_sp) # removing white spaces for the sake of correspondence
invalid_taxa <- taxa_gensp$replace[is.na(taxa_gensp$gen) | taxa_gensp$gen != "Mimosa"] # generating a list of invalid taxa
mimosa <- mimosa[!mimosa$gen_sp %in% invalid_taxa, ] # removing invalid taxa

# Correcting the dataset according to the taxa_gensp list
taxa_gensp <- taxa_gensp[!is.na(taxa_gensp$gen), ] # removing records that are NA for the gen field
for(i in 1:nrow(mimosa)){
  for(j in 1:nrow(taxa_gensp)){
    if(mimosa$gen_sp[i] == taxa_gensp$replace[j]){
      mimosa$gen_sp[i] <- paste(taxa_gensp$gen[j], taxa_gensp$sp[j], taxa_gensp$infra[j], sep = " ") # pasting all the taxa_gensp fields into one
    }
  }
}

# Removing NA strings from records of taxa without information regarding the infraspecific epithet
mimosa$gen_sp <- gsub("NA", "", mimosa$gen_sp)
mimosa$gen_sp <- trimws(mimosa$gen_sp) # removing white spaces

#======================================================================================#

#======================#
# CLEANING COORDINATES #
#======================#

# Flagging problematic record according to the 'CoordinateCleaner' package (26,396)
mimosa_coordFlagged <- mimosa %>% clean_coordinates(lon = "longitude",
                                                        lat = "latitude",
                                                        species = "gen_sp",
                                                        value = "flagged",
                                                        tests = c("equal", "gbif", 
                                                                  "institutions", 
                                                                  "outliers", "seas",
                                                                  "zeros"))

invalid_coords <- mimosa[mimosa_coordFlagged == FALSE, ] # subsetting flagged records
mimosa_clean <- mimosa[mimosa_coordFlagged  == TRUE, ] # subsetting valid records

#======================================================================================#

#============================#
# PREPARING FOR THE ANALYSES #
#============================#

# Reading the shapefile of the Brazilian terrestrial territory
br <- readOGR("shapefiles/BR/BR_UF_2020.shp")
grids_br <- readOGR("shapefiles/grids_br/grids_br.shp") 

# Projecting
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # EPSG:4326 - WGS 84
br <- spTransform(br, crswgs84)

# Intersecting coordinates with the grid cells
coords <- mimosa_clean
coordinates(coords) <- ~ longitude + latitude
proj4string(coords) <- crswgs84
coords <- over(coords, grids_br)
mimosa_clean$id_grid <- coords$id
mimosa_clean <- mimosa_clean %>% filter(!is.na(id_grid))
coordinates(mimosa_clean) <- ~ longitude + latitude
proj4string(mimosa_clean) <- crswgs84

# Writing *.csv for subsequent analyses 
write.csv(mimosa_clean, file = "datasets/mimosa-clean2.csv", row.names = F,
          fileEncoding = "UTF-8")
