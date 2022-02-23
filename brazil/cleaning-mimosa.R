#====================================================================================================#

library(tidyverse) 
library(data.table) 
library(CoordinateCleaner) 
library(flora) 
library(rgdal) 
library(raster)

#Loading functions
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
mimosa$latitude[mimosa$latitude == 0] <- NA
mimosa$longitude[mimosa$longitude == 0] <- NA

# Removing records without determiner's name (81,988)
mimosa <- mimosa %>% filter(!is.na(identifiedby))

# Removing records without identification at the species level (77,951)
mimosa <- mimosa %>% filter(!is.na(species))

# Removing records that are not based on preserved specimens and removing the 
# attribute 'basisofrecord' (77,678)
# plyr::count(mimosa$basisofrecord)
mimosa <- mimosa %>% filter(basisofrecord %in% c("Exsic", "PRESERVED_SPECIMEN",
                                                 "s", "S", "PreservedSpecimen"))
mimosa <- mimosa %>% dplyr::select(-basisofrecord)

# Removing records without coordinates (38,652)
mimosa <- mimosa %>% filter(!is.na(latitude) | !is.na(longitude))

#======================================================================================#

#=================================#
# CLEANING BY IDENTIFIER NAMES #
#=================================#

# Creating a vector to work with
identifiedby <- mimosa$identifiedby

# Extracting a vector including determiners' names. It is preferable to use this vector 
# instead of the data set itself because any changes are easily reversible. 
identifiedby <- as.character(mimosa$identifiedby)

# Ideally, it is better to start with a list of taxonomic specialists' compiled beforehand. 
# Alternatively, as experts likely indentified the majority of samples from a given taxon, 
# it is possible to infer specialists based on identification frequency. Here we looked for 
# specialists at the top of the list below. 
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# To improve accuracy, we confirmed if names in 'names.count' with at least 3 identifications
# were specialists by searching for taxonomic publications for the family of the focal group
# and authored by each name at Google Scholar. Here, we searched for: 
# allintitle: Leguminosae OR Mimosa author:"determiner".

# Next, based on the function 'replace.names', we standardized specialist's name. 
# This is done in two iterations:
# (1) The first iteration returns, for manual evaluation, the automatically replaced names
# (names above the 'top' threshold) and names that are worth to be checked 
# (names above the 'bottom' threshold but below the 'top' threshold).
# (2) In the second iteration, names that were erroneously replaced in the first iteration
# should be included in the argument 'not replace'. Likewise, names that were supposed to 
# be replaced but were below the 'top' threshold should be included in the argument 'replace'.

# Because the procedure is iterative, the user should not change the vector 'identifiedBy' 
# directly in the first iteration. For this purpose, we used a secondary vector 
# ('identifiedBy_2'). After the second iteration, when everything should be set, 
# the user must assign 'identifiedBy' as 'identifiedBy_2' before running the protocol 
# for the following name. 

# At the end of the standardizing procedure, the following vector should contain all specialists' names
specialists <- c()

# RC Barneby
replace.by <- "RC Barneby"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("R. Bameby", "R. C. Barneby 1983",
                                            "Rupert Charles Barneby; b.1911; d.2000; Barneby",
                                            "Simon, MF; Barneby, RC", "R. Barneby (NY) 1983",
                                            "Barneby, ex. num. cit.", "Vanni, R; Barneby",
                                            "Marchiori, J.N.; Barneby, R.C.", "R. BernebY"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# M Simon