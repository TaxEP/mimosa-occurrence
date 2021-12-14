#====================================================================================================#

library(tidyverse) 
library(data.table) 
library(CoordinateCleaner) 
library(flora) 
library(rgdal) 
library(raster)

library(RecordLinkage) # this package was missing from the script

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

mimosa_gbif <- fread(file = "dataset/occurrence.txt",
                     na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")


# Reducing data dimensionality by selecting only necessary columns
mimosa_gbif <- mimosa_gbif %>% dplyr::select(institutionCode,
                                             collectionCode,
                                             catalogNumber,
                                             genus,
                                             specificEpithet,
                                             infraspecificEpithet,
                                             basisOfRecord,
                                             identifiedBy,
                                             recordNumber,
                                             recordedBy,
                                             stateProvince,
                                             county,
                                             municipality,
                                             locality,
                                             decimalLongitude,
                                             decimalLatitude)

# Renaming attributes in order to match speciesLink column names
mimosa_gbif <- mimosa_gbif %>% rename("species" = specificEpithet,
                                      "institutioncode" = institutionCode,
                                      "collectioncode" = collectionCode,
                                      "catalognumber" = catalogNumber,
                                      "basisofrecord" = basisOfRecord,
                                      "identifiedby" = identifiedBy,
                                      "collector" = recordedBy,
                                      "collectornumber" = recordNumber,
                                      "stateprovince" = stateProvince,
                                      "longitude" = decimalLongitude,
                                      "latitude" = decimalLatitude,
                                      "subspecies" = infraspecificEpithet)

# Giving an unique ID number for each record
mimosa_gbif <- cbind(id = 1:nrow(mimosa_gbif), mimosa_gbif)

#=============#
# speciesLink #
#=============#

# Reading spLink (Warning message because of double quotes. Not a problem in this context)
mimosa_spLink <- fread(file = "dataset/speciesLink-20211213184635-0003085.txt", 
                       na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")



# Selecting important attributes
mimosa_spLink <- mimosa_spLink %>% dplyr::select(institutioncode,
                                                 collectioncode,
                                                 catalognumber,
                                                 genus,
                                                 species,
                                                 subspecies, 
                                                 basisofrecord,
                                                 identifiedby,
                                                 collector,
                                                 collectornumber,
                                                 stateprovince,
                                                 county,
                                                 locality,
                                                 longitude,
                                                 latitude)

# Coercing coords into numeric values (NA's are introduced by coercion in observations
# with coordinate values as 'Bloqueada')
mimosa_spLink$longitude <- as.numeric(as.character(mimosa_spLink$longitude))
mimosa_spLink$latitude <- as.numeric(as.character(mimosa_spLink$latitude))

# Giving an unique ID number for each record
mimosa_spLink <- cbind(id = (nrow(mimosa_gbif) + 1):(nrow(mimosa_gbif) + nrow(mimosa_spLink)), mimosa_spLink)

#======================================================================================#

#========================#
# CONCATENATING DATASETS #
#========================#

# Merging gbif and splink (165,758), and adding a column to define the original dataset
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

mimosa <- merge.with.source(x = mimosa_gbif,
                            y = mimosa_spLink,
                            name.x = "gbif",
                            name.y = "splink")

rm(merge.with.source, mimosa_gbif, mimosa_spLink) 

#======================================================================================#

#================#
# PRE-REFINEMENT #
#================#

# Removing duplicated registers given the institution code, the collection code and
# the catalog number (observations with NA for any of those attributes
# were not considered) (136,495)
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

# Indicating dataset specific columns
mimosa <- mimosa %>% rename("municipality_gbif" = municipality)

# Replacing 0 by NA in the coordinates 
mimosa$latitude[mimosa$latitude == 0] <- NA
mimosa$longitude[mimosa$longitude == 0] <- NA

#Removing registers without identifier name (58,930)
id_count <- data.frame(plyr::count(mimosa$identifiedby))
mimosa$identifiedby[mimosa$identifiedby %in% c("?", "-", "#", "##", "#?#", "//05/1987",
                                               "24/X/2013", 
                                               "<a href='https://bee.questagame.com/#/profile/34594?questagame_user_id=34594'>wwlearn (Low)</a>",
                                               "1/4/1983",
                                               "24/X/2013")
                    | mimosa$identifiedby == 0] <- NA
mimosa <- mimosa %>% filter(!is.na(identifiedby))

#Correcting states' names
province <- mimosa$stateprovince

lookup_states <- fread(file = "dataset/lookup_states.csv", na.strings = c("", NA),
                       encoding = "UTF-8")
#write.csv(lookup_states, file = "lists/lookup_states.csv", row.names = F)
get_states <- lookup_states$Incorreto
names(get_states) <- lookup_states$Correto
for(i in 1:length(province)){
  for(j in 1:length(get_states)){
    if(province[i] == unname(get_states[j]) & !is.na(province[i])){
      province[i] <- names(get_states[j])
    }
  }
}
mimosa$stateprovince <- province
rm(province, lookup_states, get_states, i, j, id_count)
mimosa$stateprovince[mimosa$stateprovince == "?" | mimosa$stateprovince == "-"] <- NA
plyr::count(mimosa$stateprovince) #Checking if everything has gone fine


# Filtering by Goias and Minas Gerais (26,362)
mimosa <- mimosa %>% filter(is.na(stateprovince) | stateprovince %in% c("Goias", 
                                                                        "Minas Gerais"))

# Removing records without species level identification (25,188)
mimosa <- mimosa %>% filter(!is.na(species))
mimosa <- mimosa %>% filter(!species %in% c("sp.", "sp1"))

# Removing records not based on preserved specimens and removing the 
# attribute 'basisofrecord' (30,617)
plyr::count(mimosa$basisofrecord)
mimosa <- mimosa %>% filter(basisofrecord %in% c("Exsic", "PRESERVED_SPECIMEN",
                                                 "PreservedSpecimen",
                                                 "s", "S"))
mimosa <- mimosa %>% dplyr::select(-basisofrecord)


#=====================================================================================================#

#==============================#
# CLEANING BY IDENTIFIER NAMES #
#==============================#

# Extracting a vector including determiners' names. It is preferable to use this vector instead of the data set itself because any changes are easily reversible. 
identifiedby <- as.character(mimosa$identifiedby)

# Ideally, it is better to start with a list of taxonomic specialists' compiled beforehand. Alternatively, as experts likely indentified the majority of samples from a given taxon, it is possible to infer specialists based on identification frequency. In this example we looked for specialists at the top of the list below. 
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# To improve accuracy, we confirmed if names in 'names.count' with at least 3 identifications were specialists by searching for taxonomic publications for the family of the focal group and authored by each name at Google Scholar. In this example, we searched for: allintitle: Leguminosae OR Mimosa author:"determiner".

# Next, based on the function 'replace.names', we standardized specialist's names. This is done in two iterations:
# (1) The first iteration returns, for manual evaluation, the automatically replaced names (names above the 'top' threshold) and names that are worth to be checked (names above the 'bottom' threshold but below the 'top' threshold).
# (2) In the second iteration, names that were erroneously replaced in the first iteration should be included in the argument 'not replace'. Likewise, names that were supposed to be replaced but were below the 'top' threshold should be included in the argument 'replace'.

# Because the procedure is iterative, the user should not change the vector 'identifiedBy' directly in the first iteration. For this purpose, we used a secondary vector ('identifiedBy_2'). After the second iteration, when everything should be set, the user must assign 'identifiedBy' as 'identifiedBy_2' before running the protocol for the following name. 

# At the end of the standardizing procedure, the following vector should contain all specialists' names
specialists <- c()

# RC Barneby
replace.by <- "RC Barneby"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("Rupert Charles Barneby; b.1911; d.2000; Barneby",
                                            "R. Barneby (NY) 1983", "Rupert Barneby",
                                            "confirm by Barneby, R. 1983", "R. Barneby 1986-1991",
                                            "R. BARNEBY/G. HATSCHBACH",
                                            "BArneby", "R.Barbeby", "Rupert C. Barneby",
                                            "R. Barneby (!RL2014)", "R. C. Barneby 1983",
                                            "Rupert Charles Barneby; b.1911; d.2000; Barneby",
                                            "barneby", "R. BARLEY", "R.BARNELY", "R. Barneby (!RL 2014)",
                                            "Barneby, ex. num. cit."))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# M Simon
replace.by <- "M Simon"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("M.Simon & L.M.Borges", "Marcelo Fragomeni Simon",
                                            "M. Simon, IN LIT. L. P. Queiroz", "M. Simon & L. M. Borges",
                                            "Simon M.F. & Marcelo F."))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# LM Borges
replace.by <- "LM Borges"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("I. M. Borges", "L.M.Borges (SPF)", "L.M. Borges et. al.;",
                                            "L. M. Borges & M. F. Simon", "L. M. Borges (SPF)",
                                            "L.M. Borges/13-05-2014", "L. M. Borges (SPF) 2013",
                                            "Leonardo M. Borges", "Fagg, CW; Borges, LM"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]


#Counting check
names.count <- as.data.frame(plyr::count(identifiedby))
names.count <- names.count[order(-names.count$freq), ]

#Replacing column
mimosa$identifiedby <- identifiedby

#Filtering (20,672)
specialists <- c("R. C. Barneby",
                 "L. M. Borges",
                 "L. P. de Queiroz",
                 "R. T. de Queiroz",
                 "J. S. Silva",
                 "M. S. Simon",
                 "G. P. Lewis",
                 "L. S. B. Jordão",
                 "A. P. Savassi-Coutinho",
                 "M. L. Guedes",
                 "O. S. Ribas",
                 "J. G. Nascimento",
                 "J. R. Pirani",
                 "I. B. Lima",
                 "M. Morales",
                 "C. W. Fagg",
                 "R. Grether",
                 "E. Córdula",
                 "A. Bocage",
                 "R. H. Fortunato",
                 "H. C. de Lima",
                 "A. Burkart",
                 "A. M. Miranda")
mimosa <- mimosa %>% filter(identifiedby %in% specialists)

# rm(identifiedby, specialists) ???

#===========================#
# CLEANING SCIENTIFIC NAMES #
#===========================#

#Generating a column for scientific names (without authors and including infraspecific epithet)
mimosa$gen_sp <- paste(mimosa$genus,
                       mimosa$species,
                       mimosa$subspecies,
                       sep = " ")


#Loading *.csv after manual correction 
#(???)
taxa_corrected <- read.csv("lists/Mimosa/taxa_corrected.csv", stringsAsFactors = F, 
                           na.strings = c("NA",""))

taxa_gensp <- tibble(gen = NA, sp = NA, var = NA, .rows = nrow(taxa_corrected))
for(i in 1:nrow(taxa_corrected)){
  str <- strsplit(taxa_corrected$corrected[i], split = " ")[[1]]
  if(length(str) == 2){
    taxa_gensp$gen[i] <- str[1]
    taxa_gensp$sp[i] <- str[2]
  } else if(length(str) == 4){
    taxa_gensp$gen[i] <- str[1]
    taxa_gensp$sp[i] <- str[2]
    taxa_gensp$var[i] <- paste(str[3], str[4], sep = " ")
  }
}
taxa_gensp$replace <- taxa_corrected$taxa




#Defining projection
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

pncv<-readOGR("shapefiles/PNCV/2conferido_Zoneamento_PNCV_10_03_2020.shp") #ICMBIO: https://www.gov.br/icmbio/pt-br/assuntos/biodiversidade/unidade-de-conservacao/unidades-de-biomas/cerrado/lista-de-ucs/parna-da-chapada-dos-veadeiros/parna-da-chapada-dos-veadeiros

pnsc<-readOGR("shapefiles/PNSC/parna_serra_do_cipo-polygon.shp") # ICMBIO: https://www.gov.br/icmbio/pt-br/assuntos/biodiversidade/unidade-de-conservacao/unidades-de-biomas/cerrado/lista-de-ucs/parna-da-serra-do-cipo/parna-da-serra-do-cipo
  
mun <- readOGR("shapefiles/BR/BR_UF_2020.shp") #IBGE: https://downloads.ibge.gov.br/downloads_geociencias.htm

#Projecting
proj4string(pnsc) <- crswgs84
#br <- spTransform(br, crswgs84)
mun <- spTransform(mun, crswgs84)

#plot(pncv)
#plot(pnsc)


