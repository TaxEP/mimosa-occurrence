#====================================================================================================#

library(tidyverse) 
library(data.table) 
library(CoordinateCleaner) 
library(flora) 
library(rgdal) 
library(raster)
library(flora)
library(RecordLinkage) 

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
                                            "Barneby, ex. num. cit.","Rupert Charles Barneby"))
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
                                            "Simon M.F. & Marcelo F.", "Batalha, MA; Simon, MF", "M.A. Batalha & M. Simon"))
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

# VF Dutra
replace.by <- "VF Dutra"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("V.L.Dutra","V.F.Dutra (VIC)", "VALQUIRIA F DUTRA, 21-04-2020",
                                            "V.F.Dutra (VIC),","Valquiria Ferreira Dutra",
                                            "Valquíria Ferreira Dutra"
                                            ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

#AP Savassi-Coutinho
replace.by <- "AP Savassi-Coutinho"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

#GP Lewis
replace.by <- "GP Lewis"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("g. P. Lewis","G.P. LEWIS & J.S. PAGE"))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

#OS Ribas
replace.by <- "OS Ribas"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("O.S. Ribas; J. Cordeiro","O.S. Ribas; J. Cordeiro",
                                            "\"\"Ribas, O.S.; Cordeiro, J.\"\"","O. S. Ribas, J. Cordeiro",
                                            "Cordeiro, J.; Ribas, O.S.", "Cordeiro, J; Ribas, OS",
                                            "O.S. Ribas & J. Cordeiro", "O. S. Ribas,  J. Cordeiro",
                                            "O.S. Ribas & J. Codeiro", "O. dos S. Ribas"))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#JS Silva
replace.by <- "JS Silva"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Silva, M.J.","Silva, MJ","Silva, AF","Silva, J.B.",
                                                "Silva, R.R.","Silva, J.L." , "Silva, M.G. da",
                                                "Silva, T.C.","Silva, M.A. da Silva","Silva, M.A.",
                                                "Silva, CFS","Silva, J.P.","Silva, ACB","SILVA, R.R.",
                                                "Silva, J","Silva, N.T.", "Silva, DP","Silva, RR",
                                                "SILVEIRA, J.E."),
                                replace = c("Juliana Santos Silva (UEC),","J. Santos Silva (UEC) 2010-10-27",
                                            "Juliana Santos Silva (UEC)","J. S. Silva 2009-08-14",
                                            "J. S. Silva - UEC" ,"J.S. Silva (UEC),","Santos-Silva, J",
                                            "Santos, JS","Santos-Silva, J.","J. Santos Silva","J. Santos S.",
                                            "J. Santos S.,","J. S. Silva - UFRPe"
                                            
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#LSB Jordão
replace.by <- "LSB Jordão"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c( ),
                                replace = c("L. Jordão & H.C. Lima","fide Jordão et al.",
                                            "L. Jordão; H.C. Lima","Lucas, S.B. Jordão",
                                            "Lucas Josdão","Lucas Jordão & M. Barros"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

#LP Queiroz
replace.by <- "LP Queiroz"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c( ),
                                replace = c( "L. P. Queiroz (HUEFS) 2002-11","L.P.de Queiroz & R.M.Santos",
                                             "L. P. Queiroz (HUEFS) 2002-01","L. P. de Queiroz"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# A Burkart
replace.by <- "A Burkart"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c( ),
                                replace = c("Buskart","Bukart","A. E. Burkart; L. B. Smith",
                                            "Arturo Burkart", "A. Burkart (!RL 2014)", 
                                            "A. Macedo & A. Burkart"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# AF Silva
replace.by <- "AF Silva"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Silva, M.J.","Silva, T.C.","SILVA, R.R.",
                                                "Silva, J","Silva, N.T.","Silva, CFS",
                                                "Silva, J.P.", "Silva, ACB","V.F. Silva",
                                                "Silva, RR", "Silva, DP","Silva, M.G. da",
                                                "L.A. Silva","Silva, R.R.", "SILVEIRA, J.E.",
                                                "Silva, MJ","Silva, M.A. da Silva","Silva, J.B.",
                                                "Silva, M.A.","Silva, J.L." ),
                                replace = c("Starling, M.F.V & Silva, A.F."
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# CW Fagg
replace.by <- "CW Fagg"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# RG Grether
replace.by <- "RG Grether"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("R. Grether (MEXU)","R. Gether"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#JY Tamashiro
replace.by <- "JY Tamashiro"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("Jorge Yoshio Tamashiro","J.Y. Tamashiro & L.D. Meireles",
                                            "J.Y.T"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# R	Vanni
replace.by <- "R	Vanni"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("R.O. Vanni", "R.Vanni","Vanni, RO",
                                            "Vanni, R.","R. Vanni","Vanni","Vanni, R. O.","Vanni RO",
                                            "Vanni, R","R. O. Vanni", "Vanni, R. O., //"
                                            
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# M Flores-Cruz
replace.by <- "M Flores-Cruz"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# WA Rodrigues
replace.by <- "WA Rodrigues"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Rodrigues, T; Ramos, J", "Rodrigues, T"),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# RH Fortunato
replace.by <- "RH Fortunato"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("Reenée Fortunato","Reneé H. Fortunato","R. Fortunato (BAB) 2005-07-19",
                                            "Renée H. Fortunato" 
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# M Morales
replace.by <- "M Morales"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("M. Morols" ),
                                replace = c("Matias Morales" 
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# JR Pirani
replace.by <- "JR Pirani"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Pira, J.C."),
                                replace = c( 
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# FCP Garcia
replace.by <- "FCP Garcia"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c( 
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# B Mackinder
replace.by <- "B Mackinder"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c( 
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# GPE Rocha
replace.by <- "GPE Rocha"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c( 
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#R Liesner
replace.by <- "R Liesner"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c( "Ronald Liesner", "Ronald L. Liesner"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#M Atahuachi 
replace.by <- "M Atahuachi"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("M. Atachuachi B. (FHO)" 
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#RS Cowan

replace.by <- "RS Cowan"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("rs cowan 1973" 
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#EA Ulibarri

replace.by <- "EA Ulibarri"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# W Mantovani
replace.by <- "W Mantovani"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# G Seijo
replace.by <- "G Seijo"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# DM Glazier
replace.by <- "DM Glazier"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#ML Fonseca
replace.by <- "ML Fonseca"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# PO Rosa
replace.by <- "PO Rosa"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# AV Lopes
replace.by <- "AV Lopes"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Lopes, DV","Lopes, M.M.M.","LOPES, S.M.G."),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#G Ceccantini
replace.by <- "G Ceccantini"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#VC Souza
replace.by <- "VC Souza"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Souza, E.B.","Souza, F.S.","Souza, M.","Souza","Souza, M.L." ),
                                replace = c("V. C. Souza - UEG"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#C Proen�a
replace.by <- "C Proen�a"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c( "Incógnito; Proença, C."
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#J Fontella
replace.by <- "J Fontella"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#DMT Lins
replace.by <- "DMT Lins"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#F Filardi
replace.by <- "F Filardi"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#A Fernandes
replace.by <- "A Fernandes"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Fernandes, J.M.", "Fernandes, J.M" ),
                                replace = c("A. Fernandes & P. Dezerra" 
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#HC Lima
replace.by <- "HC Lima"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Lima, A.","Lima, D.A.", "Lima, M.P.","Lima, A.D.",
                                                "LIMA, A.C.R.","Lima, M.P.M.","Lima, MPM de","Lima, MPM",
                                                "Lima, L.C.P.","Lima, V.C." ),
                                replace = c("H. C. de Lima & L.F.G da Silva","H. C. de Lima & L. F. G. da Silva",
                                            "H.C. de Lima & Marli Pires", "H.C. Lima & C.M.J. Mattos",
                                            "Haroldo Cavalcante de Lima & Robson Daumas Ribeiro",
                                            "H. C. de Lima", "Haroldo Cavalcanti de Lima"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#HS Irwin
replace.by <- "HS Irwin"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#L Rico
replace.by <- "L Rico"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#C Romero
replace.by <- "C Romero"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Romero, R."),
                                replace = c("Carolina Romero (MO)"  
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#GP Lewis
replace.by <- "GP Lewis"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("Gwilym P. Lewis", "Ratter, JA; Lewis, GP","Nascimento, M.S.B.; Lewis, G.P." 
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# E Nunes
replace.by <- "E Nunes"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Nunes, SRDFS da","Nunes, SRDFS;"),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# A Ducke
replace.by <- "A Ducke"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# MS Oliveira
replace.by <- "MS Oliveira"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Oliveira, PA","Oliveira, CC"),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#JF Pastore
replace.by <- "JF Pastore"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#JP Silva
replace.by <- "JP Silva"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Silva, M.J.", "JS Silva", "Silva, M.G. da", "Silva, MJ",
                                                "Silva, T.C.", "Silva, ACB", "J.M. Silva", "SILVA, R.R.",
                                                "Silva, J.B.", "Silva, J","Silva, RR","Silva, R.R.",
                                                "Silva, M.A.", "Silva, N.T.", "SILVEIRA, J.E.", 
                                                "Silva, J.L.", "Silva, CFS"),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#M Sobral (M Sobral is Marcelo or Marcos Sobral?)
replace.by <- "M Sobral"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("Marcelo, F.S."
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#P Taubert
replace.by <- "P Taubert"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#VF Mansano
replace.by <- "VF Mansano"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("V.F. Mansano; R.B. Pinto"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#BAS Pereira
replace.by <- "BAS Pereira"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Pereira-Silva, G.","Pereira-Silva, G","Pereira, T.O.",
                                                "Pereira, M"),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# G Bentham
replace.by <- "G Bentham"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# GF Flores
replace.by <- "GF Flores"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#M Sousa
replace.by <- "M Sousa"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("M Sobral","Souza, M.L."),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#MVB Xavier
replace.by <- "MVB Xavier"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#JPM Brenan
replace.by <- "JPM Brenan"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#H Lorenzi
replace.by <- "H Lorenzi"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#HF Leitão-Filho
replace.by <- "HF Leitão-Filho"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#EP Heringer
replace.by <- "EP Heringer"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("E. P. Heringer 1979"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#GM Barroso
replace.by <- "GM Barroso"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c( "Barreto, M","Barreto, HLM","Barreto, M."),
                                replace = c("GRAZIELA","Graziela"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#APS Coutinho
replace.by <- "APS Coutinho"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("AP Savassi-Coutinho"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#AP Fortuna-Perez
replace.by <- "AP Fortuna-Perez"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#GM Antar
replace.by <- "GM Antar"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#G Siqueira
replace.by <- "G Siqueira"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Siqueira, L.C. de","J.C. Siqueira"),
                                replace = c("Geovane S. Siqueira"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#MPM Lima
replace.by <- "MPM Lima"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Lima, V.C.","Lima, L.C.P.","Lima, A.",
                                                "Lima, D.A.","Lima, A.D."),
                                replace = c("M.P. Morim","M. P. M. lima","Marli Pires Morim de Lima",
                                            "Morin de Lina, M.P.","Amorim, M.P.","Marli Pires Morim",
                                            "Marli Lima","Morim, MP","Morim de Lima, M.P.","M.M.P.M. Lima",
                                            "Marli P.Morim de Lima","M P M de Lima","MPMORIM",
                                            "Marli P. Morim de Lima/2003","M. P. M. de Lima"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#SRDFS Nunes
replace.by <- "SRDFS Nunes"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#RR Silva
replace.by <- "RR Silva"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Silva, M.J.", "Silva, M.G. da", "Silva, MJ", "Silva, T.C.",
                                                "Silva, ACB", "Silva, J.B.", "Silva, J", "Silva, M.A.",
                                                "Silva, N.T.", "Silva, J.L.", "Silva, CFS"),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#DA Neill
replace.by <- "DA Neill"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#JM Fernandes
replace.by <- "JM Fernandes"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("A Fernandes"),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#G Gardner
replace.by <- "G Gardner"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#RM Harley
replace.by <- "RM Harley"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#MG Caxambu
replace.by <- "MG Caxambu"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#G Malme
replace.by <- "G Malme"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("G. O. A. Malme"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#LD Meireles
replace.by <- "LD Meireles"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#RT Queiroz
replace.by <- "RT Queiroz"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("LP Queiroz"),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#CT Rizzini
replace.by <- "CT Rizzini"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("Rizzini; A. Mattos"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#VR Scalon
replace.by <- "VR Scalon"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#MCTB Messias
replace.by <- "MCTB Messias"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#LH Soares-Silva
replace.by <- "LH Soares-Silva"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#MSJ Trindade
replace.by <- "MSJ Trindade"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#AMSF Vaz
replace.by <- "AMSF Vaz"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("A. Vaz","Angela Maria Studart da Fonseca Vaz"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#A Pott
replace.by <- "A Pott"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#SA Harris 
replace.by <- "SA Harris"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#M Mendoza (considering as M. Bravo-Mendoza)
replace.by <- "M Mendoza"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c( "Mendonça, RC"),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#JNC Marchiori
replace.by <- "JNC Marchiori"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#JNC Marchiori
replace.by <- "JNC Marchiori"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#JGA Nascimento
replace.by <- "JGA Nascimento"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Nascimento, M.S.B.; Alencar, M.E.",
                                                "Nascimento, M.S.B.", "Nascimento"),
                                replace = c("J. G. A. do Nascimento","J.G.A.do Nascimento",
                                            "J.G.A. do Nascimento" 
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#G Pereira-Silva
replace.by <- "G Pereira-Silva"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c( "Pereira, T.O.","Pereira, M" ),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#DM Neves
replace.by <- "DM Neves"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Nee, M."),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#RE Gereau
replace.by <- "RE Gereau"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#FS Silveira
replace.by <- "FS Silveira"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Silva, MJ", "Silva, ACB", "Silva, J",
                                                "SILVEIRA, J.E.", "Silva, CFS"),
                                replace = c("Fernanda Schmidt Silveira","Schmidt-Silveira, F",
                                            "Schmidt-Silveira, F."
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#CFP Martius
replace.by <- "CFP Martius"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Martins, R.C.","Martins, BP"  ),
                                replace = c("C. F. P. von Martius"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#DNS Machado
replace.by <- "DNS Machado"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#AM Giulietti
replace.by <- "AM Giulietti"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#P Izaguirre
replace.by <- "P Izaguirre"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#LCP Lima
replace.by <- "LCP Lima"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Lima, D.A.", "Lima, A.D.", "Lima, V.C.", 
                                                "Lima, A.", "LIMA, A.C.R."),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#MA Silva
replace.by <- "MA Silva"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Silva, M.J.", "Silva, M.G. da", "Silva, MJ",
                                                "Silva, T.C.", "Silva, ACB", "J.M. Silva",
                                                "AF Silva", "L.A. Silva", "Silva, J.B.",
                                                "Silva, J", "Silva, N.T.", "SILVEIRA, J.E.", 
                                                "Silva, J.L.", "Silva, CFS"  ),
                                replace = c("M. A. da Silva"  
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#M Brandão
replace.by <- "M Brandão"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#TP Mendes
replace.by <- "TP Mendes"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#O Poncy
replace.by <- "O Poncy"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#RD Ribeiro
replace.by <- "RD Ribeiro"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Ribeiro, O."),
                                replace = c( "R.D.Ribeiro e H.C.de Lima"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#JS Santos 
replace.by <- "JS Santos"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Santos, F.F.M.", "Santos, J.R.", "Santos, A.P.M.",
                                                "Santos R.M.", "Santos, R.M.", "Santos, R. M" ),
                                replace = c("J. Santos S. (!RL 2014)"
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#LC Siqueira
replace.by <- "LC Siqueira"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("G Siqueira", "J.C. Siqueira"),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#MMF Silva
replace.by <- "MMF Silva"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("AF Silva"),
                                replace = c("Félix-da-Silva, M.M."
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#A Gentry
replace.by <- "A Gentry"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#JW Grimes
replace.by <- "JW Grimes"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c(
                                ))

identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)

names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#Counting check
names.count <- as.data.frame(plyr::count(identifiedby))
names.count <- names.count[order(-names.count$freq), ]

#Replacing column
mimosa$identifiedby <- identifiedby

#Filtering (20,672)
specialists <- c("RC Barneby",
                 "LM Borges",
                 "M Simon",
                 "VF Dutra",
                 "APS Coutinho",
                 "JS Silva",
                 "OS Ribas",
                 "LSB Jordão",
                 "GP Lewis",
                 "RH Fortunato",
                 "A Burkart",
                 "RG Grether",
                 "AF Silva",
                 "R	Vanni",
                 "M Morales",
                 "HC Lima",
                 "R Liesner",
                 "JY Tamashiro",
                 "CW Fagg",
                 "JR Pirani",
                 "ML Fonseca",
                 "FCP Garcia",
                 "EA Ulibarri",
                 "MPM Lima",
                 "C Proença",
                 "GPE Rocha",
                 "M Flores-Cruz",
                 "A Fernandes",
                 "WA Rodrigues",
                 "M Atahuachi",
                 "B Mackinder",
                 "M Sobral",
                 "RS Cowan",
                 "W Mantovani",
                 "PO Rosa",
                 "VC Souza",
                 "HS Irwin",
                 "G Seijo",
                 "J Fontella",
                 "DMT Lins",
                 "JGA Nascimento",
                 "GM Barroso",
                 "L Rico",
                 "P Taubert",
                 "BAS Pereira",
                 "DM Glazier",
                 "JF Pastore",
                 "AV Lopes",
                 "LH Soares-Silva",
                 "G Ceccantini",
                 "VF Mansano",
                 "GM Antar",
                 "A Ducke",
                 "CT Rizzini",
                 "FS Silveira",
                 "G Siqueira",
                 "RR Silva",
                 "C Romero",
                 "F Filardi",
                 "SRDFS Nunes",
                 "E Nunes",
                 "EP Heringer",
                 "G Gardner",
                 "G Pereira-Silva",
                 "G Bentham",
                 "MS Oliveira",
                 "RE Gereau",
                 "AMSF Vaz",
                 "DA Neill",
                 "GF Flores",
                 "HF Leitão-Filho",
                 "JPM Brenan",
                 "MVB Xavier",
                 "RD Ribeiro",
                 "SA Harris",
                 "VR Scalon",
                 "G Malme",
                 "H Lorenzi",
                 "LD Meireles",
                 "AP Fortuna-Perez",
                 "DM Neves",
                 "TP Mendes",
                 "A Pott",
                 "JS Santos",
                 "JW Grimes",
                 "M Mendoza",
                 "MCTB Messias",
                 "MG Caxambu",
                 "MSJ Trindade",
                 "O Poncy",
                 "RM Harley",
                 "AM Giulietti",
                 "CFP Martius",
                 "DNS Machado",
                 "LCP Lima",
                 "M Brandão",
                 "P Izaguirre"
                 )

mimosa <- mimosa %>% filter(identifiedby %in% specialists)

rm(identifiedby, specialists) 

#===========================#
# CLEANING SCIENTIFIC NAMES #
#===========================#

library(flora)

#Generating a column for scientific names (without authors and including infraspecific epithet)

mimosa$gen_sp <- paste(mimosa$genus,
                       mimosa$species,
                       mimosa$subspecies,
                       sep = " ")

#I added the following code from cleaning-calli

#Names
taxa <- plyr::count(mimosa$gen_sp)
taxa$taxa <- as.character(taxa$x)

#Removing NA (character derived from 'subspecies' attribute)
for(i in 1:length(taxa)){
  taxa[i] <- gsub("NA", "", taxa[i])
  taxa[i] <- trimws(taxa[i])
}

#Suggesting with flora (and retrieving a few additional information that may be useful)
taxa_suggested <-get.taxa(taxa, vegetation.type = TRUE, 
                          habitat = TRUE, domain = TRUE, life.form = TRUE)

#Writing *.csv for manual checking
write.csv(taxa_suggested, file = "taxa_suggested.csv", row.names = F)

#Loading *.csv after manual correction
taxa_corrected <- read.csv("taxa_corrected_mimosa.csv", stringsAsFactors = F, 
                           na.strings = c("NA",""))

#Establishing a data set with information on genus, species and varieties
taxa_gensp <- tibble(gen = NA, sp = NA, var = NA, .rows = nrow(taxa_corrected))
for(i in 1:nrow(taxa_corrected)){
  str <- strsplit(taxa_corrected$corrected.name[i], split = " ")[[1]]
  if(length(str) == 2){
    taxa_gensp$gen[i] <- str[1]
    taxa_gensp$sp[i] <- str[2]
  } else if(length(str) == 4){
    taxa_gensp$gen[i] <- str[1]
    taxa_gensp$sp[i] <- str[2]
    taxa_gensp$var[i] <- paste(str[3], str[4], sep = " ")
  }
}

#Original names and correspondent corrected names 
taxa_gensp$replace <- taxa_corrected$original.search

#Removing invalid taxa (6,564)
mimosa$gen_sp <- gsub("NA", "", mimosa$gen_sp)
mimosa$gen_sp <- trimws(mimosa$gen_sp)
invalid_taxa <- taxa_gensp$replace[is.na(taxa_gensp$gen) | taxa_gensp$gen != "Mimosa"]
mimosa <- mimosa[!mimosa$gen_sp %in% invalid_taxa, ]

#Correcting the data set (varieties suppressed)
taxa_gensp <- taxa_gensp[!is.na(taxa_gensp$gen), ]
for(i in 1:nrow(mimosa)){
  for(j in 1:nrow(taxa_gensp)){
    if(mimosa$gen_sp[i] == taxa_gensp$replace[j]){
      mimosa$gen_sp[i] <- paste(taxa_gensp$gen[j], taxa_gensp$sp[j], sep = " ")
    }
  }
}

rm(taxa, taxa_suggested, taxa_corrected, taxa_gensp, invalid_taxa, str)

#======================================================================================#

#================#
# SPLITTING DATA #
#================#

#Coords (6,159)
mimosa_coord <- mimosa %>% filter(!is.na(latitude) & !is.na(longitude))

#No coords (6,942)
mimosa_noCoord <- mimosa %>% filter(is.na(latitude) | is.na(longitude))

#======================#
# CLEANING COORDINATES #
#======================#

#Cleaning (3,796)
mimosa_coordFlagged <- mimosa_coord %>% clean_coordinates(lon = "longitude",
                                                        lat = "latitude",
                                                        species = "gen_sp",
                                                        value = "flagged",
                                                        tests = c("equal", "gbif", 
                                                                  "institutions", 
                                                                  "outliers", "seas",
                                                                  "zeros"))

invalid_coords <- mimosa_coord[mimosa_coordFlagged == FALSE, ]
mimosa_coordClean <- mimosa_coord[mimosa_coordFlagged  == TRUE, ]

#Binding invalid_coords to mimosa_noCoord (6,943)
mimosa_noCoord <- rbind(mimosa_noCoord, invalid_coords)

rm(invalid_coords, mimosa_coordFlagged, mimosa_coord)


# Checking occurrence of  previously selected species from my project

selected_species<- read.csv("dataset/mimosa_spp_phd2.csv", header=T, sep = ';')
colnames(selected_species)<-c("genus_and_species", "clade")

selected_species<-separate(data = selected_species, col = genus_and_species, into = c("genus", "species", "variety"), sep = "\\_") 

mimosa_coordClean$selected_species<- mimosa_coordClean$species %in% selected_species$species

mimosa_noCoord$selected_species<- mimosa_noCoord$species %in% selected_species$species

# Checking occurrence of species  present in the phylogeny
phylogeny_species<- read.csv("dataset/Mimosa_data_updated_140520-VASCONCELOS2020.csv", 
                             sep = ",")

phylogeny_species<-separate(data = phylogeny_species, col = cleaned_name, 
                            into = c("genus", "species", "variety"), sep = "\\_")

mimosa_coordClean$phylogeny_species<- mimosa_coordClean$species %in% phylogeny_species$species

mimosa_noCoord$phylogeny_species<- mimosa_noCoord$species %in% phylogeny_species$species

#======================================================================================================#

#------------------------#
# Ploting occurence maps #
#------------------------#


library(raster)
library(rgdal)
library(sf)

#Reading shapefiles: Brazilian terrestrial territory, Chapada dos Veadeiros, Serra do Cipo
br <- readOGR("shapefiles/BR/BR_UF_2020.shp", encoding = "UTF-8") # IBGE: https://downloads.ibge.gov.br/downloads_geociencias.htm
pncv <- readOGR("shapefiles/PNCV/2conferido_Zoneamento_PNCV_10_03_2020.shp", encoding = "UTF-8") # ICMBIO: https://www.gov.br/icmbio/pt-br/assuntos/biodiversidade/unidade-de-conservacao/unidades-de-biomas/cerrado/lista-de-ucs/parna-da-chapada-dos-veadeiros/parna-da-chapada-dos-veadeiros
pnsc <- readOGR("shapefiles/PNSC/parna_serra_do_cipo-polygon.shp", encoding = "UTF-8") # ICMBIO: https://www.gov.br/icmbio/pt-br/assuntos/biodiversidade/unidade-de-conservacao/unidades-de-biomas/cerrado/lista-de-ucs/parna-da-serra-do-cipo/parna-da-serra-do-cipo

# Defining projection
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Reprojecting
br <- spTransform(br, crswgs84)
pncv <- spTransform(pncv, crswgs84)
pnsc<- spTransform(pnsc, crswgs84)

# Visualizing 
plot(br)
plot(pncv, add = TRUE)
plot(pnsc, add = TRUE)

# GO and MG states
br@data
go_mg <- br[br@data$SIGLA_UF %in% c("GO", "MG"), ]

# Visualizing 
plot(go_mg)
plot(pncv, add = TRUE)
plot(pnsc, add = TRUE)

#Removing information that we do not need
pncv2 <- aggregate(pncv)

#PNCV map 
extent(pncv)
coords_pncv<- subset(mimosa_coordClean, longitude >= -47.90742 & longitude <= -46.97436 
                     & latitude >= -14.21634  & latitude <= -13.61845)
coordinates(coords_pncv) <- ~ longitude + latitude
proj4string(coords_pncv) <- crswgs84

plot(pncv)
points(coords_pncv)

#PNSC map
extent(pnsc)
coords_pnsc<- subset(mimosa_coordClean, longitude > -43.63072 & longitude < -43.45546 
                & latitude > -19.54183  & latitude < -19.2201)
coordinates(coords_pnsc) <- ~ longitude + latitude
proj4string(coords_pnsc) <- crswgs84

plot(pnsc)
points(coords_pnsc)

#======================================================================================================#

#--------------------#
# Exploring datasets #
#--------------------#

# Species occurence with coords

#PNCV

#Which and how many species record in CV? 
unique(coords_pncv@data$species, incomparables = F)

#Which and how many species of CV are present on previously selected species list?
unique(subset(coords_pncv@data, selected_species == "TRUE")$species, incomparables = F)

#Which and how many specis of CV are present in the phylogeny?
unique(subset(coords_pncv@data, phylogeny_species == "TRUE")$species, incomparables = F)

#What are the most abundant species? 
coords_pncv@data %>% count(species)

#What are the most abundant place?
coords_pncv@data %>% count(municipality_gbif)

#PNSC
#Which and how many species record in SC? 
unique(coords_pnsc@data$species, incomparables = F)

#Which and how many species of SC are present on previously selected species list?
unique(subset(coords_pnsc@data, selected_species == "TRUE")$species, incomparables = F)

#Which and how many specis of CV are present in the phylogeny?
unique(subset(coords_pnsc@data, phylogeny_species == "TRUE")$species, incomparables = F)

#What are the most abundant species? 
coords_pnsc@data %>% count(species)

#What are the most abundant place?
coords_pnsc@data %>% count(municipality_gbif)


# Species occurence with nocoords

#Cleaning by locality
mimosa_noCoord2<- mimosa_noCoord %>% select (species, subspecies, stateprovince, locality,municipality_gbif,
                           selected_species, phylogeny_species) %>% 
  filter (grepl('cip�|cipo|Cip�|Cipo|CIPO|CIP�|veadeiros|Veadeiros|VEADEIROS', locality))

#Which and how many species record? 
unique(mimosa_noCoord2$species, incomparables = F)

#Which and how many species are present on previously selected species list?
unique(subset(mimosa_noCoord2, selected_species == "TRUE")$species, incomparables = F)

#Which and how many specis are present in the phylogeny?
unique(subset(mimosa_noCoord2, phylogeny_species == "TRUE")$species, incomparables = F)

#What are the most abundant species? 
mimosa_noCoord2 %>% count(species)

#What are the most abundant place?
mimosa_noCoord2 %>% count(municipality_gbif)

#write.csv2(mimosa_coordClean, file = "lists/mimosa_coords.csv")
#write.csv2(mimosa_noCoord2, file = "lists/mimosa_noCoords.csv")



