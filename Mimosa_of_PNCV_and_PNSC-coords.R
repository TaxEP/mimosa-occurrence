#====================================================================================================#

library(tidyverse) 
library(data.table) 
library(CoordinateCleaner) 
library(flora) 
library(rgdal) 
library(raster)

library(RecordLinkage) #this package was missing from the script

#Loading functions
source('functions.R') 

#====================================================================================================#

#==================#
# READING DATASETS #
#==================#

#======#
# GBIF #
#======#

#Reading GBIF data

mimosa_gbif <- fread(file = "dataset/occurrence.txt",
                     na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")


#Reducing data dimensionality by selecting only necessary columns
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

#Renaming attributes in order to match speciesLink column names
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

#Giving an unique ID number for each record
mimosa_gbif <- cbind(id = 1:nrow(mimosa_gbif), mimosa_gbif)

#=============#
# speciesLink #
#=============#

#Reading spLink (Warning message because of double quotes. Not a problem in this context)
mimosa_spLink <- fread(file = "dataset/speciesLink-20211213184635-0003085.txt", 
                       na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")



#Selecting important attributes
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

#Coercing coords into numeric values (NA's are introduced by coercion in observations
#with coordinate values as 'Bloqueada')
mimosa_spLink$longitude <- as.numeric(as.character(mimosa_spLink$longitude))
mimosa_spLink$latitude <- as.numeric(as.character(mimosa_spLink$latitude))

#Giving an unique ID number for each record
mimosa_spLink <- cbind(id = (nrow(mimosa_gbif) + 1):(nrow(mimosa_gbif) + nrow(mimosa_spLink)), mimosa_spLink)

#======================================================================================#

#========================#
# CONCATENATING DATASETS #
#========================#

#Merging gbif and splink (116,987), and adding a column to define the original dataset
#for each observation
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

# rm(merge.with.source, mimosa_gbif, mimosa_spLink) ???

#======================================================================================#

#================#
# PRE-REFINEMENT #
#================#

#Removing duplicated registers given the institution code, the collection code and
#the catalog number (observations with NA for any of those attributes
#were not considered) (89,918)
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
# rm(a, a.prime, a.na) ?????

#Indicating dataset specific columns
mimosa <- mimosa %>% rename("municipality_gbif" = municipality)

#Replacing 0 by NA in the coordinates 
mimosa$latitude[mimosa$latitude == 0] <- NA
mimosa$longitude[mimosa$longitude == 0] <- NA

#Removing registers without identifier name (58,930)
mimosa$identifiedby[mimosa$identifiedby %in% c("?", "-", "#", "##", "//05/1987",
                                               "24/X/2013") | mimosa$identifiedby == 0] <- NA
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
# rm(province, lookup_states, get_states, i, j) ???
mimosa$stateprovince[mimosa$stateprovince == "?" | mimosa$stateprovince == "-"] <- NA
plyr::count(mimosa$stateprovince) #Checking if everything has gone fine


#Filtering by states in which the campos rupestres occur 
#according to Silveira et al (2016) (31,575) 
#I'm using just GO and MG
mimosa <- mimosa %>% filter(is.na(stateprovince) | stateprovince %in% c("Goias", 
                                                                        "Minas Gerais" 
                                                                        ))

#Removing records without species level identification (30,750)
mimosa <- mimosa %>% filter(!is.na(species))
mimosa <- mimosa %>% filter(!species %in% c("sp.", "sp1"))

#Removing records not based on preserved specimens and removing the 
#attribute 'basisofrecord' (30,617)
#plyr::count(mimosa$basisofrecord)
mimosa <- mimosa %>% filter(basisofrecord %in% c("Exsic", "PRESERVED_SPECIMEN",
                                                 "s", "S"))
mimosa <- mimosa %>% dplyr::select(-basisofrecord)


#=====================================================================================================#

#=================================#
# CLEANING BY IDENTIFIER NAMES #
#=================================#

#Creating a vector to work with
identifiedby <- mimosa$identifiedby

#Generating a list of identifiers ordered by frequency of identifications (more identifications
#is supposed to be related with a better identification)
#Counting check
#names.count <- as.data.frame(plyr::count(identifiedby))
#names.count <- names.count[order(-names.count$freq), ]

#R. C. Barneby
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = c("R. C. Barneby", "R. C. BARNEBY",
                                           "Barneby, R. C", "BARNEBY, R. C."),
                              replace.by = "R. C. Barneby",
                              not.replace = c("R. Bameby", "R. M. Harley"),
                              replace = c("R. Barneby (NY)", "Rupert Barneby",
                                          "BArneby", "R. BARNEBY/G. HATSCHBACH",
                                          "R. Barneby 85", "R. Barneby (!RL 2014)",
                                          "R. Barneby, (!RL 2014)", "R. Barneby (!RL2014)",
                                          "R. Barneby 84", "R. Barneby 93", "Rupert C. Barneby",
                                          "Rupert Charles Barneby", "R Brandeby", "R. Barneby 1983",
                                          "R. Barneby (NY) 1983", "R. Barneby 1986-1991", 
                                          "Rupert Charles Barneby; b.1911; d.2000; Barneby", "R. Bameby"))

#L. M. Borges
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = c("L. M. Borges", "Borges, L. M.",
                                           "L. M. BORGES", "BORGES, L. M."),
                              replace.by = "L. M. Borges",
                              replace = c("Leonardo M. Borges", "L.M. Borges et. al.;",
                                          "L.M. Borges (SPF)", "L.M.Borges (SPF)", 
                                          "L.M. Borges/13-05-2014"),
                              not.replace = c("L. M. G. Nogueira"))

#M. F. Simon
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = c("M. F. Simon", "Simon, M. F.",
                                           "M. F. SIMON", "SIMON, M. F."),
                              replace.by = "M. F. Simon",
                              not.replace = c("M.S. Simo"),
                              replace = c("Marcelo F. Simon", "M. Simon, IN LIT. L. P. Queiroz",
                                          "M. Simon & L. M. Borges", "Simon, MF; Queiroz, LP",
                                          "Marcelo F Simon", "Marcelo Simon","Simon, MF; Barneby, RC"),
                              return = TRUE)

#L. P. de Queiroz
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = c("L. P. de Queiroz", "Queiroz, L. P.",
                                           "L. P. DE QUEIROZ", "QUEIROZ, L. P."),
                              replace.by = "L. P. de Queiroz",
                              replace = c("Luciano Paganucci de Queiroz", "L.P.de Queiroz & E.R.de Souza",
                                          "L.P.de Queiroz & D.Cardoso", "L. P. Queiroz (HUEFS) 2002-01",
                                          "L.P.de Queiroz & R.M.Santos", "L.P.de Queiroz & T.S.Nunes",
                                          "D. Queiroz", "L. P. Queiroz (HUEFS) 2002-11",
                                          "Santos, R.M. dos; Queiroz, L.P. de", "Santos, RM dos; Queiróz, LP de",
                                          "Santos, RM; Queiroz, LP de","L. Paganucci", "Queiroz, LP de; Lima, MPM de"),
                              not.replace = c("Queiroz, R.T.", "Queiroz, RT"),
                              return = T)

#R. T. de Queiroz
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = c("R. T. de Queiroz", "Queiroz, R. T",
                                           "R. T. DE QUEIROZ", "QUEIROZ, R. T."),
                              replace.by = "R. T. de Queiroz",
                              not.replace = c("L. P. de Queiroz"),
                              replace = c("R.T, Queiroz"),
                              return = T)

#V. F. Dutra
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = c("V. F. Dutra", "V. F. DUTRA",
                                           "Dutra, V. F.", "DUTRA, V. F."),
                              replace.by = "V. F. Dutra",
                              replace = c("V.F.Dutra (VIC)", "V.F.Dutra (VIC), (!RL 2014)",
                                          "Valquíria Ferreira Dutra", "Valquiria Ferreira Dutra"),
                              return = T)

#J. S. Silva
identifiedby <- replace.names(x = identifiedby, top = 0.9, bottom = 0.7, 
                              check.by = c("J. S. Silva", "Silva, J. S",
                                           "J. S. SILVA", "SILVA, J. S."),
                              replace.by = "J. S. Silva",
                              not.replace = c("Silva, J.P.", "SILVA, R.R.", "Silva, MS",
                                              "Silva, J.L.", "Silva, P.", "L. S. Silva",
                                              "Silva, J.B.", "Silva, R.", "Silva, MJ" ),
                              replace = c("J.S. Silva (2)", "Silva, JS; Nascimento, JGA",
                                          "J. S. Silva 2009-08-14", "J. Santos Silva (UEC) 2010-10-27",
                                          "J.S. Silva; M. Sales", "Silva, J.S.; Sales, M.", "Santos, J.S.",
                                          "J.S. Silva (UEC)", "J. Santos S.", "J. Santos Silva", "J. Santos S. (!RL 2014)",
                                          "Juliana Santos Silva (UEC)", "J.S.Silva/14-VIII-2009",
                                          "Juliana Santos Silva", "Juliana Santos Silva; Universidade Estadual de Campinas"),
                              return = T)

#G. P.  Lewis
identifiedby <- replace.names(x = identifiedby, top = 0.9, bottom = 0.7, 
                              check.by = c("G. P. Lewis", "Lewis, G. P.",
                                           "G. P. LEWIS", "LEWIS, G. P."),
                              replace.by = "G.  P. Lewis",
                              replace = c("Lewis, G.P.; Page, J.S.", "G.P. LEWIS & J.S. PAGE",
                                          "G,P LEWIS", "Lewis, G.P.; Nascimento, M.S.B.",
                                          "Lewis, G.P.; Clark, R.", "G.P. Lewis 78",
                                          "Lewis, G.; Rico, L.", "Gwilym P. Lewis",
                                          "G. P. Lewis & J. S. Page", "E.P. LEWIS", "G.P. Lewis e",
                                          "Lewis, G.P.;Klitgaard, B.", "Lewis, G.P.; Ratter, J.A.",
                                          "G.P. Lewis & J.S. Page", "G. P. Lewis; J. S. Page",
                                          "Lewis, GP; Ratter, JA"))

#M. L. Guedes
identifiedby <- replace.names(x = identifiedby, top = 0.80, bottom = 0.7, 
                              check.by = c("M. L. Guedes", "Guedes, M. L.",
                                           "M. L. GUEDES", "GUEDES, M. L."),
                              replace.by = "M. L. Guedes",
                              not.replace = c("M. L. Fonseca"),
                              return = TRUE)

#A. P. Savassi-Coutinho
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = c("A. P. Savassi-Coutinho", "Savassi-Coutinho, A. P.",
                                           "A. P. SAVASSI-COUTINHO", "SAVASSI-COUTINHO, A.P."),
                              replace.by = "A. P. Savassi-Coutinho",
                              replace = c("Coutinho, A.P.S.", "Coutinho, APS"),
                              return = T)

#L. S. B. Jordão
identifiedby <- replace.names(x = identifiedby, top = 0.9, bottom = 0.7, 
                              check.by = c("L. S. B. Jordão", "Jordão, L. S. B.",
                                           "L. S. B. JORDÃO", "JORDÃO, L. S. B."),
                              replace.by = "L. S. B. Jordão",
                              replace  = c("Lucas S.B. Jordão", "LUCAS S.B.JORDÃO",
                                           "Lucas S. B. Jordão", "LUCAS JORDÃO", "Lucas Josdão",
                                           "Lucas, S.B. Jordão", "L.S.B, Jordão", "Lucas Jordão",
                                           "LucasS.B. Jordão", "Jordao, Lucas", "L.S.B., Jordão",
                                           "L.S. B. Jordão", "L. Jordão","Jordão, L.","L.S.B.Jordão", 
                                           "Jordão, LSB", "L.S.B. Jordão"),
                              return = T)

#O. S. Ribas
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = c("O. S. Ribas", "Ribas, O. S.",
                                           "O. S. RIBAS", "RIBAS, O. S."),
                              replace.by = "O. S. Ribas",
                              replace = c("Ribas, OS; Cordeiro, J", "\"\"Ribas, O.S.; Cordeiro, J.\"\"",
                                          "Ribas, OS; Barbosa, E", "O.S. Ribas & J. Codeiro", 
                                          "O.S. Ribas & J. Cordeiro", "Cordeiro, J; Ribas, OS",
                                          "Ribas, OS; Larocca, P", "O.S. Ribas; J. Cordeiro"),
                              return = T)

#I. B. Lima
identifiedby <- replace.names(x = identifiedby, top = 0.93, bottom = 0.7, 
                              check.by = c("I. B. Lima", "Lima, I. B.",
                                           "I. B. LIMA", "LIMA, I. B."),
                              replace.by = "I. B. Lima",
                              return = T)

#M. Morales
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                              check.by = c("M. Morales", "Morales, M.",
                                           "M. MORALES", "MORALES, M."),
                              replace.by = "M. Morales",
                              replace = c("Matias Morales", "M. Morols", "m. Morols"),
                              not.replace = c("M.Moraes"),
                              return = T)

#C. W. Fagg
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                              check.by = c("C. W. Fagg", "Fagg, C. W.",
                                           "C. W. FAGG", "FAGG, C. W."),
                              replace.by = "C. W. Fagg",
                              replace = c("Fagg, CW; Borges, LM"),
                              return = T)

#J. G. Nascimento
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                              check.by = c("J. G. Nascimento", "Nascimento, J. G.",
                                           "J. G. NASCIMENTO", "NASCIMENTO, J. G."),
                              replace.by = "J. G. Nascimento",
                              not.replace = c("L. M. Nascimento", "Nascimento, AFS",
                                              "Nascimento, FHF", "Nascimento, L.M.",
                                              "Nascimento, M.S.B.", "J.M. Nascimento",
                                              "Nascimento"),
                              replace = c("J. G. A. do Nascimento", "J.G.A.do Nascimento",
                                          "J.G.A. do Nascimento"),
                              return = T)

#R. Grether
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                              check.by = c("R. Grether", "Grether, R.",
                                           "R. GRETHER", "GRETHER, R."),
                              replace.by = "R. Grether",
                              return = T)

#E. Córdula
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                              check.by = c("E. Córdula", "Córdula, E.",
                                           "E. CÓRDULA", "CÓRDULA, E"),
                              replace.by = "E. Córdula",
                              return = T)

#A. Bocage
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                              check.by = c("A. Bocage", "Bocage, A.",
                                           "A. BOCAGE", "BOCAGE, A."),
                              replace.by = "A. Bocage",
                              replace = c("A.L.du Bocage Lima", "Bocage, A., Marques, J.S."),
                              return = T)

#J. R. Pirani
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                              check.by = c("J. R. Pirani", "Pirani, J. R.",
                                           "J. R. PIRANI", "PIRANI, J. R."),
                              replace.by = "J. R. Pirani",
                              replace = c("Pirani, JR; Siniscalchi, CM"),
                              return = T)

#R. H. Fortunato
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                              check.by = generate.names("R. H.", "Fortunato"),
                              replace.by = "J. R. Pirani",
                              replace = c("Reenée Fortunato", "Renée H. Fortunato",
                                          "Renee H. Fortunato", "Renée H. Forunato"),
                              return = T)

#H. C. de Lima
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                              check.by = generate.names("H. C.", "de Lima"),
                              replace.by = "H. C. de Lima",
                              not.replace = c("R .C. de Lima"),
                              replace = c("H.C. Lima", "Lima, HC",
                                          "Lima, H.C.", "H. C. de Lima & L.F.G da Silva",
                                          "H C Lima", "H.C. Lima & C.M.J. Mattos",
                                          "H.C.LIMA", "H.C. Lima , C. M. B. Correia",
                                          "H.C. de LIma", "H.C. de Lima & Marli Pires",
                                          "Haroldo C. de Lima", "H.C.Lima", "H. C. de Lima & L. F. G. da Silva"),
                              return = T)

#A. Burkart
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                              check.by = generate.names("A.", "Burkart"),
                              replace.by = "A. Burkart",
                              replace = c("A. Burkart (!RL 2014)", "Buskart",
                                          "A. E. Burkart; L. B. Smith"),
                              return = T)

#A. M. Miranda
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                              check.by = generate.names("A. M.", "Miranda"),
                              replace.by = "A. M. Miranda",
                              replace = c("A.M.Miranda et H.Gomes", "A.M. Miranda; M.L. Guedes"),
                              return = T)

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


