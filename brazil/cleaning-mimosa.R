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

# Plotting a bar plot of total identifications per name (before standardization)
ggplot(mimosa %>% filter(identifiedby %in% names.count$x[names.count$freq >= 50]), 
       aes(x=reorder(identifiedby,identifiedby, 
                     function(x)-length(x)))) + 
  geom_bar()

# Where do records identified by "Administrador" come from?
administrador <- mimosa %>% filter(identifiedby %in% c("Administrador")) %>% dplyr::select(institutioncode, species, identifiedby)
ggplot(administrador, 
       aes(x=reorder(institutioncode,institutioncode, function(x)-length(x)))) + 
  geom_bar()

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
                                            "Marchiori, J.N.; Barneby, R.C.", "R. BernebY", 
                                            "Rupert Barneby" ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# M Simon
replace.by <- "M Simon"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c( 
                                  ),
                                replace = c("Marcelo Fragomeni Simon", "M. Simon, IN LIT. L. P. Queiroz",
                                            "M. Simon & L. M. Borges", "M.Simon & L.M.Borges",
                                            "M. F. Simon et. al.;","M.F.Jimon", 
                                            "Simon M.F. & Marcelo F.", "Marcelo Simon", 
                                            "Marcelo, F.S.", "Marcelo Fragomeni Simon",
                                            "Colombini, M.A.G; Starling, M.F.V confirmado por M.F.Simon, na exsicata escreveram calodrendon",
                                            "Marcelo Simon", "Marcelo F. Simon", "Marcelo, F.S.", "Marcelo Fragomeni Simon", 
                                            "R.Rodrigues da Silva - confirmado por M.F.Simon", "Proença, CEB; Simon, MF"
                                            ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# LM Borges 
replace.by <- "LM Borges"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Borges, S." 
                                ),
                                replace = c("Leonardo Maurici Borges; L.M.Borges", "I. M. Borges",
                                            "L. M. Borges (SPF) 2013", "L.M.Borges (SPF)",
                                            "Leonardo M. Borges", "L.M. Borges et. al.;",
                                            "L. M. Borges & M. F. Simon"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# JS Silva
replace.by <- "JS Silva"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Silva, J.M.", "Silva, L.A.", "Silveira, FS", "Silva, TR",
                                                "Silva, CAS", "Silva, A.F.B.", "Silva, D.F.", "Silva, A.S.",
                                                "S.M. Silva", "Silva, RR", "Silva, M.A. da Silva", 
                                                "Silva, N.T. da", "Silva, C.A.S. da", "Silveira, AM", 
                                                "Silva, W. L.", "Silva, WL", "Silva, R.R.", "Silva, JM",
                                                "Silva, G.S", "Silva, A.K.C.", "Silva, MG da", "Silveira, N.",
                                                "SILVEIRA, J.E.", "Silva, BD", "Silva, AS; Crísci, DB",
                                                "Silva, D.R.", "Silva, R.Q.", "Silva Filho, PJS", 
                                                "Silveira, F.S.", "Silva, J.P.", "Silva, ACC", "Silva, ML", 
                                                "Silvia, ACC", "Silva, MJ", "Silva, T.C.", "Silva, E.D.",
                                                "Silva, KC", "Silva,L.A.", "Silva, W",  "Silva, R. R.",
                                                "Silva, S. Y. B.", "Silva, M.L.", "SILVA R.R.", "Silva, DR",
                                                "Silva, BG", "Silva, CAS da", "SILVA, R.R.", "Silva, M.G. da",
                                                "Silva, F.M. da", "Silva, A.S.L. da", "Silva, ED", "Silva, GP",
                                                "Silva, MA", "Silva, MS", "Silva, S.M.", "Silva, M.A.", 
                                                "Silva, AG da", "Silva, M.F. da", "Silva, M.S.S.", "J. M. Silva",
                                                "Silva, R.A.", "Silva, MF da", "Silveira, F. S.", "Silva, ASL da",
                                                "Silva, WLS", "Silva, NT da", "Silva, ACB", "Silva, N.C.B.",
                                                "Silveira, F.R.", "Silva, J.L.", "Silva, A S L da",  "Silva, CA",
                                                "Silva, CFS", "Silva, DP", "Silva, A.G. da"
                                ),
                                replace = c("Sales, M.; Silva, J.S.", "J.S.Silva at M.Sales", "Santos Silva, JA",
                                            "Juliana Santos Silva; Universidade Estadual de Campinas", "Santos,J",
                                            "Santos, J.S.", "Juliana Santos Silva", "J.S.Silva; M. Sales", 
                                            "M.Sales; J.S.Silva", "Santos, JS", "J.S. Silva (1) & M.C. Abreu", 
                                            "J.S. Silva (1) et al.", "J.S. Silva; M. Sales", "J.S.Silva/14-VIII-2009",
                                            "Juliana Santos Silva (UEC)", "J. S. Silva 2009-08-14", "Santos-Silva, J",
                                            "J. S. Silva - UEC", "J. Santos S.", "J. Santos Silva", "J. S. Silva - UFRPe"
                                            ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# LP Queiroz
replace.by <- "LP Queiroz"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Queiroz, D.L.", "Queiroz, R.T.", "Queiroz", 
                                                "Queiroz, RT", "Queinz, R.", "Queiroz, R.F."
                                ),
                                replace = c( "L. P. de Queiroz (HUEFS)", "L.P.de Queiroz & T.S.Nunes", 
                                             "L.P.de Queiroz & D.Cardoso", "L.P.de Queiroz & R.M.Santos",
                                             "L.P.de Queiróz", "L. Paganucci", "Guedes, ML; Queiroz, LP de",
                                             "L. P. de Queiroz", "Santos, RM; Queiroz, LP de"
                                            ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# ML Guedes 
# Obs: She was not considered an expert according to the google scholar search. 
# However, through personal communication with other specialists, we decided to include her as an identifier.

replace.by <- "ML Guedes"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(
                                ),
                                replace = c("Silva, ML"
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# GP Lewis
replace.by <- "GP Lewis"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(
                                ),
                                replace = c("Ratter J.A. & Lewis G.P.", "G.P. Lewis & J.S. Page",
                                            "g. P. Lewis", "G.P. LEWIS & J.S. PAGE", "G. P. Lewis & J. S. Page",
                                            "G. P. Lewis; J. S. Page", "G. P. Lewis; R. Clark", "G.P.L"
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# LSB Jordão
replace.by <- "LSB Jordão"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(
                                ),
                                replace = c("L. Jordão & H.C. Lima", "Lucas Jodão & D. Machado", 
                                            "Lucas Jordão & M. Barros", "Lucas S. B. Jordão", 
                                            "Lucas Jordão (via foto)", "Lucas, S.B. Jordão",
                                            "Lucas Jordão; D. Machado", "L. Jordão; D.N.S. Machado", 
                                            "L. Jordão; H.C. Lima", "Lucas Jordão & D. Machado",
                                            "Machado & L. Jordão", "D. Machado & L. Joerdão"
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# VF Dutra
replace.by <- "VF Dutra"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(
                                ),
                                replace = c("VALQUIRIA F DUTRA, 21-04-2020", "V.F.Dutra (VIC)", "V.L.Dutra",
                                            "Valquíria Ferreira Dutra"    
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# OS Ribas
replace.by <- "OS Ribas"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(
                                ),
                                replace = c("O.S. Ribas; G. Felitto", "O.S. Ribas & J.M. Silva",
                                            "O. S. Ribas & J. Cordeiro, 2000.", "O.S.Ribas & J. Cordeiro",
                                            "Cordeiro, J.; Ribas, O.S.", "O.S.Ribas; R.R.Silva", 
                                            "O. S. Ribas,  J. Cordeiro", "O. S. Rivas", "Cordeiro, J; Ribas, OS",
                                            "Motta, J.T.; Ribas, O.S."  
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
                                not.replace = c("M. Morols", "Moraes, MD", "Moraes, L.A.", "Moraes, L.A",
                                                "Moraes, L. A", "V. M. Morales", "Moraes, JG"
                                ),
                                replace = c("Matias Morales" 
                                  
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# APS Coutinho
replace.by <- "APS Coutinho"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Couto, LA"
                                ),
                                replace = c("Savassi-Coutinho", "A.P.Savassi-Coutinho", "A.P. Savassi-Coutinho/25-09-2010",
                                            "Savassi-Coutinho, A. P. (2008)", "Savassi, AP", "Savassi-Coutinho,A.P.", 
                                            "Savassi-Coutinho, A. P.", "Savassi- Coutinho, A.P", "A. P. Savassi-Coutinho",
                                            "A.P. Savassi-Coutinho", "A. P. Savassi-Coutinho 2010-09-25", 
                                            "Savassi-Coutinho, AP", "A. Savassi-Coutinho", "Savassi-Coutinho, A.P."
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# JNC Marchiori
replace.by <- "JNC Marchiori"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Marchiori, MS"
                                ),
                                replace = c("A.A. Carpanezzi apud Marchiori"
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# A Fernandes 
replace.by <- "A Fernandes"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Fernanda Cruz", "Fernandes-Bulhão C.", "D.S. Fernandes"
                                ),
                                replace = c("A. Fernandes | E. Nunes", "A. Fernandes | P. Bezerra",
                                            "A.Fernandes | P.Bezerra", "A.Fernandes | E.Nunes",
                                            "A. Fernandes; E. P. Nunes", "A. Fernandes; P. Bezerra", 
                                            "A. G. Fernandes; P. Bezerra" 
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# MG Caxambu
replace.by <- "MG Caxambu"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# G Hatschbach (?) 
# Obs: He was not considered an expert according to the google scholar search. 
# However, through personal communication with other specialists, we decided to include him as an identifier.

replace.by <- "G Hatschbach"
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

# JGA Nascimento
replace.by <- "JGA Nascimento"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Nascimento, M.B.", "Nascimento, FHF",
                                                "J.M. Nascimento", "Nascimento, A.",
                                                "Nascimento, M.S.B.", "Nascimento, M.G.P.",
                                                "Nascimento, MP do", "Nascimento, AFS",
                                                "Nascimento, O.C. do", "Nascimento, L.M.",
                                                "Nascimento, M.S.B.; Alencar, M.E.",
                                                "Nascimento, J.C.F.", "Nascimento - Júnior, JE",
                                                "Nascimento, VT; L.G.Sousa", "Nascimento, M.S.B.; Lewis, G.P.",
                                                "J.M.Nascimento", "Nascimento, EAP","Nascimento-Junior, JE"),
                                replace = c("J.G.A. do Nascimento", "J.G.A.do Nascimento",
                                            "J. G. A. do Nascimento"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# E Nunes
replace.by <- "E Nunes"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("D. Nunes","Nunes, TS", "Nunes, SS", "Nunes, GP",
                                                "Nunes, SRDFS;", "E.M.B. Nunes"),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# RR Silva
replace.by <- "RR Silva"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Silva, J.M.", "Silva, ACC", "Silva, D.F.",
                                                "Silveira, AM", "Silva, A.G. da", "Silva, ASL da",
                                                "Silva, J.L.", "Silva, L.A.", "Silvia, ACC",
                                                "Silva, S.M.", "Silva, M.F. da", "Silva, WLS",
                                                "Silva, CFS", "Silva, CAS", "Silva, CAS da",
                                                "Silva, M.S.S.", "Silveira, N.", "Silveira, FS",
                                                "Silva, T.C.", "Silva, A.S.", "Silva, W. L.", 
                                                "Silva, G.S", "Silva, NT da", "R.Q. Silva",
                                                "Silva, A.F.B.", "Silva, M.L.",
                                                "Silva, M.G. da", "Silva, M.A.", "Silva, A.K.C.",
                                                "Silva, ACB", "Silveira, F.S.", "Silva, E.D.", 
                                                "R.H. Silva", "Silva, N.T. da", "Silva, AG da",
                                                "Silva, MF da", "Silva, N.C.B.", "Silva, J.P.",
                                                "Silva,L.A.", "Silva, F.M. da", "Silva, MG da",
                                                "Silveira, F.R."),
                                replace = c("R.R.Silva", "R. R. Silva"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# FS Silveira
replace.by <- "FS Silveira"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Oliveira, OF", "Silva, A.S.", "Silva, G.S",
                                                "Siveira, N.", "Silvia, ACC", "Silva, CAS da",
                                                "Silva, ACB", "Silva Filho, PJS", "Oliveira, MS",
                                                "Silveira, AM", "Silva, MF da", "SILVEIRA, J.E.",
                                                "Silva, CAS", "Silva, S.M.", "Silva, ASL da",
                                                "Silveira, F.R.", "Silva, A.F.B.", "Silveira A.M.",
                                                "Silva, WLS", "Oliveira, F", "Silva, ACC",
                                                "Silva, D.F.", "Silva, M.S.S.", "Silveira, N.",
                                                "S.F. Silveira"),
                                replace = c("F.Schmidt-Silveira", "Fernanda Schmidt Silveira",
                                            "Schmidt-Silveira, F", "Schmidt, F", "F.S. Schmidt",
                                            "Schmidt, FS", "Schimdt-Silveira, F", "Schimidt-Silveira, F.",
                                            "Schmit-Silveira, F.", "Schmidt Silveira, F.", "F. Schmidt-Silveira",
                                            "F.Schmidt- Silveira", "Schmidt-Silveira, F.", "Schimidt, F.",
                                            "F. Schmidt- Silveira", "Schmidt-Silveira"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# A Pott
replace.by <- "A Pott"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Pontes, RA"),
                                replace = c("A. Pott & V.J. Pott", 
                                            "A.Pott, V.J. Pott, O.S. Ribas", "A. Pott; G.A. Damasceno-Júnior; F. Macedo-Alves",
                                            "V.J.Pott & A.Pott", "A. Pott; V. J. Pott", "A. Pott; V.J. Pott", "A. Pott ;V.J. Pott"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# AA Carpanezzi
replace.by <- "AA Carpanezzi"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Carneiro, AS", "F.B. Carpanezzi"),
                                replace = c("A.A. Carpanezzi & R.R. Völtz",
                                            "Völtz, RR; Carpanezzi, AA"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# VC Souza
replace.by <- "VC Souza"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Sousa, V.F.", "Souza, E.B.", "Souza, H.F.",
                                                "Souza, D.P.", "Souza, R.T.B", "Souza, M.L.",
                                                "Souza, G.A.B.", "Souza, M.K.F.", "Souza, JRP",
                                                "Souza, J.F.O.", "Souza, F.S.", "Souza, J.P.",
                                                "Souza, RS", "Souza, I.", "Souza, R.S."),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# M Flores-Cruz
replace.by <- "M Flores-Cruz"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Flores, A", "Flores, AS"),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# A Burkart
replace.by <- "A Burkart"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("A. E. Burkart; L. B. Smith"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# R Vanni
replace.by <- "R Vanni"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# WA Rodrigues
replace.by <- "WA Rodrigues"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Rodrigues, LMO", "Rodrigues, R.R.", "Rodrigues",
                                                "Rodrigues Silva, R", "Rodrigues, RR; ; Souza, VC",
                                                "Rodrigues, M.D", "Rodrigues, RR; Souza, VC",
                                                "Rodrigues, M.L.", "Rodrigues, M. dos S."),
                                replace = c("Oliveira, E; Rodrigues, WA", "W.A. Rodrigues (conf. R.S. Cowan 1973)"
                                            ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# E Cordula
replace.by <- "E Cordula"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(  
                                ),
                                replace = c("E.Córdula", "E. CORDULA & M. ALVES", "Mattos, C; Córdula"
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]


#D Cardoso
replace.by <- "D Cardoso"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c( "Cardoso, A"  
                                ),
                                replace = c("D.Cardoso & P.W.Moonlight" 
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# VJ Pott
replace.by <- "VJ Pott"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(  
                                ),
                                replace = c("V. J. Pott; I. M. Bortolotto" 
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]


#ALB Sartori
replace.by <- "ALB Sartori"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(  
                                ),
                                replace = c("A.L.B. Sartori; F.J. Kochanovski"
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# HS  Irwin
replace.by <- "HS  Irwin"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(  
                                ),
                                replace = c(
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#EP Seleme
replace.by <- "EP Seleme"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(  
                                ),
                                replace = c(
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#RH Fortunato
replace.by <- "RH Fortunato"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(  
                                ),
                                replace = c("René Fortunato"
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]


#	HC Lima
replace.by <- "HC Lima"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Lima, L.C.P.", "Lima, M.R.", "Lima, IB", "Lima, AS", "Lima, RB",
                                                "Lima, M.P.", "Lima, A. A.", "Lima, G.A.", "Lima, J", "Lima, J.",
                                                "Lima, MPM de", "Lima, BCS", "Lima, M.E.L.", "Lima, M. de", "Lima, MPM",
                                                "Lima, I.B", "Lima", "Lima, J.R.", "Lima, M.P.M.", "Lima, MP de", "Lima H.S.",
                                                "Lima, F. .F.", "Lima, V.C.", "LIMA, A.C.R."
                                ),
                                replace = c("H. C. de Lima & L. F. G. da Silva", "Haroldo Cavalcante de Lima & Robson Daumas Ribeiro",
                                            "\"\"Lima, H.C. de; Ribeiro, R.D.\"\"", "H.C. de Lima & T.M. Rodrigues", 
                                            "H.C. de Lima & Marli Pires", "H. C. de Lima & L.F.G da Silva", "H.C. Lima; C.M.J. Mattos",
                                            "H.C. Lima & T.M. Rodrigues", "H.C. Lima & C. Mattos",  "H.C. Lima; C.M.J. Mattos; A. Santos",
                                            "H.C. Lima; C.M.J. Mattos & A. Santos", "H.C. Lima, C.M.J. Mattos & A. Santos", 
                                            "H.C. Lima & C.M.J. Mattos", "H C de Lima", "H.C.de Lima; J.E.Meireles", "H. C. de Lima"
                                            
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
                                not.replace = c(  
                                ),
                                replace = c(
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#STS Miotto
replace.by <- "STS Miotto"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(  
                                ),
                                replace = c(
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#A Ducke
replace.by <- "A Ducke"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(  
                                ),
                                replace = c(
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#FC Hoehne 
# Obs: He was not considered an expert according to the google scholar search. 
# However, through personal communication with other specialists, we decided to include him as an identifier.

replace.by <- "FC Hoehne"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("W. Hoehne"
                                ),
                                replace = c("F. C. Hochne"
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
                                not.replace = c(  
                                ),
                                replace = c("Jorge Yoshio Tamashiro", "J.Y.T", "J.Y. Tamashiro & L.D. Meireles"    
                                ))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# Replacing column
mimosa_exp <- mimosa
mimosa_exp$identifiedby <- identifiedby

# Plotting a bar plot of total identifications per name (after standardization), with a minimum of 25 identifications
#identifyby_dataframe <- as.data.frame(identifiedby)
#ggplot(identifyby_dataframe %>% filter(identifiedby %in% names.count$x[names.count$freq >= 25]), 
#       aes(x=reorder(identifiedby,identifiedby, function(x)-length(x)))) + 
#  geom_bar() +
#  theme(axis.text.x = element_text(angle = 90))

# Plotting a bar plot of total identifications per name (after standardization), with a minimum of 95 identifications
#identifiedby_dataframe <- as.data.frame(identifiedby)
#ggplot(identifiedby_dataframe %>% filter(identifiedby %in% names.count$x[names.count$freq >= 95]), 
#       aes(x=reorder(identifiedby,identifiedby, function(x)-length(x)))) + 
#  geom_bar() +
#  theme(axis.text.x = element_text(angle = 45))

# Filtering by specialists (27,324)
mimosa_exp1 <- mimosa_exp %>% filter(identifiedby %in% specialists)

# Filtering by specialists with a minimum of 26 identifications (27,324)
#mimosa_exp2 <- mimosa_exp %>% filter(identifiedby %in% 
#                                      specialists[specialists %in% 
#                                                    identifiedby_dataframe[identifiedby_dataframe$identifiedby %in%
#                                                                           names.count$x[names.count$freq > 25], ]])  
                                                            
# Filtering by specialists with a minimum of 95 identifications (26,559)
#mimosa_exp3 <- mimosa_exp %>% filter(identifiedby %in% 
#                                       specialists[specialists %in% 
#                                                     identifiedby_dataframe[identifiedby_dataframe$identifiedby %in%
#                                                                              names.count$x[names.count$freq > 94], ]])  

mimosa <- mimosa_exp1

#======================================================================================#

#===========================#
# CLEANING SCIENTIFIC NAMES #
#===========================#

# Generating a column for scientific names (without authors and including infraspecific epithet)
mimosa$gen_sp <- paste(mimosa$genus,
                      mimosa$species,
                      mimosa$subspecies,
                      sep = " ")

# Extracting scientific names
taxa <- plyr::count(mimosa$gen_sp)
taxa <- as.character(taxa$x)

# Removing NA (character derived from 'subspecies' attribute) COMENTAR LINHA A LINHA
for(i in 1:length(taxa)){
  taxa[i] <- gsub("NA", "", taxa[i])
  taxa[i] <- trimws(taxa[i])
  taxa[i] <- gsub(pattern = "  ", x = taxa[i], replacement = " ")
}

#Suggesting with flora (and retrieving a few additional information that may be useful)
taxa_suggested <-get.taxa(taxa, vegetation.type = TRUE, 
                          habitat = TRUE, domain = TRUE, life.form = TRUE)

#Writing *.csv for manual checking
write.csv(taxa_suggested, file = "taxa_suggested.csv", row.names = F)

