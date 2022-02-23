# Checking occurrence of  previously selected species from my project

selected_species<- read.csv("lists/mimosa_spp_phd2.csv", header=T, sep = ';')
colnames(selected_species)<-c("genus_and_species", "clade")

selected_species<-separate(data = selected_species, col = genus_and_species, into = c("genus", "species", "variety"), sep = "\\_") 

mimosa_coordClean$selected_species<- mimosa_coordClean$species %in% selected_species$species

mimosa_noCoord$selected_species<- mimosa_noCoord$species %in% selected_species$species

# Checking occurrence of species  present in the phylogeny
phylogeny_species<- read.csv("lists/Mimosa_data_updated_140520-VASCONCELOS2020.csv", 
                             sep = ",")

phylogeny_species<-separate(data = phylogeny_species, col = cleaned_name, 
                            into = c("genus", "species", "variety"), sep = "\\_")

mimosa_coordClean$phylogeny_species<- mimosa_coordClean$species %in% phylogeny_species$species

mimosa_noCoord$phylogeny_species<- mimosa_noCoord$species %in% phylogeny_species$species

mimosa_coordClean<-merge(mimosa_coordClean, phylogeny_species[,c("species", "clade")], by.x= "species", 
                         by.y="species", all.x=T)

mimosa_noCoord<-merge(mimosa_noCoord, phylogeny_species[,c("species", "clade")], by.x= "species", 
                      by.y="species", all.x=T)


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
  filter (grepl('cipó|cipo|Cipó|Cipo|CIPO|CIPÓ|veadeiros|Veadeiros|VEADEIROS', locality))

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