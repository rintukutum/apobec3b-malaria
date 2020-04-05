rm(list = ls())
load('./data/pop.DF.prov.RData')
#--------------
# approximate the state level information

country.LL <- plyr::dlply(
  pop.DF.prov,
  'region'
)
country.shapefiles <- list.files(
  path = './data/gadm/',
  pattern = '_1.shp',
  recursive = TRUE,
  full.names = TRUE
)
country3letter <- as.character(sapply(
  country.shapefiles,
  function(x){
    strsplit(x,split='\\_')[[1]][2]
  }
))
source('./0A-Z-function.R')



