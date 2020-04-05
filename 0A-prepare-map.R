rm(list=ls())
#------------------
# logitude and latitute to 
# country provinance
#------------------
long.lat <- read.csv(
  './pop-long-lat.csv',
  stringsAsFactors=FALSE,header=FALSE)
colnames(long.lat) <- c('subregion','lat','long')
rownames(long.lat) <- long.lat$subregion
#-----------------
pop.dd <- read.csv('./pop-dd.csv',stringsAsFactors=FALSE)
rownames(pop.dd) <- pop.dd$subregion
#-----------------
pop.dd.new <- cbind(pop.dd[,-c(1,6:7)],long.lat[rownames(pop.dd),])
rownames(pop.dd.new) <- 1:nrow(pop.dd.new)
save(pop.dd.new, file='./data/pop.dd.new.RData')
#-----------------
rm(list=ls())
library(maptools)
library(raster)
library(rgdal)

sh.files <- list.files(
  path='./data/gadm/',
  recursive = TRUE,
  pattern= '_1.shp',
  full.names = TRUE
)

load('./data/pop.dd.new.RData')
pop.country <- plyr::dlply(
  pop.dd.new,
  'region'
)
source('./0A-Z-function.R')
pop.DF <- list()
cat('---------------*---------------\n')
for(i in 1:length(sh.files)){
  #i = 19
  map.data <- readShapeFile(shapefile = sh.files[i])
  country <- levels(map.data@data$NAME_0)
  provinace <- levels(map.data@data$NAME_1)
  # present <- any(names(pop.country) %in% country)
  # if(!present){
  #   print(i)
  # }
  # cat('---------------*---------------\n')
  pop.LL <- pop.country[country][[1]]
  prov.OUT <- list()
  for(j in 1:length(provinace)){
    prov.coord <- data.frame(
      map.data@polygons[[j]]@Polygons[[1]]@coords
    )
    pop.coord <- pop.LL[,c('long','lat')]
    colnames(prov.coord) <- colnames(pop.coord)
    OUT <- list()
    for(z in 1:nrow(pop.coord)){
      long = pop.coord[z,1]
      lat = pop.coord[z,2]
      #-------------------------
      library(geosphere)
      message(i,'. j:: ',j,'/',length(provinace), ' z:: ',z)
      STdist <- apply(prov.coord,1,function(x){
        distm(x = as.numeric(x),y = c(long,lat))
      })
      idx<- order(STdist)[1]
      out <- prov.coord[idx,]
      out$dist <- as.numeric(STdist[idx])
      #----------------------
      out$idx <- idx
      out$proviance <- j
      out$z <- z
      OUT[[z]] <- out
    }
    OUT.df <- plyr::ldply(OUT)
    prov.OUT[[j]] <- OUT.df
  }
  prov.df <- plyr::ldply(prov.OUT)
  prov.closest <- plyr::dlply(
    prov.df,
    'z',
    function(x){
      idx.closest <- order(x$dist)[1]
      x$proviance[idx.closest]
    }
  )
  pop.LL$provID <- unlist(prov.closest)
  pop.LL$prov <- provinace[pop.LL$provID]
  pop.DF[[i]] <- pop.LL
}
pop.DF.prov <- plyr::ldply(pop.DF)
save(pop.DF.prov,file='./data/pop.DF.prov.RData')
