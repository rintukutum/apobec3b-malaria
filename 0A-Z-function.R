#############################
##### SECTION 0A
#############################
#--------------
# write the code to extract the regin of polygon from 
# a defined polygon coordinates
readShapeFile <- function(shape.file){
  # layer
  layer <- ogrListLayers(shape.file)
  # read the shape file
  sp.df <- readOGR(shape.file, layer=layer)
  return(sp.df)
}
#############################
##### SECTION 0B
#############################
getProv4mCountry <- function(country, province){
  # extract the province spatial polygon object 
  # china.shp.file <- './data/gadm/gadm36_CHN_shp/gadm36_CHN_1.shp'
  # country <- readShapeFile(china.shp.file)
  # province <- c('Hubei','Liaoning')
  # sp.df <- getProv4mCountry(country = country,province = province)
  idx.provice <- country$NAME_1 %in% province
  if(any(idx.provice)){
    ncount <- length(which(idx.provice))
    if(ncount==1){
      entry <- 'entry'
    }else{
      entry <- 'entries'
    }
    message('------------***------------')
    message('Found ',ncount,' ',entry)
    #---------------------------------
    INDEX.province <- which(idx.provice)
    SP.province <- country[INDEX.province,]
  }else{
    message('------------***------------')
    message('Sorry!')
    #---------------------------------
    SP.province <- NULL
  }
  return(SP.province)
}
getINDEX <- function(x){
  library(maptools)
  data("wrld_simpl")
  c <- raster(
    nrow=360, ncol=720, 
    crs=proj4string(wrld_simpl), 
    xmn=-180, 
    xmx=180, 
    ymn=-90, 
    ymx=90)
  c[] <- 1:length(c)
  #-----------------------
  # INDEX for the raster
  clip1 <- raster::crop(c, raster::extent(x)) 
  # crops the raster to the polygon boundary
  clip2 <- raster::mask(clip1,x) 
  # extracts data from the raster based on the polygon bound
  ext<- raster::extract(clip2,x) 
  #-----------------------
  # remap
  # author: RINTU KUTUM
  col.vec <- rep(NA,length(c))
  col.vec[ext[[1]]] <- 'col'
  col.mat <- matrix(data = col.vec,nrow = 720,ncol = 360)
  col.mat <- col.mat[,ncol(col.mat):1]
  idx.col.ext <- which(col.mat == 'col')
  return(idx.col.ext)
}
get.INDEX.worldmap <- function(sp.df){
  #' source: https://gis.stackexchange.com/questions/130522/increasing-speed-of-crop-mask-extract-raster-by-many-polygons-in-r
  #-----------------------
  # "SpatialPolygonsDataFrame"
  # 
  if(class(sp.df) == "SpatialPolygonsDataFrame"){
    idx.col.ext <- list()
    for(i in 1:length(sp.df)){
      idx.col.ext[[i]] <- getINDEX(x = sp.df[i,])
    }
    out <- list(
      INDEX = idx.col.ext,
      metadata = as.data.frame(sp.df@data)
    )
  }else{
    out <- list(
      INDEX = NULL,
      metadata = NULL
    )
  }
  # china.shp.file <- './data/gadm/gadm36_CHN_shp/gadm36_CHN_1.shp'
  # country <- readShapeFile(china.shp.file)
  # province <- c('Hubei','Liaoning')
  # sp.df <- getProv4mCountry(country = country,province = province)
  # out <- get.INDEX.worldmap(sp.df = sp.df)
  return(out)
}

