#------------------------
# Spatial Analysis with R book
# set.seed(0)
# xy <- cbind(x=runif(1000, 0, 100), y=runif(1000, 0, 100))
# income <- (runif(1000) * abs((xy[,1] - 50) * (xy[,2] - 50))) / 500
# 
# par(mfrow=c(1,3), las=1)
# plot(sort(income), col=rev(terrain.colors(1000)), pch=20, cex=.75, ylab='income')
# hist(income, main='', col=rev(terrain.colors(10)), xlim=c(0,5), breaks=seq(0,5,0.5))
# plot(xy, xlim=c(0,100), ylim=c(0,100), cex=income, 
#      col=rev(terrain.colors(50))[10*(income+1)])
# 
# n <- length(income)
# G <- (2 * sum(sort(income) * 1:n)/sum(income) - (n + 1)) / n
# 
# ibrary(raster)
# r1 <- raster(ncol=1, nrow=4, xmn=0, xmx=100, ymn=0, ymx=100, crs=NA)
# r1 <- rasterize(xy, r1, income, mean)
# r2 <- raster(ncol=4, nrow=1, xmn=0, xmx=100, ymn=0, ymx=100, crs=NA)
# r2 <- rasterize(xy, r2, income, mean)
# r3 <- raster(ncol=2, nrow=2, xmn=0, xmx=100, ymn=0, ymx=100, crs=NA)
# r3 <- rasterize(xy, r3, income, mean)
# r4 <- raster(ncol=3, nrow=3, xmn=0, xmx=100, ymn=0, ymx=100, crs=NA)
# r4 <- rasterize(xy, r4, income, mean)
# r5 <- raster(ncol=5, nrow=5, xmn=0, xmx=100, ymn=0, ymx=100, crs=NA)
# r5 <- rasterize(xy, r5, income, mean)
# r6 <- raster(ncol=10, nrow=10, xmn=0, xmx=100, ymn=0, ymx=100, crs=NA)
# r6 <- rasterize(xy, r6, income, mean)
# par(mfrow=c(2,3), las=1)
# plot(r1); plot(r2); plot(r3); plot(r4); plot(r5); plot(r6)
# 
# par(mfrow=c(1,3),las=1)
# hist(r4, main='',col=rev(terrain.colors(10)), xlim=c(0,5), breaks=seq(0, 5, 0.5))
# hist(r5, main='',col=rev(terrain.colors(10)), xlim=c(0,5), breaks=seq(0, 5, 0.5))
# hist(r6, main='',col=rev(terrain.colors(10)), xlim=c(0,5), breaks=seq(0, 5, 0.5))
# 
# 
# #--------------------
# # AUTOCORRELATION
# set.seed(0)
# d <- sample(100, 10)
# d
# # compute auto-correlation
# a <- d[-length(d)]
# b <- d[-1]
# plot(a, b, xlab='t', ylab='t-1')



library(raster)
p <- shapefile(system.file("external/lux.shp", package="raster"))
p <- p[p$NAME_1=="Diekirch", ]
p$value <- c(10, 6, 4, 11, 6)
par(mai=c(0,0,0,0))
plot(p, col=2:7)
xy <- coordinates(p)
points(xy, cex=6, pch=20, col='white')
text(p, 'ID_2', cex=1.5)
#------------
library(spdep)
w <- poly2nb(p, row.names=p$Id)
plot(p, col='gray', border='blue', lwd=2)
plot(w, xy, col='red', lwd=2, add=TRUE)
wm<- nb2mat(w, style='B')

n <- length(p)
y <- p$value
ybar <- mean(y)

dy <- y - ybar
g <- expand.grid(dy, dy)
yiyj <- g[,1] * g[,2]

yi <- rep(dy, each=n)
yj <- rep(dy)
yiyj <- yi * yj

pm <- matrix(yiyj, ncol=n)
pmw <- pm * wm
spmw <- sum(pmw)

smw <- sum(wm)
sw <- spmw / smw

vr <- n / sum(dy^2)
MI <- vr * sw

ww <- nb2listw(w, style='B')
moran(p$value, ww, n=length(ww$neighbours), S0=Szero(ww))
moran.test(p$value, ww, randomisation=FALSE)

moran.mc(p$value, ww, nsim=99)


n <- length(p)
ms <- cbind(id=rep(1:n, each=n), y=rep(y, each=n), value=as.vector(wm * y))
ms <- ms[ms[,3] > 0, ]

ams <- aggregate(ms[,2:3], list(ms[,1]), FUN=mean)
ams <- ams[,-1]
colnames(ams) <- c('y', 'spatially lagged y')
head(ams)

plot(ams)
reg <- lm(ams[,2] ~ ams[,1])
abline(reg, lwd=2)
abline(h=mean(ams[,2]), lt=2)
abline(v=ybar, lt=2)
#------------------------
#------------------------
if (!require("rspatial")) remotes::install_github('rspatial/rspatial')
library(rspatial)
d <- sp_data('precipitation')
d$prec <- rowSums(d[, c(6:17)])
plot(sort(d$prec), ylab='Annual precipitation (mm)', las=1, xlab='Stations')

library(sp)
dsp <- SpatialPoints(d[,4:3], proj4string=CRS("+proj=longlat +datum=NAD83"))
dsp <- SpatialPointsDataFrame(dsp, d)
CA <- sp_data("counties")
# define groups for mapping
cuts <- c(0,200,300,500,1000,3000)
# set up a palette of interpolated colors
blues <- colorRampPalette(c('yellow', 'orange', 'blue', 'dark blue'))
pols <- list("sp.polygons", CA, fill = "lightgray")
spplot(dsp, 'prec', cuts=cuts, col.regions=blues(5), sp.layout=pols, pch=20, cex=2)


TA <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000+datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0")

library(rgdal)
dta <- spTransform(dsp, TA)
cata <- spTransform(CA, TA)

# null model
RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}
null <- RMSE(mean(dsp$prec), dsp$prec)

library(dismo)
v <- voronoi(dta)
plot(v)

ca <- aggregate(cata)
vca <- intersect(v, ca)
spplot(vca, 'prec', col.regions=rev(get_col_regions()))

r <- raster(cata, res=10000)
vr <- rasterize(vca, r, 'prec')
plot(vr)

#-------------
# STEP 01: DATA FRAME to spatial obeject
#-------------
# AQ data
library(rspatial)
x <- sp_data("airqual")
x$OZDLYAV <- x$OZDLYAV * 1000

library(sp)
coordinates(x) <- ~LONGITUDE + LATITUDE
proj4string(x) <- CRS('+proj=longlat +datum=NAD83')
TA <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000+datum=NAD83 +units=km +ellps=GRS80")
library(rgdal)
aq <- spTransform(x, TA)
#-------------
# STEP 02: THE MAP aligned to the INPUT DATA
cageo <- sp_data('counties.rds')
# align to the input coordinates
ca <- spTransform(cageo, TA)
r <- raster(ca)
res(r) <- 10 # 10 km if your CRS's units are in km
g <- as(r, 'SpatialGrid')
#-------------
RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}
f1 <- function(x, test, train) {
  nmx <- x[1]
  idp <- x[2]
  if (nmx < 1) return(Inf)
  if (idp < .001) return(Inf)
  m <- gstat(formula=OZDLYAV~1, locations=train, nmax=nmx, set=list(idp=idp))
  p <- predict(m, newdata=test, debug.level=0)$var1.pred
  RMSE(test$OZDLYAV, p)
}
set.seed(20150518)
i <- sample(nrow(aq), 0.2 * nrow(aq))
tst <- aq[i,]
trn <- aq[-i,]
opt <- optim(c(8, .5), f1, test=tst, train=trn)
opt
#-------------
# Our optimal Inverse Distance Weighted (IDW) model
#-------------
# A more commonly used method is “inverse distance weighted” 
# interpolation. The only difference with the nearest
# neighbour approach is that points that are further away 
# get less weight in predicting a value a location.
library(gstat)
m <- gstat(formula=OZDLYAV~1, locations=aq, nmax=opt$par[1], set=list(idp=opt$par[2]))
idw <- interpolate(r, m)
## [inverse distance weighted interpolation]
idw <- mask(idw, ca)
plot(idw)
#-------------
# A thin plate spline model
library(fields)
m <- Tps(coordinates(aq), aq$OZDLYAV)
tps <- interpolate(r, m)
tps <- mask(tps, idw)
plot(tps)
