# INTERPOLATION TEST
## looking for a method to capture population data from hh listing datasets and allocate it over BF
## testing different methods over a piece of the dataset to make easier computation and comparison


library(terra)
library(sp)
library(sf)
library(raster)
library(dplyr)
library(rgeos)
library(tidyverse)
library(RColorBrewer)
library(spatstat)
library(rspatial)
library(rgdal)
library(spatstat)
# Set wd
wd <- "C:/GIS/UNFPA GIS/HTI PopGrid"
setwd(wd)

# Bring layers to test different methodologies

# Building footprints points Downsized for experiment
bf <- vect("C:/GIS/UNFPA GIS/HTI PopGrid/Interpolation/layers/test_bf_pts_26718.gpkg")

# downsize of the 
abs_ds <- rast("C:/GIS/UNFPA GIS/HTI PopGrid/Interpolation/layers/abs_ds.tif")

# base raster
br <- rast()
ext(br) <- ext(abs_ds)
crs(br) <- crs(abs_ds)
res(br) <- 100

# PREPARE LAYERS
## To economize computing resources group BF into raster with building count then to centroids. 
## And extract the BF that are out of the HH listing dataset coverage, only interested in HHlisting dataset gaps
bfcount <- terra::rasterize(bf,br,'count', fun=sum)

## what pixels are out abs raster?
## Inverse mask to extract pixels out of the HHlisting coverage that we will use to fill hhlisting gaps
bf_out <- terra::mask(bfcount,abs_ds, inverse=TRUE,filename="Interpolation/layers/bf_out.tif", overwrite=T)

## convert to points (and correct class) to perform "get data from nearest point 
bf_out_p <- as.points(bf_out) %>% sf::st_as_sf()
abs_ds_p <- as.points(abs_ds) %>% sf::st_as_sf()


# FIRST APPROACH: Get value from the closest pixel.----

abs_ds_out<- st_join(bf_out_p, abs_ds_p, join=st_nearest_feature)

## st_write(abs_ds_out,"Interpolation/layers/abs_ds_out.gpkg",delete_layer = T)
abs_ds_out <- vect(abs_ds_out)
abs_ds_out_r <- rasterize(abs_ds_out, br,'abs_ds',fun=sum, filename="Interpolation/layers/abs_ds_out_r.tif")

## Merge hhlisting abs with interpolated abs
abs_inter_1 <-merge(abs_ds,abs_ds_out_r,filename="Interpolation/layers/abs_inter_1.tif")

# SECOND APPROACH: Get distance to closest point, use that distance to set buffer and then calculate mean of pixels within that buffer----

## get distance of the closest point/pixel
## we are using rgeos> coerce to sp
bf_out_p_sp <- as(bf_out_p ,"Spatial")
abs_ds_p_sp <- as(abs_ds_p ,"Spatial")

## Find closest point from bf to abs
closest <- apply(gDistance(bf_out_p_sp,abs_ds_p_sp,byid = T ),MARGIN=2, FUN = which.min)
## set the value 
bf_out_p_sp$closest_value <- abs_ds_p_sp$value[closest]
## THIS WILL TAKE FOREVER, TRYING INTERPOLATION INSTEAD

# INTERPOLATION METHODS ----
# IDW 

## correct classes for spatial input
bf_out_p_sp <- as(bf_out_p ,"Spatial")
abs_ds_p_sp <- as(abs_ds_p ,"Spatial")
abs_ds_p_df <- as.data.frame(abs_ds_p_sp)
names(abs_ds_p_df)[names(abs_ds_p_df) == "coords.x1"] <- "x"
names(abs_ds_p_df)[names(abs_ds_p_df) == "coords.x2"] <- "y"

## I need to convert br into image to set observation window

## Observation window
# obs_window <- owin(abs_ds_p_sp@bbox[1,],abs_ds_p_sp@bbox[2,])
obs_window <- owin(ext(bfcount)[c(1,2)],ext(bfcount)[c(3,4)]) # we define observation window using generic raster


# mask <- as.mask(obs_window, dimyx = c(ncol(bfcount),nrow(bfcount)))
# mask <- as.mask(obs_window, eps = 99)
# mask <- as.mask(obs_window, dimyx = c(100,100))
# obs_window <- as.owin(mask)
mask <- owin2mask(mask, op="cover")
mask <- owin2mask(obs_window, op="cover",dimyx = c(ncol(bfcount),nrow(bfcount)))
mask <- spatstat.geom::as.mask(obs_window, dimyx = c(nrow(bfcount),ncol(bfcount)))


# obs_window <- owin(mask = (as.mask(as.matrix(br))))

## turn abs points into ppp object with the average population size data we want to interpolate
ppp_abs <- ppp(abs_ds_p_df$x, abs_ds_p_df$y, 
               marks = abs_ds_p_df$abs_ds, window = obs_window)
plot(ppp_abs)

## IDW interpolation (raster)
# xy<-xyFromCell(bfcount,obs_window) # coordinates of the cells

# mask <- as.mask(obs_window,dimyx =c(ncol(bfcount),nrow(bfcount)))

idw_abs_r <- idw(ppp_abs, power=0.2
               , at="pixels", as.mask = mask)

plot(idw_abs_r,
     col=(heat.colors(20)),
     main= "adb interpolated in ds region")

## IDW interpolation (points)
idw_abs_p <- idw(ppp_abs, power=0.2
                 , at="points")

plot(idw_abs_p,
     col=(heat.colors(20)),
     main= "adb interpolated in ds region")

## Find the optimal power to adjust IDW interpolation
## calculate Mean Square Error
Metrics::mse(ppp_abs$marks, idw_abs_p)

## plotting different power option IDW interpolations

par(mfrow = c(2,2))
plot(idw(ppp_abs, power=0.1, at="pixels"),col=heat.colors(20), main="power = 0.1")
plot(idw(ppp_abs, power=0.5, at="pixels"),col=heat.colors(20), main="power = 0.5")
plot(idw(ppp_abs, power=1, at="pixels"),col=heat.colors(20), main="power = 1")
plot(idw(ppp_abs, power=2, at="pixels"),col=heat.colors(20), main="power = 1.5")

## Cross Validate results and select the minimum MSE
## First attack broad ranges to save resources
powers <- seq(0.1,5,0.5)
mse_result <- NULL

for (power in powers){
  cv_idw <- idw(ppp_abs, power=power, at="points")
  mse_result <-  c(mse_result,
                   Metrics::mse(ppp_abs$marks,cv_idw))
}
optimal_power <- powers[which.min(mse_result)]
optimal_power
plot(powers, mse_result)

## Same process in a smaller range (from 1.5 up to 2.5)
powers <- seq(1.5,2.5,0.1)
mse_result <- NULL

for (power in powers){
  cv_idw <- idw(ppp_abs, power=power, at="points")
  mse_result <-  c(mse_result,
                   Metrics::mse(ppp_abs$marks,cv_idw))
}
optimal_power <- powers[which.min(mse_result)]
optimal_power
plot(powers, mse_result)

## then we run the interpolation using the optimal power parameter
idw_abs_r <- idw(ppp_abs, power=1.9
                 , at="pixels")

plot(idw_abs_r,
     col=(heat.colors(20)),
     main= paste0("adb interpolated in ds region - Optimal power = "))

## to convert result into raster and export to tif 
# idw_abs_r_sr <- rast(idw_abs_r,crs = crs(bfcount), resolution = res(bfcount))
# terra::writeRaster(idw_abs_r_sr, "Interpolation/layers/idw_abs_r_sr.tif", overwrite = T)

# NO WAY TO CONTROL THE OUTPUT RESOLUTION USING SPATSTAT PACKAGE (MASK FUNCTION), SEEMS LIKE A BUG OR SOMETHING, AT LEAST WE HAVE CALIBRATED MODEL TO FIND OPTIMAL POWER

# IDW USING TERRA
library(gstat)
idw <- gstat(id = "abs", formula = abs_ds~1, locations = ~x+y, data = abs_ds_p_df,
              nmax = 7, set=list(idp=1.9)) # it works but export into "gstat" "list" 

idw2 <- gstat::idw(formula = abs_ds~1, locations = abs_ds_p_sp, nmax =7, debug.level=0)

idw2 <- gstat::idw(formula = abs_ds_p_sp$abs_ds~1, locations = abs_ds_p_sp)
idw_abs <- interpolate(br,abs_ds_p_df)

##### test the rest of the operations in ds here -------------------------------------------------------------

## we have already the IDW tuned in spatstat, we run idw in gstat envrinonment to control output resolution

## IDW interpolation using gstat library

# PREPARE LAYERS
## I need the points with data that are going to be interpolated across an empty raster
## Points with Average Building Size data
# abs_p <- as.points(abs) %>% sf::st_as_sf()
# abs_p_sp <- as(abs_p ,"Spatial")
# abs_p_df <- as.data.frame(abs_p_sp)

# create raster using raster library from the begining otherwise IDW does not work (Raster class compatible with gstat)
# r <- 100
# e <- raster::extent(abs_p_sp )
# br<- raster::raster(e,res=c(r,r))
# crs(br) <- crs(abs_p_sp)
# br[] <- 1:raster::ncell(br)
# class(br)

br_ds_sf <- as(br_dsr,"SpatialPointsDataFrame") %>%  st_as_sf() # convert raster into simplefeature of points to include into the interpolation 

idw <- gstat::idw(layer~1, # formula that includes the column we want to interpolate
                  abs_ds_p_sp, # locations and data 
                  newdata = br_ds_sf, 
                  nmax = 25, # number of pixels around taken
                  maxdist = 5000, # max distance to get data from
                  idp=optimal_power)
idw
# rasterize in terra, easier to coerce classes

idw_r <- terra::rasterize(vect(idw),rast(br_dsr),'var1.pred',fun = max) #!!!use br_dsr that is the raster we have used for the raster operations in raster lib
plot(idw_r,
     col=(heat.colors(20)),
     main= paste0("adb interpolated in ds region - Optimal power = op"))

writeRaster(idw_r,"GIS_Model/GIS_intermediate/idw_r_ds.tif", overwrite=T)
writeRaster(bf_ds,"GIS_Model/GIS_intermediate/bf_ds.tif", overwrite=T)
## to get population we multiply ABS interpolated by the hh count from bf
pop_abs_r <- idw_r*bf_ds
writeRaster(pop_abs_r,"GIS_Model/GIS_intermediate/pop_abs_r.tif", overwrite=T)





