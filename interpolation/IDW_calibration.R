### Interpolation Method for Population - Haiti Population Grid 2022
### Luis de la Rua - November 2022

### Inverse Distance Weighted IDW interpolation method utilized to allocate population data from HH listing data to the 
### Building Footprint pixels utilized to fill the data gaps from the HH listing dataset.

## We Calibrate Interpolation finding the optimal power parameter in spatstat environment (using less resolution too to speed up the process)
## We run this in a portion of the full dataset Downsize "ds" region
## Libraries
library(terra)
library(sp)
library(sf)
library(dplyr)
library(tidyverse)
library(spatstat)
library(raster)

# Set wd
wd <- "C:/GIS/UNFPA GIS/HTI PopGrid"
setwd(wd)

## bring layers
## Admin boundaries with population (for the moment we use this until we get first Pop projections at Sub Section level)
ab_ihsi <- vect("C:/GIS/UNFPA GIS/HTI PopGrid/HHloc18_20/layers/hti_admbnda_adm3_cnigs_20181129_26718.shp")# this layer has updated sub section codes
ab <- vect("C:/GIS/UNFPA GIS/HTI PopGrid/Interpolation/layers/Haiti_adm3_2000_2015_32618.gpkg") # this one is not up to data but contains population prj to test model
ab <- ab[,-c(11:820)] # remove fields we are not using
head (ab)

## Settlement Footprint consolidated
bf <- rast("GIS_Model/GIS_Input/HTI_MS_OSM_HH.tif") # Pixel values correspond to number of buildings

## HH locations from IHSI listing cleaned in previous R process
hh <- vect("GIS_Model/GIS_Input/hhloc_clean_26718.gpkg")
## remove empty hh from hh as it will interfer into the Average Building Size calculation
hh <- subset(hh, hh$bsize != 0)

## Define Downsize area to calibrate model in order to speed up processing
ds <- vect("GIS_Model/GIS_Input/ds/ds_26718.gpkg")

## Define blank raster for zone ds
br_ds <- rast()
ext(br_ds) <- ext(ds)
crs(br_ds) <- crs(ds)
res(br_ds) <- 100

## mask datasets with ds for the calibration experiment
hh_ds <- terra::crop(hh, br_ds)
bf_ds <- terra::crop(bf, br_ds)
ext(bf_ds) <- ext(br_ds) #make both br and bf extents exactly the same

##PREPARING LAYERS
## Convert HH listing dataset to raster
### Convert vector to raster and operate rasters
hh_ds$count <- 1 # we focus on building counts not on HH
hh_ds$pop <- hh_ds$bsize # we focus on population / building not per HH
# hh_dscount <- rasterize (hh_ds, br_ds, 'count', fun=sum) # number of hh 
# hh_dspop <- rasterize(hh_ds, br_ds, 'pop', fun=sum) # total population per pixel WRONG RESULTS USING TERRA!!!

## Working on raster library - TERRA have issues with sum functions 
## coerce inputs and variables to work on raster lib
hh_dssp<- as(hh_ds,"Spatial")
hh_dssp$bsize <- as.numeric(hh_dssp$bsize) # coerce bsize to numeric
## gen br
r <- 100
e <- raster::extent(as(ds,"Spatial"))
br_dsr<- raster::raster(e,res=c(r,r))
crs(br_dsr) <- crs(hh_ds)

hh_dscount <- raster::rasterize(hh_dssp,br_dsr,'count',fun=sum)
hh_dspop <- raster::rasterize(hh_dssp,br_dsr ,'bsize', fun=sum) # lets try Raster library

abs_ds <- hh_dspop/hh_dscount # Calculate Average Building Size (average population / building) 

hh_dspop <- rast(hh_dspop) # turn back to terra
hh_dscount <- rast(hh_dscount) 
abs_ds <- rast(abs_ds) 

writeRaster(hh_dscount, "GIS_Model/GIS_intermediate/hh_dscount.tif", overwrite = T)
writeRaster(hh_dspop, "GIS_Model/GIS_intermediate/hh_dspop.tif", overwrite = T)
writeRaster(abs_ds, "GIS_Model/GIS_intermediate/abs_ds.tif", overwrite = T)

## Identify what pixels from the BF raster are out abs raster > Identify the pixels we need to assign populaiton data through interpolation
## Inverse mask to extract pixels out of the HHlisting coverage that we will use to fill hhlisting gaps
ext(bf_ds) <- ext(abs_ds) ## align rasters
bf_ds_out <- terra::mask(bf_ds,abs_ds, inverse=TRUE,filename="GIS_Model/GIS_intermediate/bf_ds_out.tif", overwrite=T)

## convert to points (and correct class) to perform "get data from nearest point 
bf_ds_out_p <- as.points(bf_ds_out) %>% sf::st_as_sf()
abs_ds_p <- as.points(abs_ds) %>% sf::st_as_sf()

## correct classes for spatial input in spatstat library
bf_ds_out_p_sp <- as(bf_ds_out_p ,"Spatial")
abs_ds_p_sp <- as(abs_ds_p ,"Spatial")
abs_ds_p_df <- as.data.frame(abs_ds_p_sp)
abs_ds_p_df$abs <- abs_ds_p_df$layer
names(abs_ds_p_df)[names(abs_ds_p_df) == "coords.x1"] <- "x"
names(abs_ds_p_df)[names(abs_ds_p_df) == "coords.x2"] <- "y"

## Set observation window
obs_window <- owin(ext(br_ds)[c(1,2)],ext(br_ds)[c(3,4)]) # we define observation window using the generic raster corners

## turn abs points into ppp object with the average population size data we want to interpolate
ppp_abs_ds <- ppp(abs_ds_p_df$x, abs_ds_p_df$y, 
               marks = abs_ds_p_df$abs, window = obs_window)
plot(ppp_abs_ds)

# IDW interpolation (raster) with optimal power calculated in IDW CALIBRATION
# Basic approach to see how this works 
idw_abs_r <- spatstat.explore::idw(ppp_abs_ds, power=0.2
                 , at="pixels")

plot(idw_abs_r,
     col=(heat.colors(20)),
     main= "adb interpolated")

## IDW interpolation (points)
idw_abs_p <- spatstat.explore::idw(ppp_abs_ds, power=0.2
                 , at="points")

plot(idw_abs_p,
     col=(heat.colors(20)),
     main= "adb interpolated")

## Find the optimal power to adjust IDW interpolation
## calculate Mean Square Error
Metrics::mse(ppp_abs_ds$marks, idw_abs_p)

## First attack broad ranges to save resources
powers <- seq(0.1,5,0.5)
mse_result <- NULL

for (power in powers){
  cv_idw <- spatstat.explore::idw(ppp_abs_ds, power=power, at="points")
  mse_result <-  c(mse_result,
                   Metrics::mse(ppp_abs_ds$marks,cv_idw))
}
optimal_power <- powers[which.min(mse_result)]
optimal_power
plot(powers, mse_result)

## Same process in a smaller range now we know where optimal could be (from 0 t0 2)
powers <- seq(0,2,0.01)
mse_result <- NULL

for (power in powers){
  cv_idw <- spatstat.explore::idw(ppp_abs_ds, power=power, at="points")
  mse_result <-  c(mse_result,
                   Metrics::mse(ppp_abs_ds$marks,cv_idw))
}
optimal_power <- powers[which.min(mse_result)]
optimal_power
plot(powers, mse_result)

## then we run the interpolation using the optimal power parameter
idw_abs_r <- spatstat.explore::idw(ppp_abs_ds, power=optimal_power
                 , at="pixels")

plot(idw_abs_r,
     col=(heat.colors(20)),
     main= paste0("adb interpolated in ds region - Optimal power = ", optimal_power))

rast(idw_abs_r)
writeRaster(rast(idw_abs_r), "GIS_Model/GIS_intermediate/idw_abs_r.tif", overwrite = T)
