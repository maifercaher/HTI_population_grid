library(gstat) # Use gstat's idw routine
library(sp)    # Used for the spsample function
library(terra)
##IDW test using other library
wd <- "C:/GIS/UNFPA GIS/HTI PopGrid"
setwd(wd)

# Bring layers to test different methodologies

# Building footprints points Downsized for experiment
bf <- vect("C:/GIS/UNFPA GIS/HTI PopGrid/Interpolation/layers/test_bf_pts_26718.gpkg")

# downsize of the 
abs_ds <- raster("C:/GIS/UNFPA GIS/HTI PopGrid/Interpolation/layers/abs_ds.tif")

# base raster In TERRA not necessary so far
# br <- rast()
# ext(br) <- ext(abs_ds)
# crs(br) <- crs(abs_ds)
# res(br) <- 100

# PREPARE LAYERS
## I need the points with data that are going to be interpolated across an empty raster
## Points with Average Building Size data
abs_ds_p <- as.points(abs_ds) %>% sf::st_as_sf()
abs_ds_p_sp <- as(abs_ds_p ,"Spatial")
abs_ds_p_df <- as.data.frame(abs_ds_p_sp)

# create raster using raster
r <- 100
e <- raster::extent(abs_ds_p_sp )
br<- raster::raster(e,res=c(r,r))
crs(br) <- crs(abs_ds_p_sp)
br[] <- 1:raster::ncell(br)
class(br)

br_sf <- as(br,"SpatialPointsDataFrame") %>%  st_as_sf() # convert raster into simplefeature of points to include into the interpolation 

idw <- gstat::idw(abs_ds~1, # Points to interpolate
                  abs_ds_p_sp, newdata = br_sf, idp=1.9)

# rasterize in terra, easier to coerce classes

# idw_sv <- vect(idw)
idw_r <- terra::rasterize(vect(idw),rast(br),'var1.pred', fun = sum)
plot(idw_r,
     col=(heat.colors(20)),
     main= paste0("adb interpolated in ds region - Optimal power = "))

writeRaster(idw_r,"C:/GIS/UNFPA GIS/HTI PopGrid/Interpolation/layers/idw_r.tif", overwrite=T)
