# HAITI POPULATION GRID 2022 --- HTI Population Grid
# MODEL V1 24/11/2022
# Luis de la Rua

# Set wd
wd <- "C:/GIS/UNFPA GIS/HTI PopGrid"
setwd(wd)

# load libraries
library(sf)
library(sp)
library(spData)
library(rgeos)
library(raster)
library(dplyr)
library(terra) # using this library for spatial operations, libraries raster and rgdal will be deprecated in 2023
library(ggplot2)
library(exactextractr)
library(gstat)

# LOADING DATA INPUT

## Admin boundaries with population (for the moment we use this until we get first Pop projections at Sub Section level)
### this layer has updated sub section codes, we use it to tabulate population output
ab_ihsi <- vect("C:/GIS/UNFPA GIS/HTI PopGrid/HHloc18_20/layers/hti_admbnda_adm3_cnigs_20181129_26718.shp") 

ab <- vect("GIS_Model/GIS_Input/Haiti_adm3_2000_2022_26718.gpkg") # this one is not up to data but contains population prj to test model
## ab <- ab[,-c(11:820)] # remove fields we are not using
head (ab)

# PREPARING RASTER FIRST DATASETS DATA 

## Settlement Footprint consolidated > Process in building integration folder in git repo 
bf <- rast("GIS_Model/GIS_Input/HTI_MS_OSM_HH.tif") # Pixel values correspond to number of buildings
## Set CRS for the entire mode
crs <- crs(ab_ihsi)

## Reproject ab (HHlisting not necessary as it comes with the right projection)
crs(ab) <- crs

## Generate base raster for Haiti how to do this Terra
br <- rast()
ext(br) <- ext(bf)
crs(br) <- crs(ab)
res(br) <- 100
## From ab we are going to create a custom boundaries with "pixelized" borders to avoid limits misalignment issues when making raster extractions

ab_px <- terra::rasterize(ab,br, field= 'OTHER_CODE', fun=max) %>%  as.polygons()
spplot(ab_px)
# writeVector(ab_px,"GIS_Model/GIS_intermediate/ab_px.gpkg", overwrite = T)

## also for the table we use for the output
ab_ihsi_px <- terra::rasterize(ab_ihsi,br, field= 'ADM3_PCODE', fun=max) %>% as.polygons()

# writeVector(ab_ihsi_px,"GIS_Model/GIS_intermediate/ab_ihsi_px.gpkg", overwrite = T)
spplot(ab_ihsi_px)

## HH locations from IHSI listing cleaned in previous R process
hh <- vect("GIS_Model/GIS_Input/hhloc_clean_26718.gpkg") 
## remove empty hh from hh as it will interfere into the Average Building Size calculation
hh <- subset(hh, hh$bsize != 0)

## Compare all input crs and set warning (To Do)

# CALCULATE POPULATION DISRTRIBUTION AND RATIOS FOR BUILDING AND POPULATION CRITERIA -------------------------------------

## BUILDING DENSITY CRITERIA ---------------------------------------------------------------------------------------------       

## Generate population distribution raster

## Calculate sum of building count per admin unit
ab_bf <- exact_extract(bf, sf::st_as_sf(ab_px),
                       fun = c('count','sum'), # Don't think we need count and also is counting pixels with 0 values.
                       full_colnames = TRUE,
                       force_df = F) %>% as.data.frame()

ab_bf <- cbind(ab_px,ab_bf)
head(ab_bf)

## Rasterize ab_bf
ab_bf_r <- rasterize(ab_bf, br, 'sum.HTI_MS_OSM_HH', fun=max )
spplot(ab_bf_r)

## Calculate BF ratio 
bf_ratio <- bf / ab_bf_r
spplot(bf_ratio)
hist(bf_ratio,breaks = 1000) # display how data looks like
writeRaster(bf_ratio,"GIS_Model/GIS_intermediate/bf_ratio.tif", overwrite = T)

## check that ratio calculation is right ## some issues with zonal statistics for this, checked in QGIS and works allright all pixels within same area sum ~1 
bf_ratiocheck <- exact_extract(bf_ratio, sf::st_as_sf(ab_px), 
                               fun = 'sum',
                               full_colnames = TRUE,
                               force_df = F) %>% as.data.frame()

bf_ratiocheck <- cbind(ab_px, bf_ratiocheck)
head(bf_ratiocheck)
View(as.data.frame(bf_ratiocheck))

## POPULATION DENSITY CRITERIA ------------------------------------------------------------------------------------------- 
## We are using population interpolated surface combining HH listing data and Building footprints

## PREPARING LAYERS
## Convert HH listing dataset to raster
### Convert vector to raster and operate rasters
hh$count <- 1 # we focus on building counts not on HH

## test rasterize using raster library 
## coerce inputs and variables to work on raster lib
hh_sp<- as(hh,"Spatial")
hh_sp$bsize <- as.numeric(hh_sp$bsize) # coerce bsize (cointains pop/building) to numeric 

## gen br_r blank raster in raster library
r <- 100
e <- raster::extent(as(ab,"Spatial"))
br_r<- raster::raster(e,res=c(r,r))
crs(br_r) <- crs(hh)

hh_count <- raster::rasterize(hh_sp,br_r,'count',fun=sum)
hh_pop <- raster::rasterize(hh_sp,br_r ,'bsize', fun=sum) # lets try Raster library

abs <- hh_pop/hh_count # Calculate Average Building Size (average population / building) 

hh_pop <- rast(hh_pop) # turn back to terra
hh_count <- rast(hh_count) 
abs <- rast(abs) 

writeRaster(hh_count, "GIS_Model/GIS_intermediate/hh_count.tif", overwrite = T)
writeRaster(hh_pop, "GIS_Model/GIS_intermediate/hh_pop.tif", overwrite = T)
writeRaster(abs, "GIS_Model/GIS_intermediate/abs.tif", overwrite = T)

# PREPARE LAYERS
## I need the points with data that are going to be interpolated across an empty raster (From HH listing dataset)
## Points with Average Building Size data
## Coercing to work in gstat (there must be a cleaner way to do this)
abs_p <- as.points(abs) %>% sf::st_as_sf()
abs_p_sp <- as(abs_p ,"Spatial")

# Create point grid to operate interpolation

br_sf <- as(raster(br),"SpatialPointsDataFrame") %>%  st_as_sf() # convert raster into simplefeature of points to include into the interpolation 

idw <- gstat::idw(layer~1, # formula that includes the column we want to interpolate
                  abs_p_sp, # locations and data 
                  newdata = br_sf, 
                  nmax = 25, # number of pixels around taken 
                  maxdist = 5000, # max distance to get data from / We need to tun this a bit 
                  idp=optimal_power) 
idw

# Interpolated ABS surface
# rasterize in terra, easier to coerce classes 

abs_inter <- terra::rasterize(vect(idw),rast(br),'var1.pred',fun = max)
plot(abs_inter,
     col=(heat.colors(20)),
     main= paste0("abs interpolated in ds region - Optimal power = ",optimal_power))

writeRaster(abs_inter,"GIS_Model/GIS_intermediate/abs_inter.tif", overwrite=T)

# Calculate population (multiply bf (building counts) * abs_inter (population per building))

pop_abs_inter <- abs_inter*bf

writeRaster(pop_abs_inter,"GIS_Model/GIS_intermediate/pop_abs_inter.tif", overwrite=T)
#curiosity how much total population we get?
global(pop_abs_inter,"sum",na.rm=T)

## Spatial distribution using Population ratio ------------------
## Rasterize projected population dataset linked to admin units layer (to generate by sex and age we would just need to loop over the population fields in ab layer (TO DO))

## Calculate total pop by admin unit using this model

ab_abs <- exact_extract(pop_abs_inter, sf::st_as_sf(ab_px),
                       fun = c('count','sum'), # Don't think we need count and also is counting pixels with 0 values.
                       full_colnames = TRUE,
                       force_df = F) %>% as.data.frame()

ab_abs <- cbind(ab_px,ab_abs)
head(ab_abs)
sum(ab_abs$sum.lyr.1) ## less population we are missing coastal pixels

## Rasterize ab_abs
ab_abs_r <- rasterize(ab_abs, br, 'sum.lyr.1', fun=max )
spplot(ab_abs_r)

## Calculate ABS ratio 
abs_ratio <- pop_abs_inter/ab_abs_r
spplot(abs_ratio)
hist(abs_ratio,breaks = 1000) # display how data looks like
writeRaster(abs_ratio,"GIS_Model/GIS_intermediate/abs_ratio.tif", overwrite = T)

## check that ratio calculation is right ## some issues with zonal statistics for this, checked in QGIS and works all right, all pixels within same area sum ~1 
abs_ratiocheck <- exact_extract(abs_ratio, sf::st_as_sf(ab_px), 
                               fun = 'sum',
                               full_colnames = TRUE,
                               force_df = F) %>% as.data.frame()

abs_ratiocheck <- cbind(ab_px, abs_ratiocheck)
head(abs_ratiocheck)



## GENERATE COMBINED ABS/BF RATIO GRID. SELECT ADMIN UNITS USING ONE OF THE METHODS BASED ON DATA COVERAGE ------------------------------
## classify admin units by comparing number of pixels from HHlisting with the pixels of the BF grid
## count hh listing pixels per admin unit
ab_px_abs <- exact_extract(abs, # Pixels with data from HH listing
                            sf::st_as_sf(ab_px),
                            fun='count', full_colnames=T,
                            force_df = F) %>% as.data.frame()

ab_px_abs <- cbind(ab_px,ab_px_abs)
ab_px_abs$abs_ount <- ab_px_abs$.
head(ab_px_abs)
## count bf pixels (>0) per admin unit
bfclamp<- clamp(bf, lower=1,values=F) # remove pixels=0
ab_px_bf <- exact_extract(bfclamp, # Pixels with data from HH listing
                           sf::st_as_sf(ab_px),
                           fun='count', full_colnames=T,
                           force_df = F) %>% as.data.frame()

ab_px_bf <- cbind(ab_px,ab_px_bf)
ab_px_bf$bf_ount <- ab_px_bf$.
head(ab_px_bf)
# Merge tables and calculate number of pixels ratio
ab_px_meth <- merge(ab_px_abs,ab_px_bf, by = 'OTHER_CODE')
ab_px_meth <-cbind(ab_px,ab_px_meth)
ab_px_meth$px_per <-  ab_px_meth$abs_ount/ab_px_meth$bf_ount
ab_px_meth <- ab_px_meth[,-c(2,3,5)] # clean columns
View(as.data.frame(ab_px_meth))

writeVector(ab_px_meth,"GIS_Model/GIS_intermediate/ab_px_meth.gpkg", overwrite = T)

## Threshold to separate methodologies
ths <- 0.4
ab_meth_abs <- subset(ab_px_meth,ab_px_meth$px_per>ths) # this is the layer we use to mask each of the ratio rasters

## Admin units selected to abs_ratio criteria
abs_ratio_sel <- mask(abs_ratio, ab_meth_abs,inverse=F)
bf_ratio_sel <- mask(bf_ratio,ab_meth_abs,inverse=T)

# Combine both pieces
pop_ratio <- terra::merge(abs_ratio_sel,bf_ratio_sel) 
plot(pop_ratio)

writeRaster(abs_ratio_sel,"GIS_Model/GIS_intermediate/abs_ratio_sel.tif", overwrite = T)
writeRaster(bf_ratio_sel,"GIS_Model/GIS_intermediate/bf_ratio_sel.tif", overwrite = T)
writeRaster(pop_ratio,"GIS_Model/GIS_intermediate/pop_ratio.tif", overwrite = T)
## multiply population fields included in ab (projected pop by age and sex) by abs ratio to get final population grid based on both abs and bf approach

# GENERATE POPULATION GRIDS BY AGE GROUP AND SEX -----------------------------

# convert population fields to numeric
ab_sf <- st_read("GIS_Model/GIS_Input/Haiti_adm3_2000_2022_26718.gpkg")

ab_sf <- ab_sf %>% mutate_at(c(6:8),as.numeric) # Convert pop fields to numeric and coerce to terra
ab <- vect(ab_sf)
head(ab)

cols <- names(ab[,6:8]) ## define the fields we are looping over 
pop_r <- NULL # Define null outputs to make the loop work
pop_abs_r <- NULL
date<-format(Sys.time(),'_%Y%m%d') # set date for file names

# For ABS criteria
for (i in cols){
  pop_r <- terra::rasterize(ab,br,   # Generate population adm unit raster for each age/sex population range 
                            field= i,fun=max)
  pop_abs_r <- pop_r*abs_ratio # Multiply pop raster by spat distribution ratio
  writeRaster(pop_abs_r,paste0("GIS_Model/GIS_output/pop_abs_",i,date,".tif"), overwrite=T)
}

# For Building Density criteria
for (i in cols){
  pop_r <- terra::rasterize(ab,br,   # Generate population adm unit raster for each age/sex population range 
                            field= i,fun=max)
  pop_bf_r <- pop_r*bf_ratio # Multiply pop raster by spat distribution ratio
  writeRaster(pop_bf_r,paste0("GIS_Model/GIS_output/pop_bf_",i,date,".tif"), overwrite=T)
}

# TRANSLATE TO RECENT ADMIN LIMITS FRAMEWORK
# Calculate Population aggregating from population grids
## for the moment we open all rasters in output
r_list <- list.files(path="GIS_Model/GIS_output/",pattern = paste0(date,".tif"), #stacking only the last rasters produced filtering by date date
                     all.files=T, full.names=T)
stck <- stack(r_list) # stack all population rasters

ab_pop <- exact_extract(stck, sf::st_as_sf(ab_ihsi_px),
                       fun = 'sum', # Don't think we need count and also is counting pixels with 0 values.
                       full_colnames = TRUE,
                       force_df = F) %>% as.data.frame()

ab_pop <- cbind(ab_ihsi_px,ab_pop) # need to do some cleaning in the future
write.csv(ab_pop, paste0("GIS_Model/GIS_output/ab_pop",date,".csv"))

a<-sum(ab_pop$sum.pop_abs_FEMMES_20221206)
b<-global(pop_abs_r,"sum",na.rm=T)
a-b
# Check pixels not captured
out <- terra::mask(pop_abs_r,ab_ihsi_px, inverse=TRUE)
global(out,"sum",na.rm=T)
## STAND BY THIS UNTIL IHSI PROVIDES THE LATEST LAYER
# Check if differences between raster and proj pop tables -----

## Total population differences caused by raster vs vector boundaries when running zonal statistics at national level and by admin unit
# Compare zonal stats results with original population figures from table
abold_r_pop <- exact_extract(stck, sf::st_as_sf(ab_px),
                           fun = 'sum', # Don't think we need count and also is counting pixels with 0 values.
                           full_colnames = TRUE,
                           force_df = F) %>% as.data.frame()

abold_r_pop <- cbind(ab_px,abold_r_pop) # need to do some cleaning in the future
abold_r_pop <- merge(abold_r_pop,ab,by = 'OTHER_CODE')
abold_r_pop$dif <- (abold_r_pop$sum.pop_abs_POPULATION_20221206 - abold_r_pop$POPULATION)*100/abold_r_pop$POPULATION
abold_r_pop$dif
write.csv(abold_r_pop, paste0("GIS_Model/GIS_output/abold_r",date,".csv"))



# Here we have to loop rasterization of admin unit layer per pop column, multiply per ratio and export to tif properly named.
# second model outcome is going to be the tables at section, Comune and Province level using the most recent admin codes and boundaries.
