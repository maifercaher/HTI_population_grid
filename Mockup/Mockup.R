# HAITI POPULATION GRID 2022 --- Mockup OPERATIONS 
# Luis de la Rua 11/2022

# Set wd
wd <- "C:/GIS/UNFPA GIS/HTI PopGrid"
setwd(wd)

# load libraries
library(sf)
library(sp)
library(maptools)
library(spData)
library(rgeos)
library(dplyr)
library(magrittr)
library(terra) # using this library for spatial operations, libraries raster and rgdal will be deprecated in 2023
library(ggplot2)
library(exactextractr)

# LOADING DATA INPUT
## Admin boundaries with population (for the moment we use this until we get first Pop projections at Sub Section level)
ab <- vect("C:/GIS/UNFPA GIS/HTI PopGrid/HHloc18_20/layers/hti_admbnda_adm3_cnigs_20181129_26718.shp") # this layer
ab <- ab[,-c(11:875)] # remove fields we are not using
head (ab)

## Settlement Footprint consolidated
bf <- vect("Interpolation/layers/MS_BF_26718.gpkg")
# bf <- project(bf,ab) # convert to UTM 18N datum NAD27

## HH locations from IHSI listing cleaned in previous R process
hh <- vect("HHloc18_20/layers/hhloc_clean_26718.gpkg") 

# PREPARING RASTER FIRST DATASETS DATA 

## Generate base raster for Haiti #how to do this Terra
br <- rast()
ext(br) <- ext(ab)
crs(br) <- crs(ab)
res(br) <- 100

## Create raster with Building Density (we are using this population distribution factor in Communes not covered in the listing, PaP and Ouest area mainly)
### Convert BF into centroids and rasterize summing the number of points
bf_pts <- centroids(bf)
bf_pts$count <- 1 # field that will simplify counts in function below
bfcount <- terra::rasterize(bf_pts, br, 'count', fun=sum)
terra::writeVector(bf_pts, "Interpolation/layers/bf_pts.gpkg", overwrite = T)
terra::writeRaster(bfcount, "Interpolation/layers/bfcount.tif", overwrite = T)

spplot(bfcount)

## Create raster with the average building size (We are using this population distribution factor in Communes well covered by the listing)
### Prepare fields
hh$pop <- as.numeric(hh$bsize)
hh <- subset(hh, hh$pop != 0) # remove empty dwellings from the dataset
hh$count <- 1 # field that will simplify counts in function below

### Convert vector to raster and operate rasters
hhcount <- rasterize (hh, br, 'count', fun=sum) # number of hh 
hhpop <- rasterize(hh, br, 'pop', fun=sum) # total population per pixel
abs <- hhpop / hhcount                    # Calculate Average Building Size (average population / building)
writeRaster(abs, "Interpolation/layers/abs.tif", overwrite = T)

spplot(abs)
# ggplot(as.data.frame(abs)) +
#   geom_histogram(aes(layer),bins=500)

# IDENTIFY WHAT COMMUNES WILL BUILDINGS OR POPULATION FACTOR

## Calculate number of pixels per commune to asses coverage of both dataset
### For BF footprints 
ab_bf <- extract(bfcount, ab,
                 fun = sum,
                 ID = T,
                 bind = T,
                 cells = T,
                 na.rm = T,
                 df = T,
                 exact = F)
ab_bf <- cbind(ab, ab_bf)
ab_bf

writeVector(ab_bf,"Interpolation/layers/ab_bf.gpkg", overwrite = T)

### Set threshold and identify what communes are not covered
thr <- 1000 # testing different thresholds
ab_nobf <- subset (ab_bf, ab_bf$lyr.1 < thr)  

## get sum per commune and then divide each pixel by this value. (Create raster SUM by Commune and operate)
ab_hhcount <- extract(hhcount, ab,
                      fun = sum,
                      ID = T,
                      bind = T,
                      cells = T,
                      na.rm = T,
                      df = T,
                      exact = F)
ab_hhcount <- cbind(ab, ab_hhcount)
ab_hhcount

writeVector(ab_hhcount,"Interpolation/layers/ab_hhcount.gpkg", overwrite = T)

### Set threshold and identify what communes are not covered
thr2 <- 1000 # testing different thresholds
ab_nohh <- subset (ab_hhcount, ab_hhcount$lyr.1 < thr2)  

# CALCULATE POPULATION DISRTRIBUTION AND COEFICIENTS - RATIOS FOR BUILDING AND HH LISTING RASTERS

## For the Building Footprint Option ----
ab_bf  # we have the total buildings / commune
## rasterize ab_bf
ab_bf_r <- rasterize(ab_bf, br, 'lyr.1', fun=max )
spplot(ab_bf_r)

## writeRaster(ab_bf_r,"Interpolation/layers/ab_br_r.tif", overwrite = T)
## Calculate ratio 
bf_ratio <- bfcount / ab_bf_r
spplot(bf_ratio)
writeRaster(bf_ratio,"Interpolation/layers/bf_ratio.tif", overwrite = T)

## check that ratio calculation is right ## some issues with zonal statistics for this, checked in QGIS and works allright all pixels within same area sum ~1 ----
ab_sf <- sf::st_as_sf(ab)

bf_ratiocheck <- exact_extract(bf_ratio, ab_sf, 
                               fun = 'sum',
                               full_colnames = TRUE,
                               force_df = F) %>% as.data.frame()

bf_ratiocheck <- cbind(ab, bf_ratiocheck)
bf_ratiocheck

# writeVector(bf_ratiocheck,"Interpolation/layers/bf_ratiocheck.gpkg", overwrite = T)

## For the HH Listing Option ----
## For this option we still need to solve misalignment, population interpolation to include pop data into the bf we use to fill listing gaps.

hhpop  # we have mockup of what the dataset would look like
## calculate population counts at commune level from listing dataset
ab_pop <- exact_extract(hhpop, ab_sf, 
                        fun = 'sum',
                        full_colnames = TRUE,
                        force_df = F) %>% as.data.frame()

ab_pop <- cbind(ab,ab_pop)
ab_pop

ab_pop_r <- rasterize(ab_pop, br, '.', fun=max ) # Covnert ab with pop totals into raster to operate and get ratios

hhpop_ratio <- hhpop / ab_pop_r # Calculate ratio poppixel / popcommune
spplot(hhpop_ratio)

writeRaster(hhpop_ratio,"Interpolation/layers/hhpop_ratio.tif", overwrite = T)

# DISTRIBUTE POPULATION
## Generate Raster Pop by Commune (with population 2015) We need to iterate this across age ranges and sex and loop operations below
ab_prpop <- ab <- vect("Interpolation/layers/Haiti_adm3_2000_2015_32618.gpkg")

ab_prpop_r <- rasterize(ab_prpop, br, 'BTOTL_2015', fun = max )
spplot(ab_prpop_r)
## writeRaster(bf_ratio,"Interpolation/layers/bf_ratio.tif", overwrite = T)

## Operate with both pop distribution coefficient grids
## HH Listing approach
pgrid_hhpop <- hhpop_ratio * ab_prpop_r
spplot(pgrid_hhpop)
writeRaster(pgrid_hhpop,"Interpolation/layers/pgrid_hhpop.tif", overwrite = T)

## BF approach
pgrid_bf <- bf_ratio * ab_prpop_r
spplot(pgrid_bf)
writeRaster(pgrid_bf,"Interpolation/layers/pgrid_bf.tif", overwrite = T)


