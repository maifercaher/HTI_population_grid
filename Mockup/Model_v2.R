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
library(openxlsx)
library(tidyterra)

# LOADING AND PREPARING DATA INPUT-------

## Import last version of the pop projection
pop_table <- read.csv("GIS_Model/GIS_Input/HTI_sections_age_population_estimate_with_id.csv")

## Admin boundaries with population (for the moment we use this until we get first Pop projections at Sub Section level)
### this layer has updated sub section codes, we use it to tabulate population output / this only contains sections
ab_ihsi <- vect("GIS_Model/GIS_Input/SECTION_IHSI_20230119_26718.gpkg") 

# this layer has been edited to ensure it captures all pixels on coastal areas. We use it to extract all raster information
ab_ihsi_zs <- vect("GIS_Model/GIS_Input/SECTION_ESPACE_URBAIN_RURAL_26718_ZStat.gpkg") # after including Urban Rural we are using this layer to calculate totals from output popgrids
ab_ihsi_zs$ADM3_PCODE_UR <- ab_ihsi_zs$sect_ur # this layer contains sections and Urban Rural partitions
ab_ihsi_zs$ADM3_PCODE <- ab_ihsi_zs$Code_Sect_1

ab <- vect("GIS_Model/GIS_Input/Haiti_adm3_2000_2022_26718.gpkg") # this one is not updated but connects with population projections to allocate pop data to
                                                                  # each of the admin units    
## ab <- ab[,-c(11:820)] # remove fields we are not using
head (ab)

## PREPARING RASTER FIRST DATASETS DATA 

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
writeVector(ab_px,"GIS_Model/GIS_intermediate/ab_px.gpkg", overwrite = T)

## also for the table we use for the output
ab_ihsi_px <- terra::rasterize(ab_ihsi,br, field= 'ADM3_PCODE', fun=max) %>% as.polygons()
writeVector(ab_ihsi_px,"GIS_Model/GIS_intermediate/ab_ihsi_px.gpkg", overwrite = T)
spplot(ab_ihsi_px)

## also for the ab_ihsi_zs that has been modified to include all pixels
ab_ihsi_zs_px <- terra::rasterize(ab_ihsi_zs,br, field= 'ADM3_PCODE_UR', fun=max) %>% as.polygons()
writeVector(ab_ihsi_px,"GIS_Model/GIS_intermediate/ab_ihsi_px_zs.gpkg", overwrite = T)
spplot(ab_ihsi_zs_px)


## HH locations from IHSI listing cleaned in previous R process
hh <- vect("GIS_Model/GIS_Input/hhloc_clean_26718.gpkg") 
## Remove empty hh from hh as it will interfere into the Average Building Size calculation
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

## Calculate population (multiply bf (building counts) * abs_inter (population per building))

pop_abs_inter <- abs_inter*bf

writeRaster(pop_abs_inter,"GIS_Model/GIS_intermediate/pop_abs_inter.tif", overwrite=T)
## curiosity how much total population we get?
global(pop_abs_inter,"sum",na.rm=T)

## Spatial distribution using Population ratio 
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

#Check that ratios are still right and sum 1 for each admin unit
pop_ratiocheck <- exact_extract(pop_ratio, sf::st_as_sf(ab_px), 
                                fun = 'sum',
                                full_colnames = TRUE,
                                force_df = F) %>% as.data.frame()

pop_ratiocheck <- cbind(ab_px, pop_ratiocheck)
head(pop_ratiocheck)
View(as.data.frame(pop_ratiocheck))
write.csv(pop_ratiocheck,"GIS_Model/GIS_intermediate/pop_ratiocheck.csv")

## CORRECTION FACTOR TO ADJUST TOTAL POPULATIONS
## From the pop_ratio checks we use this factor to erase the rounding issues that provoques differences in the totals. Areas where ratio does not add 1.0000, 
## will produce little differences when calculating the totals, we will divide the Population grids by this correction factors to adjust the totals.
pop_ratio_cf <- terra::rasterize(pop_ratiocheck,br,field='.',fun='max')
plot(pop_ratio_cf)
writeVector(pop_ratiocheck,"GIS_Model/GIS_intermediate/pop_ratiocheck.gpkg", overwrite = T)
writeRaster(pop_ratio_cf,"GIS_Model/GIS_intermediate/pop_ratio_cf.tif", overwrite = T)

## multiply population fields included in ab (projected pop by age and sex) by abs ratio to get final population grid based on both abs and bf approach

# GENERATE POPULATION GRIDS BY AGE GROUP AND SEX -----------------------------

# We used Pixelized admin boundaries and connect to population table ab_px
ab_pop <- merge(ab_px,pop_table,
      by.x='OTHER_CODE', by.y='ID_SECTION')

## Calculate modeled population grid multiplying spatial distribution ratio * Projected population in each Admin unit
## Doing this for all the sex and age groups available
### Prepare loop
cols <- names(ab_pop[,5:58]) ## define the fields we are looping over 
pop_r <- NULL # Define null outputs to make the loop work
pop_grid <- NULL
date<-format(Sys.time(),'_%Y%m%d') # set date for file names

for (i in cols){
  pop_r <- terra::rasterize(ab_pop,br,   # Generate population adm unit raster for each age/sex population range 
                            field= i,fun=max)
  pop_grid <- (pop_r*pop_ratio)/ pop_ratio_cf # Multiply pop raster by spat distribution ratio and include the CORRECTION FACTOR too
  print(global(pop_grid,"sum",na.rm=T))
  writeRaster(pop_grid,paste0("GIS_Model/GIS_output/pop_grid_",i,date,".tif"), overwrite=T)
}

# TRANSLATE TO RECENT ADMIN LIMITS FRAMEWORK
# Calculate Population aggregating from population grids
## for the moment we open all rasters in output
r_list <- list.files(path="GIS_Model/GIS_output/",pattern = paste0(date,".tif"), #stacking only the last rasters produced filtering by date date
                     all.files=T, full.names=T)
stck <- stack(r_list) # stack all population rasters

ab_pop_zs <- exact_extract(stck, sf::st_as_sf(ab_ihsi_zs_px),
                        fun = 'sum', # Don't think we need count and also is counting pixels with 0 values.
                        full_colnames = TRUE,
                        force_df = F) %>% as.data.frame()

ab_pop_zs <- cbind(ab_ihsi_zs_px,ab_pop_zs) # need to do some cleaning in the future

# Rename variables

ab_pop_zs <- as.data.frame(ab_pop_zs) %>% 
  dplyr::rename_all(
    ~stringr::str_replace_all(.,"sum.pop_grid_","")
  )  %>% 
  dplyr::rename_all(
    ~stringr::str_replace_all(.,date,"")
  ) %>% 
  dplyr::rename_all(
    ~stringr::str_replace_all(.,"5_9","05_9") # Rename 5_9 age groups for reordering
  )  %>% 
  select(sort(colnames(.)) # reorder fields alphabetically 
  )
  
ab_pop_zs <- ab_pop_zs %>% select(sort(colnames(.))) # reorder fields alphabetically 

head(ab_pop_zs)

write.csv(ab_pop_zs, paste0("GIS_Model/GIS_output/ab_pop_zs",date,".csv"))


## Check if differences between raster and proj pop tables -----

# Compare zonal stats results with original population figures from table
abold_r_pop <- exact_extract(stck, sf::st_as_sf(ab_pop),
                             fun = 'sum', # Don't think we need count and also is counting pixels with 0 values.
                             full_colnames = TRUE,
                             force_df = F) %>% as.data.frame()

abold_r_pop <- cbind(ab_pop,abold_r_pop) # need to do some cleaning in the future
write.csv(abold_r_pop, paste0("GIS_Model/GIS_output/abold_r",date,".csv"))
sum(abold_r_pop$POP)
diff <- sum(abold_r_pop$POP)-sum(ab_pop_zs$POP)
diff

#### TABULATION WORK AND PREPARING LAYERS #########
# Putting FINAL RESULTS on tabular format and the corresponding GIS layers with nice names and so on.
result <- ab_pop_zs

# Rename ADM3 (it includes urban rural strata) and create ADM1 and ADM1 Codes in table
result$UR <- substr(result$ADM3_PCODE_UR,9,9) 
result$ADM3_PCODE <- substr(result$ADM3_PCODE,1,7)
result$ADM2_PCODE <- substr(result$ADM3_PCODE,1,4)
result$ADM1_PCODE <- substr(result$ADM3_PCODE,1,2)
head(result)

# And relocate them at the beginning
result <- result %>%
  relocate(ADM1_PCODE, ADM2_PCODE, .before = ADM3_PCODE) %>% 
  relocate(ADM3_PCODE_UR, UR, .after = ADM3_PCODE)
head(result)

# Round population results
# result <- result %>%  mutate_if(is.numeric,round)

sum(result$POP) # because of the rounding we get a +11 people population difference

## Solve rounding issues in tables

##ROUND PRESERVE SUM projecting population over the point layer as source
#define formula for the round preserving sum
round_preserve_sum <- function(x, digits = 0) {
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}  

result_r <- as.data.frame(result)
result_r$POP_rps <-round_preserve_sum(result_r$POP)
dif <- sum(result_r$POP)-sum(result_r$POP_rps) 
dif # it works

result_r <- result_r %>% mutate_if(is.numeric, round_preserve_sum) 
head(result_r)
sum(result_r$POP)

result
# From here we use the result table to generate the tables for every geo level and Urban Rural Strata.
# ADM3 URBAN RURAL

# Export table connected to the IHSI layer to display results
result <- merge(ab_ihsi_zs,result, by.x='ADM3_PCODE_UR', by.y='ADM3_PCODE_UR')
sum(result$POP)
wgs84 <-"EPSG:4326"
result <- project(result,wgs84)

# clean and order fields
result <- result %>% tidyterra::rename(ADM2_PCODE_UR = com_ur,
                                       ADM1_PCODE_UR = dep_ur,
                                       Commune_adm2_UR = Com_lab_ur,
                                       Departement_adm1_UR = Dep_lab_ur,
                                       Section_adm3_UR = Sec_lab_ur,
                                       ADM3_PCODE = ADM3_PCODE.x)
result <- result[,-c(5:7,9,10,16,74,75)]
head(result)

result <- result %>%
  relocate(c(ADM2_PCODE_UR,Commune_adm2_UR), .after=ADM3_PCODE_UR) %>% 
  relocate(Commune_adm2_UR, .after = ADM2_PCODE_UR ) %>% 
  relocate(Section_adm3_UR, .after = ADM3_PCODE_UR)  %>% 
  relocate(c(ADM1_PCODE_UR,Departement_adm1_UR), .after=Commune_adm2_UR) %>% 
  relocate(ADM1_PCODE, .before = Departement) %>% 
  relocate(ADM2_PCODE, .before = Commune) %>% 
  relocate(ADM3_PCODE, .before = Section_co)

head(result)
# Vector layer ADM3 URBAN RURAL
writeVector(result,"GIS_Model/GIS_output/Population_2022_ADM3_UR.gpkg",overwrite=T)

# Excel file ADM3 URBAN RURAL
openxlsx::write.xlsx(as.data.frame(result), file = "GIS_Model/GIS_output/HTI_Population_2022.xlsx", sheetName='adm3_UR')

# ADM3

# Aggregate results by ADM3_PCODE
result_df <- as.data.frame(result)
adm3_df <- result_df %>% 
  group_by(ADM3_PCODE) %>% 
  summarize(across(c(14:66),sum)) %>% 
  mutate_if(is.numeric, round) %>% 
  as.data.frame()
sum(adm3_df$POP)

# Export table connected to the IHSI layer to display results
result <- merge(ab_ihsi,result, by.x='ADM3_PCODE', by.y='ADM3_PCODE')
sum(result$POP)
wgs84 <-"EPSG:4326"
result <- project(result,wgs84)

# clean fields
result <- result[,-c(6:9)]
head(result)

# Vector layer ADM3
writeVector(result,"GIS_Model/GIS_output/Population_2022_ADM3.gpkg",overwrite=T)

# Excel file ADM3
openxlsx::write.xlsx(as.data.frame(result), file = "GIS_Model/GIS_output/Population_2022_ADM3.xlsx", sheetName='adm3')
adm3 <- result %>% as.data.frame() %>% mutate_if(is.numeric,~round::roundX(.,digits=0,"r1.C")) 
sum(result$POP)

# Next model outcome is  the tables at Section, Commune and Province level using the most recent admin codes and boundaries.
# Import Adm2 and Adm1 layers
adm2 <- vect("GIS_Model/GIS_Input/commune_hti_ihsi_22_26718.gpkg")
adm1 <- vect("GIS_Model/GIS_Input/departement_hti_ihsi_22_26718.gpkg")

# rename admin level codes field to merge table with layerlater
adm2$ADM2_PCODE <- adm2$Code_2_ihs
adm2 <- adm2[,-c(4:7)]
adm1$ADM1_PCODE <- adm1$Code_2_ihs
adm1 <- adm1[,-c(3:5)]

# Aggregate results table at adm 2 and 1
result_df <- as.data.frame(result)
adm2_df <- result_df

adm2_df <- adm2_df %>% 
  group_by(ADM2_PCODE) %>% 
  summarize(across(c(8:60),sum)) %>% 
  mutate_if(is.numeric, round) %>% 
  as.data.frame()
sum(adm2_df$POP)

adm1_df <- result_df

adm1_df <- adm1_df %>% 
  group_by(ADM1_PCODE) %>% 
  summarize(across(c(8:60),sum)) %>%
  mutate_if(is.numeric, round) %>% 
  as.data.frame()

sum(adm1_df$POP)

# Merge with layers
result_adm1 <- merge(adm1,adm1_df, by.x='ADM1_PCODE', by.y='ADM1_PCODE')
head(result_adm1)

result_adm2 <- merge(adm2,adm2_df, by.x='ADM2_PCODE', by.y='ADM2_PCODE')
head(result_adm2)

# Export to layers
writeVector(result_adm2,"GIS_Model/GIS_output/Population_2022_ADM2.gpkg",overwrite=T)
writeVector(result_adm1,"GIS_Model/GIS_output/Population_2022_ADM1.gpkg",overwrite=T)
# Save into results excel sheet
openxlsx::write.xlsx(result_adm2, file = "GIS_Model/GIS_output/Population_2022_ADM2.xlsx", sheetName='adm2')
openxlsx::write.xlsx(result_adm1, file = "GIS_Model/GIS_output/Population_2022_ADM1.xlsx", sheetName='adm1')
