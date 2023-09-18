### Clean and process HH locations from IHSI database to create an Average HH Size grid for updated population grid for Haiti ###
### Luis de la Rua lurodriguez@unfpa.org
### October 2022

library(sf)
library(sp)
library(tidyverse)
library(readxl)
library(dplyr)
library(ggplot2)
library(plotly)
library(readr)
library(terra)
library(rgeos)
library(spData)
library(raster)


# Work directory
wd <- "C:/GIS/UNFPA GIS/HTI PopGrid/HHloc18_20"
setwd(wd)

# Load all hh locations tables (location of the original data is stored in raw data folder)
files <- list.files("raw_data", pattern = "*.csv", full.names = T)

# Merge all csv files 
merged <- files %>% 
  lapply(read_csv,  col_types = cols( .default = col_character())) %>% # use this better as solves issues with different variable formats 
  bind_rows 

str(merged) # describe result

# Convert UTM coordinates into numeric
i <- c('UtmGpsLatitudeX','UtmGpsLongitudeY') 

merged[ , i] <- apply(merged[ , i], 2,            # Specify own function within apply
                      function(x) as.numeric(as.character(x)))

# write.csv(merged,"cleaned/merged_raw.csv")
sapply(merged, class)

## COORDINATES CLEANING

# 1. Detect broadly and correct, entries with SWAPPED COORDINATES

# create variables where "corrected coordinates are going to be stored
merged$x_clean <- merged$UtmGpsLatitudeX
merged$y_clean <- merged$UtmGpsLongitudeY

# swap coordinates where coordintates seem to be swapped (x with y values)

merged$x_clean <- ifelse((merged$UtmGpsLatitudeX > 2000000), 
                         merged$UtmGpsLongitudeY, merged$x_clean)
merged$y_clean <- ifelse((merged$UtmGpsLatitudeX > 2000000), 
                         merged$UtmGpsLatitudeX, merged$y_clean)

write.csv(merged,"cleaned/merged_swap.csv")

# Delete all points located out of admin boundaries

# Load admin boundaries
ab <- st_read("C:/GIS/UNFPA GIS/HTI PopGrid/HHloc18_20/layers/hti_admbnda_adm3_cnigs_20181129_26718.shp")

utm18 <- st_crs(ab) # we use this later to define coordinate system

# Remove entries with null coordinates
merged$x_clean[is.na(merged$x_clean)] <- 0
merged$y_clean[is.na(merged$y_clean)] <- 0
# Convert merged into spatial data frame object
merged_sf <- st_as_sf(merged, 
                      coords =c("x_clean","y_clean"),
                      crs = utm18)
str(merged_sf)
crs(merged_sf)

# CORRECT PROJECTION ISSUES

# Detect broadly point with WRONG UTM ZONE (collected using UTM 19N)
# Extract points from merged points we suspect have wrong projection

merged_19 <- subset(merged, 
                    (merged$x_clean > 184455) & (merged$x_clean < 224089)) # defining longitude band where we know points with wrong UTM zone are

merged_19_sf <- st_as_sf(merged_19, 
                         coords =c("x_clean","y_clean"),
                         crs = 26719) # using 19N UTM
str(merged_19_sf)
crs(merged_19_sf)

# st_write(merged_19_sf,"layers/hhloc19.gpkg", delete_layer=T)

# Reproject into UTM18
merged_19_sf_rep <- st_transform(merged_19_sf, crs(ab))

# merge join and join_19 into one single layer and clean useless fields
hhloc_merged <- rbind(merged_sf,merged_19_sf_rep)

# REMOVE POINTS OUT ADMIN BOUNDARIES
# Spatial join ab <> merged_sf
join <- st_join(hhloc_merged ,ab,join=st_intersects)

# Remove points out of Admin Boundaries
hhloc_join <- subset(join,join$ADM0_PCODE == 'HT')

# Remove unnecessary fields
hhloc_clean <- hhloc_join[,-c(24:43)]

# Create building and hh size field 
hhloc_clean$bsize <- hhloc_clean$TotalPopulation

# REMOVE DUPLICATES (multiple hh in same building)
# for interpolating I need to remove all duplicates, Total population is recording population in the building so we are 
# having false duplicates that would affect on the interpolation. 
# Remove duplicates using gpsx/gpsy duplicates

hhloc_clean_dup <- hhloc_clean
hhloc_clean_dup <- hhloc_clean_dup %>%
  distinct(UtmGpsLatitudeX, UtmGpsLongitudeY, .keep_all = TRUE)

# Population associated to building is what we are going to use as what we are trying to assess in next stages is the Average population / building assigned to the
# building footprints that is going to be generated

# Export to GPKG
st_write(hhloc_clean_dup,"layers/hhloc_clean_26718.gpkg", delete_layer=T)

