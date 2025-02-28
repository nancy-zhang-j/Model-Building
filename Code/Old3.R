
# 1) Load libraries
########################################
# If not yet installed, uncomment:
# install.packages("VAST")
# install.packages("FishStatsUtils")
# install.packages("sf")
# install.packages("rnaturalearth")
# install.packages("rnaturalearthdata")
# install.packages("tidyr")
# install.packages("dplyr")

library(VAST)
library(dplyr)
library(FishStatsUtils)
library(tidyr)

# <<< ADDED >>>
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

########################################
# 2) Read in data
########################################

fileBIO <- read.csv('/Users/nzhang/LabFiles/FALL BOTTOM TRAWL/USEFUL FILES/22560_UNION_FSCS_SVBIO.csv')
fileCAT <- read.csv('/Users/nzhang/LabFiles/FALL BOTTOM TRAWL/USEFUL FILES/22560_UNION_FSCS_SVCAT.csv')
fileSTA <- read.csv('/Users/nzhang/LabFiles/FALL BOTTOM TRAWL/USEFUL FILES/22560_UNION_FSCS_SVSTA.csv')

# setting directory, can change based on what's needed
setwd("/Users/nzhang/LabFiles/LABIMAGES_2")

# Combine data to reduce size
fileBIO_summary <- fileBIO %>%
  group_by(ID) %>%
  summarize(
    Location  = first(STRATUM)
  )

fileCAT_filtered <- fileCAT %>%
  filter(SVSPP == 103)

fileCAT_summary <- fileCAT_filtered %>%
  group_by(ID) %>%
  summarize(
    CatchWeight = sum(EXPCATCHWT, na.rm = TRUE),
    CatchType   = first(SVSPP)
  )

fileSTA_summary <- fileSTA %>%
  group_by(ID) %>%
  summarize(
    Year       = first(EST_YEAR),
    AvgDepth   = mean(AVGDEPTH, na.rm = TRUE),
    BegLat     = first(DECDEG_BEGLAT),
    BegLon     = first(DECDEG_BEGLON),
    BottomTemp = mean(BOTTEMP, na.rm = TRUE)
  )

combined_data <- fileBIO_summary %>%
  full_join(fileCAT_summary, by = "ID") %>%
  full_join(fileSTA_summary, by = "ID")

# Set a default area-swept = 1
combined_data$AreaSwept <- 1

# Remove missing or invalid observations
combined_data <- combined_data %>%
  filter(!is.na(Year)) %>%
  filter(!is.na(BegLat), !is.na(BegLon)) %>%
  filter(!is.na(CatchType), !is.na(CatchWeight))

########################################
# 3) Add an artificial row for 2020
########################################

artificial_row <- combined_data %>%
  summarise(
    ID          = "Artificial2020",
    Location    = "Artificial",
    CatchWeight = mean(CatchWeight, na.rm = TRUE),
    CatchType   = 103,
    Year        = 2020,
    AvgDepth    = mean(AvgDepth, na.rm = TRUE),
    BegLat      = mean(BegLat, na.rm = TRUE),
    BegLon      = mean(BegLon, na.rm = TRUE),
    BottomTemp  = mean(BottomTemp, na.rm = TRUE),
    AreaSwept   = 1
  )

combined_data <- bind_rows(combined_data, artificial_row)

########################################
# 4) Create user grid (candidate lat/lon)
########################################

min_lat <- min(combined_data$BegLat, na.rm = TRUE)
max_lat <- max(combined_data$BegLat, na.rm = TRUE)
min_lon <- min(combined_data$BegLon, na.rm = TRUE)
max_lon <- max(combined_data$BegLon, na.rm = TRUE)

lat_vec <- seq(from = min_lat, to = max_lat, by = 0.05)
lon_vec <- seq(from = min_lon, to = max_lon, by = 0.05)

grid_df <- expand.grid(Lon = lon_vec, Lat = lat_vec)

grid_df$Area_km2 <- 1  
grid_df$STRATA   <- "NJ"

########################################
# 5) Remove land points using sf + rnaturalearth
########################################

# <<< ADDED >>>
# Convert grid to sf (WGS84)
grid_sf <- st_as_sf(grid_df, coords = c("Lon","Lat"), crs = 4326)

# Download land polygons
land <- rnaturalearth::ne_download(
  scale = "medium",
  type = "land",
  category = "physical",
  returnclass = "sf"
)

# Create bounding box from our grid extent
bbox <- st_bbox(c(
  xmin = min_lon, xmax = max_lon,
  ymin = min_lat, ymax = max_lat
), crs = st_crs(4326))

# Crop land to this bounding box to save time
land_cropped <- st_crop(land, bbox)

# Union into a single polygon
land_union <- st_union(land_cropped)

# Remove any grid points intersecting the land polygon
grid_sea_sf <- st_difference(grid_sf, land_union)

# Convert back to data frame (Lon/Lat columns)
grid_sea_df <- grid_sea_sf %>%
  mutate(
    Lon = st_coordinates(.)[,1],
    Lat = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry()

# We'll use grid_sea_df for VAST
input_grid_df <- grid_sea_df

########################################
# 6) Build extrapolation list
########################################

Extrapolation_List <- make_extrapolation_info(
  Region     = "User",
  input_grid = input_grid_df
)

########################################
# 7) Combine with year and fill covariates
########################################

years <- sort(unique(combined_data$Year))

# Expand to all years in your region
full_grid <- input_grid_df %>%
  crossing(Year = years)

# Join with depth/temp (if time-invariant, it will be NA for most combos)
covariate_data_full <- full_grid %>%
  left_join(combined_data, 
            by = c("Lat" = "BegLat", "Lon" = "BegLon", "Year"))

# Fill NAs with overall means if you prefer
covariate_data_full <- covariate_data_full %>%
  mutate(
    Depth      = ifelse(is.na(AvgDepth),
                        mean(combined_data$AvgDepth, na.rm = TRUE),
                        AvgDepth),
    BottomTemp = ifelse(is.na(BottomTemp),
                        mean(combined_data$BottomTemp, na.rm = TRUE),
                        BottomTemp)
  )

covariate_data_final <- covariate_data_full %>%
  select(Lat, Lon, Year, Depth, BottomTemp)

########################################
# 8) VAST settings
########################################

settings <- make_settings(
  n_x         = 100, 
  Region      = "User", 
  purpose     = "index2",
  bias.correct= FALSE
)

settings$RhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=1, "Epsilon2"=1)
settings$FieldConfig = c("Omega1" = "IID",
                         "Epsilon1" = "IID",
                         "Omega2" = "IID",
                         "Epsilon2" = "IID")

########################################
# 9) Fit the model
########################################

X1_formula = ~ Depth + BottomTemp
X2_formula = ~ Depth + BottomTemp

fit <- fit_model(
  settings        = settings,
  Lat_i           = combined_data$BegLat,
  Lon_i           = combined_data$BegLon,
  t_i             = combined_data$Year,
  b_i             = combined_data$CatchWeight,
  a_i             = combined_data$AreaSwept,
  covariate_data  = covariate_data_final,
  X1_formula      = X1_formula,
  X2_formula      = X2_formula,
  extrapolation_list = Extrapolation_List,
  working_dir     = '/Users/nzhang/LabFiles/LAB_IMAGES'
)

########################################
# 10) Inspect results
########################################

summary(fit)
print(fit$parameter_estimates)
plot(fit)
