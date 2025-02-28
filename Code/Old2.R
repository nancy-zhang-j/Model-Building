# Load required libraries
library(VAST)
library(dplyr)
library(FishStatsUtils)
library(sf)        # For spatial operations
library(zoo)       # For interpolation
library(tidyr)     # For data manipulation

#----------------------------
# 1. Read and Combine Data
#----------------------------

# Read CSV files
fileBIO <- read.csv('/Users/nzhang/LabFiles/FALL BOTTOM TRAWL/USEFUL FILES/22560_UNION_FSCS_SVBIO.csv')
fileCAT <- read.csv('/Users/nzhang/LabFiles/FALL BOTTOM TRAWL/USEFUL FILES/22560_UNION_FSCS_SVCAT.csv')
fileSTA <- read.csv('/Users/nzhang/LabFiles/FALL BOTTOM TRAWL/USEFUL FILES/22560_UNION_FSCS_SVSTA.csv')

# Summarize fileBIO
fileBIO_summary <- fileBIO %>%
  group_by(ID) %>%
  summarize(
    Location = first(STRATUM)
  )

# Filter and summarize fileCAT for species SVSPP == 103 (Golden flounder)
fileCAT_filtered <- fileCAT %>%
  filter(SVSPP == 103)

fileCAT_summary <- fileCAT_filtered %>%
  group_by(ID) %>%
  summarize(
    CatchNumber = sum(EXPCATCHWT, na.rm = TRUE),
    CatchType   = first(SVSPP)
  )

# Summarize fileSTA
fileSTA_summary <- fileSTA %>%
  group_by(ID) %>%
  summarize(
    Year        = mean(EST_YEAR),
    AvgDepth    = mean(AVGDEPTH, na.rm = TRUE),
    BegLat      = first(DECDEG_BEGLAT),
    BegLon      = first(DECDEG_BEGLON),
    BottomTemp  = mean(BOTTEMP, na.rm = TRUE)
  )

# Combine all data
combined_data <- fileBIO_summary %>%
  full_join(fileCAT_summary, by = "ID") %>%
  full_join(fileSTA_summary, by = "ID")

# Set AreaSwept (assuming a constant if actual data isn't available)
combined_data$AreaSwept <- 1

# Remove rows with missing essential information
combined_data <- combined_data %>%
  filter(
    !is.na(Year),
    !is.na(BegLat), !is.na(BegLon),
    !is.na(CatchType), !is.na(CatchNumber)
  )

#----------------------------
# 2. Create Extrapolation Grid
#----------------------------

# Define grid boundaries
min_lat <- min(combined_data$BegLat, na.rm = TRUE)
max_lat <- max(combined_data$BegLat, na.rm = TRUE)
min_lon <- min(combined_data$BegLon, na.rm = TRUE)
max_lon <- max(combined_data$BegLon, na.rm = TRUE)

# Create grid with 0.1 degree step
lat_vec <- seq(from = min_lat, to = max_lat, by = 0.1)
lon_vec <- seq(from = min_lon, to = max_lon, by = 0.1)

# Expand into every combination
grid_df <- expand.grid(Lon = lon_vec, Lat = lat_vec)

# Assign area and strata
grid_df$Area_km2 <- 1   # Adjust if necessary
grid_df$STRATA <- "NJ"

input_grid_df <- grid_df

# Create Extrapolation_List
Extrapolation_List <- make_extrapolation_info(
  Region     = "User",
  input_grid = input_grid_df
)

#----------------------------
# 3. Assign Covariate Values to Grid Points
#----------------------------

# Convert combined_data and grid_df to spatial objects
tows_sf <- st_as_sf(combined_data, coords = c("BegLon", "BegLat"), crs = 4326)
grid_sf <- st_as_sf(input_grid_df, coords = c("Lon", "Lat"), crs = 4326)

# Define unique years
years <- sort(unique(combined_data$Year))

# Create a complete grid of Lon, Lat, and Year
covariate_data_full_dynamic <- expand.grid(
  Lon = lon_vec,
  Lat = lat_vec,
  Year = years
)

# Convert to spatial objects
covariate_data_full_dynamic_sf <- st_as_sf(covariate_data_full_dynamic, coords = c("Lon", "Lat"), crs = 4326)

# Assign Depth (static) using nearest neighbor
nearest_indices_cov_depth <- st_nearest_feature(covariate_data_full_dynamic_sf, tows_sf)
covariate_data_full_dynamic$Depth <- tows_sf$AvgDepth[nearest_indices_cov_depth]

# Assign BottomTemp (dynamic) using nearest neighbor
nearest_indices_cov_bt <- st_nearest_feature(covariate_data_full_dynamic_sf, tows_sf)
covariate_data_full_dynamic$BottomTemp <- tows_sf$BottomTemp[nearest_indices_cov_bt]

# Handle missing BottomTemp if any (e.g., via interpolation)
if(any(is.na(covariate_data_full_dynamic$BottomTemp))){
  covariate_data_full_dynamic <- covariate_data_full_dynamic %>%
    group_by(Lon, Lat) %>%
    arrange(Year) %>%
    mutate(
      BottomTemp = na.approx(BottomTemp, Year, na.rm = FALSE, rule = 2)
    ) %>%
    ungroup()
}



#----------------------------
# 4. Prepare Covariate Data for VAST
#----------------------------

# Rescale covariates to have mean 0 and SD 1
covariate_data_full_dynamic <- covariate_data_full_dynamic %>%
  mutate(
    Depth = scale(Depth, center = TRUE, scale = TRUE)[,1],
    BottomTemp = scale(BottomTemp, center = TRUE, scale = TRUE)[,1]
  )

# Select relevant columns
covariate_data <- covariate_data_full_dynamic %>%
  select(Lat, Lon, Year, Depth, BottomTemp)

#----------------------------
# 5. Define VAST Settings
#----------------------------

settings <- make_settings(
  n_x            = 100,            # Number of knots
  Region         = "User",         # User-defined region
  purpose        = "index",        # Purpose of the model
  bias.correct   = FALSE,
  use_anisotropy = FALSE
)

# Configure Random Effects
settings$RhoConfig <- c(
  "Beta1"    = 3,  # Random walk in time for encounter intercept
  "Beta2"    = 0,  # No random time variation for positive catch intercept
  "Epsilon1" = 0,  # No temporal correlation for encounter
  "Epsilon2" = 0   # No temporal correlation for positive catch
)

settings$FieldConfig <- c(
  "Omega1"   = "IID",  # Spatial RE for encounter
  "Epsilon1" = "IID",  # Spatio-temporal RE for encounter (IID implies no temporal correlation)
  "Omega2"   = "IID",  # Spatial RE for positive catch
  "Epsilon2" = 0        # No spatio-temporal RE for positive catch
)

#----------------------------
# 6. Define Model Formulas
#----------------------------

# Define Model Formulas
X1_formula = ~ Depth + BottomTemp  # Encounter model
X2_formula = ~ Depth + BottomTemp  # Positive catch model

#----------------------------
# 7. Fit the VAST Model
#----------------------------

fit <- fit_model(
  settings           = settings,
  Lat_i              = combined_data$BegLat,       # Latitude of each tow
  Lon_i              = combined_data$BegLon,       # Longitude of each tow
  t_i                = combined_data$Year,         # Year of each tow
  b_i                = combined_data$CatchNumber,  # Catch quantity
  a_i                = combined_data$AreaSwept,    # Area swept
  X1_formula         = X1_formula,                 # Encounter formula
  X2_formula         = X2_formula,                 # Positive catch formula
  covariate_data     = covariate_data,             # Complete covariate data
  extrapolation_list = Extrapolation_List,         # Extrapolation grid
  working_dir        = '/Users/nzhang/LabFiles/LAB_IMAGES'  # Output directory
)

#----------------------------
# 8. Inspect and Plot Results
#----------------------------

# Summaries
summary(fit)
print(fit$parameter_estimates)

# Diagnostic Plots
plot(fit)
