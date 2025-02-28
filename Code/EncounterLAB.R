#--------------------------------------
#Load libraries/packages
library(gtools)
library(VAST)
library(FishStatsUtils)
library(dplyr)
library(tidyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata) 

#--------------------------------------
#Read files, these are taken from the FALL 2023 BOTTOM TRAWL surveys from NOAA
fileBIO <- read.csv('/Users/nzhang/LabFiles/FALL BOTTOM TRAWL/USEFUL FILES/22560_UNION_FSCS_SVBIO.csv')
fileCAT <- read.csv('/Users/nzhang/LabFiles/FALL BOTTOM TRAWL/USEFUL FILES/22560_UNION_FSCS_SVCAT.csv')
fileSTA <- read.csv('/Users/nzhang/LabFiles/FALL BOTTOM TRAWL/USEFUL FILES/22560_UNION_FSCS_SVSTA.csv')
#set working directory, where you want the info/graphs to be sent
setwd("/Users/nzhang/LabFiles/LABIMAGES_2")
#--------------------------------------
#Combine data, filtering out what is unecessary and filtering by specific id's

#Summarize fileBIO, where we group the data by its unique ID, and filter out stratum/location
fileBIO_summary <- fileBIO %>%
  group_by(ID) %>%
  summarize(
    Location  = first(STRATUM)
  )

#Filter fileCAT to keep ONLY golden flounder (SVSPP == 103)
fileCAT_filtered <- fileCAT %>%
  filter(SVSPP == 103)

#Summarize fileCAT, where we keep catchweight and catch type
fileCAT_summary <- fileCAT_filtered %>%
  group_by(ID) %>%
  summarize(
    CatchWeight = sum(EXPCATCHWT, na.rm = TRUE),
    CatchType   = first(SVSPP)
  )

#Summarize fileSTA, where we keep year, average depth, beginning latitude, beginning longitude, and bottom temperature
fileSTA_summary <- fileSTA %>%
  group_by(ID) %>%
  summarize(
    Year       = first(EST_YEAR),
    AvgDepth   = mean(AVGDEPTH, na.rm = TRUE),
    BegLat     = first(DECDEG_BEGLAT),
    BegLon     = first(DECDEG_BEGLON),
    BottomTemp = mean(BOTTEMP, na.rm = TRUE)
  )

#Combine data by 'ID', just to get rid of duplicates
combined_data <- fileBIO_summary %>%
  full_join(fileCAT_summary, by = "ID") %>%
  full_join(fileSTA_summary, by = "ID")

#Set AreaSwept = 1, as standard unit. This can be changed later.
combined_data$AreaSwept <- 1

#Remove incomplete rows, missing data, etc.
combined_data <- combined_data %>%
  filter(!is.na(Year),
         !is.na(BegLat),
         !is.na(BegLon),
         !is.na(CatchType),
         !is.na(CatchWeight))

#--------------------------------------
#Add an artificial row for 2020, VAST throws an error if not, ignore row for calculations. 
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
#--------------------------------------
#Make user defined grid, we're only going to use the latitude/longitude for the swept area, which is near NJ. 
min_lat <- min(combined_data$BegLat, na.rm = TRUE)
max_lat <- max(combined_data$BegLat, na.rm = TRUE)
min_lon <- min(combined_data$BegLon, na.rm = TRUE)
max_lon <- max(combined_data$BegLon, na.rm = TRUE)

#E.g., 0.05 degree step
lat_vec <- seq(from = min_lat, to = max_lat, by = 0.05)
lon_vec <- seq(from = min_lon, to = max_lon, by = 0.05)

#Expand into every combination
grid_df <- expand.grid(Lon = lon_vec, Lat = lat_vec)

#Set area and strata
grid_df$Area_km2 <- 1
grid_df$STRATA <- "NJ"

#--------------------------------------
#Remove landmasses from the grid
#Convert grid to sf
grid_sf <- st_as_sf(grid_df, coords = c("Lon", "Lat"), crs = 4326)

#Download global land polygons (medium resolution)
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

#Union into one big polygon, we'll subtract this from our grid
land_polygon <- st_union(world)

#Identify points on land
on_land <- st_intersects(grid_sf, land_polygon, sparse = FALSE)[,1]

#Keep only water points
grid_sf_water <- grid_sf[!on_land, ]

#Convert back to data frame
grid_df_water <- grid_sf_water %>%
  mutate(
    Lon = st_coordinates(.)[,1],
    Lat = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry()

#--------------------------------------
#Build extrapolation grid with our new water grid
input_grid_df <- grid_df_water

Extrapolation_List <- make_extrapolation_info(
  Region     = "User",
  input_grid = input_grid_df
)
#--------------------------------------
#Set up covariate data for VAST

#Find the unique years
years <- sort(unique(combined_data$Year))

#Deduplicate and sort the inputs
full_grid <- input_grid_df %>%
  crossing(Year = years)

#Combine the full grid with the data you want to use for covariates
covariate_data_full <- full_grid %>%
  left_join(
    combined_data,
    by = c("Lat" = "BegLat", "Lon" = "BegLon", "Year")
  ) %>%
  mutate( #We're checking for NA values. If there is, we replace them with an average value
    Depth = ifelse(
      is.na(AvgDepth),
      mean(combined_data$AvgDepth, na.rm = TRUE),
      AvgDepth
    ),
    BottomTemp = ifelse( #Again, checking and replacing for NA values.
      is.na(BottomTemp),
      mean(combined_data$BottomTemp, na.rm = TRUE),
      BottomTemp
    )
  )

#Make a new clean grid for covariate data once all the information has been added and verified
covariate_data_final <- covariate_data_full %>%
  select(Lat, Lon, Year, Depth, BottomTemp)
#--------------------------------------
#VAST settings
settings <- make_settings(
  n_x         = 100, #Number of knots
  Region      = "User", #Our own user grid
  purpose     = "index2", #Index of abundance calculated summing across space to get annual abundance values for each category,
  #uses gamma distribution for positive catches, restricts max_cells to 2000, 
  #and uses bias correction. This is currently recommended over "index".
  bias.correct= FALSE #Saves time
)

#Turn on random effects for encounter and positive side
settings$RhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=1, "Epsilon2"=1)
settings$FieldConfig = c("Omega1" = "IID", "Epsilon1" = "IID",
                         "Omega2" = "IID", "Epsilon2" = "IID")

#--------------------------------------
#Fit the model
X1_formula = ~ Depth + BottomTemp #Encounter probability
X2_formula = ~ Depth + BottomTemp #Positive catch/Abundance

fit <- fit_model(
  settings            = settings,
  Lat_i               = combined_data$BegLat,
  Lon_i               = combined_data$BegLon,
  t_i                 = combined_data$Year,
  b_i                 = combined_data$CatchWeight,
  a_i                 = combined_data$AreaSwept,
  covariate_data      = covariate_data_final,
  X1_formula          = X1_formula,
  X2_formula          = X2_formula,
  extrapolation_list  = Extrapolation_List,
  working_dir         = '/Users/nzhang/LabFiles/LAB_IMAGES'
)
#--------------------------------------
#results
plot(fit)
summary(fit)
print(fit$parameter_estimates)

#Optional: examine coefficients, determine others later
#fit$Report$Beta1
#fit$Report$Beta2
#----------------------------------------------------------------
#Make map_info (map_list)

map_list <- make_map_info(
  Region            = settings$Region,         #e.g., "User"
  spatial_list      = fit$spatial_list,
  Extrapolation_List= fit$extrapolation_list,
  fine_scale        = TRUE
)
#----------------------------------------------------------------
#Extract or customize PlotDF from map_list
df_transformed <- map_list$PlotDF
#Optionally, add columns or filters here:
df_transformed$Include <- TRUE  
#----------------------------------------------------------------
#Call plot_maps() for spatio-temporal variation, 
#(e.g., log-positive catch rates when using a conventional delta-model),
plot_maps(
  plot_set    = 6,                           #Epsilon1
  fit         = fit,
  PlotDF      = df_transformed,              #Must match domain used in fit
  Panel       = "Year",                      #Panels by year
  n_cells     = 1000,                        #Adjust for resolution
  working_dir = "/Users/nzhang/LabFiles/Epsilon_1"      #Where output is saved
)

plot_maps(
  plot_set    = 7,                           #Epsilon2
  fit         = fit,
  PlotDF      = df_transformed,              
  Panel       = "Year",                     
  n_cells     = 1000,                        
  working_dir = "/Users/nzhang/LabFiles/Epsilon_2"      
)

plot_maps(
  plot_set    = 16,                           #Omega_1
  fit         = fit,
  PlotDF      = df_transformed,              
  Panel       = "Year",                     
  n_cells     = 1000,                        
  working_dir = "/Users/nzhang/LabFiles/Omega_1"      
)

plot_maps(
  plot_set    = 17,                         #Omega2
  fit         = fit,
  PlotDF      = df_transformed,              
  Panel       = "Year",                      
  n_cells     = 1000,                       
  working_dir = "/Users/nzhang/LabFiles/Omega_2"     
)

#Spatial random effects for encounter (first linear predictor)
bomega1 <- fit$Report$Omega1_gc
#Spatio-temporal random effects for encounter
bepsilon1 <- fit$Report$Epsilon1_gct
#Spatial random effects for positive catch (second linear predictor)
bomega2 <- fit$Report$Omega2_gc
#Spatio-temporal random effects for positive catch
bepsilon2 <- fit$Report$Epsilon2_gct

#----------------------------------------------------------------
#Attempting predictions
#Specify knot and year, needs to be better
knot <- 19 #there are 100 total knots for this code
year_index <- 19 #1 = 1963, etc.
#----------------------------------------------------------------

#Extract latitude/longitude for the chosen knot. There are 100 knots, so we have 100 rows. 
lat_knot <- as.numeric(fit$spatial_list$latlon_x[knot, 1])
lon_knot <- as.numeric(fit$spatial_list$latlon_x[knot, 2])

#Get the actual year used
year_used <- sort(unique(covariate_data_final$Year))[year_index]

#Find depth and bottom temperature matching this lat/lon and year
row_knot_year <- covariate_data_final %>%
  filter(
    Year == year_used, #Filtering by year
    abs(Lat - lat_knot) < 0.025, #Filtering out, difference is less than 0.025
    abs(Lon - lon_knot) < 0.025  
  )

Depth_i <- row_knot_year$Depth #Knot i's depth at year r
BottomTemp_i <- row_knot_year$BottomTemp #Knot i's depth at year r
#----------------------------------------------------------------
#Section for adding/subtracing values from covariates
Depth_i2 <- Depth_i - 10
BottomTemp_i2 <- BottomTemp_i + 2
#----------------------------------------------------------------
#Fixed effects for encounter 
beta1_ft     <- 3.11871697 #density
gamma1_depth <- -0.02083465 #beta for depth
gamma1_temp  <- -0.04324361 #beta for temp

#Fixed effects for abundance
beta2_ft     <- -2.37716314 #density
gamma2_depth <- 0.01262833  #beta for depth
gamma2_temp  <- 0.02814531  #beta for temp
#----------------------------------------------------------------
#Random effects for this knot and year (check again)
Omega1_knot        <- fit$Report$Omega1_gc[knot, 1] #g = lon/latitude, c = categories/fish (only flounder), t = time/year
Epsilon1_knot_year <- fit$Report$Epsilon1_gct[knot, 1, year_index] 
Omega2_knot        <- fit$Report$Omega2_gc[knot, 1]
Epsilon2_knot_year <- fit$Report$Epsilon2_gct[knot, 1, year_index]
#----------------------------------------------------------------
#Encounter probability
lp_encounter <- beta1_ft +
  gamma1_depth * Depth_i +
  gamma1_temp  * BottomTemp_i +
  Omega1_knot +
  Epsilon1_knot_year

p_encounter <-  1 / (1 + exp(-lp_encounter))
p_encounter
#----------------------------------------------------------------
#Positive catch (log link)
lp_catch <- beta2_ft +
  gamma2_depth * Depth_i +
  gamma2_temp  * BottomTemp_i +
  Omega2_knot +
  Epsilon2_knot_year

pos_catch <- exp(lp_catch)
pos_catch
#Expected catch (delta-model)
expected_catch <- p_encounter * pos_catch
expected_catch

#----------------------------------------------------------------
#Encounter probability predicted
lp_encounter2 <- beta1_ft +
  gamma1_depth * Depth_i2 +
  gamma1_temp  * BottomTemp_i2 +
  Omega1_knot +
  Epsilon1_knot_year

p_encounter2 <-  1 / (1 + exp(-lp_encounter2))
p_encounter2
#----------------------------------------------------------------
#Positive catch (log link) predicted
lp_catch2 <- beta2_ft +
  gamma2_depth * Depth_i2 +
  gamma2_temp  * BottomTemp_i2 +
  Omega2_knot +
  Epsilon2_knot_year

pos_catch2 <- exp(lp_catch2)
pos_catch2
#Expected catch (delta-model)
expected_catch2 <- p_encounter2 * pos_catch2
expected_catch2

#Higher encounter rate with shallow + colder waters 
#Higher abundance index with deeper + warmer waters

#New stuff, February 20
#-------------------------------
#Define fixed effects
beta1_ft     <- 3.11871697   #encounter intercept
gamma1_depth <- -0.02083465  #encounter beta for depth
gamma1_temp  <- -0.04324361  #encounter beta for temperature

beta2_ft     <- -2.37716314  #positive catch intercept
gamma2_depth <- 0.01262833   #positive catch beta for depth
gamma2_temp  <- 0.02814531   #positive catch beta for temperature

#-------------------------------
#Define the adjustment scenarios
#Here, depth_changes are in meters: 0 (baseline), -0.1 (shallower by 10 cm), -1 (shallower by 100 cm)
#temp_changes are in degrees Celsius: 0 (baseline), +1, +2
depth_changes <- c(0, -0.1, -1)   
temp_changes  <- c(0, 1, 2)

#Get the number of knots from  fit; assuming fit$spatial_list$latlon_x is a matrix with rows = knots
n_knots <- nrow(fit$spatial_list$latlon_x)

#Get all unique years from covariate data
years_all <- sort(unique(covariate_data_final$Year))

#Initialize an empty data frame to store results
results <- data.frame()

#Loop over every knot and every year
for (knot in 1:n_knots) {
  #Extract knot coordinates
  lat_knot <- as.numeric(fit$spatial_list$latlon_x[knot, 1])
  lon_knot <- as.numeric(fit$spatial_list$latlon_x[knot, 2])
  
  for (yr in years_all) {
    #Find the baseline covariate values for this knot and year.
    #Using a spatial filter (difference < 0.025) similar to  example.
    row_knot_year <- covariate_data_final %>%
      filter(
        Year == yr,
        abs(Lat - lat_knot) < 0.025,
        abs(Lon - lon_knot) < 0.025
      )
    
    #If no matching record is found, skip this knot-year combination
    if(nrow(row_knot_year) == 0) next
    
    #Get baseline depth and temperature (average in case there are multiple matches)
    base_depth <- mean(row_knot_year$Depth, na.rm = TRUE)
    base_temp  <- mean(row_knot_year$BottomTemp, na.rm = TRUE)
    
    #Extract the corresponding random effects for this knot and year.
    #Note: Adjust the indexing for the time dimension as needed.
    #Here we assume the order in 'years_all' matches the time dimension of Epsilon arrays.
    time_index <- which(years_all == yr)
    Omega1_knot        <- fit$Report$Omega1_gc[knot, 1]
    Epsilon1_knot_year <- fit$Report$Epsilon1_gct[knot, 1, time_index]
    Omega2_knot        <- fit$Report$Omega2_gc[knot, 1]
    Epsilon2_knot_year <- fit$Report$Epsilon2_gct[knot, 1, time_index]
    
    #Loop over each combination of depth and temperature change scenarios
    for (d_change in depth_changes) {
      for (t_change in temp_changes) {
        
        #New covariate values after adjustment
        new_depth <- base_depth + d_change
        new_temp  <- base_temp  + t_change
        
        #Calculate encounter probability (using logistic link)
        lp_encounter <- beta1_ft + 
          gamma1_depth * new_depth + 
          gamma1_temp  * new_temp +
          Omega1_knot +
          Epsilon1_knot_year
        
        p_encounter <- 1 / (1 + exp(-lp_encounter))
        
        #Calculate positive catch (using log link)
        lp_catch <- beta2_ft + 
          gamma2_depth * new_depth + 
          gamma2_temp  * new_temp +
          Omega2_knot +
          Epsilon2_knot_year
        
        pos_catch <- exp(lp_catch)
        
        #Delta-model expected catch
        expected_catch <- p_encounter * pos_catch
        
        #Record the results
        results <- rbind(results, data.frame(
          Knot          = knot,
          Year          = yr,
          Base_Depth    = base_depth,
          Base_Temp     = base_temp,
          Depth_Change  = d_change,
          Temp_Change   = t_change,
          New_Depth     = new_depth,
          New_Temp      = new_temp,
          p_encounter   = p_encounter,
          pos_catch     = pos_catch,
          expected_catch = expected_catch
        ))
      }
    }
  }
}

#Examine the resulting predictions
head(results)

#Using plot_maps with plot_set = 3 (which typically corresponds to the abundance index)
plot_maps(
  plot_set    = 3,                   #2 usually corresponds to the abundance index in index2 settings
  fit         = fit,
  PlotDF      = df_transformed,        #Use the same PlotDF from  map info
  Panel       = "Year",
  n_cells     = 1000,
  working_dir = "/Users/nzhang/LabFiles"
)

library(ggplot2)
library(dplyr)

#---------------------------------------------------------
#Summarize abundance predictions by Year and Depth_Change 
#   (with temperature fixed at baseline)
#---------------------------------------------------------
abundance_depth <- results %>% 
  filter(Temp_Change == 0) %>% 
  group_by(Year, Depth_Change) %>%
  summarize(mean_expected = mean(expected_catch, na.rm = TRUE),
            se_expected   = sd(expected_catch, na.rm = TRUE)/sqrt(n()),
            .groups = 'drop')

#Plot: Abundance index (expected catch) versus Depth Change
ggplot(abundance_depth, aes(x = Depth_Change, y = mean_expected,
                            color = factor(Year), group = Year)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_expected - se_expected,
                    ymax = mean_expected + se_expected),
                width = 0.02) +
  labs(
    title = "Abundance Index vs. Depth Change (Temp Change = 0°C)",
    x = "Depth Change (m)",  #e.g., 0, -0.1 (10 cm shallower), -1 (100 cm shallower)
    y = "Mean Expected Abundance Index",
    color = "Year"
  ) +
  theme_minimal()

#---------------------------------------------------------
#Summarize abundance predictions by Year and Temp_Change 
#(with depth fixed at baseline)
#---------------------------------------------------------
abundance_temp <- results %>% 
  filter(Depth_Change == 0) %>% 
  group_by(Year, Temp_Change) %>%
  summarize(mean_expected = mean(expected_catch, na.rm = TRUE),
            se_expected   = sd(expected_catch, na.rm = TRUE)/sqrt(n()),
            .groups = 'drop')

#Plot: Abundance index (expected catch) versus Temperature Change
ggplot(abundance_temp, aes(x = Temp_Change, y = mean_expected,
                           color = factor(Year), group = Year)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_expected - se_expected,
                    ymax = mean_expected + se_expected),
                width = 0.1) +
  labs(
    title = "Abundance Index vs. Temperature Change (Depth Change = 0 m)",
    x = "Temperature Change (°C)",  #e.g., 0, 1, 2 °C increase
    y = "Mean Expected Abundance Index",
    color = "Year"
  ) +
  theme_minimal()

#---------------------------------------------------------
#Boxplots of abundance predictions for
#all scenario combinations to view the distribution across knots
#---------------------------------------------------------
ggplot(results, aes(x = factor(Year), y = expected_catch)) +
  geom_boxplot() +
  facet_grid(Depth_Change ~ Temp_Change,
             labeller = label_both) +
  labs(
    title = "Distribution of Abundance Index by Scenario",
    x = "Year",
    y = "Expected Abundance Index"
  ) +
  theme_minimal()

#-----------------------------------------------------
#Create a modified covariate dataset
#For example: decrease depth by 0.1 m and increase bottom temp by 1°C
covariate_data_modified <- covariate_data_final %>%
  mutate(
    Depth = Depth - 0.1,
    BottomTemp = BottomTemp + 1
  )

#-----------------------------------------------------
#Re-fit the model using the modified covariate data.
#You can re-use  original settings and data inputs.
fit_modified <- fit_model(
  settings            = settings,
  Lat_i               = combined_data$BegLat,
  Lon_i               = combined_data$BegLon,
  t_i                 = combined_data$Year,
  b_i                 = combined_data$CatchWeight,
  a_i                 = combined_data$AreaSwept,
  covariate_data      = covariate_data_modified,  #use the modified covariates
  X1_formula          = X1_formula,
  X2_formula          = X2_formula,
  extrapolation_list  = Extrapolation_List,
  working_dir         = '/Users/nzhang/LabFiles/LAB_IMAGES_modified'
)

#-----------------------------------------------------
#Create a new map info list from the updated fit.
map_list_modified <- make_map_info(
  Region            = settings$Region,
  spatial_list      = fit_modified$spatial_list,
  Extrapolation_List= fit_modified$extrapolation_list,
  fine_scale        = TRUE
)

#Optionally, adjust the PlotDF if needed
df_transformed_modified <- map_list_modified$PlotDF
df_transformed_modified$Include <- TRUE

#-----------------------------------------------------
#Replot the maps using the modified fit.
#For example, replotting the spatio-temporal variation for Epsilon1:
plot_maps(
  plot_set    = 6,  #Epsilon1
  fit         = fit_modified,
  PlotDF      = df_transformed_modified,
  Panel       = "Year",
  n_cells     = 1000,
  working_dir = "/Users/nzhang/LabFiles/Epsilon_1_modified"
)
