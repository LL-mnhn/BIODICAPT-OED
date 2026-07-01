# This script:
#   - takes pre-processed datasets.
#   - groups them together into one file.
#   - these files will be fed to the models (see next script).

##### Libraries ##### ---------------------------------------------------------
library(dplyr)
library(readr)
library(purrr)
library(cli)
library(sf)

source(here::here("R/utils_figures.R"))


##### Parameters ##### --------------------------------------------------------
source(here::here("data/config/config.R")) # Import global parameters

# biodicapt dataset
BIODICAPT_PATH_PREPROCESSED <- file.path(PROCESSED_DATA_PATH, 
        "BIODICAPT_survey_data_anonymized.csv") 
# 500eni dataset
ENI500_PATH_PREPROCESSED <- file.path(PROCESSED_DATA_PATH, 
        "ENI500_survey_data_anonymized.csv") 
# stoc dataset
STOC_PATH_PREPROCESSED <- file.path(PROCESSED_DATA_PATH, 
        "STOC_data_cleaned.csv") 

# Corine Land Cover shapefile
CORINE_SHP_PATH <- file.path(
    PROCESSED_DATA_PATH,
    paste0(
        "CLC", OBS_YEAR, 
        "_projection_france_hexagons_res", RES_KM, "km-WGS84.gpkg"
    )
)

# Chelsa shapefiles
CHELSA_SHP_PATHS <- file.path(
    PROCESSED_DATA_PATH,
    paste0(
        "CHELSA_", OBS_YEAR, "_", CHELSA_DATASETS, 
        "_projection_france_hexagons_res", RES_KM, "km-WGS84.gpkg"
    )
)


##### Helper functions ##### --------------------------------------------------
load_group_clc_chelsa_sf <- function() {

    if (!authorize_overwrite(MASTER_SF_PATH)) {
        # Check if file exists
        cli_alert_warning("Skipping pre-processing.\n\n")

    } else {
        # Load shapefiles with hexagons
        cli_alert_info("Loading shapefiles of features...")
        cli_alert_info("Loading CLC...")

        # Initialize lists
        list_shapefiles <- list(st_read(CORINE_SHP_PATH, quiet = TRUE))
        list_shp_colors <- list(read_csv(
            gsub(".gpkg", ".csv", CORINE_SHP_PATH), show_col_types = FALSE))
        list_names <- list("CLC")

        # Loop on chelsa datasets
        for (i in 1:length(CHELSA_SHP_PATHS)) {
            cli_alert_info(paste0(
                "Loading chelsa '", CHELSA_DATASETS[i], "' dataset..."))

            list_shapefiles[[length(list_shapefiles) + 1]] <- st_read(
                CHELSA_SHP_PATHS[i], quiet = TRUE)
            
            # if categorical variable, a custom color csv should exist
            if (file.exists(gsub(".gpkg", ".csv", CHELSA_SHP_PATHS[i]))) {
                list_shp_colors[[length(list_shp_colors) + 1]] <-  read_csv(
                    gsub(".gpkg", ".csv", CHELSA_SHP_PATHS[i]), 
                    show_col_types = FALSE
                )
                
            } else {
                list_shp_colors[[length(list_shp_colors) + 1]] <- NA
            }

            list_names[[length(list_names) + 1]] <- paste0(
                "CHELSA_", CHELSA_DATASETS[i])
        }

        # Quick check : compare the geometry of the cells in each shapefile (detects mismatchs)
        same_geometries <- all(
            st_geometry(list_shapefiles[[1]]) 
            == 
            st_geometry(list_shapefiles[[2]]))
        if (same_geometries) {
            cli_alert_success(paste0("Quick check: coordinates ",
            "of hexagons in shapefiles are identical.\n\n"))
        } else {
            stop(paste0("Quick check: coordinates of hexagons in shapefiles",
            " are NOT identical. Cannot proceed further."))
        }

        # Group shapefiles together (into one master_sf), avoid redundance
        cli_alert_info(paste0("Grouping shapefiles together..."))

        master_sf <- list_shapefiles[[1]] |>
        rename(!!list_names[[1]] := dominant_class)
        for (i in seq_along(list_shapefiles)[-1]) {
            df_i <- list_shapefiles[[i]] |>
                st_drop_geometry() |>          # drop geom to avoid duplication
                select(seqnum, dominant_class) |>
                rename(!!list_names[[i]] := dominant_class)

            master_sf <- master_sf |>
                left_join(df_i, by = "seqnum")
        }

        # Save result
        cli_alert_info(paste0("Saving master shapefile..."))

        st_write(master_sf, MASTER_SF_PATH, delete_dsn = TRUE, quiet = TRUE)
        cli_alert_success(paste0("Saved master shapefile!\n\n"))
    }
}

save_features_from_obs <- function(filepath, save_to) {
    if (!authorize_overwrite(save_to)) {
        # Check if file exists
        cli_alert_warning(paste0("Skipping the addition of features to ", 
        basename(filepath), ".\n\n"))

    } else {
    
        # load presence-absence observations
        cli_alert_info(paste0("Loading ", basename(filepath), "..."))
        
        df <- read_csv(filepath, show_col_types = FALSE)
        points <- st_as_sf(
            df, coords = c("LON", "LAT"), crs = 4326, remove = FALSE)
        cli_alert_info("Loading master shapefile...")


        # Load features
        master_sf <- st_read(MASTER_SF_PATH, quiet = TRUE)
        cli_alert_info("Adding features to dataset...")

        # Append table with master_sf features
        master_df <- st_join(points, master_sf, join = st_within)

        # Quick check : compare the geometry of the cells in each shapefile (detects mismatchs)
        all_points_in_hexagons <- all(!is.na(master_df$CLC))
        if (all_points_in_hexagons) {
            cli_alert_success(paste0("Sanity check: no observations ",
            "outside of shapefile hexagons.\n\n"))
            cli_alert_info("Saving dataset with features...")

        } else {
            stop(paste0("Sanity check: some observations are outside of the",
            " study area (master shapefile)."))
        }

        st_write(master_df, save_to, quiet = TRUE, delete_dsn=TRUE)
        cli_alert_success("Saved dataset with features!")
    }
}


##### Formatting shapefile data ##### -----------------------------------------
. <- load_group_clc_chelsa_sf()


##### Formatting observation data ##### ---------------------------------------
. <- save_features_from_obs(
    filepath=BIODICAPT_PATH_PREPROCESSED, save_to=BIODICAPT_OBS_FULL)
. <- save_features_from_obs(
    filepath=ENI500_PATH_PREPROCESSED, save_to=ENI500_OBS_FULL)
. <- save_features_from_obs(
    filepath=STOC_PATH_PREPROCESSED, save_to=STOC_OBS_FULL)