# This script:
#   - takes pre-processed datasets.
#   - groups them together into one file.
#   - these files will be fed to the models (see next script).

##### Libraries #####
library(dplyr)
library(readr)
library(purrr)
library(cli)
library(sf)

source(here::here("R/utils_figures.R"))


##### Parameters #####
source(here::here("data/config/config.R")) # Import global parameters

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


##### Helper functions #####
load_group_clc_chelsa_sf <- function(verbose = VERBOSE) {
    # Load shapefiles with hexagons
    if (verbose) {
        cli_alert_info("Loading shapefiles of features...")
        cli_alert_info("Loading CLC...")
    }

    # Initialize lists
    list_shapefiles <- list(st_read(CORINE_SHP_PATH, quiet = TRUE))
    list_shp_colors <- list(read_csv(
        gsub(".gpkg", ".csv", CORINE_SHP_PATH), show_col_types = FALSE))
    list_names <- list("CLC")

    # Loop on chelsa datasets
    for (i in 1:length(CORINE_SHP_PATH)) {
        if (verbose) {
            cli_alert_info(paste0("Loading chelsa '", CHELSA_DATASETS[i], "' dataset..."))
        }
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

        list_names[[length(list_names) + 1]] <- paste0("CHELSA_", CHELSA_DATASETS[i])
    }
    if (verbose) {
        cli_alert_success("Loaded shapefiles of features!\n\n")
    }

    # Quick check : compare the geometry of the cells in each shapefile (detects mismatchs)
    same_geometries <- all(st_geometry(list_shapefiles[[1]]) == st_geometry(list_shapefiles[[2]]))
    if (same_geometries) {
        if (verbose) {
            cli_alert_success("Quick check: coordinates of hexagons in shapefiles are identical.\n\n")
        }
    } else {
        stop("Quick check: coordinates of hexagons in shapefiles are NOT identical. Cannot proceed further.")
    }


    # Group shapefiles together (into one master_sf), avoid redundance of cells
    if (verbose) {
        cli_alert_info(paste0("Grouping shapefiles together..."))
    }
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
    if (!authorize_overwrite(MASTER_SF_PATH, verbose)) {
        # Check if file exists
        if (verbose) {
            cli_alert_warning("Skipping pre-processing.\n\n")
        }
    } else {
        if (verbose) {
            cli_alert_info(paste0("Saving master shapefile..."))
        }
        st_write(master_sf, MASTER_SF_PATH, delete_dsn = TRUE, quiet = TRUE)
        if (verbose) {
            cli_alert_success(paste0("Saved master shapefile!\n\n"))
        }
    }

}

format_save_stoc_featured_dataset <- function(verbose = VERBOSE) {
    # load presence-absence observations
    if (verbose) {
        cli_alert_info("Loading stoc dataset...")
    }
    stoc_df <- read_csv(STOC_PATH_PREPROCESSED, show_col_types = FALSE)
    points_stoc <- st_as_sf(stoc_df, coords = c("LON", "LAT"), crs = 4326)
    if (verbose) {
        cli_alert_success("Loaded stoc dataset!")
        cli_alert_info("Loading master shapefile...")
    }

    # Load features
    master_sf <- st_read(MASTER_SF_PATH, quiet = TRUE)
    if (verbose) {
        cli_alert_success("Loaded master shapefile!\n\n")
        cli_alert_info("Adding features to stoc dataset...\n\n")
    }

    # Append stoc table with master_sf features
    master_stoc <- st_join(points_stoc, master_sf, join = st_within)

    # Quick check : compare the geometry of the cells in each shapefile (detects mismatchs)
    all_points_in_hexagons <- all(!is.na(master_stoc$CLC))
    if (all_points_in_hexagons) {
        if (verbose) {
            cli_alert_success("Sanity check: no stoc observations outside of shapefile hexagons.\n\n")
            cli_alert_info("Saving master stoc dataset...")
        }
    } else {
        stop("Sanity check: some observations are outside of the study area (master shapefile).")
    }

    st_write(master_stoc, MASTER_OBS_PATH, quiet = TRUE)
    if (verbose) {
        cli_alert_success("Saved master stoc dataset!")
    }
}


##### Formatting shapefile data #####
. <- load_group_clc_chelsa_sf()


##### Formatting observation data #####
. <- format_save_stoc_featured_dataset()