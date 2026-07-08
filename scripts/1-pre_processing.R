# This scripts:
#   - groups datasets splitted in different files
#   - formats datasets into similar csvs/rasters
#   - anonymizes gps coordinates

##### Libraries ##### ---------------------------------------------------------
library(tidyterra)
library(terra)
library(dplyr)
library(readr)
library(tools)
library(cli)
library(sf)

source(here::here("R/utils_data.R"))
source(here::here("R/utils_figures.R"))


##### Parameters ##### --------------------------------------------------------
source(here::here("data/config/config.R")) # Import global parameters
source(here::here("data/config/seed.R")) # hidden seed (for confidentiality)

# biodicapt dataset
BIODICAPT_PATH_PREPROCESSED <- file.path(PROCESSED_DATA_PATH, 
        "BIODICAPT_survey_data_anonymized.csv") 

# 500eni dataset
ENI500_PATH_PREPROCESSED <- file.path(PROCESSED_DATA_PATH, 
        "ENI500_survey_data_anonymized.csv") 

# Corine Land Cover dataset
CORINE_BASENAME <- file.path(PROCESSED_DATA_PATH, paste0("CLC", OBS_YEAR))
SIMPLIFIED_CORINE_PATH <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("original_CLC", OBS_YEAR, "_simplified_categories.tif"))

# Chelsa dataset
CHELSA_BASENAME <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("CHELSA_", OBS_YEAR, "_"))


# stoc dataset
STOC_PATH_PREPROCESSED <- file.path(PROCESSED_DATA_PATH, 
        "STOC_data_cleaned.csv") 


##### Helper functions ##### --------------------------------------------------
preprocess_biodicapt_dataset <- function() {
    cli_alert_info("Pre-processing of BIODICAPT dataset.")
    if (!authorize_overwrite(BIODICAPT_PATH_PREPROCESSED)) {
        # Check if file exists
        cli_alert_warning("Skipping BIODICAPT pre-processing.\n\n")

    } else {
        # preprocess file if authorized
        cli_alert_info("Loading BIODICAPT files...")

        # load file and remove useless columns
        biodicapt_df <- import_biodicapt_land_surveys(
            file.path(BIODICAPT_FOLDER, BIODICAPT_FILENAMES))
    
        # Blurr GPS data (for confidentiality of land owners)
        cli_alert_info("Data anonymization...")
        anonym_biodicapt_df <- blur_coordinates(
            biodicapt_df, "LON", "LAT", RES_KM, BLUR_SEED)

        # Save resulting dataframe
        cli_alert_info("Saving file...")
        write_csv(anonym_biodicapt_df, BIODICAPT_PATH_PREPROCESSED)
        cli_alert_success("Modified BIODICAPT file saved!\n\n")
    }
}

preprocess_eni500_dataset <- function() {
    cli_alert_info("Pre-processing of 500 ENI dataset.")

    if (!authorize_overwrite(ENI500_PATH_PREPROCESSED)) {
        # Check if file exists
        cli_alert_warning("Skipping 500 ENI pre-processing.\n\n")
        
    } else {
        # preprocess file if authorized
        cli_alert_info("Loading 500 ENI files...")

        # load file and remove useless columns
        eni500_df <- read_csv(
            file.path(ENI500_FOLDER, ENI500_FILENAME),
            show_col_types = FALSE)
        cols_to_remove <- c(
            "lieu_dit", "code_postal", "pourcent_pente", "id_parcelle",
            "commentaire_parcelle", "derniere_modif_parcelle_par", "commune",
            "derniere_modif_parcelle_le", "derniere_modif_donnees_agro_par",
            "derniere_modif_donnees_agro_le", "derniere_modif_pratiques_par",
            "derniere_modif_pratiques_le", "code_parcelle", "nom_parcelle")
        eni500_df <- eni500_df |> select(-any_of(cols_to_remove))
    
        # Rename coordinates columns
        eni500_df <- eni500_df |> rename(LON = X, LAT = Y)
          
        # Blurr GPS data (for confidentiality of land owners)
        cli_alert_info("Data anonymization...")
        anonym_eni500_df <- blur_coordinates(
            eni500_df, "LON", "LAT", RES_KM, BLUR_SEED)

        # Save resulting dataframe
        cli_alert_info("Saving file...")
        write_csv(anonym_eni500_df, ENI500_PATH_PREPROCESSED)
        cli_alert_success("Modified ENI 500 file saved!\n\n")
    }
}

preprocess_raster <- function(
        raster_path,
        save_to_basename) {
    cli_alert_info(paste0("Pre-processing ", basename(raster_path), "..."))

    # Creating the paths that will be needed in the function
    raster_reprojected_path <- paste0(
        file_path_sans_ext(save_to_basename), 
        "_projection_france_res", RES_KM,"km-WGS84",
        ".tif")
        
    if (!authorize_overwrite(raster_reprojected_path)) {
        # Check if file exists
        cli_alert_warning("Skipping pre-processing of raster.\n\n")

    } else {
        cli_alert_info("Loading raster...")
        # import and reproject to a standardized grid
        raster <- rast(raster_path)
        . <- project_to_france_custom_grid(
                raster = raster, 
                save_to = raster_reprojected_path,  
                res_km = RES_KM) 
        cli_alert_success("Modified raster file saved!\n\n")
    }
}

preprocess_raster_hexagonal_sf <- function(
        raster_path,
        save_to_basename,
        extraction_mode) {
    cli_alert_info(paste0("Pre-processing of ", basename(raster_path), "..."))

    # Creating the paths that will be needed in the function
    shapefile_hexagons_path <- paste0(
        file_path_sans_ext(save_to_basename), 
        "_projection_france_hexagons_res", RES_KM,"km-WGS84",
        ".gpkg")


    if (!authorize_overwrite(shapefile_hexagons_path)) {
        # Check if file exists
        cli_alert_warning("Skipping pre-processing.\n\n")

    } else {
        cli_alert_info("Loading cropped raster...")
        # import and reproject to a standardized grid of hexagons
        raster <- rast(raster_path)
        . <- project_to_hexagons(
                raster = raster, 
                save_to = shapefile_hexagons_path, 
                res_km = RES_KM,
                method = extraction_mode) 
        cli_alert_success("Saved shapefile and CSV legend!\n\n")
    }
}

before_preprocessing_corine_raster <- function(
        raster_path,
        save_to) {
    cli_alert_info("Simplification of CLC2018 raster.")
    if (!authorize_overwrite(save_to)) {
        # Check if file exists
        cli_alert_warning("Skipping simplification of CLC2018.\n\n")

    } else {
        cli_alert_info("Loading raster...")
        raster <- rast(raster_path)

        # CLC contains (too) many categories. Subcategories can be grouped
        # together to avoid having too many modalities.
        # We keep all categories at level 1 (i.e. lowest information) ...
        # ... except crops, for which we get level 2 subcategories 
        # (level 3 is also available)
        cli_alert_info("Simplifying categories...")
        simplified_clc_raster <- simplify_CLC(
            clc_raster = raster,
            level_urban = 1,
            level_crops = 2,
            level_forests = 1,
            level_wetlands = 1,
            level_water = 1)
        
        cli_alert_info("Saving file...")
        writeRaster(simplified_clc_raster,  save_to, overwrite = TRUE)
        cli_alert_success(paste0(
            "Created simplified CORINE land cover raster file!\n\n"))
    }
}

before_preprocessing_chelsa_raster <- function(
        dataset,
        save_to) {
    cli_alert_info("Aggregation of CHELSA monthly data.")
    
    if (!authorize_overwrite(save_to)) {
        # Check if file exists
        cli_alert_warning(paste0(
            "Skipping pre-processing of monthly '", 
            dataset, "' CHELSA data.\n\n"))

    } else {
        local_chelsa_folder <- file.path(CHELSA_FOLDER, dataset)
        chelsa_12_paths <- file.path(
            local_chelsa_folder, 
            paste0(
                "CHELSA_", dataset,
                sprintf("_%02d_%s_V.2.1.tif", 1:12, OBS_YEAR)
            ))
            

        # default format is per month, we need yearly data
        cli_alert_info("Aggregating rasters together...")
        chelsa_annual_raster <- monthly_2_yearly_rasters(
            raster_paths = chelsa_12_paths, 
            buffer = RES_KM*5, 
            verbose = FALSE,
            fun = mean)
        
        # data is in kelvin, convert to celsius
        chelsa_annual_raster_celsius <- chelsa_annual_raster - 273.15
        
        cli_alert_info("Saving file...")
        writeRaster(chelsa_annual_raster_celsius,  save_to, overwrite = TRUE)
        cli_alert_success(paste0(
            "Created aggregated '", dataset, "' CHELSA raster file!\n\n"))
    }
}

preprocess_stoc_dataset <- function() {
    cli_alert_info("Pre-processing of the STOC dataset.")

    if (!authorize_overwrite(STOC_PATH_PREPROCESSED)) {
        # Check if file exists
        cli_alert_warning("Skipping STOC pre-processing.\n\n")

    } else {
        # preprocess file if authorized
        cli_alert_info("Loading STOC dataframe...")

        # loading .Rdata --> variable name: STOC
        load(STOC_FILEPATH)
        
        # select year (STOC contains years 2015-2018)
        if (OBS_YEAR %in% unique(STOC$annee.x)){
            stoc_subset <- subset(STOC, annee.x == OBS_YEAR)
        } else {
            stop(paste0(
                "Selected year '", OBS_YEAR, 
                "' is not in the loaded STOC dataset"))
        }
        
        # convert abundance to occurence
        stoc_subset <- stoc_subset |> 
            mutate(across(23:64, ~ 1L * (. > 0))) 
        # # remove useless columns (only if really necessary to reduce space)
        # cols_to_remove <- c(
        #     "NDVI", "light_pollution", "annee.x", "precip_spring", 
        #     "tmp_spring", "p_type", "p_milieu", "altitude", "urbain_p", 
        #     "departement", "habitat_principal", "habitat_secondaire", 
        #     "ouvert_p", "foret_p", "agricole_p")
        # stoc_subset <- stoc_subset |> select(-any_of(cols_to_remove))
        # Rename coordinates columns
        stoc_subset <- stoc_subset |> 
            rename(LON = longitude_wgs84, LAT = latitude_wgs84)
          
        # No blurring necessary (this is a public dataset)
        # Save resulting dataframe
        cli_alert_info("Saving file...")
        write_csv(stoc_subset, STOC_PATH_PREPROCESSED)
        cli_alert_success("Modified STOC file saved!\n\n")
    }
}

show_save_results <- function() {
    see_datasets <- readline(
        prompt = paste0(paste0(
            "Show processed datasets & Save figures ",
            "(overwrites by default)? [Y/n]: ")))
    cleaned_answer <- tolower(trimws(see_datasets))

    if (cleaned_answer %in% c("y", "yes")) {
        cli_alert_info("Showing plots to the user.")
        cli_alert_warning(paste0(
            "All figures automatically are saved as .pdf inside '", 
            FIGURES_PATH, "' ."))

        
        # 1. biodicapt dataset
        biodicapt_df <- read_csv(
            BIODICAPT_PATH_PREPROCESSED, 
            show_col_types = FALSE)
        biodicapt_plot <- ggplot_categorical_df_on_background_map(
            background_map = ggplot_get_france_base_map("national"), 
            df = biodicapt_df, 
            LON = "LON",
            LAT = "LAT",
            column = "network")
        print(biodicapt_plot)
        standardised_ggplot_save(
            figure = biodicapt_plot, 
            save_path = file.path(FIGURES_PATH, "biodicapt_networks.pdf"))

        # 2. 500 eni dataset
        eni500_df <- read_csv(
            ENI500_PATH_PREPROCESSED,
            show_col_types = FALSE)
        eni500_plot <- ggplot_categorical_df_on_background_map(
            background_map = ggplot_get_france_base_map("national"), 
            df = eni500_df, 
            LON = "LON",
            LAT = "LAT")
        print(eni500_plot)
        standardised_ggplot_save(
            figure = eni500_plot, 
            save_path = file.path(FIGURES_PATH, "500ENI_sites.pdf"))

        # # 3. Corine Land Cover Dataset
        # # raster
        # corine_raster <- rast(
        #     paste0(CORINE_BASENAME, "_projection_france_res", 
        #            RES_KM, "km-WGS84.tif")
        # )
        # corine_plot <- ggplot_categorical_raster_on_background_map(
        #     background_map = ggplot_get_france_base_map("national"), 
        #     raster = corine_raster,
        #     layer_name = "NEW_LABEL3")
        # print(corine_plot)
        # standardised_ggplot_save(
        #     figure = corine_plot, 
        #     save_path = file.path(FIGURES_PATH, "corine_raster.pdf"))
        # # shapefile  
        corine_shapefile <- st_read(
            paste0(CORINE_BASENAME, 
                "_projection_france_hexagons_res", RES_KM, "km-WGS84.gpkg"),
            quiet = TRUE)
        corine_shapefile_colors <- read_csv(
            paste0(CORINE_BASENAME, "_projection_france_hexagons_res", 
                   RES_KM, "km-WGS84.csv"),
            show_col_types = FALSE)
        corine_plot_bis <- ggplot_categorical_shapefile_on_background_map(
            background_map = ggplot_get_france_base_map("national"), 
            shapefile = corine_shapefile,
            layer_name = "dominant_class",
            color_df = corine_shapefile_colors)
        print(corine_plot_bis)
        standardised_ggplot_save(
            figure = corine_plot_bis, 
            save_path = file.path(FIGURES_PATH, "corine_hexagons.pdf"))    
        
        
        # 4. chelsa datasets
        for (i in 1:length(CHELSA_DATASETS)) {
            chelsa_dataset_basename <- paste0(
                CHELSA_BASENAME, CHELSA_DATASETS[i])
            # # raster
            # chelsa_tas_raster <- rast(
            #     paste0(chelsa_dataset_basename, 
            #         "_projection_france_res", RES_KM, "km-WGS84.tif"))
            # chelsa_tas_plot <- ggplot_quantitative_raster_on_background_map(
            #     background_map = ggplot_get_france_base_map("national"), 
            #     raster = chelsa_tas_raster,
            #     layer_name = "mean",
            #     unit="°C",
            #     limits=NULL)
            # print(chelsa_tas_plot)
            # standardised_ggplot_save(
            #     figure = chelsa_tas_plot, 
            #     save_path = file.path(
            #         FIGURES_PATH, 
            #         paste0("chelsa_", CHELSA_DATASETS[i], "_raster.pdf")))
            # # shapefile  
            chelsa_shapefile <- st_read(
                paste0(chelsa_dataset_basename, 
                    "_projection_france_hexagons_res", RES_KM, "km-WGS84.gpkg"),
                quiet = TRUE)
            chelsa_plot_bis <- ggplot_quantitative_shapefile_on_background_map(
                background_map = ggplot_get_france_base_map("national"), 
                shapefile = chelsa_shapefile,
                layer_name = "dominant_class",
                unit=CHELSA_UNITS[i],
                limits=NULL,
                precision_auto_limits = 0.1)
            print(chelsa_plot_bis)
            standardised_ggplot_save(
                figure = chelsa_plot_bis, 
                save_path = file.path(FIGURES_PATH, 
                    paste0("chelsa_", CHELSA_DATASETS[i], "_hexagons.pdf")))    
        }
        
        # 5. stoc dataset
        sp_to_show = "Pica_pica"
        stoc_df <- read_csv(STOC_PATH_PREPROCESSED, show_col_types = FALSE)
        stoc_df <- stoc_df |>
            mutate(
                across(all_of(sp_to_show), 
                ~ ifelse(. == 1, "presence", "absence"), 
                .names = "{col}"))
        cli_alert_warning(paste(
            "Stoc dataset: showing sightings of", 
            gsub("_", " ", sp_to_show), 
            "only."))

        stoc_plot <- ggplot_categorical_df_on_background_map(
            background_map = ggplot_get_france_base_map("national"), 
            df = stoc_df, 
            LON = "LON",
            LAT = "LAT",
            column = "Pica_pica",
            legend_title = paste0(
                "Sightings of ", 
                gsub("_", " ", sp_to_show), "."))
        print(stoc_plot)
        standardised_ggplot_save(
            figure = stoc_plot, 
            save_path = file.path(
                FIGURES_PATH, 
                paste0("stoc_sightings_", sp_to_show,".pdf")))

            
        cli_alert_success("Plots and PDFs are ready!")
    } else {
        cli_alert_warning("Skipping.\n\n")
    }
}


##### Transformation of raw datasets ##### ------------------------------------
if (exists("BLUR_SEED")) {
    # 1. Biodicapt dataset
    . <- preprocess_biodicapt_dataset()

    # 2. 500 eni dataset
    . <- preprocess_eni500_dataset()

    # 3. Corine Land Cover dataset
    . <- before_preprocessing_corine_raster(
        raster_path = CORINE_FILEPATH,
        save_to = SIMPLIFIED_CORINE_PATH)
    . <- preprocess_raster(
        raster_path = SIMPLIFIED_CORINE_PATH,
        save_to_basename = CORINE_BASENAME)
    . <- preprocess_raster_hexagonal_sf(
        raster_path = SIMPLIFIED_CORINE_PATH,
        save_to_basename = CORINE_BASENAME,
        extraction_mode = "categorical")

    # 4. Chelsa dataset 
    for (dataset_name in CHELSA_DATASETS) {
        chelsa_dataset_basename <- paste0(CHELSA_BASENAME, dataset_name)
        chelsa_annual_path <- paste0(
            chelsa_dataset_basename, "_annual-WGS84.tif")

        . <- before_preprocessing_chelsa_raster(
            dataset = dataset_name,
            save_to = chelsa_annual_path)
        . <- preprocess_raster(
            raster_path = chelsa_annual_path,
            save_to_basename = chelsa_dataset_basename)
        . <- preprocess_raster_hexagonal_sf(
            raster_path = chelsa_annual_path,
            save_to_basename = chelsa_dataset_basename,
            extraction_mode = "mean")
    }

    # 5. Stoc dataset
    . <- preprocess_stoc_dataset()
} else {
    cli_alert_warning(paste0(
        "Variable 'BLUR_SEED' is not available, ",
        "you are likely running this script without raw datasets installed. ",
        "Skipping pre-processing (will only be showing figures)."))
}

##### Check results ##### -----------------------------------------------------
. <- suppressWarnings(show_save_results())
