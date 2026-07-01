# This scripts checks the presence of files needed in ./data/raw_data

##### Libraries ##### ---------------------------------------------------------
library(cli)


##### Parameters ##### --------------------------------------------------------
source(here::here("data/config/config.R")) # Import global parameters


##### Helper functions #####
check_biodicapt_files <- function() {
    cli_alert_info("Checking BIODICAPT files...")

    # The folder...
    if (!file.exists(BIODICAPT_FOLDER)) {
        cli_alert_danger(
            paste0("Folder '", BIODICAPT_FOLDER, "' was not found!\n"))
        stop("Tip: did you create the necessary folders for the files?")
    }
    
    # ...must contain 5 files with the following names :
    for (filename in BIODICAPT_FILENAMES) {
        if (!file.exists(file.path(BIODICAPT_FOLDER, filename))) {
            cli_alert_danger(paste0("File '", filename, "' was not found in '", 
                BIODICAPT_FOLDER,"'!\n"))
            cli_alert_warning(paste(
                "This dataset can only be obtained from the coordinators",
                "of the BIODICAPT project. Contact them directly."))
            stop("Tip: these files must be downloaded manually.")
        }
    }
    
    cli_alert_success("Raw BIODICAPT files are available.\n\n")
}

check_eni500_files <- function() {
    cli_alert_info("Checking 500 ENI file...")

    # The folder...
    if (!file.exists(ENI500_FOLDER)) {
        cli_alert_danger(paste0("Folder '", ENI500_folder, 
            "' was not found!\n"))
        stop("Tip: did you create the necessary folders for the files?")
    }
    
    # ...must contain 1 file with the following name :
    if (!file.exists(file.path(ENI500_FOLDER, ENI500_FILENAME))) {
        cli_alert_danger(paste0("File '", filename, "' was not found in '", 
            ENI500_FOLDER,"'!\n"))
        cli_alert_warning(paste0("This dataset can only be obtained from ", 
        "the coordinators of the 500 ENI network. Contact them directly."))
        stop("Tip: this file must be downloaded manually (the version with NAs removed).")
    }
    
    cli_alert_success("Raw 500 ENI file is available.\n\n")
}

check_corine_raster <- function() {
    cli_alert_info("Checking CORINE Land Cover files...")

    # Check if file exists
    if (!file.exists(CORINE_FILEPATH)) {
        cli_alert_danger(
            paste0("File '", CORINE_FILEPATH, "' was not found!\n"))
        cli_alert_warning("If needed, download file from https://www.data.gouv.fr/datasets/corine-land-cover-edition-2018-france-metropolitaine")
        stop("Tip: did you create the necessary folders for the files?")
    }

    cli_alert_success("Raw CLC files are available.\n\n")
}

check_chelsa_rasters <- function(dataset) {
    cli_alert_info("Checking CHELSA files...")

    local_CHELSA_FOLDER <- file.path(CHELSA_FOLDER, dataset)
  
    # The folder...
    if (!file.exists(local_CHELSA_FOLDER)) {
        cli_alert_danger(paste0("Folder '", local_CHELSA_FOLDER, 
            "' was not found!\n"))
        stop("Tip: did you create the necessary folders for the files?")
    }

    # should contain 12 files (one per month of the year)
    CHELSA_12_paths <- file.path(local_CHELSA_FOLDER,
        paste0(
            "CHELSA_", dataset,
            sprintf("_%02d_%s_V.2.1.tif", 1:12, OBS_YEAR))
        )
    for (CHELSA_tiff_path in CHELSA_12_paths){
        if (!file.exists(CHELSA_tiff_path)){
            cli_alert_danger(
                paste0("File '", CHELSA_tiff_path, "' was not found!\n"))
            cli_alert_warning("If needed, download files from https://www.chelsa-climate.org/datasets/chelsa_monthly.")
            stop(paste0("Tip: download monthly CHELSA files for your dataset ",
            "and put them together in the same subfolder."))
        }
    }

    cli_alert_success("Raw CHELSA files are available.\n\n")
}

check_species_data <- function() {
    cli_alert_info("Checking STOC data file...")

    # Check if file exists
    if (!file.exists(STOC_FILEPATH)) {
        cli_alert_danger(paste0("File '", STOC_FILEPATH, "' was not found!\n"))
        cli_alert_warning("If needed, download file from https://doi.org/10.5061/dryad.bnzs7h4g3.")
        stop("Tip: did you create the necessary folders for the files?")
    }

    cli_alert_success("Raw STOC dataframe is available.\n\n")
}


##### Verify each dataset ##### -----------------------------------------------
cli_alert_warning(paste0(
    "This script can only be run on raw datasets, it should not be run ",
    "on external machines.\nIf you cloned this repository from GitHub, ",
    "ignore this script and begin usage with `1-pre_processing.R`."))

# 1. BIODICAPT dataset
. <- check_biodicapt_files()

# 2. 500 ENI dataset
. <- check_eni500_files()

# 3. CORINE Land Cover dataset
. <- check_corine_raster()

# 4. CHELSA tas dataset (tas: near-surface air temperature)
for (chelsa_dataset in CHELSA_DATASETS) {
    . <- check_chelsa_rasters(dataset = chelsa_dataset)
}

# 5. Species absence-presence dataset
. <- check_species_data()