# This script is the 'main' script of this project it:
#   - Makes a quick sanity check (verifies that necessary datasets are available)
#   - Uses several models to predict species distributions (SDM/jSDM)
#   - Finds locally optimum OEDs

##### Libraries ##### --------------------------------------------------
library(cli)
library(readr)


##### Parameters ##### -------------------------------------------------
source(here::here("data/config/config.R")) # Import global parameters

### Models
MODELS <- c("HMSC")

### OED


##### Helper functions ##### -------------------------------------------
check_dataset <- function(path, verbose = VERBOSE) {
    if (verbose) {
        cli_alert_info(paste0("Checking ", basename(path),"..."))
    }

    if (!file.exists(path)) {
        cli_alert_danger(paste0("File '", path, "' not found!\n"))
        stop("Tip: did you run `1-pre_processing.R` and `2-feature_aggregation.R`?")
    }

    if (verbose) {
        cli_alert_success("File is available.\n\n")
    }
}


##### Verify available datasets ##### ----------------------------------
. <- check_dataset(BIODICAPT_OBS_FULL)
. <- check_dataset(ENI500_OBS_FULL)
. <- check_dataset(STOC_OBS_FULL)


##### Import datasets ##### --------------------------------------------
# if datasets are available, they should load with no problem
biodicapt_csv <- read_csv(BIODICAPT_OBS_FULL, show_col_types = FALSE)
eni500_csv <- read_csv(ENI500_OBS_FULL, show_col_types = FALSE)
stoc_csv <- read_csv(STOC_OBS_FULL, show_col_types = FALSE)


##### Compare datasets ##### -------------------------------------------
# repartition of clc categories
table(biodicapt_csv$CLC)
table(eni500_csv$CLC)
table(stoc_csv$CLC)
# TODO: update column names, convert to proportions, aggregate the 3 tables in 1.
# TODO: maybe we can group them together before.

# distribution of chelsa values
# TODO

##### Select/split datasets ##### --------------------------------------
# TODO