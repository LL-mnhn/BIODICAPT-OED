##### Global Parameters #####
# These parameters are 'hidden': they can be modified, but you should not 
# change them without a *very* good reason. 
#
# The parameters grouped in this file are re-used by scripts using STOC data
suppressPackageStartupMessages(library(cli))
suppressMessages(suppressWarnings(source(here::here("R/utils_data.R"))))


### Paths ---------------------------------------------------------------------
PATH_LOCAL_RESULTS <- file.path(RESULTS_PATH, "STOC-OED")
PATH_LOCAL_BASE <- file.path(RESULTS_PATH, "STOC-OED", "base-models")
PATH_LOCAL_OED <- file.path(RESULTS_PATH, "STOC-OED", "OED")


### Dataset -------------------------------------------------------------------
X_VARIABLES <- c("p_milieu", "tmp_spring", "precip_spring") # custom selection
X_GROUPS <- c(1, 2, 2) # Group above variables in different categories
X_GROUPS_NAMES <- c("niche", "climate") # names associated with X_GROUPS categories
Y_SPECIES <- c("Sylvia_atricapilla", "Fringilla_coelebs", "Pica_pica", 
    "Periparus_ater", "Carduelis_cannabina") # can be NULL (all species)
X_FACTORS <- ("p_milieu") # just to be sure of column types


### Split ---------------------------------------------------------------------
TRAIN_SIZE <- 125 # Number of point used for training
NEW_POOL_SIZE <- 500 # Number of points from which we can pick for OED
NEW_SAMPLE_SIZE <- 50 # Number of points we can sample in new pool during OED
K_FOLD <- 5


### HMSC ----------------------------------------------------------------------
HMSC_XFORMULA <- ~ p_milieu + tmp_spring + precip_spring

SAMPLES <- 5000 # mcmc will stop after saving that much samples
THIN <- 2 # number of steps between each recording
TRANSIENT <- 0.5*SAMPLES*THIN # burn-in iterations
CHAINS <- 3


### Helpers -------------------------------------------------------------------
check_dataset <- function(path) {
    cli_alert_info(paste0("Checking ", basename(path),"..."))

    if (!file.exists(path)) {
        cli_alert_danger(paste0("File '", path, "' not found!\n"))
        stop(paste0("Tip: did you run `1-pre_processing.R` and ",
        "`2-feature_aggregation.R`?"))
    }

    cli_alert_success("File is available.\n\n")
    return(read_csv(path, show_col_types = FALSE))
}