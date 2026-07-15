##### Global Parameters #####
# These parameters are 'hidden': they can be modified, but you should not 
# change them without a *very* good reason. 
#
# The parameters grouped in this file are re-used by several scripts, having 
# them all in the same place allows for easier access, control and 
# harmonization.
library(ggplot2)
library(colorspace)


### Global variables ----------------------------------------------------------
OBS_YEAR <- "2018"  # Most recent year for CLC dataset
RES_KM <- 10        # best compromise I think
# Values for map extent (France limits with buffer, in CRS 4326)
LON_MIN <- -5
LON_MAX <- 9.55
LAT_MIN <- 41.35
LAT_MAX <- 51.05


### Global paths --------------------------------------------------------------
WORK_DIR <- file.path(".")
RAW_DATA_PATH <- file.path(WORK_DIR, "data", "raw_data")
PROCESSED_DATA_PATH <- file.path(WORK_DIR, "data", "preprocessed_data")
FIGURES_PATH <- file.path(WORK_DIR, "outputs", "figures")
RESULTS_PATH <- file.path(WORK_DIR, "outputs", "results")


### File paths ----------------------------------------------------------------
## Raw data
# BIODICAPT Dataset
BIODICAPT_FOLDER <- file.path(RAW_DATA_PATH, "BIODICAPT_pos")
BIODICAPT_FILENAMES <- c(
    "BIODICAPT_parcelles_MONTPELLIER.xlsx", 
    "Preselection_Parcelles_BIODICAPT_DYNAFOR.xlsx", 
    "Preselection_Parcelles_BIODICAPT_SCARABEE.xlsx", 
    "Preselection_Parcelles_BIODICAPT_VCG.xlsx", 
    "Preselection_Parcelles_BIODICAPT_ZAAR.xlsx"
)

# 500 ENI Dataset
ENI500_FOLDER <- file.path(RAW_DATA_PATH, "500ENI_pos")
ENI500_FILENAME <- "table_parcelle_NA_excluded.csv"

# CORINE Dataset
CORINE_FILEPATH <- file.path(RAW_DATA_PATH, 
    paste0("u", OBS_YEAR, "_clc", OBS_YEAR, "_v2020_20u1_raster100m"), "DATA", 
    paste0("U", OBS_YEAR, "_CLC", OBS_YEAR, "_V2020_20u1.tif")
)

# CHELSA Dataset
CHELSA_DATASETS <- c("tas", "spi12") # CHELSA is a database with many datasets
CHELSA_UNITS <- c("°C", "Index")
CHELSA_FOLDER <- file.path(RAW_DATA_PATH, "CHELSA-monthly")

# STOC Dataset
STOC_FILEPATH <- file.path(RAW_DATA_PATH, "STOC", "Dataframe_STOC_input.RData")

## Preprocessed data
MASTER_SF_PATH <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("MASTER_", OBS_YEAR, "_projection_france_hexagons_res", RES_KM, "km-WGS84.gpkg"))
BIODICAPT_OBS_FULL <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("BIODICAPT_", OBS_YEAR, "_obs_features_res", RES_KM, "km-WGS84.csv"))
ENI500_OBS_FULL <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("ENI500_", OBS_YEAR, "_obs_features_res", RES_KM, "km-WGS84.csv"))
STOC_OBS_FULL <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("STOC_", OBS_YEAR, "_obs_features_res", RES_KM, "km-WGS84.csv"))

# Species names
NAMES_SPECIES <- c(
    "Alauda_arvensis", "Anthus_trivialis", "Carduelis_cannabina",
    "Carduelis_carduelis", "Carduelis_chloris", "Certhia_brachydactyla",  
    "Columba_palumbus", "Corvus_corone", "Cyanistes_caeruleus",
    "Dendrocopos_major", "Emberiza_cirlus", "Emberiza_citrinella",
    "Erithacus_rubecula", "Fringilla_coelebs", "Garrulus_glandarius",
    "Hippolais_polyglotta", "Hirundo_rustica", "Luscinia_megarhynchos",
    "Motacilla_alba", "Parus_major", "Passer_domesticus", "Periparus_ater",
    "Phoenicurus_ochruros", "Phylloscopus_collybita", "Pica_pica",
    "Regulus_ignicapilla", "Saxicola_rubicola", "Serinus_serinus",
    "Sitta_europaea", "Streptopelia_decaocto", "Streptopelia_turtur",
    "Sturnus_vulgaris", "Sylvia_atricapilla", "Sylvia_communis",
    "Troglodytes_troglodytes", "Turdus_merula", "Turdus_philomelos",
    "Turdus_viscivorus"     
)


### Plot styling --------------------------------------------------------------
FONT <- "Lexend"
PALETTE <- c("#D9054E", "#28A349", "#246CBC", 
             "#5D7B84", "#C2562F", "#FFB703")
SHAPES <- c(21, 22, 24, 23, 25, 8)
SIZES <- c(1.66, 1.85, 1.5, 1.66, 1.5, 1.66)
CUSTOM_SCALES <- list(
    scale_color_manual(values = darken(PALETTE, amount = 0.66)),
    scale_fill_manual(values = PALETTE),
    scale_shape_manual(values = SHAPES),
    scale_size_manual(values = SIZES)
)
LIGHT_CUSTOM_SCALES <- list(
    scale_color_manual(values = PALETTE),
    scale_fill_manual(values = PALETTE),
    scale_shape_manual(values = SHAPES),
    scale_size_manual(values = SIZES)
)