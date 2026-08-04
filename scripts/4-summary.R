# This script is used to:
#   - Summarise scores of scripts #3 into barplots and lineplots
#   - Compare
##### Libraries ##### ---------------------------------------------------------
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(cli))

suppressPackageStartupMessages(source(here::here(file.path("R", "utils_figures.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_models.R"))))


##### Parameters ##### --------------------------------------------------------
source(here::here(file.path("data","config","config.R"))) # Global parameters

### Seed
set.seed(496) # for reproducible results

### Load variables used to generate results
SUBFOLDER <- "STOC-OED"
parameters <- readRDS(file.path(RESULTS_PATH, SUBFOLDER, "parameters.rds"))
COMBINATIONS <- parameters$COMBINATIONS
param_grid <- build_param_grid(COMBINATIONS)


##### Helper functions ##### --------------------------------------------------
cli_alert_warning("Reference values are *hard coded* in 'loop_selector'.")
loop_selector <- function(grid, loop_on) {
    # initialize base values
    combination <- list(
        train_size = c(125),
        r_effect = c("none"),
        strategy = c("none"),
        formulas = c(~ (NDVI + light_pollution + p_milieu + altitude + precip_spring + tmp_spring)),
        n_new_samples = c(0)
    )

    # replace values for loop_on by list
    combination[[loop_on]] <- unique(param_grid[[loop_on]])

    # replace formulas by number of variables
    combination$formulas <- lapply(
            lapply(combination$formulas, all.vars), 
        length)

    path_sections <- list("", "")
    path_idx <- 1
    parameter_names <- c(
        "r_effect", "strategy", "n_new_samples", "train_size", "formulas")
    path_names <- c(
        "model_random-", "_strategy-", "_new-samples-", "_training-size-", "_n-variables-")

    for (p in seq_along(parameter_names)) {
        if (length(combination[[parameter_names[p]]]) > 1) {
            path_sections[[path_idx]] <- paste0(
                path_sections[[path_idx]],
                path_names[p]
            )
            path_idx <- path_idx + 1
        } else {
            path_sections[[path_idx]] <- paste0(
                path_sections[[path_idx]],
                path_names[p], combination[[parameter_names[p]]]
            )
        }
    }
    path_sections[[2]] <- paste(path_sections[[2]], "_k")

    return(path_sections)
}


##### Main ##### --------------------------------------------------------------
cli_alert_info("------------ Score comparison ------------\n\n")
for (column in names(param_grid)) {
    filename_sections <- loop_selector(param_grid, column)
    loop_list <- unique(param_grid[[column]])

    if (column == "formulas") {
        loop_list <- lapply(lapply(loop_list, all.vars), length)
    }

    . <- barplot_raw_scores(
        parent_folder = file.path(RESULTS_PATH, SUBFOLDER), 
        save_to = paste0("barplot_", tolower(name), "_scores.pdf"),
        loop_prefix = filename_sections[[1]], 
        loop_elements = loop_list, 
        loop_suffix = filename_sections[[2]], 
        k_fold = parameters$K_FOLDS, 
        xlabel = paste0("Effect of ", tolower(name)," type on metrics"), 
        ylabel = "Average score per species",
        group_species = FALSE,
        species_names = parameters$Y_SPECIES,
        barplot = ifelse(name == "TRAIN_SIZES", FALSE, TRUE))

    . <- dotwhisker_compare_scores(
        parent_folder = file.path(RESULTS_PATH, SUBFOLDER), 
        reference_model_path = paste0(prefix, loop_list[1], suffix), 
        save_to = paste0("dotwhisker_", tolower(name), "_comparison.pdf"),
        loop_prefix = filename_sections[[1]], 
        loop_elements = loop_list[-1], 
        loop_suffix = filename_sections[[2]], 
        k_fold = parameters$K_FOLDS, 
        metric = "MSE",
        bars = "CI",
        subset_names = c("train", "val", "test"),
        xlabel = paste0(
            "Effect of ", tolower(name)," compared to '", loop_list[1], "'"), 
        group_species = TRUE,
        species_names = parameters$Y_SPECIES)

    . <- boxplot_compare_scores(
        parent_folder = file.path(RESULTS_PATH, SUBFOLDER), 
        reference_model_path = paste0(prefix, loop_list[1], suffix), 
        save_to = paste0("boxplot_", tolower(name), "_comparison.pdf"),
        loop_prefix = filename_sections[[1]], 
        loop_elements = loop_list[-1], 
        loop_suffix = filename_sections[[2]], 
        k_fold = parameters$K_FOLDS, 
        metric = "MSE",
        subset_names = c("train", "val", "test"),
        xlabel = paste0(
            "Effect of ", tolower(name)," compared to '", loop_list[1], "'"), 
        group_species = TRUE,
        species_names = parameters$Y_SPECIES)
}