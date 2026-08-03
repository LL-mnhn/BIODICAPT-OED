# This script is used to:
#   - Summarise scores of scripts #3 into barplots and lineplots
#   - Compare
##### Libraries ##### ---------------------------------------------------------
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


##### Main ##### --------------------------------------------------------------
cli_alert_info("------------ Score comparison ------------\n\n")
for (c in seq_along(COMBINATIONS)) {
    # if only one vector in the list contains several elements
    # then we loop on it to compare performances
    if (sum(lengths(COMBINATIONS[[c]]) > 1) <= 1) {
        name <- names(which(lengths(COMBINATIONS[[c]]) > 1))

        if (identical(name, character(0))) {
            prefix <- paste0("model_random-", COMBINATIONS[[c]]$R_EFFECTS,
                            "_strategy-", COMBINATIONS[[c]]$STRATEGIES,
                            "_training-size")
            loop_on <- "-" # dummy loop to show results of 1-element list
            suffix <- paste0("_k")
        } else {
            if (name == "R_EFFECTS") {
            prefix <- paste0("model_random-")
            suffix <- paste0("_strategy-", COMBINATIONS[[c]]$STRATEGIES,
                            "_training-size-", COMBINATIONS[[c]]$TRAIN_SIZES,
                            "_k")
            } else if (name == "STRATEGIES") {
                prefix <- paste0("model_random-", COMBINATIONS[[c]]$R_EFFECTS,
                                "_strategy-")
                suffix <- paste0("_training-size-", COMBINATIONS[[c]]$TRAIN_SIZES,
                                "_k")
            } else if (name == "TRAIN_SIZES") {
                prefix <- paste0("model_random-", COMBINATIONS[[c]]$R_EFFECTS,
                                "_strategy-", COMBINATIONS[[c]]$STRATEGIES,
                                "_training-size-")
                suffix <- paste0("_k")
            } else {
                stop("Unidentified error in score comparison.")
            }
            loop_on <- COMBINATIONS[[c]][name][[1]]
        }

        . <- compute_hmsc_performances(
            parent_folder = file.path(RESULTS_PATH, SUBFOLDER), 
            save_to = paste0("performance_", tolower(name), "_scores.pdf"),
            loop_prefix = prefix, 
            loop_elements = loop_on, 
            loop_suffix = suffix, 
            k_fold = parameters$K_FOLDS, 
            xlabel = paste0("Effect of ", tolower(name)," type on metrics"), 
            ylabel = "Average score per species",
            group_species = FALSE,
            species_names = parameters$Y_SPECIES,
            barplot = ifelse(name == "TRAIN_SIZES", FALSE, TRUE))
    
        . <- fine_compare_hmsc_metric(
            parent_folder = file.path(RESULTS_PATH, SUBFOLDER), 
            reference_model_path = paste0(prefix, loop_on[1], suffix), 
            save_to = paste0("performance_", tolower(name), "_comparison.pdf"),
            loop_prefix = prefix, 
            loop_elements = loop_on[-1], 
            loop_suffix = suffix, 
            k_fold = parameters$K_FOLDS, 
            metric = "MSE",
            subset_names = c("train", "val", "test"),
            xlabel = paste0("Effect of ", tolower(name)," compared to '", loop_on[1], "'"), 
            group_species = TRUE)
    }

}