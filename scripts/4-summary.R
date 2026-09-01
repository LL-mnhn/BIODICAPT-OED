# This script is used to:
#   - Summarise scores of scripts #3 into barplots and lineplots
#   - Compare
##### Libraries ##### ---------------------------------------------------------
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(cli))

suppressPackageStartupMessages(source(here::here(file.path(
    "R", "utils_figures.R"))))
suppressPackageStartupMessages(source(here::here(file.path(
    "R", "utils_models.R"))))


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
for (combination in COMBINATIONS) {
    if (sum(lapply(combination, length) > 1) > 1) {
        stop(paste(
            "Error in combination, found more than one parameter",
            "with a vector of values of length > 1."
        ))
    } else if (sum(lapply(combination, length) > 1) == 1) {
        loop_on <- names(which(lapply(combination, length) > 1))
        # if list is numerical, sort it
        if (any(sapply(combination[[loop_on]], is.numeric))) {
            combination[[loop_on]] <- sort(combination[[loop_on]])
        }
        ref_model <- list()
        loop_model <- list()

        for (name in names(combination) ) {
            if (name == loop_on) {
                ref_model[[name]] <- combination[[name]][1]
                loop_model[[name]] <- combination[[loop_on]][-1]
            } else {
                ref_model[[name]] <- combination[[name]]
                loop_model[[name]] <- combination[[name]] 
            }
        }
    } else {
        loop_on <- "single_model"
        ref_model <- NULL
        loop_model <- combination
    }

    cli_alert_info(paste0("Plotting ", loop_on, "..."))


    for (bool in c(FALSE, TRUE)) {
        if (bool == FALSE) {
            sp <- "_per_species"
        } else {
            sp <- ""
        }

        . <- suppressMessages(boxplot_compare_scores(
            parent_folder = file.path(RESULTS_PATH, SUBFOLDER), 
            reference_model_combination = ref_model, 
            loop_model_combination = loop_model, 
            loop_on = loop_on,
            k_fold = parameters$K_FOLDS, 
            metric = "MSE",
            subset_names = c("train", "val", "test"),
            xlabel = paste0(
                "Effect of ", tolower(loop_on),
                ", compared to '", ref_model[[loop_on]], "'"), 
            group_species = bool,
            species_names = parameters$Y_SPECIES,
            save_to = file.path(
                FIGURES_PATH, SUBFOLDER,
                paste0("boxplot_", tolower(loop_on), sp, "_comparison.pdf"))
        ))

        if (is.numeric(combination[[loop_on]]) | (loop_on == "HMSC_XFORMULAS")) {
            . <- suppressMessages(lineplot_model_scores(
                parent_folder = file.path(RESULTS_PATH, SUBFOLDER), 
                loop_model_combination = combination, 
                loop_on = loop_on,
                k_fold = parameters$K_FOLDS, 
                metric = "MSE",
                xlabel = paste0("Effect of ", tolower(loop_on)," type on MSE"), 
                ylabel = "Average score per species",
                group_species = bool,
                species_names = parameters$Y_SPECIES,
                save_to = file.path(
                    FIGURES_PATH, SUBFOLDER,
                    paste0("lineplot_", tolower(loop_on), sp, "_scores.pdf"))
                ))
        } else {
            . <- suppressMessages(dotwhisker_model_scores(
                parent_folder = file.path(RESULTS_PATH, SUBFOLDER), 
                loop_model_combination = combination, 
                loop_on = loop_on,
                k_fold = parameters$K_FOLDS, 
                metric = "MSE",
                xlabel = paste0(
                    "Estimated MSE coefficient relative to '", ref_model[[loop_on]], "'"),
                ylabel = paste0("Effect of ", tolower(loop_on)," type on MSE"),
                group_species = bool,
                species_names = parameters$Y_SPECIES,
                save_to = file.path(
                    FIGURES_PATH, SUBFOLDER,
                    paste0("dotwhisker_", tolower(loop_on), sp, "_scores.pdf")),
            ))
        }

    }

    . <- suppressMessages(boxplot_sp_improvements(
        parent_folder = file.path(RESULTS_PATH, SUBFOLDER),
        reference_model_combination = ref_model, 
        loop_model_combination = loop_model, 
        loop_on = loop_on,
        k_fold = parameters$K_FOLDS, 
        metric = "MSE",
        xlabel = paste0(
                "Effect of ", tolower(loop_on),
                ", compared to '", ref_model[[loop_on]], "'"), 
        subset_names = c("train", "val", "test"),
        proportion = 1/100,
        save_to = file.path(
            FIGURES_PATH, SUBFOLDER,
            paste0("improvement_", tolower(loop_on), sp, "_scores.pdf"))
    ))
}