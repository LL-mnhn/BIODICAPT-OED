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
param_grid <- build_param_grid(COMBINATIONS)


##### Helper functions ##### --------------------------------------------------
cli_alert_warning(
    "[WARNING]: Reference values are *hard coded* in 'loop_selector'.")
loop_selector <- function(
        grid, loop_on, base = parameters$BASE_COMBINATION){
    # initialize base values
    combination <- list(
        train_size = base$TRAIN_SIZES,
        r_effect = base$R_EFFECTS,
        strategy = base$STRATEGIES[1],
        formulas = base$HMSC_XFORMULAS,
        n_new_samples = base$NEW_SAMPLE_SIZE[1]
    )

    # replace formulas by number of variables
    combination$formulas <- lapply(
            lapply(combination$formulas, all.vars), 
        length)
    
    # if looping on strategies, change base new_sample_size
    if (loop_on == "strategy") {
        combination$n_new_samples <- base$NEW_SAMPLE_SIZE[2]
    }

    # if looping on new_samples, change base strategy
    if (loop_on == "n_new_samples") {
        combination$strategy <- base$STRATEGIES[2]
    }

    # replace single value by list to loop on
    combination[[loop_on]] <- unique(param_grid[[loop_on]])

    # a preffix, and a suffix in two parts
    path_sections <- list("", "")
    path_idx <- 1
    parameter_names <- c(
        "r_effect", "strategy", "n_new_samples", "train_size", "formulas")
    path_names <- c(
        "model_random-", "_strategy-", "_new-samples-", "_training-size-", 
        "_n-variables-")

    # when reaching vector of item in list, switch to next path section
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
    path_sections[[2]] <- paste0(path_sections[[2]], "_k")

    return(path_sections)
}


##### Main ##### --------------------------------------------------------------
cli_alert_info("------------ Score comparison ------------\n\n")
for (column in names(param_grid)) {
    cli_alert_info(paste0("Plotting ", column, "..."))
    filename_sections <- loop_selector(param_grid, column)
    loop_list <- unique(param_grid[[column]])

    # using names would be too long, we use number of x variables in formula
    if (column == "formulas") {
        loop_list <- sapply(lapply(loop_list, all.vars), length)
    }

    # if list is numerical, sort it
    if (any(sapply(loop_list, is.numeric))) {
        loop_list <- sort(loop_list)
    }

    for (bool in c(FALSE, TRUE)) {
        if (bool == FALSE) {
            sp <- "_per_species"
        } else {
            sp <- ""
        }

        . <- suppressMessages(boxplot_compare_scores(
            parent_folder = file.path(RESULTS_PATH, SUBFOLDER), 
            reference_model_path = paste0(
                filename_sections[[1]], loop_list[1], filename_sections[[2]]), 
            loop_prefix = filename_sections[[1]], 
            loop_elements = loop_list[-1], 
            loop_suffix = filename_sections[[2]], 
            k_fold = parameters$K_FOLDS, 
            metric = "MSE",
            subset_names = c("train", "val", "test"),
            xlabel = paste0(
                "Effect of ", tolower(column),
                ", compared to '", loop_list[1], "'"), 
            group_species = bool,
            species_names = parameters$Y_SPECIES,
            save_to = file.path(
                FIGURES_PATH, SUBFOLDER,
                paste0("boxplot_", tolower(column), sp, "_comparison.pdf"))
        ))

        if (is.numeric(loop_list)) {
            . <- suppressMessages(lineplot_model_scores(
                parent_folder = file.path(RESULTS_PATH, SUBFOLDER), 
                loop_prefix = filename_sections[[1]], 
                loop_elements = if (column != "strategy") {
                        loop_list
                    } else {
                        loop_list[-1]
                    }, 
                loop_suffix = filename_sections[[2]], 
                k_fold = parameters$K_FOLDS, 
                metric = "MSE",
                xlabel = paste0("Effect of ", tolower(column)," type on MSE"), 
                ylabel = "Average score per species",
                group_species = bool,
                species_names = parameters$Y_SPECIES,
                save_to = file.path(
                    FIGURES_PATH, SUBFOLDER,
                    paste0("lineplot_", tolower(column), sp, "_scores.pdf"))
                ))
        } else {
            . <- suppressMessages(dotwhisker_model_scores(
                parent_folder = file.path(RESULTS_PATH, SUBFOLDER), 
                loop_prefix = filename_sections[[1]], 
                loop_elements = loop_list, 
                loop_suffix = filename_sections[[2]], 
                k_fold = parameters$K_FOLDS, 
                metric = "MSE",
                xlabel = paste0(
                    "Estimated MSE coefficient relative to '", loop_list[1], "'"),
                ylabel = paste0("Effect of ", tolower(column)," type on MSE"),
                group_species = bool,
                species_names = parameters$Y_SPECIES,
                save_to = file.path(
                    FIGURES_PATH, SUBFOLDER,
                    paste0("dotwhisker_", tolower(column), sp, "_scores.pdf")),
            ))
        }

    }

    . <- suppressMessages(boxplot_sp_improvements(
        parent_folder = file.path(RESULTS_PATH, SUBFOLDER),
        reference_model_path = paste0(
                filename_sections[[1]], loop_list[1], filename_sections[[2]]), 
        loop_prefix = filename_sections[[1]], 
        loop_elements = loop_list[-1], 
        loop_suffix = filename_sections[[2]], 
        k_fold = parameters$K_FOLDS, 
        metric = "MSE",
        xlabel = paste0(
                "Effect of ", tolower(column),
                ", compared to '", loop_list[1], "'"), 
        subset_names = c("train", "val", "test"),
        proportion = 1/100,
        save_to = file.path(
            FIGURES_PATH, SUBFOLDER,
            paste0("improvement_", tolower(column), sp, "_scores.pdf"))
    ))
}