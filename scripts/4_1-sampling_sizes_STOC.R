# Using the base model obtained from script 3_1.R, test different
# training sizes to assess HMSC turning point in training 
# (the moments where it becomes "good")
##### Libraries ##### ---------------------------------------------------------
suppressPackageStartupMessages(library(geosphere))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(cli))

suppressPackageStartupMessages(source(here::here(file.path("R", "utils_models.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_figures.R"))))


##### Parameters ##### --------------------------------------------------------
source(here::here(file.path("data","config","config.R"))) # Global parameters
source(here::here(file.path("data","config","config-STOC.R"))) # Local params

### Dataset
set.seed(790231) # for reproducible results
STRATEGIES <- c("business-as-usual", "gap-filling", "simplified-uncertainty") # TODO : true-model-uncertainty
R_EFFECT <- "none"

cli_alert_warning("OVERWRITING DEFAULT PARAMETERS")
TRAIN_SIZES <- c(5, 10, 18, 25, 37, 50, 66, 75, 100, 125)


##### Load data and paths ##### -----------------------------------------------
# create folder to save results
if (authorize_overwrite(PATH_LOCAL_SIZE)) {
    unlink(PATH_LOCAL_SIZE, recursive=TRUE)
    dir.create(PATH_LOCAL_SIZE)
} else {
    stop("Cannot proceed: output folder cannot be written.")
}

# load dataset
stoc_df <- read_csv(STOC_OBS_FULL, show_col_types = FALSE)
stoc_df <- stoc_df |>
    mutate(across(all_of(X_FACTORS), as.factor))

# load train/val/test subsets
k_splits <- readRDS(file.path(PATH_LOCAL_RESULTS, "k_fold_points.rds"))


###### Helper functions ##### -------------------------------------------------
# none


##### Fitting OED ##### -------------------------------------------------------
# Goal: train a model on different training sizes, find when the accuracy goes up
cli_alert_info("------------ Fitting base models ------------\n\n")
for (k in seq(K_FOLD)) {
    cli_alert_warning(paste0("------ k-fold iteration: ", k, "/", K_FOLD, " ------\n\n"))

    # get subsets
    train_subset <- subset(
        stoc_df, 
        stoc_df$id_point_annee %in% k_splits$training_points[[k]])
    new_pool_subset <- subset(
        stoc_df, 
        stoc_df$id_point_annee %in% k_splits$new_pool_points[[k]])
    val_train_subset <- subset(
        stoc_df, 
        stoc_df$id_point_annee %in% k_splits$val_training_points[[k]])
    val_new_pool_subset <- subset(
        stoc_df, 
        stoc_df$id_point_annee %in% k_splits$val_new_pool_points[[k]])
    test_subset <- subset(
        stoc_df, 
        stoc_df$id_point_annee %in% k_splits$test_points[[k]])
    
    # Try all sampling sizes
    for (s in seq_along(TRAIN_SIZES)) {
        cli_alert_warning(paste0("------ Training size: ", s, "/", length(TRAIN_SIZES), " ------\n\n"))
        # create folder for saving results
        path_local_model_results <- file.path(
            PATH_LOCAL_SIZE,
            paste0("model_train-size-", TRAIN_SIZES[s], "_random-", R_EFFECT, "_k", k))
        dir.create(path_local_model_results, recursive = TRUE)

        # previous model fitted on this dataset
        local_training_subset <- slice_head(train_subset, n=TRAIN_SIZES[s])

        # 1. Fit model
        cli_alert_info(paste0(s, ".1. Fitting model..."))
        base_model <- prepare_hmsc_training(
            subset = local_training_subset,
            x_cols = X_VARIABLES, 
            y_cols = Y_SPECIES,
            formula = HMSC_XFORMULA,
            random_effect = R_EFFECT)
        fitted_model <- fitting_hmsc(
            hM = base_model, 
            save_to = file.path(path_local_model_results, "train_outputs.rds"),
            nchains = NCHAINS,
            thin = THIN,
            nsamples = NSAMPLES,
            ntransient = NTRANSIENT,
            freq_verbose = (NSAMPLES*2 + NTRANSIENT)/10,
            allow_parallel = TRUE)
        
        # # 2. Analysis of convergence 
        # # Convergence diagnostics are not necessary since "base-model" 
        # # converged easily. However, it can be useful to check results...
        # # Uncomment the following section for details
        # cli_alert_info(paste0(s, ".2. Convergence diagnostics..."))
        # . <- convergence_hmsc(
        #     hM = fitted_model, 
        #     nchains = NCHAINS, 
        #     thin = THIN, 
        #     save_folder = path_local_model_results)
        # cat("\n")
        cli_alert_warning(paste0(s, ".2. Skipping convergence diagnostics (see comment in code)."))

        # 3. Analysis of performance
        cli_alert_info(paste0(s, ".3. Performance evaluation..."))
        
        # Explanatory power
        cli_alert_info("Computing training scores...")
        train_scores <- evaluate_hmsc_performances(
            hM = fitted_model, subset = local_training_subset, 
            x_cols = X_VARIABLES, sp_cols = Y_SPECIES)

        # Prediction power
        cli_alert_info("Computing testing scores...")
        val_scores <- evaluate_hmsc_performances(
            hM = fitted_model, subset = rbind(val_train_subset), 
            x_cols = X_VARIABLES, sp_cols = Y_SPECIES)
        test_scores <- evaluate_hmsc_performances(
            hM = fitted_model, subset = test_subset, 
            x_cols = X_VARIABLES, sp_cols = Y_SPECIES)

        cli_alert_info("Saving scores...")
        write_csv(
            as.data.frame(train_scores), 
            file.path(path_local_model_results, "train_scores.csv"))
        write_csv(
            as.data.frame(val_scores), 
            file.path(path_local_model_results, "val_scores.csv"))
        write_csv(
            as.data.frame(test_scores), 
            file.path(path_local_model_results, "test_scores.csv"))
        # row names must be saved seperatly. (dropped)
        cli_alert_success("Scores saved!\n\n")
    }
}


##### Compare performances ##### ----------------------------------------------
cli_alert_info("------------ Results ------------\n\n")
. <- compute_hmsc_performances(
        parent_folder = PATH_LOCAL_SIZE, 
        filename = "training-size_performances.pdf",
        prefix = "model_train-size-", 
        model_types = TRAIN_SIZES, 
        sufix = paste0("_random-", R_EFFECT,"_k"), 
        k_fold = K_FOLD, 
        xlabel = "Effect of sampling size on metrics",
        ylabel = "Average score per species",
        group_species = FALSE,
        barplot = FALSE)
