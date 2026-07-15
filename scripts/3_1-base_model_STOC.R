# This script is used to find the "base" model that will be used to 
# compute OED on STOC data
##### Libraries ##### ---------------------------------------------------------
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(cli))

suppressPackageStartupMessages(source(here::here(file.path("R", "utils_models.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_data.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_figures.R"))))


##### Parameters ##### --------------------------------------------------------
source(here::here(file.path("data","config","config.R"))) # Global parameters
source(here::here(file.path("data","config","config-STOC.R"))) # Local params

### Dataset
set.seed(496) # for reproducible results


##### Load data and paths ##### -----------------------------------------------
# create folder to save results
if (authorize_overwrite(PATH_LOCAL_RESULTS)) {
    unlink(PATH_LOCAL_RESULTS, recursive=TRUE)
    dir.create(PATH_LOCAL_RESULTS)
} else {
    stop("Cannot proceed: output folder cannot be written.")
}

# load dataset
stoc_df <- check_dataset(STOC_OBS_FULL)
stoc_df <- stoc_df |>
    mutate(across(all_of(X_FACTORS), as.factor))
. <- explore_dataset(stoc_df, X_VARIABLES, NULL, PATH_LOCAL_RESULTS)

# split train/test subsets 'k_fold' times
k_fold_points <- split_stoc_points_k_fold_subsets(
    stoc_df, TRAIN_SIZE, NEW_POOL_SIZE, K_FOLD) 
saveRDS(k_fold_points, file=file.path(PATH_LOCAL_RESULTS, "k_fold_points.rds"))


##### Fitting base models ##### ------------------------------------------------
# Goal: compare models and find the one that best describes the datasets given
cli_alert_info("------------ Fitting base models ------------\n\n")
for (k in seq(K_FOLD)) {
    cli_alert_warning(paste0("---- k-fold iteration: ", k, "/", K_FOLD, "----\n\n"))

    ### prepare subsets
    train_subset <- subset(
        stoc_df, 
        stoc_df$id_point_annee %in% k_fold_points$training_points[[k]])
    val_subset <- subset(
        stoc_df, 
        stoc_df$id_point_annee %in% k_fold_points$val_training_points[[k]])
    test_subset <- subset(
        stoc_df, 
        stoc_df$id_point_annee %in% k_fold_points$test_points[[k]])

    ### test model with and without random effects
    r_effects <- c("none", "units") # "spatial" is available but way longer to run
    for (r in seq_along(r_effects)) { 
        # 0. Setup
        r_effect <- r_effects[r]
        path_local_model_results <- file.path(
            PATH_LOCAL_BASE,
            paste0("base-model_", r_effect, "-random-effect_k", k))
        dir.create(path_local_model_results, recursive = TRUE)
        cli_alert_info(paste0(r, ".0. Setting up model."))
        cli_alert_info(paste0("Random effects in model: '", r_effect, "'."))

        # 1. Fit model
        base_model <- prepare_hmsc_training(
            subset = train_subset,
            x_cols = X_VARIABLES, 
            y_cols = if (is.null(Y_SPECIES)) NAMES_SPECIES else Y_SPECIES,
            formula = HMSC_XFORMULA,
            random_effect = r_effect)
        cli_alert_info(paste0(r, ".1. Fitting model..."))
        fitted_model <- fitting_hmsc(
            hM = base_model, 
            save_to = file.path(path_local_model_results, "train_outputs.rds"),
            nchains = NCHAINS,
            thin = THIN,
            nsamples = NSAMPLES,
            ntransient = NTRANSIENT,
            freq_verbose = (NSAMPLES*2 + NTRANSIENT)/10,
            allow_parallel = TRUE)
        
        # 2. Analysis of convergence
        cli_alert_info(paste0(r, ".2. Convergence diagnostics..."))
        . <- convergence_hmsc(
            hM = fitted_model, 
            nchains = NCHAINS, 
            thin = THIN, 
            save_folder = path_local_model_results)
        cat("\n")

        # 3. Analysis of performance
        cli_alert_info(paste0(r, ".3. Performance evaluation..."))
        
        # Explanatory power
        cli_alert_info("Computing training scores...")
        train_scores <- evaluate_hmsc_performances(
            hM = fitted_model, subset = train_subset, 
            x_cols = X_VARIABLES, sp_cols = Y_SPECIES)

        # Prediction power
        cli_alert_info("Computing testing scores...")
        val_scores <- evaluate_hmsc_performances(
            hM = fitted_model, subset = val_subset, 
            x_cols = X_VARIABLES, sp_cols = Y_SPECIES)
        test_scores <- evaluate_hmsc_performances(
            hM = fitted_model, subset = test_subset, 
            x_cols = X_VARIABLES, sp_cols = Y_SPECIES)

        cli_alert_info("Saving scores...")
        write_csv(as.data.frame(train_scores), file.path(path_local_model_results, "train_scores.csv"))
        write_csv(as.data.frame(val_scores), file.path(path_local_model_results, "val_scores.csv"))
        write_csv(as.data.frame(test_scores), file.path(path_local_model_results, "test_scores.csv"))
        # row names must be saved seperatly. (dropped)
        cli_alert_success("Scores saved!\n\n")

        # 4. Result analysis
        cli_alert_info(paste0(r, ".4. Associations..."))
        . <- analyses_hmsc(
            hM = fitted_model, 
            save_folder = path_local_model_results, 
            x_groups_cats = X_GROUPS_CATS, 
            x_groups_names = X_GROUPS_NAMES, 
            supportLevel = 0.05)
        cli_alert_warning("TODO: make maps!!!")
        cli_alert_success("Model ran without errors!\n\n")
    }
}


##### Compare performances ##### ----------------------------------------------
cli_alert_info("------------ Results ------------\n\n")
agg_scores_df <- compute_hmsc_performances(
        parent_folder = PATH_LOCAL_BASE, 
        filename = "base-model_performances.pdf",
        prefix = "base-model_", 
        model_types = r_effects, 
        sufix = "-random-effect_k", 
        k_fold = K_FOLD, 
        xlabel = "Effect of model type on metrics", 
        ylabel = "Average score per species",
        group_species = FALSE)

cli_alert_info(
    "Given the results, there is little to no difference in prediction power ",
    "on new (test) data.")
cli_alert("For OED, we'll therefore use a model without random effects.")
cli_alert("Training base model on whole dataset...")


##### Predict on map of france ##### ------------------------------------------
cli_alert_info("------------ Distribution map (example) ------------\n\n")
for (k in seq(K_FOLD)) {
    a <- map_results_hmsc(
        df = stoc_df,
        x_cols = X_VARIABLES,
        parent_folder = PATH_LOCAL_BASE, 
        prefix = "base-model_", 
        model_type = "none", 
        sufix = "-random-effect_k", 
        k_fold = k, 
        sp = "Sylvia_atricapilla")
}
