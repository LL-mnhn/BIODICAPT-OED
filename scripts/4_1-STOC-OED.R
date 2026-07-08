# Using the base model obtained from script 3_1.R, finds OED for STOC dataset
##### Libraries ##### ---------------------------------------------------------
suppressPackageStartupMessages(library(geosphere))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(cli))

suppressMessages(suppressWarnings(source(here::here("R/utils_figures.R"))))
suppressMessages(suppressWarnings(source(here::here("R/utils_models.R"))))


##### Parameters ##### --------------------------------------------------------
source(here::here("data/config/config.R")) # Import global parameters
source(here::here("data/config/config-STOC.R")) # Import local parameters

### Dataset
set.seed(790231) # for reproducible results
STRATEGIES <- c("business-as-usual", "gap-filling", "simplified-uncertainty") # TODO : true-model-uncertainty
R_EFFECT <- "none"


##### Load data and paths ##### -----------------------------------------------
# create folder to save results
if (authorize_overwrite(PATH_LOCAL_OED)) {
    unlink(PATH_LOCAL_OED, recursive=TRUE)
    dir.create(PATH_LOCAL_OED)
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
extended_training_design <- function(
    strat, base_subset, extension_subset, hM) {
    # dataset changes depending on chosen strat
    if (strat == "business-as-usual") {
        # Business-as-usual: random sampling
        cli_alert_info(paste0(s, ".0. Pulling ", NEW_SAMPLE_SIZE, " new random samples..."))
        extended_set <- bind_rows(
            base_subset,
            extension_subset |> slice_sample(n = NEW_SAMPLE_SIZE))
        
    } else if (strat == "gap-filling") {
        # Gap-filling: selecting points that are the farthest from current design
        cli_alert_info(paste0(s, ".0. Searching for the most distant points..."))
        extended_set <- base_subset  # no need for cbind() here at all

        # Selecting one point at a time (re-computing distances after adding each point)
        for (i in 1:NEW_SAMPLE_SIZE) {
            dist_matrix <- distm(
                extended_set[, c("LON", "LAT")],
                extension_subset[, c("LON", "LAT")],
                fun = distGeo)
            idx_highest_min_distance <- order(
                apply(dist_matrix, 2, min), decreasing = TRUE)[1]
            extended_set <- bind_rows(
                extended_set,
                extension_subset[idx_highest_min_distance, ])
        }
    } else if (strat == "simplified-uncertainty") {
        # For what I know, this is the method used in doi.org/10.1111/2041-210X.14355
        # In short: sample where uncertainty is the highest.
        # In this version, we create one "layer" and pick the most uncertain samples at once
        cli_alert_warning("Running simplified-uncertainty, this is an Alpha (see coments)...")

        # TODO: We exclude base_subset but ACTUALLY it might be interesting
        # to resample on already sampled positions to improve performance...
        uncertainty <- get_uncertainty_hmsc(
            hM = hM, 
            df = extension_subset, 
            x_cols = X_VARIABLES) 
        idx_most_uncertain <- order(uncertainty, decreasing = TRUE)[1:NEW_SAMPLE_SIZE]
        extended_set <- bind_rows(
                base_subset,
                extension_subset[idx_most_uncertain, ])
    

    } else if (strat == "true-uncertainty") {
        # In the previous version, we picked the most uncertain points at once.
        # This in a good approximation but it lacks precision: by re-training
        # the model with some new points, the uncertainty of prediction changes
        # so, we must use optimisation algorithms. 

    } else {
        stop(paste0("Unknown strategy: '", strat, "', skipping iteration.\n\n"))
    }
    return(extended_set)
}

fine_compare_hmsc_metric <- function(
        subset_names = c("train", "val", "test"),
        metric = "MSE") {
    # load old results
    cli_alert_info("Loading base results...")
    base_mse <- list()
    k_fold_base_mse <- list()
    for (subset_type in subset_names) {
        base_mse[[subset_type]] <- list()

        for (k in seq(K_FOLD)) {
            path_base_model_results <- file.path(
                    PATH_LOCAL_BASE,
                    paste0("base-model_", R_EFFECT, "-random-effect_k", k))
            
            # load local scores
            if (metric == "MSE") {
                base_mse[[subset_type]][[k]] <- read_csv(
                    file.path(path_base_model_results, paste0(subset_type, "_scores.csv")),
                    show_col_types = FALSE)$RMSE^2
            } else if (metric %in% c("RMSE", "AUC", "TjurR2")) {
                base_mse[[subset_type]][[k]] <- read_csv(
                    file.path(path_base_model_results, paste0(subset_type, "_scores.csv")),
                    show_col_types = FALSE)[[metric]]
            } else {
                stop(paste0("Metric '", metric, "' is not handled by this function."))
            }

        }

        # Get mean per k_fold (individual mean for each subset, model_type and species)
        k_fold_base_mse[[subset_type]] <- colMeans(
            do.call(rbind, base_mse[[subset_type]]), na.rm = TRUE)
    }

    # load new results
    cli_alert_info("Loading OED results...")
    strat_mse <- list()
    k_fold_strat_mse <- list()
    for (subset_type in subset_names) {
        strat_mse[[subset_type]] <- list()
        k_fold_strat_mse[[subset_type]] <- list()

        for (s in seq_along(STRATEGIES)) {
            strat_mse[[subset_type]][[STRATEGIES[s]]] <- list()
            
            for (k in seq(K_FOLD)) {
                path_local_model_results <- file.path(
                    PATH_LOCAL_OED,
                    paste0("model_", STRATEGIES[s], "_random-", R_EFFECT, "_k", k))
                
                # load local scores
                if (metric == "MSE") {
                    strat_mse[[subset_type]][[STRATEGIES[s]]] [[k]] <- read_csv(
                        file.path(path_local_model_results, paste0(subset_type, "_scores.csv")),
                        show_col_types = FALSE)$RMSE^2
                } else if (metric %in% c("RMSE", "AUC", "TjurR2")) {
                    strat_mse[[subset_type]][[STRATEGIES[s]]] [[k]] <- read_csv(
                        file.path(path_local_model_results, paste0(subset_type, "_scores.csv")),
                        show_col_types = FALSE)[[metric]]
                } else {
                    stop(paste0("Metric '", metric, "' is not handled by this function."))
                }
            }

            # Get mean per k_fold (individual mean for each subset, model_type and species)
            k_fold_strat_mse[[subset_type]][[STRATEGIES[s]]] <- colMeans(
                do.call(rbind, strat_mse[[subset_type]][[STRATEGIES[s]]]), na.rm = TRUE)
        }
    }
        
    # format dataframe for ggplot
    results_mse <- tibble()
    for (subset_type in subset_names) {
        for (s in seq_along(STRATEGIES)) {
            if (grepl("MSE", metric)) {
                correct_diff = -1 # MSE is better when its low, so reverse order
            } else {
                correct_diff = 1
            }

            # compute average/sd accross species
            results_mse <- bind_rows(
                results_mse,
                tibble(
                    average_metric = correct_diff * mean(k_fold_strat_mse[[subset_type]][[STRATEGIES[s]]] - k_fold_base_mse[[subset_type]], na.rm = TRUE),
                    sd_metric = sd(k_fold_strat_mse[[subset_type]][[STRATEGIES[s]]] - k_fold_base_mse[[subset_type]]),
                    subset = subset_type,
                    strategy = STRATEGIES[s]))
        }
    }

    # add new fields
    results_mse <- results_mse |>
        mutate(strategy = factor(strategy, levels = STRATEGIES)) |>
        mutate(subset = factor(subset, levels = c("train", "val", "test")))

    p <- ggplot(results_mse, 
        aes(y = average_metric, x = strategy, color = subset)) +
        geom_point(
            stat = "identity", 
            position = position_dodge(width = 0.66), 
            size = 3) +
        geom_errorbar(
            aes(ymin = average_metric - sd_metric, ymax = average_metric + sd_metric),
            position = position_dodge(width = 0.66),
            width = 0.2) +
        labs(
            caption = "SD and mean computed per k_fold (over all species)",
            color = "Subset") +
        ylab(paste("Delta in average", metric)) +
        xlab("Strategy for new samples") +
        geom_hline(yintercept = 0, linetype = "dashed")

    p <- my_custom_ggplot_theme(p) + 
        scale_color_manual(values = c(PALETTE[2], PALETTE[3], PALETTE[1]))
    print(p)
    standardised_ggplot_save(
        figure = p, 
        save_path = file.path(PATH_LOCAL_OED, paste0("compared_OED_", metric,".pdf")))
    cli_alert_success("Plot of compared performances is saved!\n\n")
    
    return(results_mse)
}


##### Fitting OED ##### -------------------------------------------------------
# Goal: train a model on extended subsets (compared to base models)
# find the subset that maximises certainty of predictions
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
    
    # Try 3 sampling strategies:
    for (s in seq_along(STRATEGIES)) {
        # create folder for saving results
        path_local_model_results <- file.path(
            PATH_LOCAL_OED,
            paste0("model_", STRATEGIES[s], "_random-", R_EFFECT, "_k", k))
        dir.create(path_local_model_results, recursive = TRUE)

        # previous model fitted on this dataset
        path_previous_model <- file.path(
            PATH_LOCAL_BASE,
            paste0("base-model_", R_EFFECT, "-random-effect_k", k),
            "train_outputs.rds")
        # TODO - 3. Model uncertainty: sample locations that have the model is most uncertain about

        cli_alert_warning(paste0(
            "---- ", s, ". Strategy: '", STRATEGIES[s], "' ----\n\n"))
        extended_training_set <- extended_training_design(
            strat = STRATEGIES[s],
            base_subset = train_subset, 
            extension_subset = new_pool_subset,
            hM = readRDS(path_previous_model))
        # save names of training points (faster than re-running 
        # `extended_training_design` if investigation of results is needed)
        saveRDS(
            extended_training_set$id_point_annee, 
            file = file.path(path_local_model_results, "training_points.rds"))

        # 1. Fit model
        cli_alert_info(paste0(s, ".1. Fitting model..."))
        base_model <- prepare_hmsc_training(
            subset = extended_training_set,
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
            hM = fitted_model, subset = extended_training_set, 
            x_cols = X_VARIABLES, sp_cols = Y_SPECIES)

        # Prediction power
        cli_alert_info("Computing testing scores...")
        val_scores <- evaluate_hmsc_performances(
            hM = fitted_model, subset = rbind(val_train_subset, val_new_pool_subset), 
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
        parent_folder = PATH_LOCAL_OED, 
        prefix = "model_", 
        model_types = STRATEGIES, 
        sufix = paste0("_random-", R_EFFECT,"_k"), 
        k_fold = K_FOLD, 
        xlabel = "Sampling strategy")
. <- fine_compare_hmsc_metric(metric = "MSE")
