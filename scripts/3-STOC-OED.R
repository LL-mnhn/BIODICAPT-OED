# This script is used to:
#   - find optimal model to model STOC data
#   - find optimal strategy for model training
##### Libraries ##### ---------------------------------------------------------
suppressPackageStartupMessages(library(geosphere))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(cli))

suppressPackageStartupMessages(source(here::here(file.path("R", "utils_models.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_data.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_figures.R"))))


##### Parameters ##### --------------------------------------------------------
source(here::here(file.path("data","config","config.R"))) # Global parameters

### Seed
set.seed(496) # for reproducible results

### Variables to use in dataset
X_VARIABLES <- c(
    "NDVI", "light_pollution", "p_milieu", "altitude", 
    "precip_spring", "tmp_spring")
X_FACTORS <- c("p_milieu") 
Y_SPECIES <- c("Sylvia_atricapilla", "Parus_major", "Pica_pica", 
    "Carduelis_cannabina", "Periparus_ater")

X_GROUPS_CATS <- c(1, 2, 3, 4, 5, 6) # Group above variables in different categories
X_GROUPS_NAMES <- X_VARIABLES # names associated with X_GROUPS categories

### Data splits
MAX_TRAIN_SIZE <- 125 # Number of point used for training
NEW_POOL_SIZE <- 500 # Number of points from which we can pick for OED
NEW_SAMPLE_SIZE <- 50 # Number of points we can sample in new pool during OED
K_FOLDS <- 5

### Model
HMSC_XFORMULA <- ~ (NDVI + light_pollution + p_milieu + altitude + 
    precip_spring + tmp_spring)

### MCMC
NSAMPLES <- 5000 # mcmc will stop after saving that much samples
THIN <- 2 # number of steps between each recording
NTRANSIENT <- 0.5*NSAMPLES*THIN # burn-in iterations
NCHAINS <- 3


##### Parameters: loops ##### -------------------------------------------------
# When running this scripts, the following list sets a combination of 
# parameters to loop on. This allows to run the model on different combination
# of parameters
COMBINATIONS <- list(
    list(
        TRAIN_SIZES = c(15, 25, 50, 75, 100, 125), # must contain values <= MAX_TRAIN_SIZE
        R_EFFECTS = c("none"), # c("none", "units", "spatial")
        STRATEGIES = c("none") # c("none", "business-as-usual", "gap-filling", "simplified-uncertainty")
    ),
    list(
        TRAIN_SIZES = c(125), # must contain values <= MAX_TRAIN_SIZE
        R_EFFECTS = c("none", "units"), # c("none", "units", "spatial")
        STRATEGIES = c("none") # c("none", "business-as-usual", "gap-filling", "simplified-uncertainty")
    ),
    list(
        TRAIN_SIZES = c(125), # must contain values <= MAX_TRAIN_SIZE
        R_EFFECTS = c("none"), # c("none", "units", "spatial")
        STRATEGIES =  c("business-as-usual", "gap-filling", "simplified-uncertainty") # c("none", "business-as-usual", "gap-filling", "simplified-uncertainty")
    )
)


##### Helper functions #####
PATH_STOC_RESULTS <- file.path(RESULTS_PATH, "STOC-OED")

prepare_necessities <- function() {
    df <- read_csv(STOC_OBS_FULL, show_col_types = FALSE)
    df <- df |>
        mutate(across(all_of(X_FACTORS), as.factor))
    
    current_params <- list(
        X_VARIABLES = X_VARIABLES,
        X_FACTORS = X_FACTORS,
        Y_SPECIES = Y_SPECIES,
        X_GROUPS_CATS = X_GROUPS_CATS,
        X_GROUPS_NAMES = X_GROUPS_NAMES,
        MAX_TRAIN_SIZE = MAX_TRAIN_SIZE,
        NEW_POOL_SIZE = NEW_POOL_SIZE,
        NEW_SAMPLE_SIZE = NEW_SAMPLE_SIZE,
        K_FOLDS = K_FOLDS,
        HMSC_XFORMULA = HMSC_XFORMULA,
        NSAMPLES = NSAMPLES,
        THIN = THIN,
        NCHAINS = NCHAINS,
        NTRANSIENT = NTRANSIENT,
        COMBINATIONS = COMBINATIONS
    )


    remake_files <- function() {
        cli_alert_info("Creating base files...")
        # save parameters and dataset plots
        saveRDS(current_params, file.path(PATH_STOC_RESULTS, "parameters.rds"))
        . <- explore_dataset(df, X_VARIABLES, NULL, PATH_STOC_RESULTS)

        # save data splits for reproducability
        k_fold_points <- split_stoc_points_k_fold_subsets(
            df, MAX_TRAIN_SIZE, NEW_POOL_SIZE, K_FOLDS) 
        saveRDS(k_fold_points, 
                file = file.path(PATH_STOC_RESULTS, "k_fold_points.rds")) 
        return(k_fold_points)
    }

    if (!file.exists(PATH_STOC_RESULTS)) {
        cli_alert_info("Specified output folder not found, creating it...")
        dir.create(PATH_STOC_RESULTS, recursive = TRUE)
        k_fold_points <- remake_files()
    } else if (
        !file.exists(file.path(PATH_STOC_RESULTS, "parameters.rds")) ||
        !file.exists(file.path(PATH_STOC_RESULTS, "k_fold_points.rds")) ) {
        cli_alert_info("Output folder exists but does not contain base files.")
        k_fold_points <- remake_files()
    } else {
        cli_alert_info("Output folder exists, reading parameters...")
        previous_params <- readRDS(
                file.path(PATH_STOC_RESULTS, "parameters.rds"))
        
        if (identical(current_params, previous_params)) {
            # load file
            cli_alert_info("Loading data splits...")
            k_fold_points <- readRDS(
                file.path(PATH_STOC_RESULTS, "k_fold_points.rds"))
        } else if (authorize_overwrite(PATH_STOC_RESULTS)) {
                k_fold_points <- remake_files()
        } else {
            stop("User refused overwriting of local files.")
        }      
    }

    return(list(df, k_fold_points))
}

prepare_subsets <- function(df, current_k, splits) {
    list(
        train = subset(
            df, 
            df$id_point_annee %in% splits$training_points[[current_k]]),
        new_pool = subset(
            df, 
            df$id_point_annee %in% splits$new_pool_points[[current_k]]),
        val_train = subset(
            df, 
            df$id_point_annee %in% splits$val_training_points[[current_k]]),
        val_new_pool = subset(
            df, 
            df$id_point_annee %in% splits$val_new_pool_points[[current_k]]),
        test = subset(
            df, 
            df$id_point_annee %in% splits$test_points[[current_k]]) 
    )
}

extended_training_design <- function(
    strat, base_subset, extension_subset, hM) {
    # dataset changes depending on chosen strat
    if (strat == "none") {
        # No addition of data to base_subet
        extended_set <- base_subset
    } else if (strat == "business-as-usual") {
        # Business-as-usual: random sampling
        cli_alert_info(paste0("Pulling ", NEW_SAMPLE_SIZE, " new random samples..."))
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


##### Main ##### --------------------------------------------------------------
necessities <- prepare_necessities()
stoc_df <- necessities[1][[1]]
k_splits <- necessities[2][[1]]
total_loops <- sum(apply(sapply(COMBINATIONS, lengths), 2, prod)) * K_FOLDS

cli_alert_info("------------ Fitting models ------------\n\n")
id_loop <- 0
for (c in seq_along(COMBINATIONS)) {
    parameters <- COMBINATIONS[[c]]
    
    for (k in seq(K_FOLDS)) {
        subsets <- prepare_subsets(stoc_df, k, k_splits)

        for (m in seq_along(parameters$TRAIN_SIZES)) {
            train_size <- parameters$TRAIN_SIZES[m]
            training_set <- slice_head(subsets$train, n=train_size)

            for (r in seq_along(parameters$R_EFFECTS)) { 
                r_effect <- parameters$R_EFFECTS[r]

                for (s in seq_along(parameters$STRATEGIES)) {
                    strategy <- parameters$STRATEGIES[s]

                    # Display current setup
                    id_loop <- id_loop + 1
                    cli_alert_warning(paste0("----- Running loop: ", id_loop, 
                        "/", total_loops, " -----\n"))    
                    cli_alert_info(paste0("- Random effect: ", r_effect, "\n"))
                    cli_alert_info(paste0("- Strategy: ", strategy, "\n"))
                    cli_alert_info(paste0("- Train size: ", train_size, "\n"))
                    cli_alert_info(paste0("- k-fold: ", k, "\n\n"))

                    # 0. Setup
                    cli_alert_info(("0. Setting up path..."))
                    path_local_results <- file.path(
                        PATH_STOC_RESULTS, paste0(
                            "model_random-", r_effect,
                            "_strategy-", strategy,
                            "_training-size-", train_size,
                            "_k", k
                        )
                    )

                    # pre-check : verify that this combination of variables
                    # does not already exists (quick optimisation for
                    # potential dulicates in COMBINATIONS)
                    if (file.exists(path_local_results)) {
                        cli_alert_info("Model already fitted! Skipping...\n\n")
                        next
                    }
                    
                    # else, continue as always
                    dir.create(path_local_results, recursive = FALSE)
                    path_base_model <- file.path(
                        PATH_STOC_RESULTS, 
                        paste0( "model_random-", r_effect,
                                "_strategy-", "none",
                                "_training-size-", train_size,
                                "_k", k),
                        "chains.rds"
                    ) # only used when strategy is "simplified-uncertainty"
                    training_set <- extended_training_design(
                        strat = parameters$STRATEGIES[s],
                        base_subset = training_set, 
                        extension_subset = subsets$new_pool,
                        hM = readRDS(path_base_model))

                    # 1. Fit model
                    cat("\n")
                    cli_alert_info("1. Fitting model...")
                    base_model <- prepare_hmsc_training(
                        subset = training_set,
                        x_cols = X_VARIABLES, 
                        y_cols = if (is.null(Y_SPECIES)) NAMES_SPECIES else Y_SPECIES,
                        formula = HMSC_XFORMULA,
                        random_effect = r_effect)
                    fitted_model <- fitting_hmsc(
                        hM = base_model, 
                        save_to = file.path(path_local_results, "chains.rds"),
                        nchains = NCHAINS,
                        thin = THIN,
                        nsamples = NSAMPLES,
                        ntransient = NTRANSIENT,
                        freq_verbose = (NSAMPLES*2 + NTRANSIENT)/10,
                        allow_parallel = TRUE)

                    # 2. Analysis of convergence
                    cli_alert_info("2. Convergence diagnostics...")
                    . <- convergence_hmsc(
                        hM = fitted_model, 
                        nchains = NCHAINS, 
                        thin = THIN, 
                        save_folder = path_local_results)
                    
                    # 3. Analysis of performance
                    cli_alert_info("3. Performance evaluation...")
                    # Explanatory power
                    cli_alert_info("Computing training scores...")
                    train_scores <- evaluate_hmsc_performances(
                        hM = fitted_model, subset = training_set, 
                        x_cols = X_VARIABLES, sp_cols = Y_SPECIES)
                    # Prediction power
                    cli_alert_info("Computing testing scores...")
                    if (strategy == "none") {
                        val_scores <- evaluate_hmsc_performances(
                            hM = fitted_model, subset = subsets$val_train, 
                            x_cols = X_VARIABLES, sp_cols = Y_SPECIES)
                    } else {
                        val_scores <- evaluate_hmsc_performances(
                            hM = fitted_model, 
                            subset = bind_rows(
                                subsets$val_train, subsets$val_new_pool),
                            x_cols = X_VARIABLES, sp_cols = Y_SPECIES)
                    }
                    test_scores <- evaluate_hmsc_performances(
                        hM = fitted_model, subset = subsets$test, 
                        x_cols = X_VARIABLES, sp_cols = Y_SPECIES)

                    cli_alert_info("Saving scores...")
                    write_csv(
                        as.data.frame(train_scores), 
                        file.path(path_local_results, "train_scores.csv"))
                    write_csv(
                        as.data.frame(val_scores), 
                        file.path(path_local_results, "val_scores.csv"))
                    write_csv(
                        as.data.frame(test_scores), 
                        file.path(path_local_results, "test_scores.csv"))
                    # row names must be saved seperatly. (dropped for now)
                    cli_alert_success("Scores saved!\n\n")

                    # 4. Result analysis
                    cli_alert_info(".4. Associations...")
                    . <- analyses_hmsc(
                        hM = fitted_model, 
                        save_folder = path_local_results, 
                        x_groups_cats = X_GROUPS_CATS, 
                        x_groups_names = X_GROUPS_NAMES, 
                        supportLevel = 0.05)
                    cli_alert_success("Model ran without errors!\n\n")
                }
            }
        }
    }
}