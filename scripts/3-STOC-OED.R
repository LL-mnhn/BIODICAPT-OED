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

### Paths
PATH_STOC_RESULTS <- file.path(RESULTS_PATH, "STOC-OED")
PATH_STOC_FIGURES <- file.path(FIGURES_PATH, "STOC-OED")

### Variables
X_AVAILABLES <- c(
    "NDVI", "light_pollution", "p_milieu", 
    "altitude", "precip_spring", "tmp_spring")
X_FACTORS <- c("p_milieu") # ensures the right format for X_AVAILABLES
Y_SPECIES <- c(
    # 5 most sighted
    "Sylvia_atricapilla",
    "Fringilla_coelebs",
    "Turdus_merula",
    "Parus_major",
    "Columba_palumbus",
    # 5 midly sigthed
    "Dendrocopos_major",
    "Carduelis_chloris",
    "Pica_pica",
    "Phoenicurus_ochruros",
    "Sylvia_communis",
    # 5 least sighted
    "Carduelis_cannabina",
    "Serinus_serinus",
    "Turdus_viscivorus",
    "Anthus_trivialis",
    "Periparus_ater")

### Data
MAX_TRAIN_SIZE <- 125 # Number of point used for training
NEW_POOL_SIZE <- 500 # Number of points from which we can pick for OED
K_FOLDS <- 10

### MCMC
NSAMPLES <- 1500 # mcmc will stop after saving that much samples
THIN <- 1 # number of steps between each recording
NTRANSIENT <- 0.5*NSAMPLES*THIN # burn-in iterations
NCHAINS <- 3


##### Parameters: loops ##### -------------------------------------------------
# Each list is a combination of parameters to use to train the model
COMBINATIONS <- list(
    # Change the number of training samples
    list(
        TRAIN_SIZES = c(15, 25, 50, 75, 100, 125), # values <= MAX_TRAIN_SIZE
        R_EFFECTS = c("none"), # c("none", "units", "spatial")
        STRATEGIES = c("none"), # c("none", "business-as-usual", "gap-filling", "simplified-uncertainty")
        NEW_SAMPLE_SIZES = c(0),
        HMSC_XFORMULAS = c(
            ~ (NDVI + light_pollution + p_milieu + altitude + precip_spring + 
                tmp_spring)
        )
    ),
    # Change the number of explicative variables
    list(
        TRAIN_SIZES = c(125), 
        R_EFFECTS = c("none"), 
        STRATEGIES = c("none"), 
        NEW_SAMPLE_SIZES = c(0),
        HMSC_XFORMULAS = c(
            ~ (NDVI + light_pollution + p_milieu + altitude + precip_spring + 
                tmp_spring),
            ~ (p_milieu + light_pollution + altitude),
            ~ (p_milieu + light_pollution))
    ),
    # Add random effects
    list(
        TRAIN_SIZES = c(125), 
        R_EFFECTS = c("none", "units", "spatial"), 
        STRATEGIES = c("none"),
        NEW_SAMPLE_SIZES = c(0),
        HMSC_XFORMULAS = c(
            ~ (NDVI + light_pollution + p_milieu + altitude + precip_spring + 
                tmp_spring)
        )
    ),
    # Add new samples through different strategies
    list(
        TRAIN_SIZES = c(125), 
        R_EFFECTS = c("none"), 
        NEW_SAMPLE_SIZES = c(50),
        STRATEGIES = c("none", "business-as-usual", "gap-filling", 
            "simplified-uncertainty"),
        HMSC_XFORMULAS = c(
            ~ (NDVI + light_pollution + p_milieu + altitude + precip_spring + 
                tmp_spring)
        )
    ),
    # Add a different number of new samples
    list(
        TRAIN_SIZES = c(125), 
        R_EFFECTS = c("none"), 
        NEW_SAMPLE_SIZES = c(0, 10, 25, 50, 100, 150),
        STRATEGIES =  c("simplified-uncertainty"),
        HMSC_XFORMULAS = c(
            ~ (NDVI + light_pollution + p_milieu + altitude + precip_spring + 
                tmp_spring)
        )
    )
)

BASE_COMBINATION <- list( # must be one of COMBINATIONS
    # Change the number of training samples
    TRAIN_SIZES = c(125), 
    R_EFFECTS = c("none"), 
    STRATEGIES = c("none", "simplified-uncertainty"), # "none" is base, "simplified-uncertainty" is when looping on NEW_SAMPLE_SIZES
    NEW_SAMPLE_SIZES = c(0, 50), # 0 is base, 50 is when looping on STRATEGIES
    HMSC_XFORMULAS = c(
        ~ (NDVI + light_pollution + p_milieu + altitude + precip_spring + 
            tmp_spring)
    )
)

##### Helper functions ##### --------------------------------------------------
prepare_necessities <- function() {
    df <- read_csv(STOC_OBS_FULL, show_col_types = FALSE)
    df <- df |>
        mutate(across(all_of(X_FACTORS), as.factor))
    
    current_params <- list(
        X_AVAILABLES = X_AVAILABLES,
        X_FACTORS = X_FACTORS,
        Y_SPECIES = Y_SPECIES,
        MAX_TRAIN_SIZE = MAX_TRAIN_SIZE,
        NEW_POOL_SIZE = NEW_POOL_SIZE,
        K_FOLDS = K_FOLDS,
        NSAMPLES = NSAMPLES,
        THIN = THIN,
        NTRANSIENT = NTRANSIENT,
        NCHAINS = NCHAINS,
        COMBINATIONS = COMBINATIONS,
        BASE_COMBINATION = BASE_COMBINATION
    )

    remake_files <- function() {
        cli_alert_info("Creating base files...")
        # save parameters and dataset plots
        saveRDS(current_params, file.path(PATH_STOC_RESULTS, "parameters.rds"))
        . <- explore_dataset(df, X_AVAILABLES, NULL, PATH_STOC_FIGURES)

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
        dir.create(PATH_STOC_FIGURES, recursive = TRUE)
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
            cli_alert_info("Parameters are identical, loading data splits...")
            k_fold_points <- readRDS(
                file.path(PATH_STOC_RESULTS, "k_fold_points.rds"))
        } else if (authorize_overwrite(PATH_STOC_RESULTS)) {
                cli_alert_info("Parameters found are different from current ones..")
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
    strat, base_subset, extension_subset, new_sample_size, variables, hM) {
    # dataset changes depending on chosen strat
    if ((strat == "none") || (new_sample_size == 0)) {
        # No addition of data to base_subet
        extended_set <- base_subset
    } else if (strat == "business-as-usual") {
        # Business-as-usual: random sampling
        cli_alert_info(
            paste0("Pulling ", new_sample_size, " new random samples..."))
        extended_set <- bind_rows(
            base_subset,
            extension_subset |> slice_sample(n = new_sample_size))
        
    } else if (strat == "gap-filling") {
        # Gap-filling: select points that are the farthest from current design
        cli_alert_info(paste0("Searching for the most distant points..."))
        extended_set <- base_subset  # no need for cbind() here at all

        # Select one point at a time 
        # (re-compute distances each time a point is added)
        for (i in 1:new_sample_size) {
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
        # This is the method used in doi.org/10.1111/2041-210X.14355
        # In short: sample where uncertainty is the highest.
        # Here, we create a "layer" and pick the most uncertain samples at once
        cli_alert_warning("Running simplified-uncertainty...")

        # TODO: We exclude base_subset but ACTUALLY it might be interesting
        # to resample on already sampled positions to improve performance...
        uncertainty <- get_uncertainty_hmsc(
            hM = hM, 
            df = extension_subset, 
            x_cols = variables) 
        idx_most_uncertain <- order(
            uncertainty, decreasing = TRUE)[1:new_sample_size]
        extended_set <- bind_rows(
                base_subset,
                extension_subset[idx_most_uncertain, ])
    
    } else if (strat == "true-uncertainty") {
        # In the previous version, we picked the most uncertain points at once.
        # This in a good approximation but it lacks precision: by re-training
        # the model with some new points, the uncertainty of prediction changes
        # so, we must use optimisation algorithms. 
        stop("'true-uncertainty' method is not available yet.")

    } else {
        stop(paste0("Unknown strategy: '", strat, "'."))
    }
    return(extended_set)
}


##### Main ##### --------------------------------------------------------------
param_grid <- build_param_grid(COMBINATIONS)
total_loops <- nrow(param_grid) * K_FOLDS

necessities <- prepare_necessities()
stoc_df <- necessities[1][[1]]
k_splits <- necessities[2][[1]]


cli_alert_info("------------ Fitting models ------------\n\n")
id_loop <- 0
for (k in seq(K_FOLDS)) {
    subsets <- prepare_subsets(stoc_df, k, k_splits)

    for (p in seq_len(nrow(param_grid))) {
        # get parameters
        train_size <- param_grid$train_size[p]
        r_effect <- param_grid$r_effect[p]
        strategy <- param_grid$strategy[p]
        n_new_samples <- param_grid$n_new_samples[p]
        formula <- param_grid$formulas[[p]]

        training_set <- slice_head(subsets$train, n = train_size)
        
        # from formula, get x_variables
        x_variables <- all.vars(formula)
        x_groups_cats <- seq(1:length(x_variables)) # For variance proportion       

        # if strategy is none, overwrite n_new_samples to avoid errors
        if (strategy == "none") {
            n_new_samples <- 0
        }
        # if n_new_samples is 0, no need to apply a strategy, use "none"
        if (n_new_samples == 0) {
            strategy <- "none"
        }

        # Display current setup
        id_loop <- id_loop + 1
        cli_alert_warning(paste0("----- Running loop: ", id_loop, 
            "/", total_loops, " -----\n"))    
        cli_alert_info(paste0("- Random effect: ", r_effect, "\n"))
        cli_alert_info(paste0("- Strategy: ", strategy, "\n"))
        cli_alert_info(paste0("- New samples: ", n_new_samples, "\n"))
        cli_alert_info(paste0("- Train size: ", train_size, "\n"))
        cli_alert_info(paste0("- N variables: ", length(x_variables), "\n"))
        cli_alert_info(paste0("- K-fold: ", k, "\n\n"))

        # 0. Setup
        cli_alert_info(("0. Setting up path..."))
        path_local_results <- file.path(
            PATH_STOC_RESULTS, paste0(
                "model_random-", r_effect,
                "_strategy-", strategy,
                "_new-samples-", n_new_samples,
                "_training-size-", train_size,
                "_n-variables-", length(x_variables),
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
        # only used when strategy is "simplified-uncertainty"
        path_base_model <- file.path(
            PATH_STOC_RESULTS, 
            paste0(
                "model_random-", BASE_COMBINATION$R_EFFECTS,
                "_strategy-", BASE_COMBINATION$STRATEGIES[1],
                "_new-samples-", BASE_COMBINATION$NEW_SAMPLE_SIZES[1],
                "_training-size-", BASE_COMBINATION$TRAIN_SIZES,
                "_n-variables-", 
                length(all.vars(BASE_COMBINATION$HMSC_XFORMULAS[[1]])),
                "_k", k),
            "chains.rds"
        )
        if ((!file.exists(path_base_model)) && (strategy != "none")) {
            stop(paste0(
                path_base_model, ", path not found.",
                " Did you call a similar model with strategy='none' first?"))
        }

        # if strategy is X_parts, divide the number of samples and fit as many
        # times as necessary
        pattern_n_parts <- "([^_]+)-part-"
        regex_parts <- regmatches(strategy, regexec(pattern_n_parts, strategy))
        n_parts <- ifelse(
            is.na(regex_parts[[1]][2]), 1, as.numeric(regex_parts[[1]][2]))
        strategy_alone <- substring(
            strategy, 
            ifelse((n_parts == 1), 1, nchar(regex_parts[[1]][1])+1), 
            nchar(strategy))
        path_previous_model <- path_base_model
        n_new_samples_per_part <- split_evenly(n_new_samples, n_parts)
        
        # If loop runs n times (n>1), the model used to compute the second part
        # of the dataset becomes the model computed on part n-1.
        # This *could* be a problem since `training_set` is updated but not
        # `subsets$new_pool`. (samples in new_pool could be sampled n times).
        for (part in 1:n_parts) {
            training_set <- extended_training_design(
                strat = strategy_alone,
                base_subset = training_set, 
                extension_subset = subsets$new_pool,
                new_sample_size = n_new_samples_per_part[part],
                variables = x_variables,
                hM = readRDS(path_previous_model))

            # 1. Fit model
            cat("\n")
            cli_alert_info("1. Fitting model...")
            base_model <- prepare_hmsc_training(
                subset = training_set,
                x_cols = x_variables, 
                y_cols = if (is.null(Y_SPECIES)) NAMES_SPECIES else Y_SPECIES,
                formula = formula,
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
            
            # previous model becomes the one we just computed
            path_previous_model <- file.path(path_local_results, "chains.rds")
        }


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
            x_cols = x_variables, sp_cols = Y_SPECIES)
        # Prediction power
        cli_alert_info("Computing testing scores...")
        if (strategy == "none") {
            val_scores <- evaluate_hmsc_performances(
                hM = fitted_model, subset = subsets$val_train, 
                x_cols = x_variables, sp_cols = Y_SPECIES)
        } else {
            val_scores <- evaluate_hmsc_performances(
                hM = fitted_model, 
                subset = bind_rows(
                    subsets$val_train, subsets$val_new_pool),
                x_cols = x_variables, sp_cols = Y_SPECIES)
        }
        test_scores <- evaluate_hmsc_performances(
            hM = fitted_model, subset = subsets$test, 
            x_cols = x_variables, sp_cols = Y_SPECIES)

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
        cli_alert_info("4. Associations...")
        . <- analyses_hmsc(
            hM = fitted_model, 
            save_folder = path_local_results, 
            x_groups_cats = x_groups_cats, 
            x_groups_names = x_variables, 
            supportLevel = 0.05)
        cli_alert_success("Model ran without errors!\n\n")
    }
}