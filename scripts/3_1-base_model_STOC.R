# This script is used to find the "base" model that will be used to 
# compute OED on STOC data
##### Libraries ##### ---------------------------------------------------------
suppressPackageStartupMessages(library(Hmsc))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(cli))
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(MCMCvis))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(abind))

suppressMessages(suppressWarnings(source(here::here("R/utils_models.R"))))
suppressMessages(suppressWarnings(source(here::here("R/utils_figures.R"))))
suppressMessages(suppressWarnings(source(here::here("R/utils_data.R"))))



##### Parameters ##### --------------------------------------------------------
source(here::here("data/config/config.R")) # Import global parameters
source(here::here("data/config/config-STOC.R")) # Import local parameters

### Dataset
set.seed(496) # for reproducible results


##### Helper functions ##### --------------------------------------------------
explore_dataset <- function(df, top = 5) {
    ### NUMERICS
    cli_alert_info("Exploration of dataset.")

    # for fastest computation of figure: average per square ("carre")
    draftman_df <- df |>
        select(c("carre", all_of(X_VARIABLES))) |>
        select(where(is.numeric)) |>
        group_by(carre) |>
        summarise_all(mean)

    draft_plot <- ggplot_custom_draftman(
        draftman_df, 
        columns = setdiff(names(draftman_df), "carre"))
    suppressMessages(print(draft_plot))
    standardised_ggplot_save(
        figure = draft_plot, 
        save_path = file.path(
            PATH_LOCAL_RESULTS, 
            paste0("draftman_plot.pdf")))
    cli_alert_success("Saved Draftman's plot.")

    ### FACTORS
    factor_df <- df |>
        select(where(is.factor))
    for (col in names(factor_df)) {
        freq <- table(factor_df[[col]])
        prop <- round(freq/sum(freq), digits=3)
        combined_table <- rbind(freq = as.vector(freq), 
                            prop = as.vector(prop))

        colnames(combined_table) <- names(freq)

        cli_alert_info(paste0("Distribution of '",  col ,"' in dataset:"))
        print(t(combined_table))
        cli_alert_info(" ")
    }

    ### PREDICTED VARIABLES
    top_sp <- df[NAMES_SPECIES] |>
        summarise(across(where(is.numeric), sum)) |>
        pivot_longer(everything(), names_to = "column", values_to = "sum")
    cli_alert_info(paste0(
        "Top-", top, " MOST sighted species (from occurences):"))
    print(as.data.frame(top_sp |> slice_max(sum, n = top)))
    cli_alert_info(" ")
    cli_alert_info(paste0(
        "Top-", top," LEAST sighted species (from occurences):"))
    print(as.data.frame(top_sp |> slice_min(sum, n = top)))
    cat("\n")

}

split_points_k_fold_subsets <- function(df) {
    # To approach independent sampling : 
    #   - training : each sample is selected in a different square
    #   - new pool (simulation of 500 ENI) : each sample is selected in a 
    #       different square
    #   - test data : remaining samples after selection of training + new pool
    cli_alert_info("Splitting datasets...")

    if ((TRAIN_SIZE + NEW_POOL_SIZE) >= length(unique(df$carre))) {
        stop(paste0(
            "Cannot split dataset, inconsistent given sizes of subsets.\n",
            "train_size + new_pool_size = ", TRAIN_SIZE + NEW_POOL_SIZE, 
            " elements.", "Should be <", length(unique(df$carre)), 
            " (number of squares in dataset)."))
    }

    # Storing IDs in list (instead of full dataframes) to save storage
    k_fold_list <- list(
        training_points = list(),                   # train set
        training_squares_unknown_points = list(),   # val set
        new_pool_points = list(),                   # train OED set
        new_pool_squares_unknown_points = list(),   # val OED set
        test_points = list()                        # test set
    )
    for (k in seq(K_FOLD)) {
        # select squares
        random_order <- sample(
            unique(df$carre), length(unique(df$carre)), replace = FALSE)
        train_carre <- random_order[1:TRAIN_SIZE]
        new_pool_carre <- random_order[(TRAIN_SIZE+1):(TRAIN_SIZE+NEW_POOL_SIZE)]

        # for train sets : select one point per square
        train_pts <- sample(
            subset(df, df$carre %in% train_carre)$id_point_annee, 
            TRAIN_SIZE, 
            replace = FALSE)
        new_pool_pts <- sample(
            subset(df, df$carre %in% new_pool_carre)$id_point_annee, 
            NEW_POOL_SIZE, 
            replace = FALSE)

        # for validation sets : select remaining points in each square
        train_square_df <- subset(df, df$carre %in% train_carre)
        val_train_pts <- subset(
            train_square_df, 
            !train_square_df$id_point_annee %in% train_pts)$id_point_annee

        pool_square_df <- subset(df, df$carre %in% new_pool_carre)
        val_pool_pts <- subset(
            pool_square_df, 
            !pool_square_df$id_point_annee %in% new_pool_pts)$id_point_annee
        

        # for test : select squares not seen before
        test_never_seen_points <- subset(
            df, !df$carre %in% c(train_carre, new_pool_carre))$id_point_annee

        # assign values to lists
        k_fold_list$training_points[[k]] <- train_pts
        k_fold_list$new_pool_points[[k]] <- new_pool_pts
        k_fold_list$val_training_points[[k]] <- val_train_pts
        k_fold_list$val_new_pool_points[[k]] <- val_pool_pts
        k_fold_list$test_points[[k]] <- test_never_seen_points
    }

    cli_alert_success("Point ID of splits are ready for each split!\n\n")
    return(k_fold_list)
}

prepare_training <- function(subset, random_effect = "none"){
    cli_alert_info("Preparing training set for HMSC...")

    Y <- as.matrix(subset[if (is.null(Y_SPECIES)) NAMES_SPECIES else Y_SPECIES]) 
    XData <- as.data.frame(setNames(
        lapply(X_VARIABLES, function(col) subset[[col]]),
        X_VARIABLES)
    )

    studyDesign <- data.frame(
        units = as.factor(subset$point), 
        spatial = as.factor(subset$id_point_annee)
    )

    
    if (random_effect == "none") {
        cli_alert_info("Creation of Hmsc object without random effects...")
        hmsc_object <- Hmsc(
            Y = Y, XData = XData, 
            XFormula = HMSC_XFORMULA, 
            distr = "probit",           # For occurence data (among "normal", "probit", "poisson", "lognormal poisson")
            studyDesign = studyDesign)
        
    } else if (random_effect == "units") {
        cli_alert_info("Creation of Hmsc object with 'units' effect...")
        # a random effect associated to each point
        rL.points = HmscRandomLevel(units = unique(studyDesign$point))  # random point effect
        rL.points = setPriors(rL.points, nfMin =1,  nfMax = 1)          # limit number of latent variables
        hmsc_object <- Hmsc(
            Y = Y, XData = XData, 
            ranLevels = list("units" = rL.points),
            XFormula = HMSC_XFORMULA, 
            distr = "probit",           # For occurence data
            studyDesign = studyDesign)

    }  else if (random_effect == "spatial") {
        cli_alert_info("Creation of Hmsc object with 'spatial' effect...")
        # convert coordinates to metric
        coords_sf <- st_as_sf(
            subset, 
            coords = c("LON", "LAT"), 
            crs = 4326)
        coords_proj <- st_transform(coords_sf, crs = 2154)

        xy <- st_coordinates(coords_proj)
        rownames(xy) <- as.character(subset$id_point_annee)
        colnames(xy) <- c("longitude_grid_2154", "latitude_grid_2154")
        xy_spatial <- xy[!duplicated(rownames(xy)), ]

        if (nrow(xy_spatial) < 999) {
            rL.spatial <- HmscRandomLevel(sData = xy_spatial)
        } else {
            # too many samples, need optimization for faster computation
            rL.spatial <- HmscRandomLevel(
                sData = xy_spatial, sMethod = "NNGP", nNeighbours = 10)
        }
        rL.spatial = setPriors(rL.spatial, nfMin =1,  nfMax = 1)    # limit number of latent variables

        # # Cover 100m to ~1000km (on a log scale)
        # alpha_values <- c(0, exp(seq(log(1e2), log(1e6), length.out = 9)))
        # alphapw <- cbind(
        #   alpha_values, rep(1/length(alpha_values), length(alpha_values)))

        hmsc_object <- Hmsc(
            Y = Y, XData = XData, 
            ranLevels = list("spatial" = rL.spatial),
            XFormula = HMSC_XFORMULA, 
            distr = "probit",                           # For occurence data
            studyDesign = studyDesign)

    } else {
        stop(paste0(paste0("Random effect must be one of 'none', 'units' or",
        " 'spatial'. Got ", random_effect)))
    }

    cli_alert_success("Created Hmsc object!\n\n")
    return(hmsc_object)
}

fitting_hmsc <- function(
    hmsc_object, 
    save_to,
    freq_verbose = 100,
    allow_parallel = TRUE) {
    # Fit model
    start <- Sys.time()
    cli_alert_info(paste0("Started fitting at: ", start))

    if (!is.null(freq_verbose) & (CHAINS > 1) & allow_parallel) {
        cli_alert_warning(paste0("Cannot display fitting progress when running",
        " chains in parallel."))
        cli_alert_warning("Set `allow_parallel` to `FALSE` to see progress.")
    }

    if (allow_parallel) {
        fitted.hmsc <-  sampleMcmc(
            hmsc_object, 
            thin = THIN, samples = SAMPLES, transient = TRANSIENT, 
            nChains = CHAINS, nParallel = CHAINS, updater=list(GammaEta=FALSE),
            verbose = freq_verbose)
    } else {
        fitted.hmsc <-  sampleMcmc(
            hmsc_object, 
            thin = THIN, samples = SAMPLES, transient = TRANSIENT, 
            nChains = CHAINS, nParallel = 1, updater=list(GammaEta=FALSE),
            verbose = freq_verbose)
    }
    stop <- Sys.time()
    cli_alert_info(paste0("Completed fitting at: ", stop))
    cli_alert_info(paste0("Time elapsed: ", round(stop-start, 2)))

    # Save run 
    cli_alert_info("Saving model...")
    save(fitted.hmsc, file=save_to)
    cli_alert_success("Model saved!\n\n")
    return(fitted.hmsc)
}

convergence_hmsc <- function(hmsc_output, save_folder) {
    coda_outputs <- convertToCodaObject(hmsc_output)
    
    # Summary plots (not utilized here but they could be!)
    # MCMCsummary(object = coda_outputs$Beta, round = 2) 
    # MCMCplot(object = coda_outputs$Beta)

    # For HMSC parameters are named, we are interested in:
    #   - Beta: fixed effects
    #   - Omega: random variation in co-occurence
    for (param in c("Beta", "Omega")) {
        if (!(param %in% names(coda_outputs))) {
            next
        } 
        cli_alert_info(paste0("*** [Parameter: ", param, "] ***"))

        ### Convergence diagnostics
        if (param == "Omega") {
            chains <- coda_outputs$Omega[[1]]  # 3D array: [iter, sp, sp]
        } else {
            chains <- coda_outputs[[param]]
        }
        
        ## Traceplot (Rhat and effective size)      
        cli_alert_info("Computation of traceplots, can take some time...")
        # faster version:
        tryCatch({ # avoids debugger mode
            MCMCtrace(
                object = chains,
                pdf = TRUE,
                filename = file.path(
                    save_folder, 
                    paste0(param, "_all_traceplots.pdf")),
                ind = TRUE,
                open_pdf = FALSE,
                plot = TRUE,
                Rhat = TRUE, # ajoute le Rhat
                n.eff = TRUE, # ajoute la taille d’échantillon effective
                type = "both"   # explicitly request trace + density
            )
        }, error = function(e) {
            message("Error on iteration ", i, ": ", e$message)
        })
        # # Slower but more beautiful version
        # traceplots <- ggplot_custom_MCMCtrace(
        #     coda_object = coda_fitted_model$Beta,
        #     show_Rhat = TRUE,
        #     show_Neff = TRUE)       
        # for (i_plot in seq_along(traceplots)) {
        #     standardised_ggplot_save(
        #         figure = (traceplots[[i_plot]]$trace + traceplots[[i_plot]]$density), 
        #         save_path = file.path(
        #             PATH_LOCAL_RESULTS, 
        #             paste0("Beta_", i_plot, "_traceplot.pdf")))
        # }
        cli_alert_info("Traceplots saved!")


        ## Effective size
        cli_alert_info("-> Effective size:")
        eff_size <- effectiveSize(chains)
        eff_size_plot <- ggplot_bars(
                as.data.frame(eff_size), "eff_size", 
                breaks = ceiling(c(0, seq(400, max(eff_size), length.out = 19))),
                underlayers = list(
                    bad = annotate("rect", 
                        xmin = -Inf, xmax = 400, ymin = -Inf, ymax = Inf,
                        fill = PALETTE[1], alpha = 0.5),
                    acceptable = annotate("rect", 
                        xmin = 400, xmax = 1000, ymin = -Inf, ymax = Inf,
                        fill = PALETTE[5], alpha = 0.5),
                    good = annotate("rect", 
                        xmin = 1000, xmax = Inf, ymin = -Inf, ymax = Inf,
                        fill = PALETTE[2], alpha = 0.5),
                    captions = labs(caption = "Green: good, Orange: acceptable, Red: bad."))
                )
        print(eff_size_plot)
        interpret_diagnostics(
            eff_size, bad = 100, good = 400, order = "high_better", mode = "quick")
        standardised_ggplot_save(
            figure = eff_size_plot, 
            save_path = file.path(
                save_folder, 
                paste0("hist_neff_", param, ".pdf")))

        ## Gelman-Rubin convergence diagnostic
        if (CHAINS > 1) {
            cli_alert_info("-> Gelman-Rubin convergence diagnostic:")
            psrf <- gelman.diag(
                chains,  multivariate = FALSE)$psrf[, "Point est."]
            psrf_plot <- ggplot_bars(
                as.data.frame(psrf), "psrf", bins = 10,                
                underlayers = list(
                    good = annotate("rect", 
                        xmin = -Inf, xmax = 1.05, ymin = -Inf, ymax = Inf,
                        fill = PALETTE[2], alpha = 0.5),
                    acceptable = annotate("rect", 
                        xmin = 1.05, xmax = 1.1, ymin = -Inf, ymax = Inf,
                        fill = PALETTE[5], alpha = 0.5),
                    good = annotate("rect", 
                        xmin = 1.1, xmax = Inf, ymin = -Inf, ymax = Inf,
                        fill = PALETTE[1], alpha = 0.5),
                    captions = labs(caption = "Green: good, Orange: acceptable, Red: bad."))
                )
            print(psrf_plot)
            interpret_diagnostics(psrf, bad=1.1, good=1.05, mode = "quick")
            standardised_ggplot_save(
                figure = psrf_plot, 
                save_path = file.path(
                    save_folder, 
                    paste0("hist_psrf_", param, ".pdf")))

            # ## Geweke diagnostic
            # cli_alert_info("-> Geweke convergence diagnostic:")
            # # Rule of thumb: <2 (no proof of non-convergence), >2 (monitor convergence closely)
            # geweke_estimates_all <- geweke.diag(chains)
            # for (i in seq(nchain(chains))) {
            #     cli_alert_info(paste0("[Chain ", i, "]"))
            #     interpret_diagnostics(geweke_estimates_all[[i]][[1]], good=2, mode = "quick")
            # }

            ## Autocorrelation
            cli_alert_info("-> Autocorrelation:")
            for (lag in c(20, 50)) {
                autocorr_estimates <- autocorr.diag(
                    chains, lags=c(lag/THIN))
                cli_alert_info(paste0("[Lag ", lag, "]:"))
                interpret_diagnostics(autocorr_estimates, good=0.1, mode = "quick")
            }
            
        }
    }
}

evaluate_hmsc_performances <- function(hM, subset) {
    local_preds_list <- predict_hmsc(
        hM = hM, 
        df = subset, 
        x_variables = X_VARIABLES)
    local_preds <- abind(local_preds_list, along = 3)

    local_Y <- as.matrix(subset[Y_SPECIES])  # actual observed values for test set
    evaluateModelFitCustom(hM = hM, Y = local_Y, predY = local_preds)
}

analyses_hmsc <- function(hmsc_output, save_folder, supportLevel = 0.95) {
    # X-Y associations
    for (param in c("Beta", "Omega")) {
        if (is.null(hmsc_output$ranLevels) & (param=="Omega")){
            next
        }
        post_association = getPostEstimate(hmsc_output, parName = param)
        XY_grid <- ggplot_custom_plotBeta(
            hmsc_output,
            post = post_association,
            supportLevel = supportLevel
            )
        standardised_ggplot_save(
            figure = XY_grid, 
            save_path = file.path(save_folder, paste0(param, "_XY_associations.pdf")))
    }

    if (!is.null(hmsc_output$ranLevels)) {
        rand_XX_grid <- ggplot_custom_random_corr_associations(
            hmsc_output, supportLevel = supportLevel)
        standardised_ggplot_save(
            figure = rand_XX_grid, 
            save_path = file.path(save_folder, "random_XX_associations.pdf"))
    }
    
    ## Variance partitionning 
    VP = computeVariancePartitioning(
        hmsc_output, 
        group = X_GROUPS, # c(1,2,2)
        groupnames = X_GROUPS_NAMES) # c("habitat","climate"))
    variance_bars <- ggplot_custom_plotVariancePartitioning(hmsc_output, VP = VP)
    standardised_ggplot_save(
        figure = variance_bars, 
        save_path = file.path(save_folder, "variance_partitioning.pdf"))
}

compare_hmsc_performances <- function(save_folder) {
    cli_alert_info("Fetching scores...")
    scores_df <- data.frame()
    for (k in seq(K_FOLD)) {
        r_effects <- c("none", "units") # "spatial" is available but way longer to run
        for (r in seq_along(r_effects)) { 
            # 0. Setup
            r_effect <- r_effects[r]
            path_local_scores <- file.path(
                PATH_LOCAL_BASE,
                paste0("base-model_", r_effect, "-random-effect_k", k))


            # load train/val/test/scores
            train_scores <- read_csv(
                file.path(path_local_scores, "train_scores.csv"), 
                show_col_types = FALSE)
            val_scores <- read_csv(
                file.path(path_local_scores, "val_scores.csv"),
                show_col_types = FALSE)
            test_scores <- read_csv(
                file.path(path_local_scores, "test_scores.csv"),
                show_col_types = FALSE)
            
            # add columns 
            train_scores <- train_scores |>
                pivot_longer(
                    cols = c(RMSE, AUC, TjurR2), 
                    names_to = "metric", values_to = "score") |>
                mutate(effect = r_effect, k_fold = k, dataset = "train")
            val_scores <- val_scores |>
                pivot_longer(
                    cols = c(RMSE, AUC, TjurR2), 
                    names_to = "metric", values_to = "score") |>
                mutate(effect = r_effect, k_fold = k, dataset = "val")  
            test_scores <- test_scores |>
                pivot_longer(
                    cols = c(RMSE, AUC, TjurR2), 
                    names_to = "metric", values_to = "score") |>
                mutate(effect = r_effect, k_fold = k, dataset = "test")  

            # concatenate
            scores_df <- rbind(train_scores, val_scores, test_scores, scores_df)
        }
    }

    ### Plot scores
    aggregated_df <- scores_df |>
        group_by(metric, effect, dataset) |>
        summarise(
            avg_score = mean(score, na.rm = TRUE),
            sd_score = sd(score, na.rm = TRUE),
            .groups = "drop_last") |>
        mutate(dataset = factor(dataset, levels = c("train", "val", "test")))

    cli_alert_info("Creating plot...")
    p <- ggplot(aggregated_df, aes(y = avg_score, x = effect, fill = dataset)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.66), width=0.66) +
        geom_errorbar(
            aes(ymin = avg_score - sd_score, ymax = avg_score + sd_score),
            position = position_dodge(width = 0.66),
            width = 0.2
        ) +
        labs(caption = "SD and mean computed per k_fold (over all species)") +
        facet_wrap(~ metric, scales = "free_y")

    p <- my_custom_ggplot_theme(p) + 
        scale_fill_manual(values = c(PALETTE[2], PALETTE[3], PALETTE[1]))
    print(p)
    standardised_ggplot_save(
        figure = p, 
        save_path = file.path(save_folder, "base_average_performances.pdf"))
    cli_alert_success("Plot of performances saved!\n\n")
    return(aggregated_df)
}

map_results_hmsc <- function(r_effect = "none", k_fold = 1, sp = "Sylvia_atricapilla") {
    cli_alert_info("Mapping results on whole dataset...")
    cli_alert_info(paste0("\tRandom effect: ", r_effect))
    cli_alert_info(paste0("\tK_fold: ", k_fold))
    local_folder <- file.path(
            PATH_LOCAL_BASE,
            paste0("base-model_", r_effect, "-random-effect_k", k_fold))
    load(file.path(local_folder, "train_outputs.Rda")) # loads `fitted.hmsc`
    load(file.path(PATH_LOCAL_RESULTS, "k_fold_points.Rda")) # loads `k_fold_points`

    cli_alert_info("Making predictions, can take a few minutes...")
    full_preds_list <- predict_hmsc(
        hM = fitted.hmsc,
        df = stoc_df,
        x_variables = X_VARIABLES)
    cli_alert_info("Predictions are ready.")

    # Predictions are DISTRIBUTIONS: let's take the average for each point
    cli_alert_info("Converting distributions to point average...")
    full_preds_avg_occurrence <- apply(
        simplify2array(full_preds_list), c(1,2), mean)
    stoc_df_preds <- stoc_df |>
        mutate(suitability = full_preds_avg_occurrence[, sp]) |>
        mutate(uncertainty = 1-abs(suitability - 0.5) * 2) |>
        mutate(subset = ifelse(
            id_point_annee %in% k_fold_points$training_points[[k_fold]], 
            "train",
            ifelse(
                id_point_annee %in% k_fold_points$val_training_points[[k_fold]], 
                "val", 
                "test")))
    # To cite 10.5281/ZENODO.11067678: "The models cannot provide a direct 
    # indication of the probability of presence – they instead give an 
    # index of habitat suitability, measured on a scale from 0 to 1."
    cli_alert_info(paste0("\tSelected species: ", sp))

    # plot training samples
    p0.1 <- ggplot_categorical_df_on_background_map(
        background_map = ggplot_get_france_base_map("national"), 
        df = stoc_df_preds |> filter(subset == "train"), 
        LON = "LON",
        LAT = "LAT",
        column = NULL,
        legend_title = "Training location")
    print(p0.1)
    standardised_ggplot_save(
        figure = p0.1, 
        save_path = file.path(local_folder, "training_samples.pdf"))
    p0.2 <- ggplot_categorical_df_on_background_map(
        background_map = ggplot_get_france_base_map("national"), 
        df = stoc_df_preds |> filter(subset == "test"), 
        LON = "LON",
        LAT = "LAT",
        column = NULL,
        legend_title = "Training location")
    print(p0.2)
    standardised_ggplot_save(
        figure = p0.2, 
        save_path = file.path(local_folder, "test_samples.pdf"))

    # show map of suitability per site
    p1.1 <- ggplot_quantitative_df_on_background_map(
            background_map = ggplot_get_france_base_map("national"), 
            df = stoc_df_preds, 
            LON = "LON",
            LAT = "LAT",
            column = "suitability",
            unit = paste0("Estimated HSI for ", sp)) +
        labs(caption = "HSI = 'Habitat Suitability Index'")
    print(p1.1)
    standardised_ggplot_save(
        figure = p1.1, 
        save_path = file.path(local_folder, "suitability_for_species.pdf"))
    shp1.2 <- interpolate_scattered_points_to_hexagons(
            df = stoc_df_preds,  
            column = "suitability",
            res_km = RES_KM, 
            LON = "LON",
            LAT = "LAT",
            idp = 2,           
            maxdist_m = 100  )
    p1.2 <- ggplot_quantitative_shapefile_on_background_map(
        background_map = ggplot_get_france_base_map("national"),
        shapefile = shp1.2,
        layer_name = "interpolated_value",
        unit = paste0("Estimated HSI for ", sp),
        limits = NULL,
        precision_auto_limits = 1) +
        labs(caption = "HSI = 'Habitat Suitability Index'")
    print(p1.2)
    standardised_ggplot_save(
        figure = p1.2, 
        save_path = file.path(local_folder, "suitability_for_species_hexagons.pdf"))    
    
    

    # show map of uncertainty per site (function of suitability)
    p2.1 <- ggplot_quantitative_df_on_background_map(
            background_map = ggplot_get_france_base_map("national"), 
            df = stoc_df_preds, 
            LON = "LON",
            LAT = "LAT",
            column = "uncertainty",
            unit = paste0("Rough uncertainty of HSI for ", sp)) +
        labs(caption = "HSI = 'Habitat Suitability Index'")
    print(p2.1)
    standardised_ggplot_save(
        figure = p2.1, 
        save_path = file.path(local_folder, "certainty_of_suitability.pdf"))
    shp2.2 <- interpolate_scattered_points_to_hexagons(
            df = stoc_df_preds,  
            column = "uncertainty",
            res_km = RES_KM, 
            LON = "LON",
            LAT = "LAT",
            idp = 2,           
            maxdist_m = 100)
    p2.2 <- ggplot_quantitative_shapefile_on_background_map(
        background_map = ggplot_get_france_base_map("national"),
        shapefile = shp2.2,
        layer_name = "interpolated_value",
        unit = paste0("Rough uncertainty of HSI for ", sp),
        limits = NULL,
        precision_auto_limits = 1) +
        labs(caption = "HSI = 'Habitat Suitability Index'")
    print(p2.2)
    standardised_ggplot_save(
        figure = p2.2, 
        save_path = file.path(local_folder, "certainty_of_suitability_hexagons.pdf"))    

    return(stoc_df_preds)
}


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
. <- explore_dataset(stoc_df)

# split train/test subsets 'k_fold' times
k_fold_points <- split_points_k_fold_subsets(stoc_df) 
save(k_fold_points, file=file.path(PATH_LOCAL_RESULTS, "k_fold_points.Rda"))


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
        base_model <- prepare_training(train_subset, random_effect = r_effect)
        cli_alert_info(paste0(r, ".1. Fitting model..."))
        fitted_model <- fitting_hmsc(
            base_model, 
            save_to =  file.path(path_local_model_results, "train_outputs.Rda"),
            freq_verbose = (SAMPLES*2+TRANSIENT)/10,
            allow_parallel = TRUE
        )
        
        # 2. Analysis of convergence
        cli_alert_info(paste0(r, ".2. Convergence diagnostics..."))
        . <- convergence_hmsc(fitted_model, save_folder = path_local_model_results)
        cat("\n")

        # 3. Analysis of performance
        cli_alert_info(paste0(r, ".3. Performance evaluation..."))
        
        # Explanatory power
        cli_alert_info("Computing training scores...")
        train_scores <- evaluate_hmsc_performances(hM = fitted_model, subset = train_subset)

        # Prediction power
        cli_alert_info("Computing testing scores...")
        val_scores <- evaluate_hmsc_performances(hM = fitted_model, subset = val_subset)
        test_scores <- evaluate_hmsc_performances(hM = fitted_model, subset = test_subset)

        cli_alert_info("Saving scores...")
        write_csv(as.data.frame(train_scores), file.path(path_local_model_results, "train_scores.csv"))
        write_csv(as.data.frame(val_scores), file.path(path_local_model_results, "val_scores.csv"))
        write_csv(as.data.frame(test_scores), file.path(path_local_model_results, "test_scores.csv"))
        # row names must be saved seperatly. (dropped)
        cli_alert_success("Scores saved!\n\n")

        # 4. Result analysis
        cli_alert_info(paste0(r, ".4. Associations..."))
        . <- analyses_hmsc(
            hmsc_output = fitted_model, 
            save_folder = path_local_model_results,
            supportLevel = 0.05
        )
        cli_alert_warning("TODO: make maps!!!")
        cli_alert_success("Model ran without errors!\n\n")
    }
}


##### Compare performances ##### ----------------------------------------------
cli_alert_info("------------ Results ------------\n\n")
agg_scores_df <- compare_hmsc_performances(PATH_LOCAL_BASE)

cli_alert_info(
    "Given the results, there is little to no difference in prediction power ",
    "on new (test) data.")
cli_alert("For OED, we'll therefore use a model without random effects.")
cli_alert("Training base model on whole dataset...")



##### Predict on map of france ##### ------------------------------------------
cli_alert_info("------------ Distribution map (example) ------------\n\n")
for (k in seq(K_FOLD)) {
    a <- map_results_hmsc(k_fold = k)
}
