# Set of functions used to create / handle / analyse models
##### Libraries ##### ---------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(abind)
library(readr)
library(Hmsc)
library(cli)
library(sf)


##### Parameters ##### --------------------------------------------------------
source(here::here("data/config/config.R")) # Import global parameters


##### Functions ##### ---------------------------------------------------------
# A function to prepare dataset for Hmsc training, outputs a model ready for training
# ARGS:
#   - combinations: a list of lists containing elements TRAIN_SIZES, R_EFFECTS, STRATEGIES, HMSC_XFORMULAS, NEW_SAMPLE_SIZE
build_param_grid <- function(combinations) {
    # create a grid, each line is a combination of variables
    # that way, we loop on this table instead of having nested loops
    grids <- lapply(combinations, function(p) {
        expand.grid(
            train_size = p$TRAIN_SIZES,
            r_effect = p$R_EFFECTS,
            strategy = p$STRATEGIES,
            formulas = p$HMSC_XFORMULAS,
            n_new_samples = p$NEW_SAMPLE_SIZES,
            stringsAsFactors = FALSE
        )
    })
    dplyr::distinct(dplyr::bind_rows(grids))
}

# A function to prepare dataset for Hmsc training, outputs a model ready for training
# ARGS:
#   - subset: a dataframe. Must contain columns listed in x_cols and y_cols.
#   - x_cols: a list of strings. The columns containing explanatory variables.
#   - y_cols: a list of strings. The columns containing species occurrences.
#   - formula: a formula for the Hmsc model. Based on names in x_cols/y_cols.
#   - random_effect: whether to add a random effect to Hmsc ("units" or "spatial") or not ("none", default).
prepare_hmsc_training <- function(
        subset, x_cols, y_cols, formula, random_effect = "none"){
    cli_alert_info("Preparing training set for HMSC...")

    # Hmsc only accepts data in a very specific format
    ydata <- as.matrix(subset[y_cols]) 
    xdata <- as.data.frame(setNames(
        lapply(x_cols, function(col) subset[[col]]),
        x_cols)
    )
    studyDesign <- data.frame(
        units = as.factor(subset$point), 
        spatial = as.factor(subset$id_point_annee)
    )

    
    if (random_effect == "none") {
        cli_alert_info("Creation of Hmsc object without random effects...")
        hmsc_object <- Hmsc(
            Y = ydata, XData = xdata, 
            XFormula = formula, 
            distr = "probit",           # For occurence data (among "normal", "probit", "poisson", "lognormal poisson")
            studyDesign = studyDesign)
        
    } else if (random_effect == "units") {
        cli_alert_info("Creation of Hmsc object with 'units' effect...")
        # a random effect associated to each point
        rL.points = HmscRandomLevel(units = unique(studyDesign$point))  # random point effect
        rL.points = setPriors(rL.points, nfMin =1,  nfMax = 1)          # limit number of latent variables
        hmsc_object <- Hmsc(
            Y = ydata, XData = xdata, 
            ranLevels = list("units" = rL.points),
            XFormula = formula, 
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
            Y = ydata, XData = xdata, 
            ranLevels = list("spatial" = rL.spatial),
            XFormula = formula, 
            distr = "probit",                           # For occurence data
            studyDesign = studyDesign)

    } else {
        stop(paste0(paste0("Random effect must be one of 'none', 'units' or",
        " 'spatial'. Got ", random_effect)))
    }

    cli_alert_success("Created Hmsc object!\n\n")
    return(hmsc_object)
}

# A function to fit a Hmsc model and save results
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - save_to: a string. The path where the file will be saved (must end with .rds).
#   - nchains: a numeric. The number of chains to run.
#   - thin: a numeric. The number of steps between each recording of a sample.
#   - nsamples: a numeric. The number of samples to collect
#   - ntransient: a numeric. The number of steps to wait for before collecting samples.
#   - freq_verbose: a numeric (default is 100). The frequency of verbose messages, you get one every "freq_verbose" step.
#   - allow_parallel: a boolean (default is TRUE). If nchains > 1, allows to compute chains in parallel for faster fitting. Removes verbose messages during fitting.
fitting_hmsc <- function(
    hM, 
    save_to,
    nchains,
    thin,
    nsamples,
    ntransient,
    freq_verbose = 100,
    allow_parallel = TRUE) {
    # Fit model
    start <- Sys.time()
    cli_alert_info(paste0("Started fitting at: ", start))

    if (!is.null(freq_verbose) & (nchains > 1) & allow_parallel) {
        cli_alert_warning(paste0("Cannot display fitting progress when running",
        " chains in parallel."))
        cli_alert_warning("Set `allow_parallel` to `FALSE` to see progress.")
    }

    if (allow_parallel) {
        fitted.hmsc <-  sampleMcmc(
            hM, 
            thin = thin, samples = nsamples, transient = ntransient, 
            nChains = nchains, nParallel = nchains, updater=list(GammaEta=FALSE),
            verbose = freq_verbose)
    } else {
        fitted.hmsc <-  sampleMcmc(
            hM, 
            thin = thin, samples = nsamples, transient = ntransient, 
            nChains = nchains, nParallel = 1, updater=list(GammaEta=FALSE),
            verbose = freq_verbose)
    }
    stop <- Sys.time()
    cli_alert_info(paste0("Completed fitting at: ", stop))
    cli_alert_info(paste0("Time elapsed: ", round(stop-start, 2)))

    # Save run 
    cli_alert_info("Saving model...")
    saveRDS(fitted.hmsc, file = save_to)
    cli_alert_success("Model saved!\n\n")
    return(fitted.hmsc)
}

# A function to display convergence diagnostics for a Hmsc model
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - nchains: a numeric. The number of chains to run.
#   - thin: a numeric. The number of steps between each recording of a sample.
#   - save_folder: a string. The path where plotted PDFs will be saved.
convergence_hmsc <- function(hM, nchains, thin, save_folder) {
    coda_outputs <- convertToCodaObject(hM)
    
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
            )}, 
        error = function(e) {
            message("Error in MCMCtrace: ", e$message)
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
        #             save_folder, 
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
        if (nchains > 1) {
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
                    chains, lags=c(lag/thin))
                cli_alert_info(paste0("[Lag ", lag, "]:"))
                interpret_diagnostics(autocorr_estimates, good=0.1, mode = "quick")
            }
            
        }
    }
    cat("\n")
}

# A function to display convergence diagnostics for a Hmsc model
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - subset: a dataframe. Must contain columns listed in x_cols and y_cols.
#   - x_cols: a list of strings. The columns containing explanatory variables.
#   - sp_cols: a list of strings. The columns containing species occurrences.
evaluate_hmsc_performances <- function(hM, subset, x_cols, sp_cols) {
    local_preds_list <- predict_hmsc(
        hM = hM, 
        df = subset, 
        x_variables = x_cols)
    local_preds <- abind(local_preds_list, along = 3)

    local_Y <- as.matrix(subset[sp_cols])  # actual observed values for test set
    evaluateModelFitCustom(hM = hM, Y = local_Y, predY = local_preds)
}

# A function to display XX and XY associations after Hmsc fitting
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - x_groups_cats: a list of integers. For each x variable, a number assigning it to a group (for variance partitioning)
#   - x_groups_names: a list of strings. A label for each number in x_groups_cat.
#   - save_folder: a string. The path where plotted PDFs will be saved.
#   - supportLevel: a numeric between 0 and 1. The minimum confidence to display results (default is 0.95)
analyses_hmsc <- function(
        hM, save_folder, x_groups_cats, x_groups_names, supportLevel = 0.95) {
    # X-Y associations
    for (param in c("Beta", "Omega")) {
        if (is.null(hM$ranLevels) & (param=="Omega")){
            next
        }
        post_association = getPostEstimate(hM, parName = param)
        XY_grid <- ggplot_custom_plotBeta(
            hM, post = post_association, supportLevel = supportLevel)
        standardised_ggplot_save(
            figure = XY_grid, 
            save_path = file.path(save_folder, paste0(param, "_XY_associations.pdf")))
    }

    if (!is.null(hM$ranLevels)) {
        rand_XX_grid <- ggplot_custom_random_corr_associations(
            hM, supportLevel = supportLevel)
        standardised_ggplot_save(
            figure = rand_XX_grid, 
            save_path = file.path(save_folder, "random_XX_associations.pdf"))
    }
    
    # Variance partitionning 
    vp = computeVariancePartitioning(
        hM, 
        group = x_groups_cats, # c(1,2,2)
        groupnames = x_groups_names) # c("habitat","climate"))
    variance_bars <- ggplot_custom_plotVariancePartitioning(hM, VP = vp)
    standardised_ggplot_save(
        figure = variance_bars, 
        save_path = file.path(save_folder, "variance_partitioning.pdf"))
}

# A function to compare training scores between k_folds, subset and model type.
# ARGS:
#   - parent_folder: a string. The path where subfolders of each model are located.
#   - save_to: a string. Name of the file for saving resulting plot (should end with .pdf). If NULL, does not save pdf.
#   - loop_prefix: a string. The prefix of the subfolders names.
#   - loop_elements: a list of strings. Middle elements for subfolders names.
#   - loop_suffix: a string. The suffix of the subfolders names.
#   - k_fold: a numeric. The number of cross-validation subsets to make. Goes after suffix.
#   - xlabel: a string. The label for the x-axis of the plot (default is "Model type").
#   - ylabel: a string. The label for the y-axis of the plot (default is "Model type").
#   - group_species: whether to take the mean accross all k_folds and species (TRUE, default) or only accross k_folds (FALSE).
#   - species_names: names of species in CSV (rownames are not available from csvs).
#   - barplot: whether to make a barplot (TRUE, default) or a pointplot (FALSE).
barplot_raw_scores <- function(
        parent_folder, 
        save_to,
        loop_prefix, 
        loop_elements, 
        loop_suffix, 
        k_fold, 
        xlabel = "Model type", 
        ylabel = "Average Score",
        group_species = TRUE,
        species_names = NULL,
        barplot = TRUE) {
    cli_alert_info("Fetching scores...")
    scores_df <- data.frame()
    for (k in seq(k_fold)) {
        for (i in seq_along(loop_elements)) { 
            # 0. Setup
            loop_element <- loop_elements[i]
            path_local_scores <- file.path(
                parent_folder, paste0(loop_prefix, loop_element, loop_suffix, k))


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
                mutate(species = species_names) |>
                pivot_longer(
                    cols = c(RMSE, AUC, TjurR2), 
                    names_to = "metric", values_to = "score") |>
                mutate(loop_element = loop_element, k_fold = k, dataset = "train")
            val_scores <- val_scores |>
                mutate(species = species_names) |>
                pivot_longer(
                    cols = c(RMSE, AUC, TjurR2), 
                    names_to = "metric", values_to = "score") |>
                mutate(loop_element = loop_element, k_fold = k, dataset = "val")  
            test_scores <- test_scores |>
                mutate(species = species_names) |>
                pivot_longer(
                    cols = c(RMSE, AUC, TjurR2), 
                    names_to = "metric", values_to = "score") |>
                mutate(loop_element = loop_element, k_fold = k, dataset = "test")  

            # concatenate
            scores_df <- rbind(train_scores, val_scores, test_scores, scores_df)
        }
    }

    ### Plot scores
    if (group_species) {
        aggregated_df <- scores_df |>
            group_by(metric, loop_element, dataset) |>
            summarise(
                avg_score = mean(score, na.rm = TRUE),
                sd_score = sd(score, na.rm = TRUE),
                .groups = "drop_last") |>
            mutate(dataset = factor(dataset, levels = c("train", "val", "test")))
    } else {
        aggregated_df <- scores_df |>
            group_by(metric, loop_element, dataset, species) |>
            summarise(
                avg_score = mean(score, na.rm = TRUE),
                sd_score = sd(score, na.rm = TRUE),
                .groups = "drop_last") |>
            mutate(dataset = factor(dataset, levels = c("train", "val", "test")))
    }


    cli_alert_info("Creating plot...")

    if (barplot) {
        aggregated_df <- aggregated_df |>
            mutate(loop_element = factor(loop_element, levels = loop_elements))
        p <- ggplot(
                aggregated_df, 
                aes(y = avg_score, x = loop_element, fill = dataset, 
                    ymin = avg_score - sd_score, ymax = avg_score + sd_score)) +
            geom_bar(stat = "identity", position = position_dodge(width = 0.66), width=0.66) +
            geom_errorbar(
                aes(),
                position = position_dodge(width = 0.66),
                width = 0.2) +
            scale_x_discrete(labels = function(x) {
                sapply(strsplit(x, "-"), function(words) {
                    if (length(words) == 1) {
                        words
                    } else {
                        paste0(toupper(substr(words, 1, 1)), collapse = "")
                    }
                })
            })
    } else {
        p <- ggplot(
                aggregated_df, 
                aes(y = avg_score, x = loop_element,
                    ymin = avg_score - sd_score, ymax = avg_score + sd_score)) +
            geom_ribbon(alpha = 0.33, aes(fill = dataset)) +
            geom_point(size = 0.33, aes(color = dataset)) +
            geom_line(aes(color = dataset))
    }

    p <- p +
        labs(caption = "SD and mean computed per k_fold (over all species)")

    if (group_species) {
        p <- p + facet_wrap(~ metric)
    } else {
        p <- p + facet_grid(metric ~ species, scales = "free_x")
    }

    p <- my_custom_ggplot_theme(p) + 
        xlab(xlabel) + 
        ylab(ylabel)

    
    p <- p + scale_fill_manual(values = c(PALETTE[2], PALETTE[3], PALETTE[1]))
    if (!barplot) {
        p <- p + scale_color_manual(values = c(PALETTE[2], PALETTE[3], PALETTE[1]))
    }
    
    print(p)
    if (!is.null(save_to)) {
        standardised_ggplot_save(
            figure = p, 
            save_path = file.path(parent_folder, save_to))
        cli_alert_success("Plot of performances saved!")
    }
    cli_alert_success("Plot of performances ready!\n\n")

    return(aggregated_df)
}

# A function to compare training scores between k_folds, subset and model type.
# ARGS:
#   - parent_folder: a string. The path where subfolders of each model to compare against the reference are located.
#   - reference_model_path: a string. The template path to a reference folder of the model (parent_folder goes before, k_fold number goes at the end).
#   - save_to: a string. Name of the file for saving resulting plot (should end with .pdf). If NULL, does not save pdf.
#   - loop_prefix: a string. The prefix of the subfolders names.
#   - loop_elements: a list of strings. Middle elements for subfolders names.
#   - loop_suffix: a string. The suffix of the subfolders names.
#   - k_fold: a numeric. The number of cross-validation subsets to make. Goes after suffix and reference_model_path .
#   - bars: a string. Either "CI" (default) or "SD", the type of error bars to plot.
#   - metric: a string. The metric to plot, either "MSE" (default), "RMSE", "AUC" or "TjurR2".
#   - subset_names: a list of strings. Score values are fetch from files with subset_name + "_scores.csv". Default is c("train", "val", "test").
#   - xlabel: a string. The label for the x-axis of the plot (default is "Effect").
#   - group_species: whether to take the mean accross all k_folds and species (TRUE, default) or only accross k_folds (FALSE).
#   - species_names: a list of strings. Optional vector giving the species label for each row. Default is NULL ("species_1", "species_2", ...).
dotwhisker_compare_scores <- function(
        parent_folder, 
        reference_model_path, 
        save_to,
        loop_prefix, 
        loop_elements, 
        loop_suffix, 
        k_fold, 
        bars = "CI",
        metric = "MSE",
        subset_names = c("train", "val", "test"),
        xlabel = "Effect", 
        group_species = TRUE,
        species_names = NULL) {
    
    # load reference results
    cli_alert_info("Loading base results...")  
    ref_scores <- list()
    ref_scores_k_means <- list()
    for (subset_name in subset_names) {
        ref_scores[[subset_name]] <- list()

        for (k in seq(k_fold)) {
            # setup path
            reference_model_scores_path <- file.path(
                parent_folder,
                paste0(reference_model_path, k), 
                paste0(subset_name, "_scores.csv")
            )   
            
            # load scores
            if (metric == "MSE") {
                ref_scores[[subset_name]][[k]] <- read_csv(
                    reference_model_scores_path,
                    show_col_types = FALSE)$RMSE^2
            } else if (metric %in% c("RMSE", "AUC", "TjurR2")) {
                ref_scores[[subset_name]][[k]] <- read_csv(
                    reference_model_scores_path,
                    show_col_types = FALSE)[[metric]]
            } else {
                stop(paste0("Metric '", metric, "' is not handled by this function."))
            }

        }

        # Get mean per k_fold (individual mean for each subset, model_type and species)
        ref_scores_k_means[[subset_name]] <- colMeans(
            do.call(rbind, ref_scores[[subset_name]]), na.rm = TRUE)
    }

    # load new results
    cli_alert_info("Loading other results...")
    other_scores <- list()
    other_scores_k_means <- list()
    for (subset_name in subset_names) {
        other_scores[[subset_name]] <- list()
        other_scores_k_means[[subset_name]] <- list()

        for (i in seq_along(loop_elements)) { 
            loop_element <- loop_elements[i]
            other_scores[[subset_name]][[loop_element]] <- list()
            
            for (k in seq(k_fold)) {
                other_model_scores_path <- file.path(
                    parent_folder, 
                    paste0(loop_prefix, loop_element, loop_suffix, k),
                    paste0(subset_name, "_scores.csv")
                )
                
                # load scores
                if (metric == "MSE") {
                    other_scores[[subset_name]][[loop_element]][[k]] <- read_csv(
                        other_model_scores_path,
                        show_col_types = FALSE)$RMSE^2
                } else if (metric %in% c("RMSE", "AUC", "TjurR2")) {
                    other_scores[[subset_name]][[loop_element]][[k]] <- read_csv(
                        other_model_scores_path,
                        show_col_types = FALSE)[[metric]]
                } else {
                    stop(paste0("Metric '", metric, "' is not handled by this function."))
                }
            }

            # mean per k_fold (over each subset, model_type and species)
            other_scores_k_means[[subset_name]][[loop_element]] <- colMeans(
                do.call(rbind, other_scores[[subset_name]][[loop_element]]), 
                na.rm = TRUE)
        }
    }

    # Number of species = number of rows, assign names
    if (!group_species) {
        n_species <- length(ref_scores[[subset_names[1]]][[1]])
        if (is.null(species_names)) {
            species_names <- paste0("species_", seq_len(n_species))
        } else if (length(species_names) != n_species) {
            stop(paste0(
                "`species_names` has length ", length(species_names),
                " but the score files have ", n_species, " rows (species)."
            ))
        }
    }

    # format dataframe for ggplot
    scores_df <- tibble()
    for (subset_name in subset_names) {
        for (i in seq_along(loop_elements)) {
            loop_element <- loop_elements[i]

            if (grepl("MSE", metric)) {
                score_order <- -1 # MSE is better when its low, so reverse order
            } else {
                score_order <- 1
            }

            # Build a (k_fold x species) matrix of *raw* differences before
            # any averaging, so both fold- and species-level variability are
            # still available when we compute the CI/SD below.
            diff_matrix <- do.call(rbind, lapply(seq(k_fold), function(k) {
                score_order * (
                    other_scores[[subset_name]][[loop_element]][[k]] -
                    ref_scores[[subset_name]][[k]]
                )
            }))
            # rows = folds, columns = species

            if (group_species) {
                # Average across species within each fold first, so each
                # fold contributes exactly one observation. The interval
                # then reflects fold-to-fold (across k_fold) variability.
                diff_obs_list <- list(all = rowMeans(diff_matrix, na.rm = TRUE))
                row_species <- list(all = NA_character_)
            } else {
                # Keep species separate: for each species, the observations
                # are that species' differences across the k folds. This
                # gives one CI/SD per species, computed across k_fold.
                n_species <- ncol(diff_matrix)
                diff_obs_list <- setNames(
                    lapply(seq_len(n_species), function(s) diff_matrix[, s]),
                    species_names
                )
            }

            for (obs_name in names(diff_obs_list)) {
                diff_obs <- diff_obs_list[[obs_name]]

                n <- sum(!is.na(diff_obs))
                avg_metric <- mean(diff_obs, na.rm = TRUE)
                sd_metric <- sd(diff_obs, na.rm = TRUE)

                # t critical value is more appropriate than a normal
                # approximation when n (e.g. k_fold) is small; falls back
                # to NA when there's fewer than 2 observations
                # (which should not happen in reality)
                crit_value <- if (n > 1) qt(0.975, df = n - 1) else NA_real_
                conf_margin <- crit_value * sd_metric / sqrt(n)

                scores_df <- bind_rows(
                    scores_df,
                    tibble(
                        average_metric = avg_metric,
                        sd_lower = avg_metric - sd_metric,
                        sd_upper = avg_metric + sd_metric,
                        ci_lower = avg_metric - conf_margin,
                        ci_upper = avg_metric + conf_margin,
                        n_obs = n,
                        conf_margin = conf_margin,
                        sd_metric = sd_metric,
                        subset = subset_name,
                        loop_element = loop_element,
                        species = if (group_species) NA_character_ else obs_name))
            }
        }
    }

    # add new fields
    scores_df <- scores_df |>
        mutate(loop_element = factor(loop_element, levels = loop_elements)) |>
        mutate(subset = factor(subset, levels = subset_names))

    if (!group_species) {
        scores_df <- scores_df |>
            mutate(species = factor(species, levels = species_names))
    }

    # create captions
    if (bars == "CI") {
        bar_name <- "95% CI"
    } else if (bars == "SD") {
        bar_name <- "SD"
    } else {
        stop("Unkown value for 'bars'")
    }
    if (group_species) {
        bottom_caption <- paste(bar_name, 
            "computed across k_fold (species-averaged within each fold)")
    } else {
        bottom_caption <- paste(bar_name, 
            "computed across k_fold, per species")
    }

    p <- ggplot(scores_df, 
        aes(y = average_metric, x = loop_element, color = subset)) +
        geom_point(
            stat = "identity", 
            position = position_dodge(width = 0.66), 
            size = 3) +
        labs(
            caption = bottom_caption,
            color = "Subset") +
        ylab(paste("Delta in average", metric)) +
        xlab(xlabel) +
        geom_hline(yintercept = 0, linetype = "dashed")

    if (bars == "CI") {
        p <- p + geom_errorbar(
            aes(
                ymin = ci_lower, 
                ymax = ci_upper
            ),
            position = position_dodge(width = 0.66),
            width = 0.2)
    } else if (bars == "SD") {
        p <- p + geom_errorbar(
            aes(
                ymin = sd_lower, 
                ymax = sd_upper
            ),
            position = position_dodge(width = 0.66),
            width = 0.2)        
    } else {
        stop("Unknown value for 'bars'")
    }

    # When no grouping, one subplot per species laid out horizontally.
    if (!group_species) {
        p <- p + facet_wrap(~species, nrow = 1)
    }

    p <- my_custom_ggplot_theme(p) + 
        scale_color_manual(values = c(PALETTE[2], PALETTE[3], PALETTE[1])) +
        scale_x_discrete(labels = function(x) {
            sapply(strsplit(x, "-"), function(words) {
                if (length(words) == 1) {
                    words
                } else {
                    paste0(toupper(substr(words, 1, 1)), collapse = "")
                }
            })
        })
    print(p)

    if (!is.null(save_to)) {
        standardised_ggplot_save(
            figure = p, 
            save_path = file.path(parent_folder, save_to))
        cli_alert_success("Plot of compared performances saved!")
    }
    
    cli_alert_success("Plot of compared performances ready!\n\n")

    return(scores_df)
}

# A function to compare training scores between k_folds, subset and model type.
# ARGS:
#   - parent_folder: a string. The path where subfolders of each model to compare against the reference are located.
#   - reference_model_path: a string. The template path to a reference folder of the model (parent_folder goes before, k_fold number goes at the end).
#   - save_to: a string. Name of the file for saving resulting plot (should end with .pdf). If NULL, does not save pdf.
#   - loop_prefix: a string. The prefix of the subfolders names.
#   - loop_elements: a list of strings. Middle elements for subfolders names.
#   - loop_suffix: a string. The suffix of the subfolders names.
#   - k_fold: a numeric. The number of cross-validation subsets to make. Goes after suffix and reference_model_path .
#   - metric: a string. The metric to plot, either "MSE" (default), "RMSE", "AUC" or "TjurR2".
#   - subset_names: a list of strings. Score values are fetch from files with subset_name + "_scores.csv". Default is c("train", "val", "test").
#   - xlabel: a string. The label for the x-axis of the plot (default is "Effect").
#   - group_species: whether to take the mean accross all k_folds and species (TRUE, default) or only accross k_folds (FALSE).
#   - species_names: a list of strings. Optional vector giving the species label for each row. Default is NULL ("species_1", "species_2", ...).
boxplot_compare_scores <- function(
        parent_folder, 
        reference_model_path, 
        save_to,
        loop_prefix, 
        loop_elements, 
        loop_suffix, 
        k_fold, 
        metric = "MSE",
        subset_names = c("train", "val", "test"),
        xlabel = "Effect", 
        group_species = TRUE,
        species_names = NULL) {
    
    # load reference results
    cli_alert_info("Loading base results...")  
    ref_scores <- list()
    ref_scores_k_means <- list()
    for (subset_name in subset_names) {
        ref_scores[[subset_name]] <- list()

        for (k in seq(k_fold)) {
            # setup path
            reference_model_scores_path <- file.path(
                parent_folder,
                paste0(reference_model_path, k), 
                paste0(subset_name, "_scores.csv")
            )   
            
            # load scores
            if (metric == "MSE") {
                ref_scores[[subset_name]][[k]] <- read_csv(
                    reference_model_scores_path,
                    show_col_types = FALSE)$RMSE^2
            } else if (metric %in% c("RMSE", "AUC", "TjurR2")) {
                ref_scores[[subset_name]][[k]] <- read_csv(
                    reference_model_scores_path,
                    show_col_types = FALSE)[[metric]]
            } else {
                stop(paste0("Metric '", metric, "' is not handled by this function."))
            }

        }

        # Get mean per k_fold (individual mean for each subset, model_type and species)
        ref_scores_k_means[[subset_name]] <- colMeans(
            do.call(rbind, ref_scores[[subset_name]]), na.rm = TRUE)
    }

    # load new results
    cli_alert_info("Loading other results...")
    other_scores <- list()
    other_scores_k_means <- list()
    for (subset_name in subset_names) {
        other_scores[[subset_name]] <- list()
        other_scores_k_means[[subset_name]] <- list()

        for (i in seq_along(loop_elements)) { 
            loop_element <- loop_elements[i]
            other_scores[[subset_name]][[loop_element]] <- list()
            
            for (k in seq(k_fold)) {
                other_model_scores_path <- file.path(
                    parent_folder, 
                    paste0(loop_prefix, loop_element, loop_suffix, k),
                    paste0(subset_name, "_scores.csv")
                )
                
                # load scores
                if (metric == "MSE") {
                    other_scores[[subset_name]][[loop_element]][[k]] <- read_csv(
                        other_model_scores_path,
                        show_col_types = FALSE)$RMSE^2
                } else if (metric %in% c("RMSE", "AUC", "TjurR2")) {
                    other_scores[[subset_name]][[loop_element]][[k]] <- read_csv(
                        other_model_scores_path,
                        show_col_types = FALSE)[[metric]]
                } else {
                    stop(paste0("Metric '", metric, "' is not handled by this function."))
                }
            }

            # Get mean per k_fold (individual mean for each subset, model_type and species)
            other_scores_k_means[[subset_name]][[loop_element]] <- colMeans(
                do.call(rbind, other_scores[[subset_name]][[loop_element]]), 
                na.rm = TRUE)
        }
    }

    # Now that we've loaded actual score vectors, we know how many species
    # there are (= vector length). Build/validate species_names for use
    # when group_species = FALSE.
    if (!group_species) {
        n_species <- length(ref_scores[[subset_names[1]]][[1]])
        if (is.null(species_names)) {
            species_names <- paste0("species_", seq_len(n_species))
        } else if (length(species_names) != n_species) {
            stop(paste0(
                "`species_names` has length ", length(species_names),
                " but the score files have ", n_species, " rows (species)."
            ))
        }
    }

    # format dataframe for ggplot
    scores_df <- tibble()
    raw_diffs_df <- tibble()
    for (subset_name in subset_names) {
        for (i in seq_along(loop_elements)) {
            loop_element <- loop_elements[i]

            if (grepl("MSE", metric)) {
                score_order <- -1 # MSE is better when its low, so reverse order
            } else {
                score_order <- 1
            }

            # Build a (k_fold x species) matrix of *raw* differences, before
            # any averaging, so both fold- and species-level variability are
            # still available when we compute the CI/SD below.
            diff_matrix <- do.call(rbind, lapply(seq(k_fold), function(k) {
                score_order * (
                    other_scores[[subset_name]][[loop_element]][[k]] -
                    ref_scores[[subset_name]][[k]]
                )
            }))
            # rows = folds, columns = species

            if (group_species) {
                # Average across species within each fold first, so each
                # fold contributes exactly one observation. The interval
                # then reflects fold-to-fold (across k_fold) variability.
                diff_obs_list <- list(all = rowMeans(diff_matrix, na.rm = TRUE))
            } else {
                # Keep species separate: for each species, the observations
                # are that species' differences across the k folds. This
                # gives one CI/SD (and one box) per species, computed
                # across k_fold.
                n_species <- ncol(diff_matrix)
                diff_obs_list <- setNames(
                    lapply(seq_len(n_species), function(s) diff_matrix[, s]),
                    species_names
                )
            }

            for (obs_name in names(diff_obs_list)) {
                diff_obs <- diff_obs_list[[obs_name]]

                # --- raw observations, kept for the boxplot ---
                raw_diffs_df <- bind_rows(
                    raw_diffs_df,
                    tibble(
                        diff_value = diff_obs,
                        subset = subset_name,
                        loop_element = loop_element,
                        species = if (group_species) NA_character_ else obs_name))

                # --- summary stats, kept for reference / the return value ---
                n <- sum(!is.na(diff_obs))
                avg_metric <- mean(diff_obs, na.rm = TRUE)
                sd_metric <- sd(diff_obs, na.rm = TRUE)

                # t critical value is more appropriate than a normal
                # approximation when n (e.g. k_fold) is small; falls back
                # to NA when there's fewer than 2 observations.
                crit_value <- if (n > 1) qt(0.975, df = n - 1) else NA_real_
                conf_margin <- crit_value * sd_metric / sqrt(n)

                scores_df <- bind_rows(
                    scores_df,
                    tibble(
                        average_metric = avg_metric,
                        n_obs = n,
                        subset = subset_name,
                        loop_element = loop_element,
                        species = if (group_species) NA_character_ else obs_name))
            }
        }
    }

    # add new fields
    scores_df <- scores_df |>
        mutate(loop_element = factor(loop_element, levels = loop_elements)) |>
        mutate(subset = factor(subset, levels = subset_names))
    raw_diffs_df <- raw_diffs_df |>
        mutate(loop_element = factor(loop_element, levels = loop_elements)) |>
        mutate(subset = factor(subset, levels = subset_names))

    if (!group_species) {
        scores_df <- scores_df |>
            mutate(species = factor(species, levels = species_names))
        raw_diffs_df <- raw_diffs_df |>
            mutate(species = factor(species, levels = species_names))
    }

    # create caption
    if (group_species) {
        bottom_caption <- paste(
            "Distribution of per-fold differences",
            "(species-averaged within each fold, across k_fold")
    } else {
        bottom_caption <- paste(
            "Distribution of per-fold differences across k_fold,",
            " per species")
    }

    p <- ggplot(raw_diffs_df, 
        aes(y = diff_value, x = loop_element, fill = subset)) +
        geom_boxplot(
            position = position_dodge(width = 0.75),
            width = 0.6,
            outlier.shape = 16) +
        labs(
            caption = bottom_caption,
            fill = "Subset") +
        ylab(paste("Delta in average", metric)) +
        xlab(xlabel) +
        geom_hline(yintercept = 0, linetype = "dashed")

    # One subplot per species, laid out horizontally, when species aren't
    # pooled together.
    if (!group_species) {
        p <- p + facet_wrap(~species, nrow = 1)
    }

    p <- my_custom_ggplot_theme(p) + 
        scale_fill_manual(values = c(PALETTE[2], PALETTE[3], PALETTE[1])) +
        scale_x_discrete(labels = function(x) {
            sapply(strsplit(x, "-"), function(words) {
                if (length(words) == 1) {
                    words
                } else {
                    paste0(toupper(substr(words, 1, 1)), collapse = "")
                }
            })
        })
    print(p)

    if (!is.null(save_to)) {
        standardised_ggplot_save(
            figure = p, 
            save_path = file.path(parent_folder, save_to))
        cli_alert_success("Plot of compared performances saved!")
    }
    cli_alert_success("Plot of compared performances ready!\n\n")

    return(list(summary = scores_df, raw_diffs = raw_diffs_df))
}

# A function to compute a habitat suitability map for a species based on Hmsc predictions.
# ARGS:
#   - df: a dataframe. Must contain columns listed in x_cols and y_cols.
#   - x_cols: a list of strings. The columns containing explanatory variables.
#   - parent_folder: a string. The path where subfolders of each model are located.
#   - loop_prefix: a string. The prefix of the subfolder name.
#   - loop_element: a string. Middle element for subfolder name.
#   - loop_suffix: a string. The suffix of the subfolder name.
#   - k_fold: a numeric. The number of cross-validation subsets to make. Goes after suffix.
#   - xlabel: a string. The label for the x-axis of the plot (default is "Model type").
#   - sp: a string. The name of a species (in df) to plot.
map_results_hmsc <- function(
        df,
        x_cols,
        parent_folder, 
        loop_prefix, 
        loop_element, 
        loop_suffix, 
        k_fold, 
        sp) {

    local_folder <- file.path(
            parent_folder,
            paste0(loop_prefix, loop_element, loop_suffix, k_fold))
    local_model <- readRDS(file.path(local_folder, "train_outputs.rds")) 

    cli_alert_info("Making predictions, can take a few minutes...")
    full_preds_list <- predict_hmsc(
        hM = local_model, df = df, x_variables = x_cols)
    cli_alert_info("Predictions are ready.")

    # Predictions are DISTRIBUTIONS: let's take the average for each point
    cli_alert_info("Converting distributions to point average...")
    full_preds_avg <- apply(
        simplify2array(full_preds_list), c(1,2), mean)
    full_preds_sd <- apply(
        simplify2array(full_preds_list), c(1,2), sd)
    df_preds <- df |>
        mutate(suitability = full_preds_avg[, sp]) |>
        mutate(uncertainty = full_preds_sd[, sp]) |>
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

    # # plot training samples
    # p0.1 <- ggplot_categorical_df_on_background_map(
    #     background_map = ggplot_get_france_base_map("national"), 
    #     df = df_preds |> filter(subset == "train"), 
    #     LON = "LON",
    #     LAT = "LAT",
    #     column = NULL,
    #     legend_title = "Training location")
    # print(p0.1)
    # standardised_ggplot_save(
    #     figure = p0.1, 
    #     save_path = file.path(local_folder, "training_samples.pdf"))
    # p0.2 <- ggplot_categorical_df_on_background_map(
    #     background_map = ggplot_get_france_base_map("national"), 
    #     df = df_preds |> filter(subset == "test"), 
    #     LON = "LON",
    #     LAT = "LAT",
    #     column = NULL,
    #     legend_title = "Training location")
    # print(p0.2)
    # standardised_ggplot_save(
    #     figure = p0.2, 
    #     save_path = file.path(local_folder, "test_samples.pdf"))

    # show map of suitability per site
    p1.1 <- ggplot_quantitative_df_on_background_map(
            background_map = ggplot_get_france_base_map("national"), 
            df = df_preds, 
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
            df = df_preds,  
            column = "suitability",
            res_km = RES_KM, 
            LON = "LON",
            LAT = "LAT",
            idp = 2,           
            maxdist_m = 100)
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
            df = df_preds, 
            LON = "LON",
            LAT = "LAT",
            column = "uncertainty",
            unit = paste0("Standard deviation of HSI for ", sp)) +
        labs(caption = "HSI = 'Habitat Suitability Index'")
    print(p2.1)
    standardised_ggplot_save(
        figure = p2.1, 
        save_path = file.path(local_folder, "certainty_of_suitability.pdf"))
    shp2.2 <- interpolate_scattered_points_to_hexagons(
            df = df_preds,  
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
        unit = paste0("Standard deviation of HSI for ", sp),
        limits = NULL,
        precision_auto_limits = 1e-5) +
        labs(caption = "HSI = 'Habitat Suitability Index'")
    print(p2.2)
    standardised_ggplot_save(
        figure = p2.2, 
        save_path = file.path(local_folder, "certainty_of_suitability_hexagons.pdf"))    

    return(df_preds)
}

# A function to automatise verbose interpretation of diagnostic vectors
# ARGS:
#   - vector: a vector of values to analyse.
#   - bad: a numeric. The threshold under which values reveal a bad fit.
#   - good: a numeric. The threshold over which values reveal a good fit.
#   - order: a string. Indicates if lower is better ("low_better", default) or higer is bette (high_better).
#   - mode: a sting. Indicates if showing full analysis ("full", default) or a quick 1-line summary ("quick").
interpret_diagnostics <- function(
    vector, 
    good, 
    bad = NULL, 
    mode = "full",
    order = "low_better") {
    
    if (order == "low_better") {
        n_good <- sum(vector < good)    

        if (is.null(bad)) {
            n_bad <- sum(vector > good)
            if (mode == "full") {
                cli_alert_info(paste0(
                    "Rule of thumb: <", good," (good)"))
            }

            n_acceptable <- 0
        } else {
            n_acceptable <- sum((vector > good) & (vector < bad))
            n_bad <- sum(vector > bad)    
            if (mode == "full") {
                cli_alert_info(paste0(
                    "Rule of thumb: <", good," (good), ",
                    good, "-", bad, ", (acceptable), >",
                    bad, " (bad)"))
            }
        }
        
    } else if (order == "high_better") {
        n_good <- sum(vector > good)    

        if (is.null(bad)) {
            n_bad <- sum(vector < good)  
            if (mode == "full") {
                cli_alert_info(paste0(
                    "Rule of thumb: >", good," (good)"))
            }

            n_acceptable <- 0
        } else {
            n_acceptable <- sum((vector < good) & (vector > bad))
            n_bad <- sum(vector < bad)
            if (mode == "full") {
                cli_alert_info(paste0(
                    "Rule of thumb: >", good," (good), ",
                    good, "-", bad, ", (acceptable), <",
                    bad, " (bad)"))   
            }
        }
     
    }
    else {
        stop(paste0("Mode should be one of 'low_better' or 'high_better'. ",
        "Got '", mode, "' ."))
    }
    
    if (mode == "full") {
        if (n_good > 0) {
            cli_alert_success(paste0(
            "- Number of 'good' estimates: ", n_good, 
            " (", round(100*n_good/length(vector), 2), "% of given values)"
            ))
        }
        if ((n_acceptable > 0) & !(is.null(bad))) {
            cli_alert_info(paste0(
            "- Number of 'acceptable' estimates: ", n_acceptable, 
            " (", round(100*n_acceptable/length(vector), 2), "% of given values)"
            ))  
        }
        if (n_bad > 0) {
            cli_alert_info(paste0(
            "- Number of 'bad' estimates: ", n_bad, 
            " (", round(100*n_bad/length(vector), 2), "% of given values)"
            ))
        }
    } else if (mode == "quick") {
        if (is.null(bad)) {
            cli_alert_info(paste0(
                n_good, " (", round(100*n_good/length(vector), 2), "%) are good ",
                "and ", n_bad ," (", round(100*n_bad/length(vector), 2),"%) are bad."))
        } else {
            cli_alert_info(paste0(
                n_good, " (", round(100*n_good/length(vector), 2), "%) are good, ",
                n_acceptable, " (", round(100*n_acceptable/length(vector), 2), "%) are acceptable ",
                "and ", n_bad ," (", round(100*n_bad/length(vector), 2),"%) are bad."))
        }
        
    }

}

# A function that wraps the process to make predictions with a hmsc model.
# ARGS:
#   - hM: a Hmsc fitted model object
#   - df: a dataframe with columns "point", "id_point_annee" and names in x_variables
#   - x_variables: a list of strings. The names of columns to keep in data.
predict_hmsc <- function(hM, df, x_variables) {
    # Format explanatory variables
    XData <- as.data.frame(setNames(
        lapply(x_variables, function(col) df[[col]]),
        x_variables))
    
    # Format study design (grouping of samples together)
    studyDesign <- data.frame(
        units = as.factor(df$point), 
        spatial = as.factor(df$id_point_annee)
    )

    # Make prediction on new dataset
    full_preds_list <- predict(
        hM, XData = XData, studyDesign = studyDesign, expected = TRUE)
    
    return(full_preds_list)
}

# A function that mimciks Hmsc::evaluateModelFit but can also work on non-training data.
# ARGS :
#   - hM: a Hmsc fitted model object
#   - y: a matrix of species observation ("ground truth")
#   - predY: the predictions made by the model
evaluateModelFitCustom <- function(hM, Y, predY) {

    ns <- ncol(Y)
    mPredY <- apply(predY, c(1, 2), mean)  # posterior mean prediction per obs/species

    RMSE <- rep(NA, ns)     # RMSE, the lower the better (corr with AUC)
    AUC <- rep(NA, ns)      # AUC (the closer to 1, the better)
    TjurR2 <- rep(NA, ns)   # Tjur R² (% of variance explained)

    for (j in seq_len(ns)) {
        sel <- !is.na(Y[, j])
        obs <- Y[sel, j]
        pred <- mPredY[sel, j]

        # RMSE
        RMSE[j] <- sqrt(mean((obs - pred)^2))

        # AUC (only meaningful if both 0s and 1s present)
        if (length(unique(obs)) == 2) {
            AUC[j] <- as.numeric(pROC::auc(obs, pred, quiet = TRUE))
        }

        # Tjur R2: difference in mean predicted probability between
        # presences and absences
        if (length(unique(obs)) == 2) {
            TjurR2[j] <- mean(pred[obs == 1]) - mean(pred[obs == 0])
        }
    }

    names(RMSE) <- names(AUC) <- names(TjurR2) <- colnames(Y)

    return(list(RMSE = RMSE, AUC = AUC, TjurR2 = TjurR2))
}

# A function that computes the uncertainty of a Hmsc model on its predictions
#   - hM: a Hmsc fitted model object
#   - df: a dataframe with columns "point", "id_point_annee" and names in x_variables
#   - x_variables: a list of strings. The names of columns to keep in data.
get_uncertainty_hmsc <- function(hM, df, x_cols) {
    predicted_occurrences <- predict_hmsc(
        hM = hM, df = df, x_variables = x_cols)
    # predicted_lists contains a list of length = number of samples.
    # for each sample, we get a matrix of n_obs x n_species
    
    # Here, we define uncertainty as the average of the standard deviation 
    # accross observed species:
    sd_point_sp <- apply(simplify2array(predicted_occurrences), c(1,2), sd)
    uncertainty_per_point <- as_tibble(sd_point_sp) |> 
        mutate(average_sd = rowMeans(across(everything())))
    
    return(uncertainty_per_point$average_sd)
}