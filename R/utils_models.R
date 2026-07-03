# Set of functions used to create / handle / analyse models
##### Libraries ##### ---------------------------------------------------------
library(cli)

##### Functions ##### ---------------------------------------------------------
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
    XData <- as.data.frame(setNames(
        lapply(x_variables, function(col) df[[col]]),
        x_variables))
    studyDesign <- data.frame(
        units = as.factor(df$point), 
        spatial = as.factor(df$id_point_annee)
    )
    full_preds_list <- predict(
        hM,
        XData = XData,
        studyDesign = studyDesign,
        expected = TRUE
    )
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