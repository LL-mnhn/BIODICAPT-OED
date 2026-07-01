# Set of functions used to create / handle / analyse models

##### Functions #####
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