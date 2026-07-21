# Set of functions used to plot harmonised figures

##### Liraries #####
library(colorspace)
library(tidyterra)
library(patchwork)
library(reshape2)
library(MCMCvis)
library(ggplot2)
library(GGally)
library(dplyr)
library(terra)
library(coda)
library(cli)
library(sf)

source(here::here("R/utils_data.R")) 


##### Parameters #####
source(here::here("data/config/config.R")) # all parameters are grouped together


##### Functions #####
# A wrapper to create a custom ggplot theme 
# ARGS:
#   - figure: a ggplot object.
#   - with_palette: a boolean. If TRUE (default) uses custom color/fill/shape/sizes.
my_custom_ggplot_theme <- function(figure, with_palette=TRUE, LIGHT=FALSE){
    customised_fig <- figure +
        theme_linedraw(
            base_family = FONT
        ) +
        theme(
            # title and subtitle styling
            plot.title.position = "plot",
            plot.title = element_text(
                size = 18,
                face = "bold",
                color = "#000000",  
                margin = margin(b = 10)
            ),
            plot.subtitle = element_text(
                size = 14,
                color = "#777777", 
                margin = margin(b = 10)
            ),
            
            # plot styling
            plot.caption.position = "plot",
            plot.caption = element_text(
                size = 9,
                color = "#999999", 
                margin = margin(t = 15),
                hjust = 0
            ),
            axis.text = element_text(
                size = 11,
                color = "#000000"
            ),
            aspect.ratio = 1,
            
            # external grid
            axis.ticks = element_line(
                linetype = "solid",
                linewidth = 0.50,
                color = "#000000"
            ),
            panel.border = element_rect(
                colour = "#000000",
                linewidth = 1,
                fill = NA
            ),
            
            # internal grid
            panel.grid.major = element_line(
                linetype = "solid",
                linewidth = 0.15,
                color = "#999999"
            ),
            panel.grid.minor = element_blank(),
        )

        if (with_palette){
            if (LIGHT) {
                return(customised_fig + LIGHT_CUSTOM_SCALES)
            } else {
                return(customised_fig + CUSTOM_SCALES)
            }
        } else {
            return(customised_fig)
        }
}

# A function to save a ggplot figure to pdf
# ARGS:
#   - figure: a ggplot object.
#   - save_path: a filepath to create a pdf file.
standardised_ggplot_save <- function(figure, save_path, .width = 18, .height = 6){
    # check if string ends with ".pdf"
    if (!endsWith(save_path, ".pdf")){
        stop(paste("Provided save_path must end with '.pdf', got", save_path))
    }
  
    ggsave(
        filename = save_path,
        plot = figure,
        dpi = 300,
        width = .width,             # large width to account for plots with very wide legends
        height = .height,
        device = cairo_pdf)
}

# A function that creates a simple background map of france in ggplot2
# ARGS:
#   - borders_type: a string. Either "national" (default) or "regional". if "regional", draws highest level inner borders ofthe country.
ggplot_get_france_base_map <- function(borders_type="national"){
    europe_shp = rnaturalearth::ne_countries(continent = "Europe", scale = "large", returnclass = "sf")
    france_shp = get_metropolitan_france_shapefile(borders_type)
    
    base_map <- ggplot(europe_shp) +
        geom_sf(fill = "grey80", color = "white") +                     # color of countries
        geom_sf(data = france_shp, fill = "white", color = "black") +   # color of France
        theme_minimal() +
        theme(
            panel.background = element_rect(fill = "lightcyan1", color = NA),
            panel.grid.major = element_line(color = "lightcyan1"),
            panel.grid.minor = element_line(color = "lightcyan1")
        ) +
        labs(
            x = "longitude",
            y = "latitude"
        ) +
        coord_sf(
            xlim = c(LON_MIN, LON_MAX), 
            ylim = c(LAT_MIN, LAT_MAX))
    
    return(my_custom_ggplot_theme(base_map, with_palette=FALSE))
}

# A function to create a ggplot that shows scattered locations on a map
# ARGS:
#   - background_map: a ggplot object. The background map that will be used.
#   - df: a dataframe.
#   - LON: a string. The name of the column with longitude values.
#   - LAT: a string. The name of the column with latitude values.
#   - column: a string (optional). The name of a column with categorical values.
ggplot_categorical_df_on_background_map <- function(
        background_map, 
        df, 
        LON = "LON",
        LAT = "LAT",
        column = NULL,
        legend_title = "Sampling location") {
       
    # make sure that coordinates are in the right coordinates system
    data <- st_as_sf(df, coords = c(LON, LAT), crs = 4326)
    
    # shuffle to avoid biased overlaps
    data <- slice_sample(data, prop = 1)
    
    if (!is.null(column)) {
        uniques_vals <- unique(data[[column]])
        if (length(uniques_vals) > 6){
            stop(paste(c(
                "This function can handle up to 6 unique values, found",
                length(uniques_vals),
                "in column",
                column,
                ". (tbh: sorry I was just being lazy at the time. If needed, append the CUSTOM_SCALES in data/config/config.R)."
            )))
        } else {
            # make plot
            map_obs <- suppressMessages(background_map +
                geom_sf(
                    data = data,
                    stroke = 0.8,
                    aes(
                        color = .data[[column]],
                        fill = .data[[column]],
                        shape = .data[[column]],
                        size = .data[[column]],
                    )
                ) +
                labs(
                    x = "longitude",
                    y = "latitude",
                    color = legend_title,
                    fill = legend_title,
                    shape = legend_title,
                    size = legend_title
                ) +
                coord_sf(
                    xlim = c(LON_MIN, LON_MAX),
                    ylim = c(LAT_MIN, LAT_MAX))
            )
        }
        return(my_custom_ggplot_theme(map_obs, with_palette = TRUE))
        
    } else {
        # make plot
        map_obs <- suppressMessages(background_map +
            geom_sf(
                data = data,
                size = SIZES[1],
                shape = SHAPES[1],        
                fill = PALETTE[1], 
                color = darken(PALETTE[1], amount = 0.66), 
                stroke = 0.8
            ) +
            labs(
                x = "longitude",
                y = "latitude"
            ) +
            coord_sf(
                xlim = c(LON_MIN, LON_MAX),
                ylim = c(LAT_MIN, LAT_MAX))
        )
        return(my_custom_ggplot_theme(map_obs, with_palette = FALSE))
    }
}

# A function to create a ggplot that shows scattered locations on a map
# ARGS:
#   - background_map: a ggplot object. The background map that will be used.
#   - df: a dataframe.
#   - LON: a string. The name of the column with longitude values.
#   - LAT: a string. The name of the column with latitude values.
#   - column: a string (optional). The name of a column with quantitative values.
#   - unit: a string. A label that will be shown along the palette displayed.
#   - limits: a vector of 2 values (optional). Imposes hard limits on the values considered by the palette.
#   - precision_auto_limits: when limits is NULL, precision of color scale (values are rounded to closest precision_auto_limits)
ggplot_quantitative_df_on_background_map <- function(
        background_map, 
        df, 
        LON = "LON",
        LAT = "LAT",
        column = NULL,
        unit = NULL,
        limits = NULL,
        precision_auto_limits = 1e-5) {
       
    # make sure that coordinates are in the right coordinates system
    data <- st_as_sf(df, coords = c(LON, LAT), crs = 4326)
    
    # shuffle to avoid biased overlaps
    data <- slice_sample(data, prop = 1)

   
    if (!is.null(column)) {
        if (is.vector(limits) && length(limits) == 2){
            low_limit <- limits[1]
            high_limit <- limits[2]
        } else if (is.null(limits)) {
            # round palette scale to the bottom and top nearest multiple of 5
            low_limit <- floor(min(df[[column]], na.rm = TRUE) / precision_auto_limits) * precision_auto_limits
            high_limit <- ceiling(max(df[[column]], na.rm = TRUE) / precision_auto_limits) * precision_auto_limits
        } else {
            stop(paste("'limits' is not recognised. Expected vector of length 2 or NULL, got", limits))
        }

        
        # make plot
        map_obs <- suppressMessages(background_map +
            geom_sf(
                data = data,
                stroke = 0.8,
                aes( 
                    fill = .data[[column]], 
                    color = .data[[column]]
                )
            ) +
            scale_fill_continuous(
                na.value = "transparent", 
                palette = "turbo",
                limits = c(low = low_limit, high = high_limit)) +
            scale_color_continuous(
                na.value = "transparent", 
                palette = "turbo",
                limits = c(low = low_limit, high = high_limit)) +
            labs(
                x = "longitude",
                y = "latitude",
                fill = unit,
                color = unit) +
            coord_sf(
                xlim = c(LON_MIN, LON_MAX),
                ylim = c(LAT_MIN, LAT_MAX))
        )
        
        return(my_custom_ggplot_theme(map_obs, with_palette = FALSE))
        
    } else {
        # make plot
        map_obs <- suppressMessages(background_map +
            geom_sf(
                data = data,
                size = SIZES[1],
                shape = SHAPES[1],        
                fill = PALETTE[1], 
                color = darken(PALETTE[1], amount = 0.66), 
                stroke = 0.8
            ) +
            labs(
                x = "longitude",
                y = "latitude"
            ) +
            coord_sf(
                xlim = c(LON_MIN, LON_MAX),
                ylim = c(LAT_MIN, LAT_MAX))
        )
        return(my_custom_ggplot_theme(map_obs, with_palette = FALSE))
    }
}

# A function to create a ggplot that shows a raster of continuous values on a map
# ARGS:
#   - background_map: a ggplot object. The background map that will be used.
#   - raster: a raster or the path to a raster.
#   - layer_name: a string. The name of the layer containing the values to show.
#   - unit: a string. A label that will be shown along the palette displayed.
#   - limits: a vector of 2 values (optional). Imposes hard limits on the values considered by the palette.
#   - precision_auto_limits: when limits is NULL, precision of color scale (values are rounded to closest precision_auto_limits)
ggplot_quantitative_raster_on_background_map <- function(
        background_map, 
        raster,
        layer_name,
        unit = "°C",
        limits = NULL,
        precision_auto_limits = 1e-5){
    if (class(raster)[1] == "character") {
        # convert to dataframe for ggplot2
        raster <- rast(raster)  
    } else if (class(raster)[1] != "SpatRaster") {
        stop(paste("Was expecting a string or SpatRaster object, got", class(raster)))
    }

    raw_df <- as.data.frame(raster, xy = TRUE)
    
    if (is.vector(limits) && length(limits) == 2){
        low_limit <- limits[1]
        high_limit <- limits[2]
    } else if (is.null(limits)) {
        # round palette scale to the bottom and top nearest multiple of 5
        low_limit <- floor(min(raw_df[[layer_name]], na.rm = TRUE) / precision_auto_limits) * precision_auto_limits
        high_limit <- ceiling(max(raw_df[[layer_name]], na.rm = TRUE) / precision_auto_limits) * precision_auto_limits
    } else {
        stop(paste("'limits' is not recognised. Expected vector of length 2 or NULL, got", limits))
    }

    map_quantity_grid <- suppressMessages(background_map +
        geom_raster(
            data = raw_df, 
            alpha = 0.8,
            aes(x = x, y = y, fill = .data[[layer_name]])) +
        scale_fill_continuous(
            na.value = "transparent", 
            palette = "turbo",
            limits = c(low = low_limit, high = high_limit)) +
        labs(x = "longitude", y = "latitude", fill = unit) +
        coord_sf(
            xlim = c(LON_MIN, LON_MAX), 
            ylim = c(LAT_MIN, LAT_MAX))
    )
    
    return(my_custom_ggplot_theme(map_quantity_grid, with_palette = FALSE))
}

# A function to create a ggplot that shows a shapefile of continuous values on a map
# ARGS:
#   - background_map: a ggplot object. The background map that will be used.
#   - shapefile: a shapefile.
#   - layer_name: a string. The name of the layer containing the values to show.
#   - unit: a string. A label that will be shown along the palette displayed.
#   - limits: a vector of 2 values (optional). Imposes hard limits on the values considered by the palette.
#   - precision_auto_limits: when limits is NULL, precision of color scale (values are rounded to closest precision_auto_limits)
ggplot_quantitative_shapefile_on_background_map <- function(
    background_map,
    shapefile,
    layer_name,
    unit="°C",
    limits=NULL,
    precision_auto_limits = 1e-5) {

    if (is.vector(limits) && length(limits) == 2){
        low_limit <- limits[1]
        high_limit <- limits[2]
    } else if (is.null(limits)) {
        # round palette scale to the bottom and top nearest multiple of 5
        low_limit <- floor(min(shapefile[[layer_name]], na.rm = TRUE) / precision_auto_limits) * precision_auto_limits
        high_limit <- ceiling(max(shapefile[[layer_name]], na.rm = TRUE) / precision_auto_limits) * precision_auto_limits
    } else {
        stop(paste("'limits' is not recognised. Expected vector of length 2 or NULL, got", limits))
    }

    map_quantity_grid <- suppressMessages(background_map +
        geom_sf(
            data = shapefile,
            alpha = 0.8,
            color = NA,
            aes(fill = .data[[layer_name]])) +
        scale_fill_continuous(
            na.value = "transparent", 
            palette = "turbo",
            limits = c(low = low_limit, high = high_limit)) +
        scale_color_continuous(
            na.value = "transparent", 
            palette = "turbo",
            limits = c(low = low_limit, high = high_limit)) +
        labs(x = "longitude", y = "latitude", fill = unit) +
        coord_sf(
            xlim = c(LON_MIN, LON_MAX), 
            ylim = c(LAT_MIN, LAT_MAX))
    )

    return(my_custom_ggplot_theme(map_quantity_grid, with_palette = FALSE))
}

# A function to create a ggplot that shows a raster of categorical values on a map
# ARGS:
#   - background_map: a ggplot object. The background map that will be used.
#   - raster: a raster or the path to a raster.
#   - layer_name: a string. The name of the layer containing the values to show.
ggplot_categorical_raster_on_background_map <- function(
        background_map, 
        raster, 
        layer_name) {
    
    if (class(raster)[1] == "character") {
        # get raster data
        raster <- rast(raster)  
    } else if (class(raster)[1] != "SpatRaster") {
        stop(paste("Was expecting a string or SpatRaster object, got", class(raster)))
    }
        
    # Convert to dataframe for ggplot2
    df <- as.data.frame(raster, xy = TRUE)

    # build colors from color table
    coltab <- coltab(raster)[[1]]
    colors_hex <- rgb(
        coltab$red,
        coltab$green,
        coltab$blue,
        alpha = coltab$alpha,
        maxColorValue = 255)
    lvls <- levels(raster)[[1]] 

    # Match colors to labels
    label_colors <- setNames(
        colors_hex[match(lvls$Value, coltab$value)], # one is maj, not the other, idk why
        lvls[[layer_name]]
    )

    # make ggplot
    map_category_grid <- suppressMessages(background_map +
        geom_raster(
            data = df, 
            alpha = 0.8,
            aes(x = x, y = y, fill = .data[[layer_name]])) +
        scale_fill_manual(
            values = label_colors,
            na.value = "transparent"
        ) +
        labs(
            x = "longitude", 
            y = "latitude", 
            fill = "Land Cover") +
        coord_sf(
            xlim = c(LON_MIN, LON_MAX), 
            ylim = c(LAT_MIN, LAT_MAX))
    )

    return(my_custom_ggplot_theme(map_category_grid, with_palette = FALSE))
}

# A function to create a ggplot that shows a shapefile of categorical values on a map
# ARGS:
#   - background_map: a ggplot object. The background map that will be used.
#   - shapefile: a shapefile.
#   - layer_name: a string. The name of the layer containing the values to show.
#   - color_df: a dataframe. Contains columns "Value", "hex" 
ggplot_categorical_shapefile_on_background_map <- function(
    background_map,
    shapefile,
    layer_name,
    color_df,
    label_layer_name = "NEW_LABEL3") {

    shapefile[[layer_name]] <- as.factor(shapefile[[layer_name]])
    color_key <- setNames(color_df$hex, color_df$Value)
    label_key <- setNames(color_df[[label_layer_name]], color_df$Value)
    
    # make ggplot
    map_category_grid <- suppressMessages(background_map +
        geom_sf(
            data = shapefile, 
            alpha = 0.8,
            color = NA,
            aes(fill = .data[[layer_name]])) +
        scale_fill_manual(
            values = color_key,
            labels = label_key,
            na.value = "transparent") +
        scale_color_manual(
            values = color_key,
            labels = label_key,
            na.value = "transparent") +
        labs(
            x = "longitude", 
            y = "latitude", 
            fill = "Land Cover") +
        coord_sf(
            xlim = c(LON_MIN, LON_MAX), 
            ylim = c(LAT_MIN, LAT_MAX))
    )
    
    return(my_custom_ggplot_theme(map_category_grid, with_palette = FALSE))
}

# A function to plot violin plots in a standardised way
# ARGS:
#   - df: a dataframe.
#   - values: a string. The name of the column with values for y-axis
#   - categories: a string. The name of the column with factors for x-axis
#   - xlab: a string (optional, default is NULL). Label for x-axis.
#   - ylab: a string (optional, default is NULL). Label for y-axis.
#   - order: a vector. Contains all unique levels of df[[categories]] in a custom order.
#   - violin_fill: a color (default is PALETTE[4]). The background color for the violins
#   - box_width: a float (default is 0.05). Width of the boxplot inside the violin.
ggplot_violin_box_plot <- function(
    df,
    values,
    categories,
    xlab = NULL,
    ylab = NULL,
    order = NULL,
    violin_fill = PALETTE[4],
    box_width = 0.05) {
    if (is.null(order)) {
        plot <- ggplot(
            df, 
            aes(y=.data[[values]], x=.data[[categories]]))
    } else {
        plot <- ggplot(
            df, 
            aes(y=.data[[values]], x=factor(.data[[categories]], level = order)))  
    }

    plot <- plot + 
        geom_violin(trim = FALSE, fill = violin_fill) +
        geom_boxplot(width = box_width, color = "black", fill = "white") + 
        labs(x = xlab, y = ylab)
    plot <- my_custom_ggplot_theme(plot, with_palette = TRUE)

    return(plot)
}

# A function to summarise the significance of a p-value
# ARGS
#   - p: a numeric, the p-value to evaluate.
sig_stars <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) "***"
    else if (p < 0.01) "**"
    else if (p < 0.05) "*"
    else if (p < 0.1)  "·"
    else ""
}

# A function to choose which correlation to apply to a association of variables
# when using GGally::ggpairs
# ARGS: what's expected by GGpairs
assoc_fun <- function(data, mapping, ...) {
    x <- eval_data_col(data, mapping$x)
    y <- eval_data_col(data, mapping$y) 

    if (is.numeric(x) && is.numeric(y)) {
        # both continuous -> Pearson correlation
        test <- cor.test(x, y, method = "pearson")
        r <- test$estimate
        stars <- sig_stars(test$p.value)
        lbl <- paste0("Corr:\nr = ", round(r, 2), stars)

    } else if (is.numeric(x) != is.numeric(y)) {
        # one continuous, one categorical -> correlation ratio (eta)
        if (is.numeric(x)) { 
            num_var <- x; cat_var <- factor(y) 
        } else { 
            num_var <- y; cat_var <- factor(x) 
        }

        fit <- aov(num_var ~ cat_var)
        fit_summary <- summary(fit)[[1]]
        eta2 <- fit_summary[1, "Sum Sq"] / sum(fit_summary[, "Sum Sq"])
        p_val <- fit_summary[1, "Pr(>F)"]
        stars <- sig_stars(p_val)
        lbl <- paste0("Corr ratio:\nη = ", round(sqrt(eta2), 2), stars)

    } else {
        # both categorical -> Cramer's V
        tbl <- table(x, y)
        chi <- suppressWarnings(chisq.test(tbl))
        n   <- sum(tbl)
        V   <- sqrt((chi$statistic / n) / min(nrow(tbl) - 1, ncol(tbl) - 1))
        stars <- sig_stars(chi$p.value)
        lbl <- paste0("Cramer:\nV = ", round(V, 2), stars)
    }

    ggally_text(label = lbl, mapping = aes(), color = "black", ...) +
    theme_void()
}

# A function to add a ellipse and loess to plotted points in ggpair diag
# ARGS: what's expected by GGpairs
lower_cont_fun <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
        geom_point(alpha = 0.4, color = "grey20") +
        stat_ellipse(color = PALETTE[2], type = "norm", level = 0.95, linewidth = 0.6) +
        geom_smooth(method = "loess", color = PALETTE[1], se = FALSE, linewidth = 0.6, ...)
}

# A function that create's a draftman's plot.
# ARGS:
#   - df: a dataframe.
#   - columns: a vector. Contains the names of the columns to plot against each other.
ggplot_custom_draftman <- function(df, columns) {
    plot <- ggpairs(df,
        columns = columns,
        upper = list(
            continuous = assoc_fun, 
            combo = assoc_fun, 
            discrete = assoc_fun),
        lower = list(
            continuous = lower_cont_fun, 
            combo = "box_no_facet", 
            discrete = "count")
    )

    plot <- my_custom_ggplot_theme(plot) +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank())
    return(plot)
}

# A function that creates basic histograms with ggplots
# ARGS:
#   - df: a data.frame with columns x and category.
#   - x: a string. The name of a column to compute histograms on.
#   - category: a string. The name of a column to color the histograms with (default is NULL, no color added).
#   - bins: A numeric. Controls the number of bins (default is 1).
#   - breaks: A numeric. Default is NULL (not toggled) if given, overwrites bins.
ggplot_bars <- function(df, x, category = NULL, bins = 10, breaks = NULL, underlayers = NULL) {
  if (is.null(category)) {
    graph <- ggplot(data = df, aes(x = .data[[x]])) +
      underlayers +                          # drawn first, behind bars
      geom_histogram(
        color = darken(PALETTE[4], amount = 0.5),
        fill = PALETTE[4],
        bins = if (is.null(breaks)) bins else NULL,
        breaks = breaks)
  } else {
    graph <- ggplot(data = df, aes(x = .data[[x]], color = .data[[category]], fill = .data[[category]])) +
      underlayers +                          # drawn first, behind bars
      geom_histogram(binwidth = binwidth, position = "identity", alpha = 0.5)
  }
  return(my_custom_ggplot_theme(graph, with_palette = TRUE))
}

# A function that mimicks MCMCvis::MCMCtrace
# ARGS:
#   - mcmc_chains: self-explanatory. (e.g. with hmsc: transform to conda then run `ggplot_custom_MCMCtrace(coda_hmsc$Beta)`).
#   - show_Rhat: a boolean. If TRUE (default) computes and displays Rhat at the bottom of the plot. It is a value used for convergence diagnostic (usually, we consider that the chains converged if Rhat < 1.05, a value of Rhat > 1.1 should be concerning).
#   - show_Neff: a boolean. If TRUE (default) computes and displays Neff at the bottom of the plot. It is the "effective sample size", or the number of independent samples (in practice, Neff > 400 is considered acceptable).
ggplot_custom_MCMCtrace <- function(
        coda_object,
        show_Rhat = TRUE,
        show_Neff = TRUE) {
    
    if ((show_Rhat) ||(show_Neff)) {
        model_summary <- MCMCvis::MCMCsummary(
            coda_object, 
            round = 2,
            ISB = FALSE
        )
    }
    
    # add iteration number and chain ID
    chain_list <- list()
    for (i in 1:nchain(coda_object)) {
        chain_df <- as.data.frame(coda_object[[i]])
        chain_df$Iteration <- 1:nrow(chain_df)
        chain_df$Chains <- paste("chain", i)
        chain_list[[i]] <- chain_df
    }
    params <- names(as.data.frame(coda_object[[i]]))

    # bind chains into a single table
    chains_combined <- bind_rows(chain_list)

    list_plots <- list()
    for (i in 1:length(params)) {
        local_param <- params[i]
        local_df <- chains_combined[c(local_param, "Iteration", "Chains")]
        
        trace_plot <- ggplot(
            data = local_df,
            aes(
                x = .data[["Iteration"]], 
                y = .data[[local_param]], 
                color = .data[["Chains"]])) +
            geom_line(size = 0.25, alpha = 0.5)
        trace_plot <- my_custom_ggplot_theme(trace_plot, with_palette = TRUE, LIGHT = TRUE)
        trace_plot <- trace_plot + 
            theme(legend.position = "none")
        
        if ((show_Rhat) ||(show_Neff)) {
            caption_info <- ""
            
            if (show_Neff) {
                caption_info <- paste0(caption_info, "Neff: ", model_summary[local_param, "n.eff"], ". ")
            }
            if (show_Rhat) {
                caption_info <- paste0(caption_info, "Rhat: ", model_summary[local_param, "Rhat"], ". ")
            }
            trace_plot <- trace_plot + labs(subtitle = caption_info)
        }
        
        # Create density plot
        density_plot <- ggplot(
            data = local_df,
            aes(
                x = .data[[local_param]], 
                color = .data[["Chains"]])) +
            geom_density(fill = NA, size = 0.5) +
            coord_flip() 
        density_plot <- my_custom_ggplot_theme(density_plot, with_palette = TRUE, LIGHT = TRUE)
        density_plot <- density_plot + labs(x = NULL, y = "Density")
        
        if ((show_Rhat) ||(show_Neff)) {
            density_plot <- density_plot + labs(caption = "")
        }
        
        
        # Combine plots side by side
        list_plots[[i]] <- list(
            trace = trace_plot + theme(aspect.ratio = 0.75), 
            density = density_plot + theme(aspect.ratio = 0.75))
    }
    
    return(list_plots)
}

# A function that mimicks Hmsc::plotBeta
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - post: post posterior summary of Beta parameters obtained from getPostEstimate().
#   - supportLevel: a numeric threshold for plotting, values between 0.5 and 1 (default 0.95).
ggplot_custom_plotBeta <- function(hM, post, supportLevel = 0.95) {
  
    # Reproduce the Support calculation from source
    betaP  <- post$support
    toPlot <- 2 * betaP - 1
    toPlot <- toPlot * ((betaP > supportLevel) + (betaP < (1 - supportLevel)) > 0)
    betaMat <- matrix(toPlot, nrow = hM$nc, ncol = ncol(hM$Y))

    rownames(betaMat) <- hM$covNames
    colnames(betaMat) <- hM$spNames

    # Long format for ggplot
    df <- as.data.frame(as.table(betaMat))
    colnames(df) <- c("Covariate", "Species", "value")

    plot <- ggplot(df, aes(x = Covariate, y = Species, fill = value)) +
        geom_tile(color = "grey90") +
        scale_fill_gradient2(
            low     = PALETTE[3],
            mid     = "white",
            high    = PALETTE[1],
            midpoint = 0,
            limits  = c(-1, 1),
            name    = "Support"
        ) +
        theme_minimal() +
        theme(
            axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
            axis.text.y  = element_text(face = "italic"),
            panel.grid   = element_blank()
        ) +
        labs(
            x = NULL, y = NULL, 
            subtitle = paste0("Showing support levels >=", supportLevel, ".")) 

    return(my_custom_ggplot_theme(plot, with_palette = FALSE) +
        theme(axis.text.x = element_text(angle = 45, hjust=1)))
}

# A function that mimicks Hmsc::computeAssociations + Corplot
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - supportLevel: a numeric threshold for plotting, values between 0.5 and 1 (default 0.95).
ggplot_custom_random_corr_associations <- function(hM, supportLevel = 0.95) {
    OmegaCor = computeAssociations(hM)
    toPlot = ((OmegaCor[[1]]$support>supportLevel)
        + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean

    # Convert matrix to long format for ggplot
    toPlot_df <- melt(toPlot)
    colnames(toPlot_df) <- c("Var1", "Var2", "value")

    plot <- ggplot(toPlot_df, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile(color = "white", linewidth = 0.3) +
        scale_fill_gradient2(
            low  = "blue",
            mid  = "white",
            high = "red",
            midpoint = 0,
            limits = c(-1, 1),
            name = "Correlation"
        ) +
        scale_y_discrete(limits = rev(levels(factor(toPlot_df$Var2)))) + 
        labs(subtitle = paste("random effect level:", fitted_model$rLNames[1])) +
        theme_minimal() +
        theme(
            axis.text.x  = element_text(angle = 45, hjust = 1),
            axis.text.y  = element_text(size = 8),
            axis.title   = element_blank(),
            plot.title   = element_text(hjust = 0.5),
            panel.grid   = element_blank()
        ) +
        coord_fixed()

    return(my_custom_ggplot_theme(plot, with_palette = FALSE) +
        theme(axis.text.x  = element_text(angle = 45, hjust = 1)))
}

# A function that mimicks Hmsc::plotVariancePartitioning
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - VP: a matrix obtained from Hmsc::computeVariancePartitioning.
ggplot_custom_plotVariancePartitioning <- function(hM, VP) {

    # Build labels with means
    if (!(length(fitted_model$rLNames) == 0)) {
        labels <- c(VP$groupnames, paste0("Random: ", hM$rLNames))
    } else {
        labels <- c(VP$groupnames)
    }
    
    means <- round(100 * rowMeans(VP$vals), 1)
    labels <- paste0(labels, " (mean = ", means, ")")

    # Long format
    df <- as.data.frame(VP$vals)
    df$Group <- factor(labels, levels = rev(labels))
    df_long <- pivot_longer(df, -Group, names_to = "Species", values_to = "Proportion")

    plot <- ggplot(df_long, aes(x = Species, y = Proportion, fill = Group)) +
        geom_col() +
        labs(x = "Species", y = "Variance proportion") 

    return(my_custom_ggplot_theme(plot, with_palette = TRUE) +
        theme(axis.text.x  = element_text(angle = 45, hjust = 1)))
}

# A function that summarises the contents of a dataset with X variables and 
# Y species. Saves a dataframe of occurences and a draftman's plot in PDF.
#   - df: a data.frame 
#   - x_cols: a list of strings. The columns containing explanatory variables.
#   - sp_cols: a list of strings. The columns containing species occurrences.
#   - save_folder: a string. Path to a folder where "draftman_plot.pdf" will be saved.
#   - top: a numeric (default is 5). Controls the number of species to show as most and least represented in dataset
explore_dataset <- function(df, x_cols, sp_cols, save_folder, top = 5) {
    ### PLOT
    cli_alert_info("Exploration of dataset.")

    # # for fastest computation of figure: average per square ("carre")
    # draftman_df <- df |>
    #     select(c("carre", all_of(x_cols))) |>
    #     select(where(is.numeric)) |>
    #     group_by(carre) |>
    #     summarise_all(mean)
    # draft_plot <- ggplot_custom_draftman(
    #     draftman_df, 
    #     columns = setdiff(names(draftman_df), "carre"))
    
    # All points individually (takes longer to run)
    draft_plot <- ggplot_custom_draftman(df, columns = x_cols)
    suppressMessages(print(draft_plot))
    standardised_ggplot_save(
        figure = draft_plot, 
        save_path = file.path(save_folder, "draftman_plot.pdf"),
        .width = 36, .height = 12)
    cli_alert_success("Saved Draftman's plot.")

    # ### FACTORS
    # factor_df <- df |>
    #     select(where(is.factor))

    # for (col in names(factor_df)) {
    #     # for each column, compute frequence and proportion of each variable
    #     freq <- table(factor_df[[col]])
    #     prop <- round(freq/sum(freq), digits=3)
    #     combined_table <- rbind(freq = as.vector(freq), 
    #                         prop = as.vector(prop))

    #     colnames(combined_table) <- names(freq)

    #     cli_alert_info(paste0("Distribution of '",  col ,"' in dataset:"))
    #     print(t(combined_table))
    #     cli_alert_info(" ")
    # }

    ### PREDICTED VARIABLES
    # Make a table with the number of occurences
    top_sp <- df[if (!is.null(sp_cols)) sp_cols else NAMES_SPECIES] |>
        summarise(across(where(is.numeric), sum)) |>
        pivot_longer(everything(), names_to = "column", values_to = "sum")

    # print top-x species and bottom-x species
    cli_alert_info(paste0(
        "Top-", top, " MOST sighted species (from occurences):"))
    print(as.data.frame(top_sp |> slice_max(sum, n = top)))
    cli_alert_info(" ")
    cli_alert_info(paste0(
        "Top-", top," LEAST sighted species (from occurences):"))
    print(as.data.frame(top_sp |> slice_min(sum, n = top)))
    write_csv(
            as.data.frame(top_sp), 
            file.path(save_folder, "occurences_in_full_dataset.csv"))
    cli_alert_success("Saved table of occurrences in dataset.")
    cat("\n")

}