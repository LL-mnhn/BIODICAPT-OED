# Set of functions used to plot harmonised figures

##### Liraries #####
library(tidyterra)
library(terra)
library(sf)
library(ggplot2)
library(colorspace)
source(here::here("R/utils_data.R")) 

##### Parameters #####
source(here::here("data/config/config.R")) # all parameters are grouped together


##### Functions #####
# A wrapper to create a custom ggplot theme 
# ARGS:
#   - figure: a ggplot object.
#   - with_palette: a boolean. If TRUE (default) uses custom color/fill/shape/sizes.
my_custom_ggplot_theme <- function(figure, with_palette=TRUE){
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
            return(customised_fig + CUSTOM_SCALES)
        } else {
            return(customised_fig)
        }
}

# A function to save a ggplot figure to pdf
# ARGS:
#   - figure: a ggplot object.
#   - save_path: a filepath to create a pdf file.
standardised_ggplot_save <- function(figure, save_path){
    # check if string ends with ".pdf"
    if (!endsWith(save_path, ".pdf")){
        stop(paste("Provided save_path must end with '.pdf', got", save_path))
    }
  
    ggsave(
        filename = save_path,
        plot = figure,
        dpi = 300,
        width = 18,             # large width to account for plots with very wide legends
        height = 6,
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
            map_obs <- background_map +
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
                    ylim = c(LAT_MIN, LAT_MAX)
                )
        }
        return(my_custom_ggplot_theme(map_obs, with_palette = TRUE))
        
    } else {
        # make plot
        map_obs <- background_map +
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
                ylim = c(LAT_MIN, LAT_MAX)
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
ggplot_quantitative_df_on_background_map <- function(
        background_map, 
        df, 
        LON = "LON",
        LAT = "LAT",
        column = NULL,
        unit = NULL) {
       
    # make sure that coordinates are in the right coordinates system
    data <- st_as_sf(df, coords = c(LON, LAT), crs = 4326)
    
    # shuffle to avoid biased overlaps
    data <- slice_sample(data, prop = 1)

   
    if (!is.null(column)) {
        # round palette scale to the bottom and top nearest multiple of 5
        low_limit <- floor(min(data[[column]], na.rm = TRUE) / 5) * 5
        high_limit <- ceiling(max(data[[column]], na.rm = TRUE) / 5) * 5
        
        # make plot
        map_obs <- background_map +
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
        
        return(my_custom_ggplot_theme(map_obs, with_palette = FALSE))
        
    } else {
        # make plot
        map_obs <- background_map +
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
                ylim = c(LAT_MIN, LAT_MAX)
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
#   - limits: a vector of 2 values (optional). Imposes hard limits on the values considered by the palette. If not given, takes values that are multiples of 5.
ggplot_quantitative_raster_on_background_map <- function(
        background_map, 
        raster,
        layer_name,
        unit="°C",
        limits=NULL,
        precision_auto_limits=1){
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

    map_quantity_grid <- background_map +
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
    
    return(my_custom_ggplot_theme(map_quantity_grid, with_palette = FALSE))
}

# A function to create a ggplot that shows a shapefile of continuous values on a map
# ARGS:
#   - background_map: a ggplot object. The background map that will be used.
#   - shapefile: a shapefile.
#   - layer_name: a string. The name of the layer containing the values to show.
#   - color_df: a dataframe. Contains columns "Value", "hex" 
ggplot_quantitative_shapefile_on_background_map <- function(
    background_map,
    shapefile,
    layer_name,
    unit="°C",
    limits=NULL,
    precision_auto_limits = 1) {

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

    map_quantity_grid <- background_map +
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
    map_category_grid <- background_map +
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
    map_category_grid <- background_map +
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

    return(my_custom_ggplot_theme(map_category_grid, with_palette = FALSE))
}
