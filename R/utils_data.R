# Set of utility functions used accross scripts.

##### Liraries #####
library(rnaturalearth)
library(exactextractr)
library(tidyterra)
library(stringr)
library(dggridR)
library(readxl)
library(readr)
library(dplyr)
library(gstat)
library(terra)
library(purrr)
library(tools)
library(cli)
library(sf)


##### Parameters #####
source(here::here("data/config/config.R")) # all parameters are grouped together


##### Functions #####
# A function to ask for a inputs by a user, which works in interactive and batch modes
# ARGS:
#   - prompt: the message to display before asking for input
typeline <- function(prompt = "Enter text: ") {
    if (interactive() ) {
        txt <- readline(prompt)
    } else {
        cat(prompt)
        txt <- readLines("stdin", n=1)
    }
    return(txt)
}

# A function to ask the user if he wants to overwrite a file.
# ARGS:
#   - path: a path (folder or file) to check.
#   - verbose: a boolean. If TRUE, shows info messages.
authorize_overwrite <- function(path) {
    # if file does not exist, no overwrite needed, return TRUE
    if (!file.exists(path)) {
        return(TRUE)
    } else { 
        # if file exist, ask user for what needs to be done
        user_input <- typeline(prompt = paste0(" Overwrite `", path, 
            "`? [Y/n]: "))
        cleaned_input <- tolower(trimws(user_input))

        if (cleaned_input %in% c("y", "yes")) {
            cli_alert_info("User authorized process to overwrite file.")
            return(TRUE)

        } else if (cleaned_input %in% c("n", "no")) {
            cli_alert_info("User refused to allow process to overwrite.")
            return(FALSE)

        } else {
            # if answer cannot be identified, return FALSE
            cli_alert_warning(paste0("Answer was '", user_input,  
                "', expected Y(es)/N(o)."))
            cli_alert_warning("\t↳ Defaults to `FALSE` (no overwrite).")
            return(FALSE)
        }
    }
}

# A function to import all xlsx files from BIODICAPT land surveys into one df.
# ARGS:
#   - xlsx_paths: a list of .xlsx files.
# WARNING: this function will only work correctly for biodicapt files.
import_biodicapt_land_surveys <- function(xlsx_paths) {
    
    # loop on all files and store in list
    map_dfr(xlsx_paths, function(path) {

        # get main sheet for each file.
        sheet_index <- if (grepl("MONTPELLIER", path, fixed = TRUE)) 1 else 2
        df <- read_xlsx(path, sheet = sheet_index, .name_repair = "unique_quiet")

        # remove useless columns
        cols_to_remove <- c("code_local", "Contact_Viti", "Referent_Projet", 
            "IFT total", "IFT_data_year", "lineaire(m)_haies_buffer_1km",
            "%surface_haies_buffer_1km", "%surface_ESN", "Lot", "ocsol2026"
        )
        df <- df |> 
          select(-any_of(cols_to_remove)) |>    # remove named columns
          select(-matches("^\\.\\.\\.\\d+$"))   # remove unnamed columns
      
        # Add coordinates (LON/LAT) depending on site
        if (grepl("SCARABEE", path, fixed = TRUE)) {
            # SCARABEE: axes are swapped but in the right coordinate system, 
            # rename them directly
            df <- df |> rename(LON = Y_L93, LAT = X_L93)
        } else {
            # Others: reproject from Lambert93 to WGS84
            coords <- df |>
                st_as_sf(coofs = c("X_L93", "Y_L93"), crs = 2154) |>
                st_transform(crs = 4326) |>
                st_coordinates()

            # add new columns
            df <- df |> mutate(LON = coords[, 1], LAT = coords[, 2])
          
            # remove old colums
            df <- df |> select(-X_L93, -Y_L93)
        }

        # Drop Lambert93 columns, fix types, add network label
        df |>
            mutate(
                network = path |>
                    basename() |>
                    file_path_sans_ext() |>
                    str_split("_") |>
                    pluck(1) |>
                    last()
            )
    })
}

# A function that "blurs" coordinate within a dataframe by randomly shifting
#   each coordinate on its longitutde and latitude
# ARGS:
#   - df: a dataframe with gps coordinates.
#   - x_lon: a string. The name of the column with longitude values.
#   - y_lat: a string. The name of the colum with latitude values.
#   - res_km: an integer/a float. The range within which the blurring occurs (in km).
#   - seed: an integer. When given, sets a seed at function level.
blur_coordinates <- function(df, x_lon, y_lat, res_km, seed = NULL) {
    # for reproducible results  
    if (!is.null(seed)){
        set.seed(seed) 
    }
  
    # Manual check of coordinate system (CRS 4326)
    if (    !between(min(df[[x_lon]]), LON_MIN, LON_MAX) ||
            !between(max(df[[x_lon]]), LON_MIN, LON_MAX) ||
            !between(min(df[[y_lat]]), LAT_MIN, LAT_MAX) ||
            !between(max(df[[y_lat]]), LAT_MIN, LAT_MAX)) {
        cli_alert_danger("Dataset failed basic extent test")
        stop(paste0(
            "Wrong coordinate system and/or borders are outside", 
            " France's metropolitan area."))
    }

    # get mean latitude and longitude resolutions
    lat_mean  <- (LAT_MIN + LAT_MAX) / 2
    res_lat   <- res_km / 111.0
    res_lon   <- res_km / (111.0 * cos(lat_mean * pi / 180))

    # Dividing by 1.96 so ~95% of points stay within a circle of res_km 
    df |>
        mutate(
            "{x_lon}" := .data[[x_lon]] + rnorm(dplyr::n(), mean = 0, sd = res_lon / 1.96),
            "{y_lat}" := .data[[y_lat]] + rnorm(dplyr::n(), mean = 0, sd = res_lat / 1.96)
        )
}

# A function that creates an empty raster that matches France's extent (WGS 84)
# ARGS:
#   - res_km: an integer/a float. The size of each pixel of the raster (in km, approximative).
get_france_raster_template <- function(res_km) {
    # get mean latitude and longitude resolutions
    lat_mean <- (LAT_MIN + LAT_MAX) / 2
    res_lat <- res_km / 111.0
    res_lon <- res_km / (111.0 * cos(lat_mean * pi / 180))
    
    # create empty raster with template grid
    template <- rast(
        extent = ext(LON_MIN, LON_MAX, LAT_MIN, LAT_MAX),
        resolution = c(res_lon, res_lat),
        crs = "EPSG:4326"
    )
    return(template)
}

# A function that imports a shapefile of france's borders
# ARGS:
#   - borders: a string. Either "national" (default) or "regional". if "regional", draws highest level inner borders ofthe country.
get_metropolitan_france_shapefile <- function(borders = "national") {
    # import from different functions depending on borders
    if (borders == "regional"){
        france_sf <- ne_states(country = "France", returnclass = "sf")
        
    } else if (borders == "national") {
        france_sf <- ne_countries(
            country = "France", scale = "large", returnclass = "sf")    
        
    } else {
        stop(paste0(borders, " is not recognised (should be one of 'regional' or 'national')"))
    }
    
    # Use hard coded extent of france to exclude overseas territories from shp
    st_crop(
        france_sf, 
        xmin = -5.5, 
        xmax = 9.7, 
        ymin = 41.2, 
        ymax = 51.2)
}

# A function that clips values of a raster outside of a shapefile
# ARGS:
#   - raster: a raster.
#   - shapefile: a shapefile.
clip_raster_from_shapefile <- function(raster, shapefile) {
    # verify raster and shapefile crs are marching
    if (!same.crs(raster, shapefile)) {
        shapefile <- st_transform(shapefile, crs(raster))  # fix silently
    }
    # Clip the raster using the shapefile
    mask(raster, shapefile, touches = TRUE)
}

# A function to clip a raster to france metropolitan area and ensure crs WGS84
# ARGS:
#   - raster: a raster.
#   - save_to: a string or NULL. If NULL, returns raster, else the filepath where the raster will be saved. Must end in ".tif".
#   - verbose: a boolean. If TRUE, shows info messages.
#   - res_km: an integer/a float. Is used to define a buffer around the raster (in km, approximative).
clip_raster_france_wgs84_crs <- function(
        raster,  
        buffer, 
        verbose = TRUE,
        save_to = NULL) {
    
    # import france shapefile with buffer to avoid clipping important data
    france_sf <- get_metropolitan_france_shapefile()
    france_sf_buffered <- france_sf |>
        st_transform(2154) |>              # EPSG:2154 = RGF93, better conservation of distances
        st_buffer(dist = buffer*1000) |> # buffer distance in meters
        st_transform(crs(raster))   # transform to match raster's crs
    
    # reduce size of raster
    raster_cropped <- crop(raster, ext(france_sf_buffered))
    raster_clipped <- clip_raster_from_shapefile(
        raster_cropped, france_sf_buffered)
    
    if (verbose) {
        cli_alert_info("Updating CRS, this might take some time...\n")
    }
    raster_wgs84 <- project(raster_clipped, "EPSG:4326")

    if (verbose) {
        cli_alert_info("Saving file...\n")
    }

    if (is.null(save_to)) {
        return(raster_wgs84)
    } else {
        writeRaster(raster_wgs84, save_to, overwrite = TRUE)  
    }   
}

# A function to transform a raster into a (raster) grid layered on france
# ARGS:
#   - raster: a raster.
#   - save_to: a string. The filepath where the raster will be saved. Must end in ".tif".
#   - verbose: a boolean. If TRUE, shows info messages.
#   - res_km: an integer/a float. Is used to define a buffer around the raster (in km, approximative).
project_to_france_custom_grid <- function(raster, save_to, res_km, verbose = TRUE) {   
    # fit CORINE raster to new grid
    if (verbose) {
        cli_alert_info("Fitting new grid...\n")
    }
    france_grid_template <- get_france_raster_template(res_km)
    raster_new_grid <- project(
        raster, france_grid_template, method = "near")

    if (verbose) {
        cli_alert_info("Masking...\n")
    }
    # remove values outside of france borders
    france_sf <- get_metropolitan_france_shapefile()
    france_sf_buffered <- france_sf |>
        st_transform(2154) |>              # EPSG:2154 = RGF93, better conservation of distances
        st_buffer(dist = 2*res_km*1000) |> # buffer distance in meters
        st_transform(crs(raster_new_grid)) # transform to match raster's crs
    raster_new_grid_clean <- clip_raster_from_shapefile(
        raster_new_grid, france_sf_buffered)

    if (verbose) {
        cli_alert_info("Saving file...\n")
    }
    writeRaster(raster_new_grid_clean, save_to, overwrite = TRUE)       
}

# A function to transform a raster into a (shapefile) grid of hexagons layered on france
# ARGS:
#   - raster: a raster.
#   - save_to: a string. The filepath where the shapefile will be saved. Must end in ".shp".
#   - verbose: a boolean. If TRUE, shows info messages.
#   - res_km: an integer/a float. Is used to define a buffer around the raster (in km, approximative).
#   - method: a string. If raster contains categorical data, must be "categorical". If raster contains numerical values, eithr "mean" or "median".
project_to_hexagons <- function(raster, save_to, res_km, method, verbose = TRUE) {
    # 1. Define grid of hexagons
    # get grid of hexagons over france extent
    if (verbose) {
        cli_alert_info("Creating grid of hexagons...\n")
    }
    dggs <- dgconstruct(
        spacing = res_km, metric = TRUE, resround = 'nearest')
    france_sf <- get_metropolitan_france_shapefile()
    france_sf_buffered <- france_sf |>
        st_transform(2154) |>              # EPSG:2154 = RGF93, better conservation of distances
        st_buffer(dist = 2*res_km*1000) |> # buffer distance in meters
        st_transform(crs("EPSG:4326"))     # transform to match raster's crs

    # get mean latitude and longitude resolutions
    lat_mean <- (LAT_MIN + LAT_MAX) / 2
    res_lat <- res_km / 111.0
    res_lon <- res_km / (111.0 * cos(lat_mean * pi / 180))

    full_grid <- dgshptogrid(dggs, france_sf_buffered, 
        cellsize = min(res_lat, res_lon)/2) # ensures no hexagon is forgotten

    # compute areas
    full_grid$hex_area <- st_area(full_grid)
    hex_clipped <- st_intersection(full_grid, france_sf_buffered)
    hex_clipped$clipped_area <- st_area(hex_clipped)

    # filter (keep hexagons with at least 50% area within shape)
    hex_filtered <- hex_clipped[as.numeric(hex_clipped$clipped_area / hex_clipped$hex_area) >= 0.5, ]
    hex_grid <- full_grid[full_grid$seqnum %in% hex_filtered$seqnum, ]

    # 2. Associate each hexagon to the most represented value within
    # its area in the raster
    if (verbose) {
        cli_alert_info("Extracting values, this might take some time...\n")
    }
    extracted <- exact_extract(
        raster, 
        hex_grid, 
        function(values, coverage_fracs) {
            # Keep only pixels where >50% of their surface is within the hexagon
            values <- values[coverage_fracs >= 0.5]
            if (length(values) == 0) return(NA)
            
            if (method == "categorical") {
                # Return most represented value
                as.numeric(names(sort(table(values), decreasing = TRUE))[1])
            } else if (method == "mean") {
                mean(values, na.rm = TRUE)
            } else if (method == "median") {
                median(values, na.rm = TRUE)
            } else {
                stop(paste0("Unknown method: '", method, "'. Use 'categorical', 'mean', or 'median'."))
            }
        }
    )
    hex_grid$dominant_class <- extracted

    # 3. Polish and save
    if (verbose) {
        cli_alert_info("Saving files...\n")
    }
    st_write(hex_grid, save_to, delete_dsn = TRUE, quiet = TRUE)
    if (method == "categorical") {
        # extract color and labels from raster
        category_table <- cats(raster)[[1]]
        corine_legend <- category_table %>%
            mutate(hex = rgb(Red, Green, Blue, maxColorValue = 1))
        write_csv(corine_legend, paste0(file_path_sans_ext(save_to), ".csv"))
    }
    
}

# A function to transform a dataframe of points to a shapefile of hexagons, through interpolation of missing values.
# ARGS:
#   - df: a dataframe. Must contain "LON" (longitude, WGS84) and "LAT" (latitude, WGS84) and "column".
#   - column: a string. The name of the column containing the values for interpolation.
#   - res_km: an integer/a float. The size of each hexagon on the grid (in km, approximative).
#   - LON: a string. The name of the column with longitude values.
#   - LAT: a string. The name of the column with latitude values. 
#   - idp: a numeric. Specifies the inverse distance weighting power. (lower = finer influence, higher = smoother)
#   - maxdist_m: a numeric. Only observations within a distance of maxdist from the prediction location are used for prediction or simulation.
interpolate_scattered_points_to_hexagons <- function(
        df,  
        column,
        res_km = RES_KM, 
        LON = "LON",
        LAT = "LAT",   
        idp = 2,             
        maxdist_m = Inf) {
    method <- "idw" # other methods exist, but they are not really better

    # 1. Build hexagon grid (same as project_to_hexagons)
    invisible(capture.output(dggs <- dgconstruct(
        spacing = res_km, metric = TRUE, resround = 'nearest')))
    france_sf <- get_metropolitan_france_shapefile()
    france_sf_buffered <- france_sf |>
        st_transform(2154) |>
        st_buffer(dist = 2*res_km*1000) |>
        st_transform(crs("EPSG:4326"))

    # get mean latitude and longitude resolutions
    lat_mean <- (LAT_MIN + LAT_MAX) / 2
    res_lat <- res_km / 111.0
    res_lon <- res_km / (111.0 * cos(lat_mean * pi / 180))

    # convert shapefile to grid of hexagons
    full_grid <- dgshptogrid(
        dggs, france_sf_buffered, cellsize = min(res_lat, res_lon)/2)

    full_grid$hex_area <- st_area(full_grid)
    hex_clipped <- st_intersection(full_grid, france_sf_buffered)
    hex_clipped$clipped_area <- st_area(hex_clipped)
    hex_filtered <- hex_clipped[
        as.numeric(hex_clipped$clipped_area / hex_clipped$hex_area) >= 0.5, ]
    hex_grid <- full_grid[full_grid$seqnum %in% hex_filtered$seqnum, ]

    # 2. Points as sp object (gstat predates sf and still expects Spatial*)
    points_sf <- st_as_sf(df, coords = c(LON, LAT), crs = 4326)
    points_sp <- as(points_sf, "Spatial")

    # 3. Interpolate onto hexagon centroids
    hex_centroids <- st_centroid(hex_grid)
    centroids_sp <- as(hex_centroids, "Spatial")

    formula <- as.formula(paste(column, "~ 1"))

    if (method == "idw") {
        invisible(capture.output(interp <- idw(
            formula, points_sp, 
            newdata = centroids_sp, idp = idp, maxdist = maxdist_m)))
        hex_grid$interpolated_value <- interp$var1.pred
    } else {
        stop(paste0("Unknown method: '", method, "'. Use 'idw' or 'kriging'."))
    }

    return(hex_grid)
}

# A function to aggregate rasters together (months -> year) with buffer around borders
# ARGS:
#   - raster_paths: filepaths to 12 *similar* rasters.
#   - buffer: an integer/a float. Is used to define a buffer around the raster (in km, approximative).
#   - verbose: a boolean. If TRUE, shows info messages.
#   - fun: a function. Tells how the data should be aggregated together.
monthly_2_yearly_rasters <- function(
        raster_paths, 
        buffer, 
        verbose = TRUE, 
        fun = mean){
    
    # Verify we have 12 paths (for 12 months)
    if (length(raster_paths) != 12){
        stop(paste("List of strings recieved contains", length(raster_paths), "elements, expected 12."))
    }

    # Stack rasters together
    raster_list <- c()
    for (n_path in 1:12){
        # import raster
        raster <- rast(raster_paths[n_path])
        
        if (!same.crs(raster, "EPSG:4326")){
            stop(paste0("Raster (at ", raster_paths[n_path], ") has crs ", 
                crs(raster), ". Expected 'EPSG:4326'."))
        }
        
        cropped_raster <- clip_raster_france_wgs84_crs(
            raster, verbose = verbose, buffer = buffer)     

        if (n_path == 1){
            rasters_together <- cropped_raster
        } else {
            raster_list <- c(rasters_together, cropped_raster)
        }
    }
    
    # Compute and return average/median/other on list of cropped raster
    annual_raster <- app(rasters_together, fun = fun)

    return(annual_raster)
}

# A function that collapses categories of land cover for CORINE
# ARGS:
#   - clc_raster The CLC raster to modify.
#   - level_*: an integer. For each categories, decides if the category should be collapsed together at level 1, 2 or 3 (default is 3: all subcategories are shown).
#   - verbose: a boolean. If TRUE (default), shows info messages.
simplify_CLC <- function(
        clc_raster,
        level_urban = 1,
        level_crops = 1,
        level_forests = 1,
        level_wetlands = 1,
        level_water = 1,
        verbose = TRUE){
    if (verbose) {
        cli_alert_warning("This function should only be used with CORINE Land Cover rasters!")
        cli_alert_info("Loading categories...")
    }

    # Get dataframe
    cat_table <- cats(clc_raster)[[1]]   

    # Only subcategories are named in raster. Adding other names.
    all_names <- setNames(cat_table$LABEL3, cat_table$CODE_18)
    all_names <- c(
        all_names,
        "1" = "Artificial surfaces",
        "11" = "Urban fabric",
        "12" = "Industrial, commercial and transport units",
        "13" = "Mine, dump and construction sites",
        "14" = "Artificial, non-agricultural vegetated areas",
        "2" = "Agricultural areas",
        "21" = "Arable land",
        "22" = "Permanent crops",
        "23" = "Pastures",
        "24" = "Heterogeneous agricultural areas",
        "3" = "Forest and semi-natural areas",
        "31" = "Forests",
        "32" = "Shrub and/or herbaceous vegetation associations",
        "33" = "Open spaces with little or no vegetation",
        "4" = "Wetlants",
        "41" = "Inland wetlands",
        "42" = "Coastal wetlans",
        "5" = "Water bodies",
        "51" = "Inland waters",
        "52" = "Marine waters",
        "999" = "NO_DATA"
    )

    # List for levels of collapse
    levels <- list(
        as.integer(substr(as.character(cat_table$CODE_18), 1, 1)),
        as.integer(substr(as.character(cat_table$CODE_18), 1, 2)),
        as.integer(substr(as.character(cat_table$CODE_18), 1, 3))
    )
    groups <- c(level_urban, level_crops, level_forests, 
                level_wetlands, level_water)

    if (verbose) {
        cli_alert_info("Creating new table...")
    }
    # Create new table
    cat_table <- cat_table |>
        mutate(NEW_CODE_18 = case_when(
            grepl("^1", cat_table$CODE_18) ~ as.character(levels[[ groups[[ 1 ]] ]]),
            grepl("^2", cat_table$CODE_18) ~ as.character(levels[[ groups[[ 2 ]] ]]),
            grepl("^3", cat_table$CODE_18) ~ as.character(levels[[ groups[[ 3 ]] ]]),
            grepl("^4", cat_table$CODE_18) ~ as.character(levels[[ groups[[ 4 ]] ]]),
            grepl("^5", cat_table$CODE_18) ~ as.character(levels[[ groups[[ 5 ]] ]]),
            # Default case (if none of the above match)
            TRUE ~ "999")) |>
        group_by(NEW_CODE_18) |>
        mutate(NEW_LABEL3 = all_names[NEW_CODE_18]) |>
        mutate(
            NEW_Red = mean(Red, na.rm = TRUE),
            NEW_Green = mean(Green, na.rm = TRUE),
            NEW_Blue = mean(Blue, na.rm = TRUE)) |>
        ungroup()    
    
    # Build a reclassification matrix: from Value → new integer code
    rcl <- cat_table |>
        select(Value, NEW_CODE_18) |>
        distinct() |>
        as.matrix()

    clc_raster_collapsed <- classify(clc_raster, rcl)

    # An assign new table table keyed on the new codes
    new_levels_simple <- cat_table |>
        select(NEW_CODE_18, NEW_LABEL3, NEW_Red, NEW_Green, NEW_Blue) |>
        distinct() |>
        mutate(
            Value = as.integer(NEW_CODE_18),
            Red = NEW_Red,
            Green = NEW_Green,
            Blue = NEW_Blue) |>
        select(Value, NEW_LABEL3, NEW_CODE_18, Red, Green, Blue)

    levels(clc_raster_collapsed) <- new_levels_simple 

    if (verbose) {
        cli_alert_info("Assigning new colors...")
    }
    # Build back the color table
    new_coltab <- cat_table |>
        select(NEW_CODE_18, NEW_Red, NEW_Green, NEW_Blue) |>
        distinct() |>
        mutate(
            Value = as.integer(NEW_CODE_18),
            R = round(NEW_Red * 255),
            G = round(NEW_Green * 255),
            B = round(NEW_Blue * 255),
            A = 255) |>
        select(Value, R, G, B, A) |>
        as.data.frame()
    coltab(clc_raster_collapsed) <- new_coltab

    # Show layer of names first
    activeCat(clc_raster_collapsed) <- "NEW_LABEL3" 

    if (verbose) {
        cli_alert_success("Re-classified CLC2018 is ready!")
    }
    return(clc_raster_collapsed)
}

# A function that splits STOC dataset into different subsets
#   - df: a dataframe. Must contain "carre" and "id_point_annee" columns.
#   - train_size: a numeric. The number of samples for the training samples.
#   - new_pool_size: a numeric. The number of samples for the new pool of samples.
#   - k_fold: a numeric. The number of cross-validation subsets to make.
split_stoc_points_k_fold_subsets <- function(
    df, train_size, new_pool_size, k_fold) {
    # To approach independent sampling : 
    #   - training : each sample is selected in a different square
    #   - new pool (simulation of 500 ENI) : each sample is selected in a 
    #       different square
    #   - test data : remaining samples after selection of training + new pool
    cli_alert_info("Splitting datasets...")

    if ((train_size + new_pool_size) >= length(unique(df$carre))) {
        stop(paste0(
            "Cannot split dataset, inconsistent given sizes of subsets.\n",
            "train_size + new_pool_size = ", train_size + new_pool_size, 
            " elements.", "Should be <", length(unique(df$carre)), 
            " (number of squares in dataset)."))
    }

    # Storing IDs in list (instead of full dataframes) to save storage
    k_fold_list <- list()
    for (k in seq(k_fold)) {
        # select squares
        random_order <- sample(
            unique(df$carre), length(unique(df$carre)), replace = FALSE)
        train_carre <- random_order[1:train_size]
        new_pool_carre <- random_order[(train_size+1):(train_size+new_pool_size)]

        # for train sets : select one point per square
        train_pts <- sample(
            subset(df, df$carre %in% train_carre)$id_point_annee, 
            train_size, 
            replace = FALSE)
        new_pool_pts <- sample(
            subset(df, df$carre %in% new_pool_carre)$id_point_annee, 
            new_pool_size, 
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