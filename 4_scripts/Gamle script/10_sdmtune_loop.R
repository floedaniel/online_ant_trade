library(spocc)
library(rgbif)
library(tidyverse)
library(taxize)
library(dplyr)
library(CoordinateCleaner)
library(SDMtune)
library(sp)
library(rgeos)
library(dismo)
library(terra)
library(ggplot2)
library(ggspatial)
library(viridis)
library(metR)
library(tidyterra)
library(rJava)
library(zeallot)
library(readxl)
library(rio)

# Set working directories ------------------------------------------------------
# Use your species folder (each species will be stored in its own speciesKey folder)
base_dir <- "C:/Users/dafl/OneDrive - Folkehelseinstituttet/VKM Data/27.02.2025_maur_forprosjekt_biologisk_mangfold/data/species"
master_file <- "./2_processed_data/master_data.xlsx"
master_data <- rio::import(master_file) %>% 
  as_tibble() %>% 
  filter(!is.na(speciesKey) & !is.na(species))  # assuming columns 'speciesKey' and 'species'

# Load current_predictors ------------------------------------------------------------
# Load the base current and future rasters
current <- terra::rast("C:/Users/dafl/OneDrive - Folkehelseinstituttet/VKM Data/27.02.2025_maur_forprosjekt_biologisk_mangfold/data/1_raw_data/bioclim/current/current.tif")
future  <- terra::rast("C:/Users/dafl/OneDrive - Folkehelseinstituttet/VKM Data/27.02.2025_maur_forprosjekt_biologisk_mangfold/data/1_raw_data/bioclim/future/future.tif")
# Define the directory containing additional environmental .tif files
env_dir <- "C:/Users/dafl/OneDrive - Folkehelseinstituttet/VKM Data/27.02.2025_maur_forprosjekt_biologisk_mangfold/data/1_raw_data/bioclim/env"

# List all .tif files in the directory
env_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)

# Load the additional environmental layers as a SpatRaster stack
env_stack <- terra::rast(env_files)

# Combine the base and additional layers
current_predictors <- c(current, env_stack)
future_predictors  <- c(future, env_stack)

# Validate current_predictors and crop to occurrence data extent ---------------------
raster_extent <- terra::ext(-180, 180, -60, 90)
current_predictors <- terra::crop(current_predictors, raster_extent)

# Download Natural Earth Data --------------------------------------------------
world <- tryCatch({
  rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
}, error = function(e) {
  message("Error downloading Natural Earth data: ", e$message)
  stop("Script halted due to missing geographic data.")
})

# Process one species at a time ------------------------------------------------
for (i in seq_len(nrow(master_data))) {
  speciesKey <- as.character(master_data$speciesKey[i])
  species_name <- master_data$species[i]
  message("Processing species: ", species_name, " (SpeciesKey: ", speciesKey, ")")
  
  # Create folder structure: species folder and an SDMtune subfolder
  species_folder <- file.path(base_dir, speciesKey)
  sdm_folder <- file.path(species_folder, "SDMtune")
  
  if (!dir.exists(species_folder)) dir.create(species_folder, recursive = TRUE)
  if (!dir.exists(sdm_folder)) dir.create(sdm_folder, recursive = TRUE)
  

  # Download occurrence data from GBIF ----------------------------------------
  occ_data <- tryCatch({
    occ_data(taxonKey = speciesKey, limit = 100000)$data
  }, error = function(e) NULL)
  
  if (is.null(occ_data) || nrow(occ_data) < 100) {
    message("No occurrence data found for species: ", species_name)
    writeLines("No occurrence data found in GBIF for this species.",
               file.path(sdm_folder, "no_occurrence_data.txt"))
    next
  }
  
  occ_df <- occ_data %>%
    as_tibble() %>%
    select(scientificName, decimalLongitude, decimalLatitude) %>%
    filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
    rename(species = scientificName)
  
  # Check if required columns exist -------------------------------------------
  if (!all(c("decimalLongitude", "decimalLatitude", "species") %in% colnames(occ_df))) {
    message("Required columns missing for species: ", species_name)
    writeLines("Occurrence data is missing required coordinates.",
               file.path(sdm_folder, "no_occurrence_data.txt"))
    next
  }
  
  # Replace species column with the original species name (from master file)
  occ_df <- occ_df %>% mutate(species = species_name)
  
  # Coordinate cleaning --------------------------------------------------------
  occ_df_clean <- tryCatch({
    occ_df %>%
      cc_val() %>%
      cc_equ() %>%
      cc_cap() %>%
      cc_cen() %>%
      cc_sea() %>%
      cc_zero() %>%
      cc_dupl() %>% 
      cc_outl(
        method = "quantile",   # quantile-based method for stricter detection
        mltpl = 3,             # using a lower multiplier
        verbose = TRUE,        # print details
        value = "clean"
      )
  }, error = function(e) {
    message("Error during coordinate cleaning for species: ", species_name)
    NULL
  })
  
  if (is.null(occ_df_clean) || nrow(occ_df_clean) < 80) {
    message("Insufficient data for species: ", species_name)
    writeLines("Insufficient occurrence data in GBIF for modeling this species.",
               file.path(sdm_folder, "no_occurrence_data.txt"))
    next
  }
  
  # Save cleaned occurrence data ----------------------------------------------
  cleaned_data_path <- file.path(sdm_folder, paste0("cleaned_occurrence_data_", speciesKey, ".csv"))
  rio::export(occ_df_clean, cleaned_data_path)
  
  # -------------------------------------------------------------------------
  # Define a modeling loop function using SDMtune for the species occurrence data
  modeling_loop <- function(df, current_predictors) {
    species_list <- unique(df$species)
    
    if (!dir.exists("output")) dir.create("output")
    
    for (sp in species_list) {
      message("Modeling species: ", sp)
      
      sp_data <- occ_df_clean
      
      sp_data <- subset(df, species == sp)
      coordinates(sp_data) <- ~decimalLongitude + decimalLatitude
      projection(sp_data) <- CRS("+proj=longlat +datum=WGS84")
      
      circles <- dismo::circles(sp_data, d = 2000000, lonlat = TRUE)
      merged_polygons <- rgeos::gUnaryUnion(circles@polygons)
      
      sampled_points <- sp::spsample(merged_polygons, 20000, type = "random", iter = 50)
      sampled_cells <- terra::cellFromXY(raster(current_predictors[[1]]), sampled_points)
      bg_coords <- unique(terra::xyFromCell(raster(current_predictors[[1]]), sampled_cells))
      
      world <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
      
      bg_coords <- sf::st_as_sf(x = na.omit(data.frame(bg_coords)),
                                coords = c("x", "y"),
                                crs = "+proj=longlat +datum=WGS84 +no_defs +type=crs")
      
      world <- sf::st_make_valid(world)
      bg_coords <- sf::st_make_valid(bg_coords)
      bg_coords <- sf::st_filter(bg_coords, world)
      bg_coords <- sfheaders::sf_to_df(bg_coords, fill = FALSE, unlist = NULL)
      bg_coords <- bg_coords[,3:4]
      
      p_coords <- as.data.frame(coordinates(sp_data))
      valid_coords <- terra::extract(current_predictors, as.matrix(p_coords), cells = TRUE)
      p_coords <- p_coords[!is.na(valid_coords[, 1]), ]
      
      if (nrow(p_coords) > 0) {
        bg_swd <- prepareSWD(species = sp, p = p_coords, a = bg_coords, env = current_predictors)
      } else {
        message("No valid presence points for species: ", sp)
        next
      }
      
      c(train, test) %<-% trainValTest(bg_swd, test = 0.2, only_presence = TRUE, seed = 25)
      
      model <- train(method = "Maxnet", data = train, fc = "lqpth", reg=1)
      
      vi <- varImp(model, permut = 5)
      
      top_vars <- vi %>% 
        arrange(desc(Permutation_importance)) %>% 
        slice_head(n = 10) %>% 
        pull(Variable)
      
      # Subset current_predictors to only the top 10 variables
      current_predictors <- current_predictors[[top_vars]]
      bg_swd <- prepareSWD(species = sp, p = p_coords, a = bg_coords, env = current_predictors)
      c(train, test) %<-% trainValTest(bg_swd, test = 0.2, only_presence = TRUE, seed = 25)
      
      background <- prepareSWD(species = sp,
                               a = bg_coords,
                               env = current_predictors)
      
      selected_variables_model <- varSel(
        model = train(method = "Maxnet", data = train),
        metric = "auc",
        test = test,
        bg4cor = background,
        method = "spearman",
        cor_th = 0.7,
        env = current_predictors,
        use_pc = FALSE,
        progress = TRUE,
        permut = 5
      )
      
      selected_vars <- names(selected_variables_model@data@data)
      reduced_current_predictors <- current_predictors[[selected_vars]]
      
      reduced_swd <- prepareSWD(species = sp, p = p_coords, a = bg_coords, env = reduced_current_predictors)
      c(train, test) %<-% trainValTest(reduced_swd, test = 0.2, only_presence = TRUE, seed = 25)
      
      h <- list(reg = seq(0.5, 10, 0.5),
                fc = c("lq", "lp", "lt", "lh", "qp", "qt", "qh", "pt", "ph", "th",
                       "lqp", "lqt", "lqh", "lpt", "lph", "lth", "qpt", "qph", "qth", "pth",
                       "lqpt", "lqph", "lqth", "lpth", "qpth"))
      
      optimized_model <- optimizeModel(
        model = train(method = "Maxnet", data = train),
        hypers = h,
        metric = "auc",
        test = test,
        keep_best = 0.5,
        seed = 124,
        interactive = TRUE,
        progress = TRUE
      )
      
      final_model <- optimized_model@models[[which.max(optimized_model@results$test_AUC)]]
      
      # Create a folder name (relative) for the model report
      report_folder <- paste0("modelReport_", speciesKey)
      
      # Save current working directory
      old_wd <- getwd()
      # Change working directory to the SDMtune folder (sdm_folder is already absolute)
      setwd(sdm_folder)
      
      # Now call modelReport using a relative folder name
      modelReport(final_model,
                  type = "cloglog",
                  folder = report_folder,
                  test = test,
                  response_curves = TRUE,
                  only_presence = TRUE,
                  jk = TRUE,
                  env = reduced_current_predictors,
                  permut = 2)
      
      # Revert to original working directory
      setwd(old_wd)
      
      message("Performing five-fold cross-validation for: ", sp)
      folds <- randomFolds(reduced_swd, 5, only_presence = TRUE, seed = 321)
      
      cv_model <- train("Maxent",
                        data = reduced_swd,
                        fc = final_model@model@fc,
                        reg = final_model@model@reg,
                        folds = folds)

      
      pred <- SDMtune::predict(cv_model, reduced_current_predictors, type = "cloglog", folds = folds, fun = c("mean", "max"),)
      prediction_path <- file.path(sdm_folder, paste0("predicted_distribution_", speciesKey, ".tif"))
      terra::writeRaster(pred$mean, prediction_path, overwrite = TRUE)
      
      plot_title <- paste("Maxnet Prediction for", sp)
      
      p1 <- ggplot() +
        tidyterra::geom_spatraster(data = pred$mean) +
        scale_fill_viridis_c(name = "Occurrence Probability", option = "plasma", limits = c(0, 1)) +
        labs(title = plot_title, x = "Longitude", y = "Latitude", subtitle = "Five-fold crossvalidation") +
        theme_minimal()
      
      plot_path <- file.path(sdm_folder, paste0("prediction_plot_", speciesKey, ".png"))
      ggplot2::ggsave(plot_path, p1, width = 10, height = 7)
      
      p2 <- ggplot() +
        tidyterra::geom_spatraster(data = current_predictors[[1]]) +
        geom_point(data = bg_coords, aes(x = x, y = y), col = "gray") +
        geom_point(data = p_coords, aes(x = decimalLongitude, y = decimalLatitude)) +
        scale_fill_viridis_c(name = "Occurrence records", option = "plasma") +
        labs(title = "Points", x = "Longitude", y = "Latitude") +
        theme_minimal() +
        theme(legend.position = "none")
      
      plot_path <- file.path(sdm_folder, paste0("Occurrence_plot_", speciesKey, ".png"))
      ggplot2::ggsave(plot_path, p2, width = 10, height = 7)
      
      ths_current <- thresholds(cv_model@models[[1]], type = "cloglog")
      pa_current_path <- file.path(sdm_folder, paste0("pa_map_current_", speciesKey, ".tif"))
      p_ths_current <- plotPA(pred$mean, th = ths_current[3, 2])
      ggsave(filename = file.path(sdm_folder, "ths_current.png"), plot = p_ths_current, , width = 40, height = 30, units = "cm")
      
    # future prediction -------------------------------------------------------------------------
      future_predictors_subset <- future_predictors[[ names(reduced_current_predictors) ]]
      
      future_swd <- prepareSWD(species = sp, p = p_coords, a = bg_coords, env = future_predictors_subset)
      c(train, test) %<-% trainValTest(future_swd, test = 0.2, only_presence = TRUE, seed = 25)
      
      folds <- randomFolds(future_swd, 5, only_presence = TRUE, seed = 321)
      
      future_cv_model <- train("Maxent",
                               data = future_swd,
                               fc = final_model@model@fc,
                               reg = final_model@model@reg,
                               folds = folds)
      
      future_pred <- SDMtune::predict(future_cv_model, future_predictors_subset, type = "cloglog", folds = folds, fun = c("mean", "max"))
      
      p3 <- ggplot() +
        tidyterra::geom_spatraster(data = future_pred$mean) +
        scale_fill_viridis_c(name = "Occurrence Probability", option = "plasma", limits = c(0, 1)) +
        labs(title = plot_title, x = "Longitude", y = "Latitude", subtitle = "Five-fold crossvalidation") +
        theme_minimal()
      
      future_plot_path <- file.path(sdm_folder, paste0("prediction_plot_future_", speciesKey, ".png"))
      ggplot2::ggsave(future_plot_path, p3, width = 10, height = 7)
      
      ths_future <- thresholds(future_cv_model@models[[1]], type = "cloglog")
      pa_future_path <- file.path(sdm_folder, paste0("pa_map_future_", speciesKey))
      p_ths_future <- plotPA(future_pred$mean, th = ths_future[3, 2])
      ggsave(filename = file.path(sdm_folder, "ths_future.png"), plot = p_ths_future, , width = 40, height = 30, units = "cm")
      
      
      message("Modeling completed for species: ", sp)
      flush.console()
    }
    flush.console()
    gc()  # Clear memory after processing species
  }
  
  tryCatch({
    modeling_loop(occ_df_clean, current_predictors)
    
    output_files <- list.files("output", full.names = TRUE)
    file.copy(output_files, sdm_folder, overwrite = TRUE)
    unlink("output", recursive = TRUE)
    
    message("Modeling completed for species: ", species_name)
  }, error = function(e) {
    message("Error during modeling for species: ", species_name, " - ", e$message)
    writeLines(paste("Error during modeling:", e$message),
               file.path(sdm_folder, "error_log.txt"))
  })
  
  flush.console()
  gc()  # Clear memory after processing each species
}
# END
