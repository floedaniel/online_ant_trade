# ============================================================================
# SDMtune ENSEMBLE LOOP - Random Forest + BRT
# ============================================================================
# Loops through processed species folders, creates RF + BRT ensemble models
# Reads occurrence data from species folders (same as MaxEnt pipeline)
# Performs variable selection, trains tuned RF and BRT, combines using AUC weights
#
# Why RF + BRT only?
# - MaxEnt handled separately in main pipeline (10_updated_sdmtune_loop.R)
# - RF and BRT are complementary: RF handles interactions, BRT handles gradients
# - Both methods robust to overfitting with proper tuning
# - Ensemble diversity achieved through different algorithms
# ============================================================================

# ---------- Load Libraries ----------
library(SDMtune)
library(terra)
library(tidyverse)
library(sf)
library(zeallot)
library(CoordinateCleaner)
library(spocc)
library(rgbif)
library(rio)
library(rnaturalearth)
library(ggplot2)
library(tidyterra)

# ---------- Configuration ----------
MIN_RECORDS <- 50           # Minimum records needed for modeling
CV_FOLDS <- 5               # Cross-validation folds
BUFFER_DISTANCE <- 3000000  # 3000 km background buffer
base_dir <- "./species"

# ---------- Load Environmental Data ----------
message("\n", paste(rep("=", 80), collapse = ""))
message("LOADING ENVIRONMENTAL DATA")
message(paste(rep("=", 80), collapse = ""))

current_predictors_dir <- "./1_raw_data/bioclim/current"
current_predictors_files <- list.files(current_predictors_dir, pattern = "\\.tif$", full.names = TRUE)
current_predictors <- terra::rast(current_predictors_files)
message("Loaded current bioclim: ", terra::nlyr(current_predictors), " layers")

future_predictors_dir <- "./1_raw_data/bioclim/future"
future_predictors_files <- list.files(future_predictors_dir, pattern = "\\.tif$", full.names = TRUE)
future_predictors <- terra::rast(future_predictors_files)
message("Loaded future bioclim: ", terra::nlyr(future_predictors), " layers")

# Load world map for background filtering
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# ---------- Load Master Data and Prescreening Results ----------
message("\nLoading species list...")

master_file <- "./2_processed_data/complete_ant_data.xlsx"
master_data <- rio::import(master_file)

prescreening_file <- "./5_outputs/screen_output/prescreening_bio6_spocc.xlsx"
prescreening_data <- rio::import(prescreening_file)

# Filter prescreening data for INCLUDE decisions
included_species <- prescreening_data %>%
  filter(decision == "INCLUDE") %>%
  pull(species)

message("Found ", length(included_species), " species marked as INCLUDE")

# ---------- Norwegian Native Species Filter ----------
message("Fetching GBIF keys for Norwegian native species...")

in_norway <- rio::import("./2_processed_data/Arter_Norge.xlsx")
norwegian_ant_species <- in_norway %>% dplyr::select(Species, Norway) %>% filter(Norway=="Ja") %>% pull(Species)

norwegian_keys <- map_dfr(norwegian_ant_species, function(sp) {
  match <- tryCatch({
    name_backbone(name = sp, rank = "SPECIES")
  }, error = function(e) NULL)

  if (!is.null(match) && !is.null(match$usageKey)) {
    tibble(
      input_name = sp,
      acceptedKey = match$acceptedUsageKey %||% match$usageKey
    )
  } else {
    tibble(input_name = sp, acceptedKey = NA_integer_)
  }
})

norwegian_accepted_keys <- norwegian_keys %>%
  filter(!is.na(acceptedKey)) %>%
  pull(acceptedKey) %>%
  unique()

# ---------- Filter Species List ----------
test_species <- master_data %>%
  filter(status == "ACCEPTED") %>%
  filter(!acceptedKey %in% norwegian_accepted_keys) %>%
  filter(name %in% included_species) %>%
  dplyr::select(name, acceptedKey)

message("Final species list: ", nrow(test_species), " species to process")

# ---------------------------------------------------------------------------------------------
# For testing - uncomment to run on subset
# test_species <- test_species[1:3, ]
# ---------------------------------------------------------------------------------------------

overall_start_time <- Sys.time()
success_count <- 0
skip_count <- 0
fail_count <- 0

# ============================================================================
# MAIN SPECIES LOOP
# ============================================================================

for (i in 1:nrow(test_species)) {
  species_start_time <- Sys.time()

  species_name <- test_species$name[i]
  species_key <- test_species$acceptedKey[i]

  message("\n", paste(rep("=", 80), collapse = ""))
  message("Processing species ", i, " of ", nrow(test_species))
  message("Species: ", species_name, " | Key: ", species_key)
  message(paste(rep("=", 80), collapse = ""))

  # Create folder structure
  safe_species_name <- gsub(" ", "_", species_name)
  species_folder_name <- paste0(species_key, "_", safe_species_name)
  species_folder <- file.path(base_dir, species_folder_name)
  sdm_folder <- file.path(species_folder, "SDM_maxnet")
  ensemble_folder <- file.path(sdm_folder, "ensemble")

  # Skip if ensemble folder already exists
  if (dir.exists(ensemble_folder)) {
    message("Ensemble folder already exists - SKIPPING")
    skip_count <- skip_count + 1
    next
  }

  # Check if species folder exists
  if (!dir.exists(species_folder)) {
    message("Species folder does not exist - SKIPPING")
    fail_count <- fail_count + 1
    next
  }

  # Create ensemble folder
  dir.create(ensemble_folder, recursive = TRUE, showWarnings = FALSE)

  # ---------- Download Occurrence Data ----------
  message("Downloading occurrence data...")

  occ_data <- tryCatch({
    spocc::occ(
      query = species_name,
      from = c("gbif", "bison", "ala", "inat", "idigbio"),
      limit = 100000
    )
  }, error = function(e) {
    message("Error downloading occurrence data: ", e$message)
    NULL
  })

  if (is.null(occ_data)) {
    message("Failed to download occurrence data - SKIPPING")
    fail_count <- fail_count + 1
    next
  }

  occ_df <- occ2df(occ_data)
  occ_df <- dplyr::select(occ_df, name, longitude, latitude, date, prov)
  names(occ_df) <- c("species", "decimalLongitude", "decimalLatitude", "date", "source")
  occ_df$species <- species_name

  message("Raw occurrence records: ", nrow(occ_df))

  if (nrow(occ_df) < MIN_RECORDS) {
    message("Insufficient occurrence data (< ", MIN_RECORDS, " records) - SKIPPING")
    fail_count <- fail_count + 1
    next
  }

  # ---------- Data Cleaning ----------
  message("Cleaning occurrence data...")

  # Convert to numeric
  occ_df$decimalLongitude <- as.numeric(occ_df$decimalLongitude)
  occ_df$decimalLatitude <- as.numeric(occ_df$decimalLatitude)

  # Remove NA coordinates
  occ_df <- occ_df %>% filter(!is.na(decimalLongitude) & !is.na(decimalLatitude))

  # CoordinateCleaner
  occ_df_clean <- occ_df %>%
    cc_val(verbose = FALSE) %>%
    cc_equ(verbose = FALSE) %>%
    cc_cap(verbose = FALSE) %>%
    cc_cen(verbose = FALSE) %>%
    cc_sea(verbose = FALSE) %>%
    cc_zero(verbose = FALSE) %>%
    cc_dupl(verbose = FALSE) %>%
    cc_outl(method = "quantile", mltpl = 3, value = "clean", verbose = FALSE)

  message("After cleaning: ", nrow(occ_df_clean), " records")

  if (nrow(occ_df_clean) < MIN_RECORDS) {
    message("Insufficient data after cleaning - SKIPPING")
    fail_count <- fail_count + 1
    next
  }

  # Spatial thinning
  message("Applying spatial thinning...")
  occ_thinned <- tryCatch({
    thinData(coords = as.data.frame(occ_df_clean),
             env = current_predictors[[1]],
             x = "decimalLongitude",
             y = "decimalLatitude")
  }, error = function(e) {
    message("Thinning failed: ", e$message)
    occ_df_clean
  })

  message("After thinning: ", nrow(occ_thinned), " records")

  if (nrow(occ_thinned) < MIN_RECORDS) {
    message("Insufficient data after thinning - SKIPPING")
    fail_count <- fail_count + 1
    next
  }

  p_coords <- occ_thinned %>%
    dplyr::select(decimalLongitude, decimalLatitude)

  # ---------- Generate Background Points ----------
  message("Generating background points...")

  presence_sf <- sf::st_as_sf(p_coords, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
  presence_projected <- sf::st_transform(presence_sf, crs = 3857)

  circles <- sf::st_buffer(presence_projected, dist = BUFFER_DISTANCE)
  merged_polygon <- sf::st_union(circles)

  # Fix invalid geometries (prevents "Edge X crosses edge Y" errors)
  merged_polygon <- sf::st_make_valid(merged_polygon)

  merged_polygon_geo <- sf::st_transform(merged_polygon, crs = 4326)

  bg_points <- sf::st_sample(merged_polygon_geo, size = 20000, type = "random")
  bg_coords_sf <- sf::st_sf(geometry = bg_points)

  # Filter to land
  bg_coords_land <- sf::st_filter(bg_coords_sf, world)
  bg_coords <- sf::st_coordinates(bg_coords_land)
  bg_coords <- as.data.frame(bg_coords[, c("X", "Y")])
  names(bg_coords) <- c("decimalLongitude", "decimalLatitude")

  message("Generated ", nrow(bg_coords), " background points")

  # ---------- Variable Selection ----------
  message("Performing variable selection...")

  # Prepare initial SWD with ALL variables
  initial_swd <- tryCatch({
    prepareSWD(species = species_name, p = p_coords, a = bg_coords, env = current_predictors)
  }, error = function(e) {
    message("prepareSWD failed: ", e$message)
    NULL
  })

  if (is.null(initial_swd)) {
    message("Failed to prepare SWD - SKIPPING")
    fail_count <- fail_count + 1
    next
  }

  c(train_initial, test_initial) %<-% trainValTest(initial_swd, test = 0.2, only_presence = TRUE, seed = 25)

  # Prepare background SWD for correlation analysis
  bg_swd <- prepareSWD(species = "background", p = bg_coords, a = bg_coords, env = current_predictors)

  # Variable selection with correlation threshold
  selected_vars_model <- tryCatch({
    varSel(
      model = train(method = "Maxnet", data = train_initial),
      metric = "auc",
      test = test_initial,
      bg4cor = bg_swd,
      method = "spearman",
      cor_th = 0.7,
      env = current_predictors,
      use_pc = FALSE,
      progress = FALSE,
      permut = 10
    )
  }, error = function(e) {
    message("Variable selection failed (cor_th=0.7): ", e$message)
    message("Retrying with relaxed threshold (cor_th=0.85)...")

    tryCatch({
      varSel(
        model = train(method = "Maxnet", data = train_initial),
        metric = "auc",
        test = test_initial,
        bg4cor = bg_swd,
        method = "spearman",
        cor_th = 0.85,
        env = current_predictors,
        use_pc = FALSE,
        progress = FALSE,
        permut = 5
      )
    }, error = function(e2) {
      message("Variable selection failed completely: ", e2$message)
      NULL
    })
  })

  if (is.null(selected_vars_model)) {
    message("Variable selection failed - SKIPPING")
    fail_count <- fail_count + 1
    next
  }

  selected_vars <- selected_vars_model@data@data %>% 
    names()

  message("Selected ", length(selected_vars), " variables")

  # Subset predictors to selected variables
  final_predictors <- current_predictors[[selected_vars]]
  future_predictors_subset <- future_predictors[[selected_vars]]

  # Prepare final SWD with selected variables
  final_swd <- prepareSWD(species = species_name, p = p_coords, a = bg_coords, env = final_predictors)
  c(train_data, test_data) %<-% trainValTest(final_swd, test = 0.2, only_presence = TRUE, seed = 25)

  message("Train: ", nrow(train_data@data), " | Test: ", nrow(test_data@data))

  # ============================================================================
  # TRAIN AND TUNE MODELS
  # ============================================================================

  message("\n========== TRAINING & TUNING MODELS ==========")

  # ---------- Model 1: Random Forest (Tuned) ----------
  message("Training & tuning Random Forest...")

  rf_hyperparams <- list(
    mtry = c(3, 4, 5, 6, 7),
    ntree = c(500, 750, 1000, 1200)
  )

  rf_base <- train(method = "RF", data = train_data)

  rf_tuned <- tryCatch({
    optimizeModel(
      model = rf_base,
      hypers = rf_hyperparams,
      metric = "auc",
      test = test_data,
      pop = 10,
      gen = 10,
      seed = 123,
      progress = FALSE
    )
  }, error = function(e) {
    message("RF tuning failed: ", e$message, " - using default parameters")
    rf_base
  })

  # Get best RF model and AUC
  if ("SDMmodelCV" %in% class(rf_tuned)) {
    best_idx <- which.max(rf_tuned@results$test_AUC)
    rf_model <- rf_tuned@models[[best_idx]]
    rf_auc <- rf_tuned@results$test_AUC[best_idx]
  } else {
    rf_model <- rf_tuned
    rf_auc <- SDMtune::auc(rf_model, test = test_data)
  }

  message("RF AUC: ", round(rf_auc, 3))

  # Extract RF hyperparameters
  rf_ranger <- rf_model@models[[1]]
  message("  Best hyperparameters: mtry=", rf_ranger@model@mtry, ", ntree=", rf_ranger@model@ntree)

  # ---------- Model 2: BRT (Tuned) ----------
  message("Training & tuning BRT...")

  brt_hyperparams <- list(
    n.trees = c(500, 750, 1000, 1500),
    interaction.depth = c(2, 3, 4, 5),
    shrinkage = c(0.001, 0.005, 0.01, 0.05)
  )

  brt_base <- train(method = "BRT", data = train_data,
                    n.trees = 500, interaction.depth = 3, shrinkage = 0.01)

  brt_tuned <- tryCatch({
    optimizeModel(
      model = brt_base,
      hypers = brt_hyperparams,
      metric = "auc",
      test = test_data,
      pop = 10,
      gen = 10,
      seed = 123,
      progress = FALSE
    )
  }, error = function(e) {
    message("BRT tuning failed: ", e$message, " - using default parameters")
    brt_base
  })

  # Get best BRT model and AUC
  if ("SDMmodelCV" %in% class(brt_tuned)) {
    best_idx <- which.max(brt_tuned@results$test_AUC)
    brt_model <- brt_tuned@models[[best_idx]]
    brt_auc <- brt_tuned@results$test_AUC[best_idx]
  } else {
    brt_model <- brt_tuned
    brt_auc <- SDMtune::auc(brt_model, test = test_data)
  }

  message("BRT AUC: ", round(brt_auc, 3))

  # Extract BRT hyperparameters
  brt_gbm <- brt_model@models[[1]]
  
  message("  Best hyperparameters: n.trees=", brt_gbm@model@n.trees,
          ", depth=", brt_gbm@model@interaction.depth,
          ", shrinkage=", brt_gbm@model@shrinkage)

  # ============================================================================
  # CREATE ENSEMBLE (RF + BRT only)
  # ============================================================================

  message("\n========== CREATING ENSEMBLE ==========")

  # Predict with each model - CURRENT climate
  message("Generating current climate predictions...")
  pred_rf_current <- predict(rf_ranger, data = final_predictors, type = "cloglog")
  pred_brt_current <- predict(brt_gbm, data = final_predictors, type = "cloglog")

  # Predict with each model - FUTURE climate
  message("Generating future climate predictions...")
  pred_rf_future <- predict(rf_ranger, data = future_predictors_subset, type = "cloglog")
  pred_brt_future <- predict(brt_gbm, data = future_predictors_subset, type = "cloglog")

  # Calculate AUC-based weights
  aucs <- c(rf_auc, brt_auc)
  weights <- aucs / sum(aucs)

  message("AUC-based weights:")
  message("  RF  = ", round(weights[1], 3), " (AUC: ", round(rf_auc, 3), ")")
  message("  BRT = ", round(weights[2], 3), " (AUC: ", round(brt_auc, 3), ")")

  # Calculate weighted ensemble predictions
  ensemble_current <- weights[1] * pred_rf_current + weights[2] * pred_brt_current
  ensemble_future <- weights[1] * pred_rf_future + weights[2] * pred_brt_future

  # Calculate climate change
  climate_change <- ensemble_future - ensemble_current

  # ============================================================================
  # SAVE OUTPUTS
  # ============================================================================

  message("\n========== SAVING OUTPUTS ==========")

  # Save rasters
  terra::writeRaster(ensemble_current,
                     file.path(ensemble_folder, paste0("ensemble_current_", species_key, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(ensemble_future,
                     file.path(ensemble_folder, paste0("ensemble_future_", species_key, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(climate_change,
                     file.path(ensemble_folder, paste0("ensemble_change_", species_key, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(pred_rf_current,
                     file.path(ensemble_folder, paste0("rf_current_", species_key, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(pred_brt_current,
                     file.path(ensemble_folder, paste0("brt_current_", species_key, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(pred_rf_future,
                     file.path(ensemble_folder, paste0("rf_future_", species_key, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(pred_brt_future,
                     file.path(ensemble_folder, paste0("brt_future_", species_key, ".tif")),
                     overwrite = TRUE)

  # Save evaluation table
  eval_df <- data.frame(
    Model = c("RF", "BRT"),
    AUC = round(c(rf_auc, brt_auc), 3),
    Weight = round(weights, 3),
    N_variables = length(selected_vars),
    N_presence = nrow(p_coords),
    N_background = nrow(bg_coords)
  )

  rio::export(eval_df, file.path(ensemble_folder, paste0("evaluation_", species_key, ".xlsx")))

  # Save selected variables
  var_df <- data.frame(Variable = selected_vars)
  rio::export(var_df, file.path(ensemble_folder, paste0("selected_variables_", species_key, ".xlsx")))

  # ============================================================================
  # VISUALIZE
  # ============================================================================

  message("Creating visualizations...")

  map_theme <- function() {
    theme_minimal() + theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_line(color = "white", size = 0.3),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey30"),
      legend.position = "bottom",
      legend.key.width = unit(2, "cm")
    )
  }

  # Current ensemble
  p_current <- ggplot() +
    tidyterra::geom_spatraster(data = ensemble_current, maxcell = 5e+07) +
    scale_fill_viridis_c(name = "Suitability", option = "plasma", limits = c(0, 1), na.value = "white") +
    labs(title = paste("Ensemble Current -", species_name),
         subtitle = paste("RF+BRT weighted | RF AUC:", round(rf_auc, 3), "| BRT AUC:", round(brt_auc, 3))) +
    map_theme()

  ggsave(file.path(ensemble_folder, paste0("ensemble_current_", species_key, ".png")),
         p_current, width = 12, height = 9, dpi = 300, bg = "white")

  # Future ensemble
  p_future <- ggplot() +
    tidyterra::geom_spatraster(data = ensemble_future, maxcell = 5e+07) +
    scale_fill_viridis_c(name = "Suitability", option = "plasma", limits = c(0, 1), na.value = "white") +
    labs(title = paste("Ensemble Future -", species_name),
         subtitle = "SSP585 2021-2040 | RF+BRT weighted") +
    map_theme()

  ggsave(file.path(ensemble_folder, paste0("ensemble_future_", species_key, ".png")),
         p_future, width = 12, height = 9, dpi = 300, bg = "white")

  # Climate change
  change_vals <- terra::values(climate_change, na.rm = TRUE)
  mean_change <- mean(change_vals, na.rm = TRUE)
  pos_change <- sum(change_vals > 0, na.rm = TRUE) / length(change_vals) * 100
  neg_change <- sum(change_vals < 0, na.rm = TRUE) / length(change_vals) * 100

  p_change <- ggplot() +
    tidyterra::geom_spatraster(data = climate_change, maxcell = 5e+07) +
    scale_fill_gradient2(name = "Change", low = "#d73027", mid = "#ffffbf", high = "#1a9850",
                         midpoint = 0, limits = c(-1, 1), na.value = "transparent",
                         oob = scales::squish) +
    labs(title = paste("Ensemble Climate Change -", species_name),
         subtitle = paste("Mean:", round(mean_change, 3), "| Gain:", round(pos_change, 1),
                         "% | Loss:", round(neg_change, 1), "%")) +
    map_theme()

  ggsave(file.path(ensemble_folder, paste0("ensemble_change_", species_key, ".png")),
         p_change, width = 12, height = 9, dpi = 300, bg = "white")

  # ---------- Summary ----------
  species_elapsed <- difftime(Sys.time(), species_start_time, units = "mins")
  message("\nSUCCESS - Completed in ", round(species_elapsed, 2), " minutes")
  message("Outputs saved to: ", ensemble_folder)

  success_count <- success_count + 1
}

# ============================================================================
# FINAL SUMMARY
# ============================================================================

overall_end_time <- Sys.time()
overall_elapsed <- difftime(overall_end_time, overall_start_time, units = "mins")

message("\n", paste(rep("=", 80), collapse = ""))
message("ENSEMBLE MODELING COMPLETE")
message(paste(rep("=", 80), collapse = ""))
message("Total species processed: ", nrow(test_species))
message("  Success: ", success_count)
message("  Skipped (already exists): ", skip_count)
message("  Failed: ", fail_count)
message("Total time: ", round(overall_elapsed, 2), " minutes")
if (success_count > 0) {
  message("Average per species: ", round(overall_elapsed / success_count, 2), " minutes")
}
message(paste(rep("=", 80), collapse = ""))
