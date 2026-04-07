# ============================================================================
# SIMPLIFIED SDMtune Loop - MaxEnt Best Practices
# Based on: Merow et al. 2013, Elith et al. 2011, Radosavljevic & Anderson 2014
# ============================================================================

library(spocc)
library(rgbif)
library(tidyverse)
library(CoordinateCleaner)
library(SDMtune)
library(terra)
library(sf)
library(ggplot2)
library(rio)
library(rnaturalearth)

# Configuration
MIN_RECORDS <- 100
DOWNLOAD_LIMIT <- 500
CV_FOLDS <- 4
RECENT_YEARS <- 1950
BUFFER_DISTANCE <- 2000000  # 2km

# Load data
message("Loading data...")
master_data <- rio::import("./2_processed_data/complete_ant_data.xlsx")
prescreening_data <- rio::import("./5_outputs/screen_output/prescreening_bio6_spocc.xlsx")
gabi_data <- tryCatch(rio::import("./2_processed_data/gabi_antmaps_data_clean.csv"), error = function(e) NULL)

# Load environmental data
message("Loading environmental layers...")
current_bioclim <- terra::rast(list.files("./1_raw_data/bioclim/current", pattern = "\\.tif$", full.names = TRUE))
future_bioclim <- terra::rast(list.files("./1_raw_data/bioclim/future", pattern = "\\.tif$", full.names = TRUE))
sbio_stack <- terra::rast(list.files("./1_raw_data/bioclim/sbio", pattern = "\\.tif$", full.names = TRUE))
env_stack <- terra::rast(list.files("./1_raw_data/bioclim/env", pattern = "\\.tif$", full.names = TRUE))

current_predictors <- c(current_bioclim, sbio_stack, env_stack)
future_predictors <- c(future_bioclim, sbio_stack, env_stack)

message("Layers loaded - Current: ", terra::nlyr(current_predictors), " | Future: ", terra::nlyr(future_predictors))

# Load world map
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world <- sf::st_make_valid(world)

# Filter species
included_species <- prescreening_data %>% filter(decision == "INCLUDE") %>% pull(species)
norwegian_keys <- c(...)  # Add Norwegian keys here if needed

test_species <- master_data %>%
  filter(status == "ACCEPTED") %>%
  filter(name %in% included_species) %>%
  select(name, acceptedKey)

message("Processing ", nrow(test_species), " species")

# TEST SUBSET
test_species <- test_species[1:3, ]

# =============================================================================
# MAIN LOOP
# =============================================================================

for (i in seq_len(nrow(test_species))) {

  species_name <- test_species$name[i]
  species_key <- test_species$acceptedKey[i]

  message("\n", paste(rep("=", 80), collapse = ""))
  message("Processing species ", i, " of ", nrow(test_species))
  message("Species: ", species_name, " | Key: ", species_key)

  # Create folders
  species_folder <- file.path("./species", paste0(species_key, "_", gsub("[^A-Za-z0-9_(),]", "_", species_name)))
  sdm_folder <- file.path(species_folder, "SDM_maxnet")
  raster_folder <- file.path(sdm_folder, "rasters")
  tables_folder <- file.path(sdm_folder, "tables")

  dir.create(species_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(sdm_folder, showWarnings = FALSE)
  dir.create(raster_folder, showWarnings = FALSE)
  dir.create(tables_folder, showWarnings = FALSE)

  # Download occurrence data
  message("Downloading occurrence data...")
  occ_df <- tryCatch({
    spocc::occ(query = species_name, from = c("gbif", "inat", "ala"),
               limit = DOWNLOAD_LIMIT, has_coords = TRUE) %>%
      spocc::occ2df()
  }, error = function(e) {
    message("Error downloading data: ", e$message)
    data.frame()
  })

  if (nrow(occ_df) == 0) {
    message("No occurrence data found")
    next
  }

  # Add GABI data if available
  if (!is.null(gabi_data)) {
    gabi_species <- gabi_data %>% filter(grepl(gsub(" .*", "", species_name), valid_species_name, ignore.case = TRUE))
    if (nrow(gabi_species) > 0) {
      gabi_occ <- gabi_species %>%
        select(species = valid_species_name, decimalLongitude = dec_long, decimalLatitude = dec_lat) %>%
        mutate(date = NA, source = "gabi", year = NA) %>%
        filter(!is.na(decimalLongitude) & !is.na(decimalLatitude))
      occ_df <- bind_rows(occ_df, gabi_occ)
    }
  }

  message("Raw records: ", nrow(occ_df))

  # Temporal filter
  if (sum(!is.na(occ_df$year)) > 0.5 * nrow(occ_df)) {
    occ_df <- occ_df %>% filter(is.na(year) | year >= RECENT_YEARS)
  }

  # Clean coordinates
  occ_df_clean <- tryCatch({
    occ_df %>%
      cc_val(verbose = TRUE) %>%
      cc_equ(verbose = TRUE) %>%
      cc_cap(verbose = TRUE) %>%
      cc_cen(verbose = TRUE) %>%
      cc_sea(verbose = TRUE, ref = NULL) %>%
      cc_zero(verbose = TRUE) %>%
      cc_dupl(verbose = TRUE) %>%
      cc_outl(method = "quantile", mltpl = 3, verbose = TRUE, value = "clean")
  }, error = function(e) {
    message("Cleaning failed, using minimal cleaning")
    occ_df %>% cc_zero(verbose = TRUE) %>% cc_dupl(verbose = TRUE)
  })

  if (nrow(occ_df_clean) < MIN_RECORDS) {
    message("Insufficient data after cleaning: ", nrow(occ_df_clean))
    next
  }

  # Spatial thinning
  occ_thinned <- tryCatch({
    thinData(coords = as.data.frame(occ_df_clean), env = current_predictors[[1]],
             x = "decimalLongitude", y = "decimalLatitude", verbose = FALSE)
  }, error = function(e) occ_df_clean)

  if (nrow(occ_thinned) < MIN_RECORDS) {
    message("Insufficient data after thinning: ", nrow(occ_thinned))
    next
  }

  message("Clean records: ", nrow(occ_thinned))

  # Background selection (circles method with fallback to random global)
  p_coords <- occ_thinned[, c("decimalLongitude", "decimalLatitude")]

  bg_coords <- tryCatch({
    message("Attempting circle-based background selection...")

    presence_sf <- sf::st_as_sf(p_coords, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
    presence_projected <- sf::st_transform(presence_sf, crs = 3857)
    circles <- sf::st_buffer(presence_projected, dist = BUFFER_DISTANCE)

    # Turn off S2 spherical geometry to avoid edge crossing errors with complex polygons
    sf::sf_use_s2(FALSE)

    merged_polygon <- sf::st_union(circles)
    merged_polygon <- sf::st_simplify(merged_polygon, dTolerance = 1000)
    merged_polygon_geo <- sf::st_transform(merged_polygon, crs = 4326)

    # Re-enable S2
    sf::sf_use_s2(TRUE)

    bg_points <- sf::st_sample(merged_polygon_geo, size = 20000, type = "random")
    bg_coords_land <- sf::st_filter(sf::st_sf(geometry = bg_points), world)
    bg_coords_temp <- sf::st_coordinates(bg_coords_land)
    bg_coords_temp <- as.data.frame(bg_coords_temp[, c("X", "Y")])
    colnames(bg_coords_temp) <- c("x", "y")

    message("Circle-based background selection successful: ", nrow(bg_coords_temp), " points")
    bg_coords_temp

  }, error = function(e) {
    message("Circle-based background selection FAILED: ", e$message)
    message("Falling back to random global background points...")

    # Fallback: Random points from entire land surface
    sf::sf_use_s2(TRUE)  # Ensure S2 is back on

    bg_points_global <- sf::st_sample(world, size = 20000, type = "random")
    bg_coords_global <- sf::st_coordinates(bg_points_global)
    bg_coords_global <- as.data.frame(bg_coords_global[, c("X", "Y")])
    colnames(bg_coords_global) <- c("x", "y")

    message("Random global background selection successful: ", nrow(bg_coords_global), " points")
    bg_coords_global
  })

  # Validate coordinates
  valid_p <- terra::extract(current_predictors, as.matrix(p_coords), cells = TRUE)
  p_coords <- p_coords[!is.na(valid_p[, 1]), ]

  valid_bg <- terra::extract(current_predictors, as.matrix(bg_coords), cells = TRUE)
  bg_coords <- bg_coords[!is.na(valid_bg[, 1]), ]

  message("Presence: ", nrow(p_coords), " | Background: ", nrow(bg_coords))

  if (nrow(p_coords) < MIN_RECORDS || nrow(bg_coords) < 100) {
    message("Insufficient valid coordinates")
    next
  }

  # Prepare SWD
  data_swd <- tryCatch({
    prepareSWD(species = species_name, p = p_coords, a = bg_coords, env = current_predictors)
  }, error = function(e) {
    message("SWD preparation failed: ", e$message)
    NULL
  })

  if (is.null(data_swd)) next

  # Split data
  c(train_data, test_data) %<-% trainValTest(data_swd, test = 0.2, only_presence = FALSE, seed = 123)

  # Variable selection
  message("Variable selection...")
  background <- addSamplesToBg(train_data, all = TRUE)

  selected_vars_model <- tryCatch({
    varSel(model = train(method = "Maxnet", data = train_data), metric = "auc", test = test_data,
           bg4cor = background, method = "spearman", cor_th = 0.7,
           env = current_predictors, use_pc = FALSE, progress = TRUE, permut = 10)
  }, error = function(e) {
    message("Variable selection failed, trying relaxed threshold...")
    tryCatch({
      varSel(model = train(method = "Maxnet", data = train_data), metric = "auc", test = test_data,
             bg4cor = background, method = "spearman", cor_th = 0.85,
             env = current_predictors, use_pc = FALSE, progress = TRUE, permut = 5)
    }, error = function(e2) NULL)
  })

  if (is.null(selected_vars_model)) {
    message("Variable selection failed completely")
    next
  }

  # Hyperparameter optimization - SIMPLIFIED (MaxEnt best practices)
  message("Optimizing hyperparameters...")

  h <- list(
    reg = seq(1, 5, 0.5),  # Conservative regularization (1-5)
    fc = "lqh"             # Linear + Quadratic + Hinge (sufficient for most distributions)
  )

  optimized_model <- tryCatch({
    optimizeModel(model = selected_vars_model@models[[1]], hypers = h, metric = "auc",
                  test = test_data, pop = 15, gen = 15, seed = 123)
  }, error = function(e) {
    message("Optimization failed: ", e$message)
    NULL
  })

  if (is.null(optimized_model)) next

  final_model <- optimized_model@models[[1]]

  # Predictions
  message("Making predictions...")

  pred_current <- predict(final_model, data = current_predictors, type = "cloglog")
  pred_future <- predict(final_model, data = future_predictors, type = "cloglog")
  change <- pred_future - pred_current

  # Save rasters
  terra::writeRaster(pred_current, file.path(raster_folder, paste0("current_", species_key, ".tif")), overwrite = TRUE)
  terra::writeRaster(pred_future, file.path(raster_folder, paste0("future_", species_key, ".tif")), overwrite = TRUE)
  terra::writeRaster(change, file.path(raster_folder, paste0("change_", species_key, ".tif")), overwrite = TRUE)

  # Evaluation
  auc_test <- auc(final_model, test = test_data)
  tss_test <- tss(final_model, test = test_data)

  eval_df <- data.frame(
    Species = species_name,
    AUC = auc_test,
    TSS = tss_test,
    N_presence = nrow(p_coords),
    N_background = nrow(bg_coords),
    N_variables = length(final_model@presence@pa)
  )

  rio::export(eval_df, file.path(tables_folder, paste0("evaluation_", species_key, ".xlsx")))

  # Variable importance
  vi <- varImp(final_model)
  rio::export(vi, file.path(tables_folder, paste0("variable_importance_", species_key, ".xlsx")))

  # Jackknife
  jk <- doJk(final_model, metric = "auc", test = test_data, with_only = TRUE)
  rio::export(jk, file.path(tables_folder, paste0("jackknife_", species_key, ".xlsx")))

  # Model report
  modelReport(final_model, type = "cloglog",
              folder = file.path(sdm_folder, paste0("report_", species_key)),
              test = test_data, response_curves = TRUE, only_presence = TRUE,
              jk = TRUE, env = current_predictors)

  message("Species ", species_name, " completed successfully!")

  gc()
}

message("\n", paste(rep("=", 80), collapse = ""))
message("ALL SPECIES COMPLETED!")
