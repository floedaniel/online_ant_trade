# ================================================================================
# ================================================================================
#
#                    SDM PRODUCTION PIPELINE - MaxEnt/SDMtune
#
#                    Norwegian Ant Species Distribution Modeling
#
# ================================================================================
# ================================================================================
#
# Script: 10_updated_sdmtune_loop.R
# Purpose: Production SDM pipeline with full features and comprehensive outputs
# Author: SDM Pipeline Project
# Last Updated: 2026-01-16
#
# Description:
# Complete species distribution modeling pipeline using MaxEnt (via SDMtune).
# Predicts current and future distributions under climate change (SSP585 2021-2040).
# Includes robust error handling, automatic fallbacks, and comprehensive outputs.
#
# Key Features:
# - Circle-based background selection with automatic fallback
# - Two-stage variable selection (importance + correlation filtering)
# - Conservative hyperparameter optimization (H, Q, HQ features; reg 2-10)
# - Cross-validation with uncertainty quantification
# - Clamped and unclamped predictions
# - MESS analysis for extrapolation detection
# - 20+ visualizations per species
#
# ================================================================================

# ================================================================================
# CHAPTER 1: INITIALIZATION
# ================================================================================

# ---------- 1.1 Load Required Libraries ----------
# Core SDM and data manipulation packages
library(spocc)              # Multi-source occurrence data download
library(rgbif)              # GBIF taxonomic backbone
library(tidyverse)          # Data manipulation
library(CoordinateCleaner)  # Geographic data quality control
library(SDMtune)            # MaxEnt modeling and optimization
library(terra)              # Raster processing (modern replacement for raster)
library(raster)             # Legacy raster (required by dismo::mess)
library(sf)                 # Spatial vector data
library(ggplot2)            # Visualization
library(tidyterra)          # terra integration with ggplot2
library(viridis)            # Color scales
library(zeallot)            # Multiple assignment operator %<-%
library(rio)                # Flexible data import/export
library(rnaturalearth)      # World country polygons
library(writexl)            # Excel export
library(lwgeom)             # Advanced geometry operations
library(rnaturalearthdata)  # Natural Earth data
library(dismo)              # MESS analysis
library(patchwork)          # Multi-panel plot composition

# ---------- 1.2 Configuration Parameters ----------
MIN_RECORDS <- 50           # Minimum presence records required after cleaning
DOWNLOAD_LIMIT <- 100000    # Maximum records per data source (spocc)
CV_FOLDS <- 5               # Cross-validation folds (increased from 2 for robustness)
RECENT_YEARS <- 1960        # Temporal filter: only records from 1960 onwards

# Background sampling buffer distance (meters, projected EPSG:3857)
# Rationale: 3000 km balances biogeographic realism with computational efficiency
# 200 km too restrictive (causes variable selection/optimization failures)
BUFFER_DISTANCE <- 3000000  # 3000 km radius around occurrence points

# Variable pre-filtering: Removes correlated SBIO/ENV layers before modeling
# Rationale: Keeps all BIO variables (climate), removes redundant soil/environmental layers
USE_VARIABLE_FILTER <- TRUE  # Set FALSE to disable pre-filtering


# ================================================================================
# CHAPTER 2: DATA LOADING AND SPECIES FILTERING
# ================================================================================

# ---------- 2.1 Load Master Data and Prescreening Results ----------
message("Loading master data and prescreening results...")

master_file <- "./2_processed_data/complete_ant_data.xlsx"
master_data <- rio::import(master_file)

prescreening_file <- "./5_outputs/screen_output/prescreening_bio6_spocc.xlsx"
prescreening_data <- rio::import(prescreening_file)

# ---------- 2.2 Load Local Occurrence Data (GABI AntMaps) ----------
# Supplements online data sources with curated local museum/collection records
message("Loading local occurrence data (GABI AntMaps)...")

gabi_file <- "./2_processed_data/gabi_antmaps_data_clean.csv"

gabi_data <- tryCatch({
  rio::import(gabi_file)
}, error = function(e) {
  message("Warning: Could not load GABI data: ", e$message)
  NULL
})

if (!is.null(gabi_data)) {
  message("Loaded ", nrow(gabi_data), " records from GABI AntMaps database")
} else {
  message("No local GABI data available - will use only online sources")
}

# Filter prescreening data for INCLUDE decisions
included_species <- prescreening_data %>%
  filter(decision == "INCLUDE") %>%
  pull(species)  # Assuming the column with species names is 'name'

message(sprintf("Found %d species marked as INCLUDE in prescreening", length(included_species)))

# ---------- 2.3 Norwegian Native Species Filter ----------
# Excludes Norwegian native species from modeling (focus on non-native species)
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
      gbif_name = match$canonicalName %||% NA_character_,
      usageKey = match$usageKey,
      acceptedKey = match$acceptedUsageKey %||% match$usageKey,
      status = match$status
    )
  } else {
    tibble(
      input_name = sp,
      gbif_name = NA_character_,
      usageKey = NA_integer_,
      acceptedKey = NA_integer_,
      status = "NOT_FOUND"
    )
  }
})

message(sprintf("Found GBIF keys for %d / %d Norwegian species",
                sum(!is.na(norwegian_keys$usageKey)), length(norwegian_ant_species)))

norwegian_accepted_keys <- norwegian_keys %>%
  filter(!is.na(acceptedKey)) %>%
  pull(acceptedKey) %>%
  unique()

# ---------- 2.4 Final Species List Assembly ----------
# Filter criteria:
# 1. ACCEPTED taxonomic status (removes synonyms, invalid names)
# 2. NOT in Norwegian native species list
# 3. INCLUDED in prescreening (passed data quality checks)
test_species <- master_data %>%
  filter(status == "ACCEPTED") %>%
  filter(!acceptedKey %in% norwegian_accepted_keys) %>%
  filter(name %in% included_species) %>%
  dplyr::select(name, acceptedKey)

message(sprintf("Final species list: %d species to process", nrow(test_species)))


# ================================================================================
# TESTING SUBSET (OPTIONAL)
# ================================================================================
# Uncomment one of the lines below to test pipeline on a small subset
# test_species <- test_species[1:5, ]           # First 5 species
# test_species <- test_species[20:25, ]         # Species 20-25
# test_species <- test_species %>% sample_n(10) # Random 10 species
# ================================================================================


# ================================================================================
# CHAPTER 3: ENVIRONMENTAL DATA PREPARATION
# ================================================================================

# ---------- 3.1 Load Environmental Layers ----------
# Loads bioclim (BIO1-BIO19), soil bioclim (SBIO), and environmental (ENV) layers
message("Loading environmental data...")

# Load current bioclim
current_predictors_dir <- "./1_raw_data/bioclim/current"
current_predictors_files <- list.files(current_predictors_dir, pattern = "\\.tif$", full.names = TRUE)
message("Found ", length(current_predictors_files), " current bioclim files")
if (length(current_predictors_files) == 0) {
  stop("ERROR: No current bioclim files found in ", current_predictors_dir)
}
current_predictors <- terra::rast(current_predictors_files)
message("Loaded current bioclim: ", terra::nlyr(current_predictors), " layers")

# Load future bioclim
future_predictors_dir <- "./1_raw_data/bioclim/future"
future_predictors_files <- list.files(future_predictors_dir, pattern = "\\.tif$", full.names = TRUE)
message("Found ", length(future_predictors_files), " future bioclim files")
if (length(future_predictors_files) == 0) {
  stop("ERROR: No future bioclim files found in ", future_predictors_dir)
}
future_predictors <- terra::rast(future_predictors_files)
message("Loaded future bioclim: ", terra::nlyr(future_predictors), " layers")

# Load SBIO layers
sbio_dir <- "./1_raw_data/bioclim/sbio"
sbio_files <- list.files(sbio_dir, pattern = "\\.tif$", full.names = TRUE)
message("Found ", length(sbio_files), " SBIO files")
if (length(sbio_files) == 0) {
  warning("WARNING: No SBIO files found in ", sbio_dir)
  sbio_stack <- NULL
} else {
  sbio_stack <- terra::rast(sbio_files)
  message("Loaded SBIO layers: ", terra::nlyr(sbio_stack), " layers")
}

# Load environmental layers
env_dir <- "./1_raw_data/bioclim/env"
env_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)
message("Found ", length(env_files), " environmental files")
if (length(env_files) == 0) {
  warning("WARNING: No environmental files found in ", env_dir)
  env_stack <- NULL
} else {
  env_stack <- terra::rast(env_files)
  message("Loaded environmental layers: ", terra::nlyr(env_stack), " layers")
}

# Combine all layers
if (!is.null(sbio_stack) && !is.null(env_stack)) {
  current_predictors <- c(current_predictors, sbio_stack, env_stack)
  future_predictors <- c(future_predictors, sbio_stack, env_stack)
} else if (!is.null(sbio_stack)) {
  current_predictors <- c(current_predictors, sbio_stack)
  future_predictors <- c(future_predictors, sbio_stack)
} else if (!is.null(env_stack)) {
  current_predictors <- c(current_predictors, env_stack)
  future_predictors <- c(future_predictors, env_stack)
}

message("Total layers - Current: ", terra::nlyr(current_predictors), ", Future: ", terra::nlyr(future_predictors))

# ---------- 3.2 Variable Pre-Filtering (Optional) ----------
# Reduces computational load by removing correlated SBIO/ENV variables
# Keeps all BIO variables (core climate data)
if (USE_VARIABLE_FILTER) {
  message("\nApplying variable filter (Keep ALL BIO, remove correlated SBIO/ENV)...")

  # Check if filter file exists
  filter_file <- "./5_outputs/selected_variable_names.txt"

  if (!file.exists(filter_file)) {
    message("  Filter file not found. Running quick_variable_filter.R...")
    tryCatch({
      source("./4_scripts/quick_variable_filter.R")
    }, error = function(e) {
      warning("Could not run variable filter: ", e$message)
      warning("Continuing with all variables...")
      USE_VARIABLE_FILTER <<- FALSE
    })
  }

  if (USE_VARIABLE_FILTER && file.exists(filter_file)) {
    selected_vars <- readLines(filter_file)
    message("  Loaded ", length(selected_vars), " filtered variables")

    # Check which variables are available
    available_vars <- names(current_predictors)
    missing_vars <- setdiff(selected_vars, available_vars)

    if (length(missing_vars) > 0) {
      warning("  ", length(missing_vars), " selected variables not found, using available ones")
      selected_vars <- intersect(selected_vars, available_vars)
    }

    # Filter to selected variables
    current_predictors <- current_predictors[[selected_vars]]
    future_predictors <- future_predictors[[selected_vars]]

    message("  Using ", length(selected_vars), " filtered variables")
  }
} else {
  message("\nUsing all variables (no filtering)")
}

# ---------- 3.3 Load Reference Spatial Data ----------
# Natural Earth country polygons for background sampling and visualization
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, crs = 4326)
world <- subset(world, sovereignt != "Antarctica")


# ================================================================================
# CHAPTER 4: MAIN SPECIES PROCESSING LOOP
# ================================================================================

# ---------- 4.1 Initialize Output Directory and Timing ----------
base_dir <- "./species"
if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)
overall_start_time <- Sys.time()
message("\n", paste(rep("=", 80), collapse = ""))
message("STARTING SDM PIPELINE")
message("Start time: ", format(overall_start_time, "%Y-%m-%d %H:%M:%S"))
message(paste(rep("=", 80), collapse = ""))

for (i in seq_len(nrow(test_species))) {

  # ================================================================================
  # CHAPTER 4.2: SPECIES SETUP
  # ================================================================================

  species_name <- test_species$name[i]
  species_key <- test_species$acceptedKey[i]
  species_start_time <- Sys.time()

  message("\n", paste(rep("=", 80), collapse = ""))
  message("Processing species ", i, " of ", nrow(test_species))
  message("Species: ", species_name, " | Key: ", species_key)
  message("Started: ", format(species_start_time, "%H:%M:%S"))

  # ---------- 4.2.1 Create Output Directory Structure ----------
  # Format: species/[acceptedKey]_[species_name]/SDM_maxnet/
  safe_species_name <- gsub(" ", "_", species_name)
  species_folder_name <- paste0(species_key, "_", safe_species_name)
  species_folder <- file.path(base_dir, species_folder_name)
  sdm_folder <- file.path(species_folder, "SDM_maxnet")

  # Create subfolders for organized output
  raster_folder <- file.path(sdm_folder, "rasters")
  tables_folder <- file.path(sdm_folder, "tables")

  # Skip if SDM folder already exists and contains output files
  if (dir.exists(sdm_folder)) {
    message("SDM folder already exists - SKIPPING")
    next
  }

  if (!dir.exists(species_folder)) dir.create(species_folder, recursive = TRUE)
  if (!dir.exists(sdm_folder)) dir.create(sdm_folder, recursive = TRUE)
  if (!dir.exists(raster_folder)) dir.create(raster_folder, recursive = TRUE)
  if (!dir.exists(tables_folder)) dir.create(tables_folder, recursive = TRUE)


  # ================================================================================
  # CHAPTER 4.3: OCCURRENCE DATA ACQUISITION
  # ================================================================================

  # ---------- 4.3.1 Download Online Occurrence Data ----------
  # Sources: GBIF, BISON, ALA (Australia), iNaturalist, iDigBio
  message("Downloading occurrence data for: ", species_name)

  occ_data <- tryCatch({
    spocc::occ(
      query = species_name,
      from = c("gbif", "bison", "ala", "inat", "idigbio"),
      limit = DOWNLOAD_LIMIT
    )
  }, error = function(e) {
    message("Error downloading occurrence data: ", e$message)
    NULL
  })

  if (is.null(occ_data)) {
    message("Failed to download occurrence data for: ", species_name)
    writeLines(paste("Failed to download occurrence data:", Sys.time()),
               file.path(sdm_folder, "error_log.txt"))
    next
  }

  occ_df <- occ2df(occ_data)
  occ_df <- dplyr::select(occ_df, name, longitude, latitude, date, prov)
  names(occ_df) <- c("species", "decimalLongitude", "decimalLatitude", "date", "source")
  occ_df$species <- species_name

  occ_df$year <- NA

  if (!all(is.na(occ_df$date))) {
    occ_df$year <- as.numeric(substr(as.character(occ_df$date), 1, 4))
  }

  message("Raw occurrence records from online sources: ", nrow(occ_df))

  # ---------- 4.3.2 Add Local GABI AntMaps Data ----------
  # Supplements online records with curated museum/collection data
  gabi_records <- 0

  if (!is.null(gabi_data)) {
    message("Searching for species in local GABI database...")

    # Clean species name - remove author citation (text in parentheses)
    species_name_clean <- gsub("\\s*\\([^)]+\\)\\s*$", "", species_name)
    species_name_clean <- trimws(species_name_clean)

    message("Matching species name: '", species_name_clean, "'")

    # Filter GABI data for current species (case-insensitive matching)
    gabi_species <- gabi_data %>%
      filter(tolower(trimws(valid_species_name)) == tolower(species_name_clean))

    if (nrow(gabi_species) > 0) {
      message("Found ", nrow(gabi_species), " records in GABI database")

      # Select and rename columns to match occ_df structure
      gabi_occ <- gabi_species %>%
        dplyr::select(
          species = valid_species_name,
          decimalLongitude = dec_long,
          decimalLatitude = dec_lat
        ) %>%
        mutate(
          date = NA,
          source = "gabi_antmaps",
          year = NA
        )

      # Remove records with missing coordinates
      gabi_occ <- gabi_occ %>%
        filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) 
    
      message("GABI records with valid coordinates: ", nrow(gabi_occ))
      
      # make sure all are as.character before bind_rows
      char_cols <- c("species", "decimalLongitude", "decimalLatitude")
      occ_df <- occ_df %>%  mutate(across(all_of(char_cols), as.character))
      gabi_occ <- gabi_occ %>%  mutate(across(all_of(char_cols), as.character))

      # Combine online and local data
      occ_df <- bind_rows(occ_df, gabi_occ)
      gabi_records <- nrow(gabi_occ)

      message("Total records after adding GABI data: ", nrow(occ_df))
    } else {
      message("No matching records found in GABI database")
    }
  }

  message("Raw occurrence records (total): ", nrow(occ_df),
          " (online: ", nrow(occ_df) - gabi_records, ", GABI: ", gabi_records, ")")

  # Check if we have exactly the download limit - this indicates incomplete data
  if (nrow(occ_df) <= DOWNLOAD_LIMIT && nrow(occ_df) < 200) {
    message("WARNING: Only ", nrow(occ_df), " records downloaded (limit: ", DOWNLOAD_LIMIT, ")")
    message("This may indicate insufficient data or download issues. Skipping species.")
    writeLines(paste("Insufficient occurrence data:", nrow(occ_df), "records (at or below download limit)"),
               file.path(sdm_folder, "insufficient_data.txt"))
    next
  }

  if (nrow(occ_df) < MIN_RECORDS) {
    message("Insufficient occurrence data (< ", MIN_RECORDS, " records).")
    writeLines(paste("Insufficient occurrence data:", nrow(occ_df), "records"),
               file.path(sdm_folder, "insufficient_data.txt"))
    next
  }


  # ================================================================================
  # CHAPTER 4.4: DATA QUALITY CONTROL
  # ================================================================================

  # ---------- 4.4.1 Temporal Filtering ----------
  # Rationale: Reduces sampling bias from historical vs. modern collection efforts
  # Reference: Fithian et al. 2015 (bias reduction in presence-only models)
  message("Applying temporal filtering...")

  if (sum(!is.na(occ_df$year)) > 0.5 * nrow(occ_df)) {
    occ_df_temporal <- occ_df %>% filter(is.na(year) | year >= RECENT_YEARS)

    recent_with_year <- occ_df_temporal %>% filter(!is.na(year))
    if (nrow(recent_with_year) > 20) {
      # Safe temporal bias detection with error handling
      tryCatch({
        temporal_bias_lat <- lm(decimalLatitude ~ year, data = recent_with_year)
        temporal_bias_lon <- lm(decimalLongitude ~ year, data = recent_with_year)

        # Check if coefficients exist before accessing
        if (!is.null(temporal_bias_lat) && !is.null(temporal_bias_lon)) {
          lat_coef <- summary(temporal_bias_lat)$coefficients
          lon_coef <- summary(temporal_bias_lon)$coefficients

          if (nrow(lat_coef) >= 2 && nrow(lon_coef) >= 2) {
            if (lat_coef[2, 4] < 0.05 | lon_coef[2, 4] < 0.05) {
              warning("Temporal bias detected - interpret with caution")
            }
          }
        }
      }, error = function(e) {
        message("Could not assess temporal bias: ", e$message)
      })
    }
  } else {
    occ_df_temporal <- occ_df
  }

  message("Records after temporal filtering: ", nrow(occ_df_temporal))

  # ---------- 4.4.2 Coordinate Cleaning ----------
  # Uses CoordinateCleaner package (Zizka et al. 2019)
  # Sequential filters: invalid coords, equal lat/lon, capitals, centroids, sea, zeros, duplicates, outliers
  message("Cleaning occurrence data...")

  occ_df_clean <- tryCatch({
    occ_df_temporal %>%
      cc_val(verbose = TRUE) %>%
      cc_equ(verbose = TRUE) %>%
      cc_cap(verbose = TRUE) %>%
      cc_cen(verbose = TRUE) %>%
      cc_sea(verbose = TRUE) %>%  # ref = NULL uses default reference
      cc_zero(verbose = TRUE) %>%
      cc_dupl(verbose = TRUE) %>%
      cc_outl(method = "quantile", mltpl = 3, verbose = TRUE, value = "clean")
  }, error = function(e) {
    message("Error during coordinate cleaning: ", e$message)
    message("Attempting simplified cleaning without cc_sea...")
    # Fallback: skip cc_sea if it fails
    tryCatch({
      occ_df_temporal %>%
        cc_val(verbose = TRUE) %>%
        cc_equ(verbose = TRUE) %>%
        cc_cap(verbose = TRUE) %>%
        cc_cen(verbose = TRUE) %>%
        cc_zero(verbose = TRUE) %>%
        cc_dupl(verbose = TRUE) %>%
        cc_outl(method = "quantile", mltpl = 3, verbose = TRUE, value = "clean")
    }, error = function(e2) {
      message("Simplified cleaning also failed: ", e2$message)
      message("Using minimal cleaning (duplicates and zeros only)...")
      # Minimal fallback
      occ_df_temporal %>%
        cc_zero(verbose = TRUE) %>%
        cc_dupl(verbose = TRUE)
    })
  })

  if (is.null(occ_df_clean)) {
    message("Failed to clean occurrence data for: ", species_name)
    writeLines(paste("Failed during coordinate cleaning:", Sys.time()),
               file.path(sdm_folder, "error_log.txt"))
    next
  }

  message("Records after cleaning: ", nrow(occ_df_clean))

  if (nrow(occ_df_clean) < MIN_RECORDS) {
    message("Insufficient data after cleaning (< ", MIN_RECORDS, " records).")
    writeLines(paste("Insufficient data after cleaning:", nrow(occ_df_clean), "records"),
               file.path(sdm_folder, "insufficient_data.txt"))
    next
  }

  # ---------- 4.4.3 Spatial Thinning ----------
  # Removes duplicate records within same raster cell
  # Rationale: Reduces spatial autocorrelation (violates MaxEnt independence assumption)
  message("Applying spatial thinning...")

  occ_coords_before_thin <- occ_df_clean[, c("decimalLongitude", "decimalLatitude")]

  occ_thinned <- tryCatch({
    # Use first layer for spatial thinning (just needs cell alignment, not specific variable)
    # Ensure coords is a clean data.frame (not tibble or other class)
    coords_clean <- as.data.frame(occ_df_clean)

    # Ensure we have a single-layer raster
    env_layer <- current_predictors[[1]]
    if (terra::nlyr(env_layer) > 1) {
      env_layer <- env_layer[[1]]
    }

    thinData(
      coords = coords_clean,
      env = env_layer,
      x = "decimalLongitude",
      y = "decimalLatitude",
      verbose = FALSE
    )
  }, error = function(e) {
    message("Error during spatial thinning: ", e$message)
    message("Skipping spatial thinning and continuing with unthinned data")
    occ_df_clean
  })

  if (!is.null(occ_thinned) && nrow(occ_thinned) > 0) {
    occ_df_clean <- as.data.frame(occ_thinned)
    message("Records after spatial thinning: ", nrow(occ_df_clean),
            " (removed ", nrow(occ_coords_before_thin) - nrow(occ_df_clean), " duplicates)")
  } else {
    message("Spatial thinning produced no results, keeping original data")
  }

  if (nrow(occ_df_clean) < MIN_RECORDS) {
    message("Insufficient data after spatial thinning (< ", MIN_RECORDS, " records).")
    writeLines(paste("Insufficient data after spatial thinning:", nrow(occ_df_clean), "records"),
               file.path(sdm_folder, "insufficient_data.txt"))
    next
  }


  # ================================================================================
  # CHAPTER 4.5: BACKGROUND POINT SELECTION
  # ================================================================================

  # ---------- 4.5.1 Circle-Based Background Selection with Automatic Fallback ----------
  # Primary method: 3000 km circular buffers around occurrence points
  # Rationale: Biogeographically-informed background (avoids convex hull overfitting)
  # Fallback: Random global background if circle method fails (geometry errors, S2 issues)
  #
  # S2 Spherical Geometry Handling:
  # - S2 uses strict validity rules for spherical coordinates (sf >= 1.0)
  # - Overlapping buffers create edge intersections violating S2 validity
  # - Solution: Disable S2 during union/simplify/transform, re-enable after
  # - Reference: r-spatial/sf Issue #1710
  message("Generating background points...")

  p_coords <- occ_df_clean[, c("decimalLongitude", "decimalLatitude")]

  bg_coords <- tryCatch({
    message("Attempting circle-based background selection...")

    # Convert to sf and create buffers
    presence_sf <- sf::st_as_sf(p_coords, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
    presence_projected <- sf::st_transform(presence_sf, crs = 3857)
    circles <- sf::st_buffer(presence_projected, dist = BUFFER_DISTANCE)

    # Disable S2 spherical geometry temporarily
    sf::sf_use_s2(FALSE)

    # Merge and simplify
    merged_polygon <- sf::st_union(circles)
    merged_polygon <- sf::st_simplify(merged_polygon, dTolerance = 1000)
    merged_polygon_geo <- sf::st_transform(merged_polygon, crs = 4326)

    # Re-enable S2
    sf::sf_use_s2(TRUE)

    # Sample points and filter to land
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

  # Validate presence coordinates
  valid_coords <- terra::extract(current_predictors, as.matrix(p_coords), cells = TRUE)
  p_coords <- p_coords[!is.na(valid_coords[, 1]), ]

  # Validate background coordinates
  valid_bg <- terra::extract(current_predictors, as.matrix(bg_coords), cells = TRUE)
  bg_coords <- bg_coords[!is.na(valid_bg[, 1]), ]

  message("Presence: ", nrow(p_coords), " | Background: ", nrow(bg_coords))

  if (nrow(p_coords) < MIN_RECORDS || nrow(bg_coords) < 100) {
    message("Insufficient valid coordinates for modeling")
    writeLines(paste("Insufficient valid coordinates - Presence:", nrow(p_coords), "Background:", nrow(bg_coords)),
               file.path(sdm_folder, "insufficient_data.txt"))
    next
  }


  # ================================================================================
  # CHAPTER 4.6: DATA VISUALIZATION (OCCURRENCE & BACKGROUND)
  # ================================================================================

  # ---------- 4.6.1 Data Inspection Plot ----------
  # Visualizes occurrence points and background sampling strategy
  background_plot <- tryCatch({
    ggplot(world) +
      geom_sf(fill = "gray90", color = "gray70", linewidth = 0.2) +
      geom_point(data = bg_coords, aes(x = x, y = y), color = "black", alpha = 0.5, size = 0.5) +
      geom_point(data = p_coords, aes(x = decimalLongitude, y = decimalLatitude), color = "#cc4778", size = 1) +
      coord_sf(crs = 4326, expand = FALSE) +
      labs(title = paste(species_name, "- Occurrence & Background"),
           x = "Longitude", y = "Latitude") +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "aliceblue", color = NA),
            panel.grid = element_line(color = "white", linewidth = 0.3))
  }, error = function(e) {
    message("Warning: Could not create background plot with sf: ", e$message)
    message("Creating simplified version without world map...")
    # Fallback without world polygon
    ggplot() +
      geom_point(data = bg_coords, aes(x = x, y = y), color = "black", alpha = 0.5, size = 0.5) +
      geom_point(data = p_coords, aes(x = decimalLongitude, y = decimalLatitude), color = "#cc4778", size = 1) +
      labs(title = paste(species_name, "- Occurrence & Background"),
           x = "Longitude", y = "Latitude") +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "aliceblue", color = NA))
  })

background_plot

ggsave(file.path(sdm_folder, paste0("data_inspection_", species_key, ".png")), background_plot, width = 10, height = 8, dpi = 300)

  # ---------- 4.6.2 Distribution Map ----------
  # Shows occurrence points only (without background)
distribution_plot <-  ggplot(world) +
  geom_sf(fill = "gray90", color = "gray70", linewidth = 0.2) +
    geom_point(data = p_coords, aes(x = decimalLongitude, y = decimalLatitude), color = "#bc3754", size = 1) +
    coord_sf(crs = 4326, expand = FALSE) +
    labs(title = paste(species_name, "distribution"), x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "aliceblue", color = NA),
          panel.grid = element_line(color = "white", linewidth = 0.3))

distribution_plot

ggsave(file.path(sdm_folder, paste0("distribution_", species_key, ".png")), distribution_plot, width = 10, height = 8, dpi = 300)


  # ================================================================================
  # CHAPTER 4.7: MODEL DATA PREPARATION
  # ================================================================================

  # ---------- 4.7.1 Create SWD Object and Train/Test Split ----------
  # SWD = Samples With Data (SDMtune format combining coordinates + environmental values)
  # Split: 80% training, 20% testing (presence-only split)
  message("Preparing data...")

  data_swd <- tryCatch({
    prepareSWD(species = species_name, p = p_coords, a = bg_coords, env = current_predictors)
  }, error = function(e) {
    message("Error preparing SWD: ", e$message)
    NULL
  })

  if (is.null(data_swd)) {
    message("Failed to prepare SWD for: ", species_name)
    writeLines(paste("Failed to prepare SWD:", Sys.time()),
               file.path(sdm_folder, "error_log.txt"))
    next
  }

  # Add presence to background (optional - mirrors Maxent behavior)
  data_swd <- addSamplesToBg(data_swd, all = FALSE)
  message("Background locations after adding unique presence: ", sum(data_swd@pa == 0))

  c(train_data, test_data) %<-% trainValTest(data_swd, test = 0.2, only_presence = TRUE, seed = 25)


  # ================================================================================
  # CHAPTER 4.8: VARIABLE SELECTION (TWO-STAGE)
  # ================================================================================

  # ---------- 4.8.1 Stage 1: Permutation Importance (Top 15 Variables) ----------
  # Reduces computational load by pre-filtering to most important variables
  # Uses initial HQ model (feature class hinge+quadratic, reg=1)
  message("Variable selection...")

  # Step 1: Variable importance
  initial_model <- tryCatch({
    train(method = "Maxnet", data = train_data, fc = "hq", reg = 1)
  }, error = function(e) {
    message("Error training initial model: ", e$message)
    NULL
  })

  if (is.null(initial_model)) {
    message("Failed to train initial model for: ", species_name)
    writeLines(paste("Failed to train initial model:", Sys.time()),
               file.path(sdm_folder, "error_log.txt"))
    next
  }

vi <- varImp(initial_model, permut = 5)

# Save variable importance
rio::export(vi, file.path(tables_folder, paste0("initial_variable_importance_", species_key, ".xlsx")))

  # Plot variable importance
  vi_plot <- plotVarImp(vi[, c("Variable", "Permutation_importance")], color = "steelblue")

  ggsave(file.path(sdm_folder, paste0("initial_variable_importance_", species_key, ".png")),
         vi_plot, width = 10, height = max(6, nrow(vi) * 0.3), dpi = 300)

  top_vars <- vi %>%
    arrange(desc(Permutation_importance)) %>%
    slice_head(n = min(15, nrow(vi))) %>%
    pull(Variable)

current_predictors_subset <- current_predictors[[top_vars]]

  # ---------- 4.8.2 Stage 2: Correlation-Based Removal ----------
  # Removes redundant correlated variables (Spearman ρ > 0.7)
  # Rationale: Reduces multicollinearity, prevents model over-complexity
  # Fallback: Relaxes threshold to 0.85 if primary filter fails
background <- prepareSWD(species = species_name, a = bg_coords, env = current_predictors_subset)
background <- addSamplesToBg(background, all = TRUE)

  # Save correlation matrix
  cor_matrix <- tryCatch({
    corVar(background, method = "spearman", cor_th = 0.7)
  }, error = function(e) {
    message("Could not compute correlation matrix: ", e$message)
    NULL
  })

  if (!is.null(cor_matrix) && nrow(cor_matrix) > 0) {
    rio::export(cor_matrix, file.path(tables_folder, paste0("correlated_variables_", species_key, ".xlsx")))
    message("Found ", nrow(cor_matrix), " pairs of correlated variables (>0.7)")
  }

  swd_subset <- prepareSWD(species = species_name, p = p_coords, a = bg_coords, env = current_predictors_subset)
  c(train_subset, test_subset) %<-% trainValTest(swd_subset, test = 0.2, only_presence = TRUE, seed = 25)

  selected_variables_model <- tryCatch({
    varSel(
      model = train(method = "Maxnet", data = train_subset),
      metric = "auc",
      test = test_subset,
      bg4cor = background,
      method = "spearman",
      cor_th = 0.7,
      env = current_predictors_subset,
      use_pc = FALSE,
      progress = TRUE,
      permut = 10
    )
  }, error = function(e) {
    message("Error during variable selection: ", e$message)
    message("Attempting variable selection with relaxed correlation threshold...")
    # Try with higher correlation threshold as fallback
    tryCatch({
      varSel(
        model = train(method = "Maxnet", data = train_subset),
        metric = "auc",
        test = test_subset,
        bg4cor = background,
        method = "spearman",
        cor_th = 0.85,  # More permissive
        env = current_predictors_subset,
        use_pc = FALSE,
        progress = TRUE,
        permut = 5  # Fewer permutations to save time
      )
    }, error = function(e2) {
      message("Variable selection failed even with relaxed threshold: ", e2$message)
      NULL
    })
  })

  if (is.null(selected_variables_model)) {
    message("Failed during variable selection for: ", species_name)
    writeLines(paste("Failed during variable selection:", Sys.time()),
               file.path(sdm_folder, "error_log.txt"))
    next
  }

  selected_vars <- names(selected_variables_model@data@data)
  final_predictors <- current_predictors_subset[[selected_vars]]

  message("Final variables (", length(selected_vars), "): ", paste(selected_vars, collapse = ", "))

  # Save selected variables
  selected_vars_df <- data.frame(
    variable = selected_vars,
    rank = seq_along(selected_vars)
  )
  rio::export(selected_vars_df, file.path(tables_folder, paste0("selected_variables_", species_key, ".xlsx")))

  # Get final variable importance for selected variables
  final_vi <- varImp(selected_variables_model, permut = 10)
  rio::export(final_vi, file.path(tables_folder, paste0("final_variable_importance_", species_key, ".xlsx")))

  # Plot final variable importance
  final_vi_plot <- plotVarImp(final_vi[, c("Variable", "Permutation_importance")], color = "darkgreen") 

  ggsave(file.path(sdm_folder, paste0("final_variable_importance_", species_key, ".png")),
         final_vi_plot, width = 10, height = max(6, nrow(final_vi) * 0.3), dpi = 300)

  # Prepare final SWD with selected variables
  final_swd <- prepareSWD(species = species_name, p = p_coords, a = bg_coords, env = final_predictors)
  c(train_final, test_final) %<-% trainValTest(final_swd, test = 0.2, only_presence = TRUE, seed = 25)


  # ================================================================================
  # CHAPTER 4.9: HYPERPARAMETER OPTIMIZATION
  # ================================================================================

  # ---------- 4.9.1 Define Hyperparameter Grid ----------
  # EPPO Best Practices & MaxEnt 3.4.0+ Guidelines:
  # - Threshold features omitted (not ecologically interpretable)
  # - Hinge (h): Smooth response curves, ecologically realistic
  # - Quadratic (q): Hump-shaped responses (classic ecological theory)
  # - Product (p) avoided: Risk of overfitting with small sample sizes
  #
  # Regularization: 2-10 by steps of 2 (conservative)
  # Rationale: Higher regularization prevents overfitting (Radosavljevic & Anderson 2014)
  # Appropriate for 50-1000 occurrence records
message("Optimizing hyperparameters...")

h <- list(reg = seq(2, 10, 2),  # 5 values: 2, 4, 6, 8, 10
          fc = c("h", "q", "hq"))  # 3 feature classes = 15 total combinations

  # ---------- 4.9.2 Genetic Algorithm Optimization with Fallbacks ----------
  # Primary: Genetic algorithm (pop=10, gen=10) searches 15 model configurations
  # Fallback 1: Expanded grid (all features, reg 0.5-5) if primary fails
  # Fallback 2: gridSearch (exhaustive) if fallback 1 fails
  optimized_model <- tryCatch({
    message("PRIMARY OPTIMIZATION: Testing FC=h,q,hq with reg=2,4,6,8,10...")
    result <- optimizeModel(
      model = train(method = "Maxnet", data = train_final),
      hypers = h,
      metric = "auc",
      test = test_final,
      pop = 10,  # Reduced population size to avoid "lower than population" error
      gen = 10,
      seed = 124,
      interactive = TRUE,
      progress = TRUE
    )
    message("PRIMARY OPTIMIZATION SUCCEEDED")
    result
  }, error = function(e) {
    message("PRIMARY OPTIMIZATION FAILED: ", e$message)
    message("FALLBACK 1: Retrying with expanded hyperparameter grid...")
    message("Testing all features (l,q,h,lq,lh,qh,lqp,lqh,lph,qph,lqph) with reg=0.5-5 for n=", nrow(p_coords))

      # Fallback: Full grid with all feature classes except threshold
      h_retry <- list(
        reg = seq(0.5, 5, 0.5),  # 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5 (10 values)
        fc = c(
          "l", "q", "h",
          "lq", "lh",  "qh", 
          "lqp", "lqh", "lph", "qph",
          "lqph"
        )
      )

      tryCatch({
        result <- optimizeModel(
          model = train(method = "Maxnet", data = train_final),
          hypers = h_retry,
          metric = "auc",
          test = test_final,
          pop = 15,
          gen = 15,
          seed = 124,
          interactive = TRUE,
          progress = TRUE
        )
        message("FALLBACK 1 SUCCEEDED")
        result
      }, error = function(e2) {
        message("FALLBACK 1 FAILED: ", e2$message)
        message("FALLBACK 2: Using gridSearch (exhaustive testing)...")
        message("Testing FC=h,lh,lqh,lqph with reg=1,2,3,4,5...")

        # Ultimate fallback: use gridSearch (tests all combinations)
        tryCatch({
          result <- gridSearch(
            model = train(method = "Maxnet", data = train_final),
            hypers = list(reg = seq(1, 5, 1), fc = c("h", "lh", "lqh", "lqph")),
            metric = "auc",
            test = test_final,
            interactive = TRUE,
            progress = TRUE
          )
          message("FALLBACK 2 (gridSearch) SUCCEEDED")
          result
        }, error = function(e3) {
          message("FALLBACK 2 FAILED: ", e3$message)
          message("ALL OPTIMIZATION METHODS FAILED - skipping species")
          NULL
        })
      })
  })
  
  if (is.null(optimized_model)) {
    message("Failed during hyperparameter optimization for: ", species_name)
    writeLines(paste("Failed during hyperparameter optimization:", Sys.time()),
               file.path(sdm_folder, "error_log.txt"))
    next
  }

  final_model <- optimized_model@models[[which.max(optimized_model@results$test_AUC)]]

  message("=========================================================")
  message("SELECTED BEST MODEL:")
  message("  Feature Class: ", final_model@model@fc)
  message("  Regularization: ", final_model@model@reg)
  message("  Test AUC: ", round(max(optimized_model@results$test_AUC), 4))
  message("=========================================================")

  # Save optimization results
  # Extract FC and reg parameters from each model
  fc_params <- sapply(optimized_model@models, function(m) m@model@fc)
  reg_params <- sapply(optimized_model@models, function(m) m@model@reg)

  opt_summary <- data.frame(
    model_rank = seq_along(optimized_model@results$train_AUC),
    fc = fc_params,
    reg = reg_params,
    train_AUC = optimized_model@results$train_AUC,
    test_AUC = optimized_model@results$test_AUC,
    diff_AUC = optimized_model@results$diff_AUC
  )

  # Sort by test_AUC descending so best model is at top
  opt_summary <- opt_summary[order(-opt_summary$test_AUC), ]

  rio::export(opt_summary, file.path(tables_folder, paste0("optimization_results_", species_key, ".xlsx")))
  message("Saved optimization results for ", nrow(opt_summary), " model configurations")
  message("Top 3 models by test_AUC:")
  for (i in 1:min(3, nrow(opt_summary))) {
    message("  ", i, ". FC=", opt_summary$fc[i], " reg=", opt_summary$reg[i],
            " | test_AUC=", round(opt_summary$test_AUC[i], 4))
  }


  # ================================================================================
  # CHAPTER 4.10: MODEL EVALUATION
  # ================================================================================

  # ---------- 4.10.1 Calculate Performance Metrics ----------
  # AUC: Area Under ROC Curve (0-1; >0.7 acceptable, >0.8 good)
  # TSS: True Skill Statistic (-1 to +1; >0.4 acceptable, >0.6 good)
  # Reference: Lobo et al. 2008 (AUC limitations), Allouche et al. 2006 (TSS)
  message("Calculating evaluation metrics...")

  auc_val <- SDMtune::auc(final_model, test = test_final)
  tss_val <- SDMtune::tss(final_model, test = test_final)

  ths <- thresholds(final_model, type = "cloglog")
  optimal_threshold <- ths[3, 2]

  message("AUC: ", round(auc_val, 3), " | TSS: ", round(tss_val, 3),
          " | Threshold: ", round(optimal_threshold, 3))

  # Calculate confusion matrix at optimal threshold
  cm <- confMatrix(final_model, test = test_final, th = optimal_threshold, type = "cloglog")
  cm_stats <- cm[1, ]  # Stats at optimal threshold

  # Calculate sensitivity and specificity
  sensitivity <- cm_stats$tp / (cm_stats$tp + cm_stats$fn)
  specificity <- cm_stats$tn / (cm_stats$tn + cm_stats$fp)

  message("Confusion Matrix - TP: ", cm_stats$tp, " | FP: ", cm_stats$fp,
          " | TN: ", cm_stats$tn, " | FN: ", cm_stats$fn)
  message("Sensitivity: ", round(sensitivity, 3), " | Specificity: ", round(specificity, 3))

  # ---------- 4.10.2 Cross-Validation ----------
  # k-fold cross-validation using optimized hyperparameters
  # Provides robust performance estimate across independent data subsets
  message("Performing ", CV_FOLDS, "-fold cross-validation...")

  folds <- randomFolds(final_swd, k = CV_FOLDS, only_presence = TRUE, seed = 321)

  cv_model <- tryCatch({
    train("Maxnet",
          data = final_swd,
          fc = final_model@model@fc,
          reg = final_model@model@reg,
          folds = folds)
  }, error = function(e) {
    message("Error during cross-validation: ", e$message)
    NULL
  })

  if (is.null(cv_model)) {
    message("Failed during cross-validation for: ", species_name)
    writeLines(paste("Failed during cross-validation:", Sys.time()),
               file.path(sdm_folder, "error_log.txt"))
    next
  }

  cv_auc_test <- SDMtune::auc(cv_model, test = TRUE)
  cv_tss_test <- SDMtune::tss(cv_model, test = TRUE)

  message("CV Testing - AUC: ", round(cv_auc_test, 3), " | TSS: ", round(cv_tss_test, 3))

  # ---------- 4.10.3 ROC Curve Visualization ----------
  tryCatch({
    roc_plot <- plotROC(final_model, test = test_final)
    ggsave(file.path(sdm_folder, paste0("roc_curve_", species_key, ".png")),
           roc_plot, width = 8, height = 6, dpi = 300)
  }, error = function(e) {
    message("Could not create ROC plot: ", e$message)
  })

  # ---------- 4.10.4 Jackknife Variable Importance ----------
  # Tests each variable independently (with_only) and by exclusion (without)
  # Identifies dominant and redundant variables
  message("Performing Jackknife test...")
  jk_results <- tryCatch({
    doJk(final_model, metric = "auc", test = test_final, with_only = TRUE)
  }, error = function(e) {
    message("Jackknife test failed: ", e$message)
    NULL
  })

  if (!is.null(jk_results)) {
    rio::export(jk_results, file.path(tables_folder, paste0("jackknife_results_", species_key, ".xlsx")))

    # Plot jackknife results
    jk_plot <- plotJk(jk_results, type = "test", ref = auc_val) 

    ggsave(file.path(sdm_folder, paste0("jackknife_plot_", species_key, ".png")),
           jk_plot, width = 10, height = max(6, nrow(jk_results) * 0.3), dpi = 300)
  }


  # ================================================================================
  # CHAPTER 4.11: MODEL PREDICTIONS
  # ================================================================================

  # ---------- 4.11.1 Current Climate Predictions (Clamped & Unclamped) ----------
  # CLAMPED: Restricts environmental values to training range (conservative, default)
  # UNCLAMPED: Allows extrapolation beyond training data (exploratory)
  # Rationale: Clamped predictions prevent unrealistic extrapolation, but unclamped
  #            useful for assessing model behavior under novel climates
  message("Generating predictions (clamped and unclamped)...")

  # CLAMPED predictions (default, conservative)
  # Clamps environmental values to training range (prevents unrealistic extrapolation)
  pred_current_clamped <- predict(cv_model, data = final_predictors, type = "cloglog",
                                  fun = "mean", progress = TRUE, clamp = TRUE)

  # UNCLAMPED predictions (allows extrapolation)
  # Model extrapolates beyond training data range (may produce unrealistic values)
  pred_current_unclamped <- predict(cv_model, data = final_predictors, type = "cloglog",
                                    fun = "mean", progress = TRUE, clamp = FALSE)

  # Additional uncertainty measures (using clamped version)
  pred_uncertainty <- predict(cv_model, data = final_predictors, type = "cloglog",
                              fun = c("max", "sd"), progress = TRUE, clamp = TRUE)

  pred_max <- pred_uncertainty$max
  pred_sd <- pred_uncertainty$sd  # Standard deviation across CV folds

  # Calculate clamping effect (difference map)
  clamp_diff_current <- pred_current_unclamped - pred_current_clamped


  # ================================================================================
  # CHAPTER 4.12: MESS ANALYSIS (EXTRAPOLATION DETECTION)
  # ================================================================================

  # ---------- 4.12.1 Multivariate Environmental Similarity Surfaces ----------
  # MESS identifies areas where predictions extrapolate beyond training data
  # Negative MESS = extrapolation (novel climate, less reliable predictions)
  # Positive MESS = interpolation (within training range, reliable predictions)
  # Reference: Elith et al. 2010 (The art of modelling range-shifting species)
  message("Computing MESS analysis...")

  tryCatch({
    # Extract training environment from presence-only data (from final model)
    train_env <- final_model@data@data[final_model@data@pa == 1, ]

    # Convert terra raster to raster::stack (required by dismo::mess)
    if (!require("raster", quietly = TRUE)) {
      message("Installing raster package for MESS analysis...")
      install.packages("raster")
      library(raster)
    }

    final_predictors_raster <- raster::stack(final_predictors)

    # Compute MESS using raster::stack
    mess_result <- dismo::mess(final_predictors_raster, train_env, full = FALSE)

    # Convert back to terra for consistency
    mess_raster <- terra::rast(mess_result)
    names(mess_raster) <- "MESS"

    # Save MESS raster
    terra::writeRaster(mess_raster,
                      file.path(raster_folder, paste0("mess_", species_key, ".tif")),
                      overwrite = TRUE)

    # Create MESS plot with appropriate color scale
    mess_plot <- ggplot() +
      tidyterra::geom_spatraster(data = mess_raster, maxcell = 5e+07) +
      scale_fill_gradient2(
        low = "#d73027",      # Red for extrapolation
        mid = "#ffffbf",      # Yellow for transition
        high = "#1a9850",     # Green for interpolation
        midpoint = 0,
        na.value = "transparent",
        name = "MESS",
        limits = c(-100, 100),
        oob = scales::squish
      ) +
      labs(title = paste("MESS Analysis (Extrapolation Risk) -", species_name),
           subtitle = "Red = novel climate (extrapolation, less reliable); Green = within training range",
           x = "Longitude", y = "Latitude") +
      map_theme()

    ggsave(file.path(sdm_folder, paste0("mess_analysis_", species_key, ".png")),
           mess_plot, width = 14, height = 10, dpi = 300, bg = "white")

    # Calculate extrapolation statistics
    mess_vals <- terra::values(mess_raster, na.rm = TRUE)
    pct_extrapolation <- sum(mess_vals < 0, na.rm = TRUE) / length(mess_vals) * 100
    mean_mess <- mean(mess_vals, na.rm = TRUE)

    message("MESS - Mean: ", round(mean_mess, 2), " | Extrapolation: ", round(pct_extrapolation, 1), "%")

    # Save MESS summary
    mess_summary <- data.frame(
      mean_mess = round(mean_mess, 3),
      pct_extrapolation = round(pct_extrapolation, 1),
      pct_interpolation = round(100 - pct_extrapolation, 1)
    )
    rio::export(mess_summary, file.path(tables_folder, paste0("mess_summary_", species_key, ".xlsx")))

  }, error = function(e) {
    message("MESS analysis failed: ", e$message)
    message("Continuing without MESS analysis...")
  })


  # ================================================================================
  # CHAPTER 4.13: FUTURE CLIMATE PROJECTIONS
  # ================================================================================

  # ---------- 4.13.1 Future Predictions (SSP585 2021-2040) ----------
  # Projects model to future climate scenario (high emissions pathway)
  # Both clamped and unclamped versions generated
  message("Future projections (clamped and unclamped)...")

  # Subset future predictors to match final_predictors
  # Note: Layer names are already aligned at the beginning of the script
  future_predictors_subset <- future_predictors[[names(final_predictors)]]

  # CLAMPED future predictions (conservative, default for publication)
  future_pred_clamped <- predict(cv_model, data = future_predictors_subset, type = "cloglog",
                                 fun = "mean", progress = TRUE, clamp = TRUE)

  # UNCLAMPED future predictions (allows extrapolation to novel climates)
  future_pred_unclamped <- predict(cv_model, data = future_predictors_subset, type = "cloglog",
                                   fun = "mean", progress = TRUE, clamp = FALSE)

  # Ensure same geometry for all rasters
  if (!terra::compareGeom(pred_current_clamped, future_pred_clamped, stopOnError = FALSE)) {
    future_pred_clamped <- terra::resample(future_pred_clamped, pred_current_clamped, method = "bilinear")
    future_pred_unclamped <- terra::resample(future_pred_unclamped, pred_current_clamped, method = "bilinear")
  }

  # Calculate clamping effect for future predictions
  clamp_diff_future <- future_pred_unclamped - future_pred_clamped

  # Calculate climate change differences (both clamped and unclamped versions)
  climate_change_clamped <- future_pred_clamped - pred_current_clamped
  climate_change_unclamped <- future_pred_unclamped - pred_current_unclamped

  # Calculate clamping statistics
  clamp_stats_current <- terra::values(clamp_diff_current, na.rm = TRUE)
  clamp_stats_future <- terra::values(clamp_diff_future, na.rm = TRUE)

  mean_clamp_current <- mean(abs(clamp_stats_current), na.rm = TRUE)
  mean_clamp_future <- mean(abs(clamp_stats_future), na.rm = TRUE)
  pct_clamped_current <- sum(abs(clamp_stats_current) > 0.01, na.rm = TRUE) / length(clamp_stats_current) * 100
  pct_clamped_future <- sum(abs(clamp_stats_future) > 0.01, na.rm = TRUE) / length(clamp_stats_future) * 100

  message("Clamping effect - Current: ", round(mean_clamp_current, 4), " (", round(pct_clamped_current, 1), "% affected)")
  message("Clamping effect - Future: ", round(mean_clamp_future, 4), " (", round(pct_clamped_future, 1), "% affected)")

  # Climate change statistics (using clamped version as default)
  diff_values <- terra::values(climate_change_clamped, na.rm = TRUE)
  mean_change <- mean(diff_values, na.rm = TRUE)
  pos_change_area <- sum(diff_values > 0, na.rm = TRUE) / length(diff_values) * 100
  neg_change_area <- sum(diff_values < 0, na.rm = TRUE) / length(diff_values) * 100

  message("Climate change (clamped) - Mean: ", round(mean_change, 4), " | Gain: ", round(pos_change_area, 1),
          "% | Loss: ", round(neg_change_area, 1), "%")


  # ================================================================================
  # CHAPTER 4.14: OUTPUT RASTERS
  # ================================================================================

  # ---------- 4.14.1 Save All Prediction Rasters ----------
  # Saves 9 raster outputs: current (2), future (2), change (2), clamp effect (2), uncertainty (1), MESS (1)
  message("Saving rasters (clamped and unclamped versions)...")

  # Current climate predictions
  terra::writeRaster(pred_current_clamped,
                     file.path(raster_folder, paste0("current_clamped_", species_key, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(pred_current_unclamped,
                     file.path(raster_folder, paste0("current_unclamped_", species_key, ".tif")),
                     overwrite = TRUE)

  # Future climate predictions
  terra::writeRaster(future_pred_clamped,
                     file.path(raster_folder, paste0("future_clamped_", species_key, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(future_pred_unclamped,
                     file.path(raster_folder, paste0("future_unclamped_", species_key, ".tif")),
                     overwrite = TRUE)

  # Climate change difference maps
  terra::writeRaster(climate_change_clamped,
                     file.path(raster_folder, paste0("change_clamped_", species_key, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(climate_change_unclamped,
                     file.path(raster_folder, paste0("change_unclamped_", species_key, ".tif")),
                     overwrite = TRUE)

  # Clamping effect difference maps
  terra::writeRaster(clamp_diff_current,
                     file.path(raster_folder, paste0("clamp_effect_current_", species_key, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(clamp_diff_future,
                     file.path(raster_folder, paste0("clamp_effect_future_", species_key, ".tif")),
                     overwrite = TRUE)

  # Uncertainty map (SD across CV folds)
  terra::writeRaster(pred_sd,
                     file.path(raster_folder, paste0("uncertainty_", species_key, ".tif")),
                     overwrite = TRUE)


  # ================================================================================
  # CHAPTER 4.15: VISUALIZATIONS
  # ================================================================================

  # ---------- 4.15.1 Define Map Theme ----------
  message("Creating visualizations...")

  # Define map theme
  map_theme <- function() {
    theme_minimal() + theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_line(color = "white", size = 0.3),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
      axis.title = element_text(lineheight = 12),
      legend.position = "bottom",
      legend.key.width = unit(2, "cm")
    )
  }

  # Current clamped - Global
  message("  - Current clamped (global)")
  current_clamped_plot <- ggplot() +
    tidyterra::geom_spatraster(data = pred_current_clamped, maxcell = 5e+07) +
    scale_fill_viridis_c(name = "Suitability", option = "plasma", limits = c(0, 1), na.value = "white") +
    labs(title = paste("Current Clamped -", species_name),
         subtitle = paste("Test AUC:", round(auc_val, 3), "| CV AUC:", round(cv_auc_test, 3), "| TSS:", round(tss_val, 3)),
         x = "Longitude", y = "Latitude") +
    map_theme()
  
  ggsave(file.path(sdm_folder, paste0("current_clamped_", species_key, ".png")),
         current_clamped_plot, width = 14, height = 10, dpi = 300, bg = "white")

  # Current clamped - Europe
  message("  - Current clamped (Europe)")
  current_europe_clamped_plot <- current_clamped_plot +
    coord_sf(xlim = c(-25, 45), ylim = c(34, 72), expand = FALSE) +
    labs(title = paste("Current Clamped (Europe) -", species_name))
  
  ggsave(file.path(sdm_folder, paste0("current_europe_clamped_", species_key, ".png")),
         current_europe_clamped_plot, width = 14, height = 10, dpi = 300, bg = "white")

  # Current unclamped - Global
  message("  - Current unclamped (global)")
  current_unclamped_plot <- ggplot() +
    tidyterra::geom_spatraster(data = pred_current_unclamped, maxcell = 5e+07) +
    scale_fill_viridis_c(name = "Suitability", option = "plasma", limits = c(0, 1), na.value = "white") +
    labs(title = paste("Current Unclamped -", species_name),
         subtitle = paste("Test AUC:", round(auc_val, 3), "| CV AUC:", round(cv_auc_test, 3), "| TSS:", round(tss_val, 3)),
         x = "Longitude", y = "Latitude") +
    map_theme()
  
  ggsave(file.path(sdm_folder, paste0("current_unclamped_", species_key, ".png")),
         current_unclamped_plot, width = 14, height = 10, dpi = 300, bg = "white")

  # Current unclamped - Europe
  message("  - Current unclamped (Europe)")
  current_europe_unclamped_plot <- current_unclamped_plot +
    coord_sf(xlim = c(-25, 45), ylim = c(34, 72), expand = FALSE) +
    labs(title = paste("Current Unclamped (Europe) -", species_name))
  
  ggsave(file.path(sdm_folder, paste0("current_europe_unclamped_", species_key, ".png")),
         current_europe_unclamped_plot, width = 14, height = 10, dpi = 300, bg = "white")

  # Future clamped - Global
  message("  - Future clamped (global)")
  future_clamped_plot <- ggplot() +
    tidyterra::geom_spatraster(data = future_pred_clamped, maxcell = 5e+07) +
    scale_fill_viridis_c(name = "Suitability", option = "plasma", limits = c(0, 1), na.value = "white") +
    labs(title = paste("Future Clamped -", species_name),
         subtitle = paste("Test AUC:", round(auc_val, 3), "| CV AUC:", round(cv_auc_test, 3), "| TSS:", round(tss_val, 3)),
         x = "Longitude", y = "Latitude") +
    map_theme()
  
  ggsave(file.path(sdm_folder, paste0("future_clamped_", species_key, ".png")),
         future_clamped_plot, width = 14, height = 10, dpi = 300, bg = "white")

  # Future clamped - Europe
  message("  - Future clamped (Europe)")
  future_europe_clamped_plot <- future_clamped_plot +
    coord_sf(xlim = c(-25, 45), ylim = c(34, 72), expand = FALSE) +
    labs(title = paste("Future Clamped (Europe) -", species_name))
  
  ggsave(file.path(sdm_folder, paste0("future_europe_clamped_", species_key, ".png")),
         future_europe_clamped_plot, width = 14, height = 10, dpi = 300, bg = "white")

  # Future unclamped - Global
  message("  - Future unclamped (global)")
  future_unclamped_plot <- ggplot() +
    tidyterra::geom_spatraster(data = future_pred_unclamped, maxcell = 5e+07) +
    scale_fill_viridis_c(name = "Suitability", option = "plasma", limits = c(0, 1), na.value = "white") +
    labs(title = paste("Future Unclamped -", species_name),
         subtitle = paste("Test AUC:", round(auc_val, 3), "| CV AUC:", round(cv_auc_test, 3), "| TSS:", round(tss_val, 3)),
         x = "Longitude", y = "Latitude") +
    map_theme()
  
  ggsave(file.path(sdm_folder, paste0("future_unclamped_", species_key, ".png")),
         future_unclamped_plot, width = 14, height = 10, dpi = 300, bg = "white")

  # Future unclamped - Europe
  message("  - Future unclamped (Europe)")
  future_europe_unclamped_plot <- future_unclamped_plot +
    coord_sf(xlim = c(-25, 45), ylim = c(34, 72), expand = FALSE) +
    labs(title = paste("Future Unclamped (Europe) -", species_name))
  
  ggsave(file.path(sdm_folder, paste0("future_europe_unclamped_", species_key, ".png")),
         future_europe_unclamped_plot, width = 14, height = 10, dpi = 300, bg = "white")

  # Uncertainty map
  message("  - Uncertainty map")
  uncertainty_plot <- ggplot() +
    tidyterra::geom_spatraster(data = pred_sd, maxcell = 5e+07) +
    scale_fill_viridis_c(name = "SD", option = "magma", na.value = "white") +
    labs(title = paste("Prediction Uncertainty -", species_name),
         subtitle = "Standard deviation across CV folds",
         x = "Longitude", y = "Latitude") +
    map_theme()
  
  ggsave(file.path(sdm_folder, paste0("uncertainty_", species_key, ".png")),
         uncertainty_plot, width = 14, height = 10, dpi = 300, bg = "white")

  # Climate change impact
  message("  - Climate change impact")
  climate_change_plot <- ggplot() +
    tidyterra::geom_spatraster(data = climate_change_clamped, maxcell = 5e+07) +
    scale_fill_gradient2(name = "Change\n(Future - Current)",
                         low = "#d73027", mid = "#ffffbf", high = "#1a9850",
                         midpoint = 0, limits = c(-1, 1), na.value = "transparent",
                         oob = scales::squish) +
    labs(title = paste("Climate Change Impact (CLAMPED) -", species_name),
         subtitle = paste("Mean change:", round(mean_change, 3), "| Gain:", round(pos_change_area, 1),
                          "% | Loss:", round(neg_change_area, 1), "%"),
         x = "Longitude", y = "Latitude") +
    map_theme()
  
  ggsave(file.path(sdm_folder, paste0("climate_change_impact_clamped_", species_key, ".png")),
         climate_change_plot, width = 14, height = 10, dpi = 300, bg = "white")

  # Clamp effect current
  message("  - Clamping effect (current)")
  clamp_effect_current_plot <- ggplot() +
    tidyterra::geom_spatraster(data = clamp_diff_current, maxcell = 5e+07) +
    scale_fill_gradient2(name = "Difference\n(Unclamped - Clamped)",
                         low = "blue", mid = "white", high = "red",
                         midpoint = 0, na.value = "transparent") +
    labs(title = paste("Clamping Effect - Current -", species_name),
         subtitle = paste("Mean effect:", round(mean_clamp_current, 4), "|", round(pct_clamped_current, 1), "% affected"),
         x = "Longitude", y = "Latitude") +
    map_theme()
  
  ggsave(file.path(sdm_folder, paste0("clamp_effect_current_", species_key, ".png")),
         clamp_effect_current_plot, width = 14, height = 10, dpi = 300, bg = "white")

  # Clamp effect future
  message("  - Clamping effect (future)")
  clamp_effect_future_plot <- ggplot() +
    tidyterra::geom_spatraster(data = clamp_diff_future, maxcell = 5e+07) +
    scale_fill_gradient2(name = "Difference\n(Unclamped - Clamped)",
                         low = "blue", mid = "white", high = "red",
                         midpoint = 0, na.value = "transparent") +
    labs(title = paste("Clamping Effect - Future -", species_name),
         subtitle = paste("Mean effect:", round(mean_clamp_future, 4), "|", round(pct_clamped_future, 1), "% affected"),
         x = "Longitude", y = "Latitude") +
    map_theme()
  
  ggsave(file.path(sdm_folder, paste0("clamp_effect_future_", species_key, ".png")),
         clamp_effect_future_plot, width = 14, height = 10, dpi = 300, bg = "white")

  # MESS analysis
  message("  - MESS analysis")
  mess_plot <- ggplot() +
    tidyterra::geom_spatraster(data = mess_raster, maxcell = 5e+07) +
    scale_fill_gradient2(
      low = "#d73027",      # Red for extrapolation
      mid = "#ffffbf",      # Yellow for transition
      high = "#1a9850",     # Green for interpolation
      midpoint = 0,
      na.value = "transparent",
      name = "MESS",
      limits = c(-100, 100),
      oob = scales::squish
    ) +
    labs(title = paste("MESS Analysis (Extrapolation Risk) -", species_name),
         subtitle = "Red = novel climate (extrapolation, less reliable); Green = within training range",
         x = "Longitude", y = "Latitude") +
    map_theme()
  
  ggsave(file.path(sdm_folder, paste0("mess_analysis_", species_key, ".png")),
         mess_plot, width = 14, height = 10, dpi = 300, bg = "white")

  message("Visualization plots complete")

  # ---------- 4.15.2 Clamping Comparison Plots (3-Panel) ----------
  # Side-by-side comparison: Clamped | Unclamped | Difference
message("Creating clamping comparison plots...")

current_clamped_plt <- ggplot() +
  tidyterra::geom_spatraster(data = pred_current_clamped, maxcell=5e+07) +
  scale_fill_viridis_c(name = "Suitability", option = "plasma", limits = c(0, 1), na.value = "white") +
  labs(title = "Current (Clamped)",
       subtitle = paste(species_name)) +
  map_theme() + theme(legend.position = "right")

current_unclamped_plt <- ggplot() +
  tidyterra::geom_spatraster(data = pred_current_unclamped, maxcell=5e+07) +
  scale_fill_viridis_c(name = "Suitability", option = "plasma", limits = c(0, 1), na.value = "white") +
  labs(title = "Current (Unclamped)") +
  map_theme() + theme(legend.position = "right")

current_diff_plt <- ggplot() +
  tidyterra::geom_spatraster(data = clamp_diff_current, maxcell=5e+07) +
  scale_fill_gradient2(name = "Difference", low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                       midpoint = 0, limits = c(-0.5, 0.5), oob = scales::squish) +
  labs(title = "Difference") +
  map_theme() + theme(legend.position = "right")

#current_comparison <- current_clamped_plt + current_unclamped_plt + current_diff_plt +
#  plot_layout(ncol = 3) +
#  plot_annotation(title = paste("Clamping Comparison - Current -", species_name),
#                  subtitle = "Left: Clamped (conservative) | Middle: Unclamped (extrapolates) | Right: Difference",
#                  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
#
#ggsave(file.path(sdm_folder, paste0("current_clamping_comparison_", species_key, ".png")),
#       current_comparison, width = 20, height = 7, dpi = 300, bg = "white")

# Side-by-side comparison: Future CLAMPED vs UNCLAMPED
future_clamped_plt <- ggplot() +
  tidyterra::geom_spatraster(data = future_pred_clamped, maxcell=5e+07) +
  scale_fill_viridis_c(name = "Suitability", option = "plasma", limits = c(0, 1), na.value = "white") +
  labs(title = "Future (Clamped)",
       subtitle = paste(species_name, "- SSP585 2021-2040")) +
  map_theme() + theme(legend.position = "right")

future_unclamped_plt <- ggplot() +
  tidyterra::geom_spatraster(data = future_pred_unclamped, maxcell=5e+07) +
  scale_fill_viridis_c(name = "Suitability", option = "plasma", limits = c(0, 1), na.value = "white") +
  labs(title = "Future (Unclamped)") +
  map_theme() + theme(legend.position = "right")

future_diff_plt <- ggplot() +
  tidyterra::geom_spatraster(data = clamp_diff_future, maxcell=5e+07) +
  scale_fill_gradient2(name = "Difference", low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                       midpoint = 0, limits = c(-0.5, 0.5), oob = scales::squish) +
  labs(title = "Difference") +
  map_theme() + theme(legend.position = "right")

#future_comparison <- future_clamped_plt + future_unclamped_plt + future_diff_plt +
#  plot_layout(ncol = 3) +
#  plot_annotation(title = paste("Clamping Comparison - Future -", species_name),
#                  subtitle = "Left: Clamped (conservative) | Middle: Unclamped (extrapolates) | Right: Difference",
#                  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
#
#ggsave(file.path(sdm_folder, paste0("future_clamping_comparison_", species_key, ".png")),
#       future_comparison, width = 20, height = 7, dpi = 300, bg = "white")

  # ---------- 4.15.3 Binary Threshold Maps ----------
  # Converts continuous suitability to binary presence/absence using max TSS threshold
message("Creating threshold maps...")

# Get thresholds
ths <- thresholds(final_model, type = "cloglog")

# Display threshold information
message("Thresholds:")
message("  Minimum training presence: ", round(ths[1, 2], 4))
message("  Equal sens/spec: ", round(ths[2, 2], 4))
message("  Max TSS: ", round(ths[3, 2], 4))

# Create binary maps using plotPA (using CLAMPED predictions)
# Max TSS threshold (recommended)
current_binary_maxTSS <-  plotPA(pred_current_clamped, th = ths[3, 2])
ggsave(filename = file.path(sdm_folder, "th_current_binary_maxTSS_clamped.png"), plot = current_binary_maxTSS, width = 40, height = 30, units = "cm")

future_binary_maxTSS <-   plotPA(future_pred_clamped, th = ths[3, 2])
ggsave(filename = file.path(sdm_folder, "th_future_binary_maxTSS_clamped.png"), plot = future_binary_maxTSS, width = 40, height = 30, units = "cm")

  # ---------- 4.15.4 Response Curves ----------
  # Shows predicted suitability as function of each environmental variable
  # Standard response: using all background data
  # Marginal response: using presence data only
message("Creating response curves...")

# Create directory for response curves
response_dir <- file.path(sdm_folder, "response_curves")
if (!dir.exists(response_dir)) dir.create(response_dir, recursive = TRUE)

# Generate response curves for each variable
for (var_name in names(final_predictors)) {
    tryCatch({
      # Standard response curve
      p1 <- plotResponse(
        model = final_model,
        var = var_name,
        type = "cloglog",
        only_presence = FALSE,
        marginal = FALSE,
        rug = TRUE,
        color = "blue"
      )


ggsave(filename = file.path(response_dir, paste0("response_", var_name, ".png")),
        plot = p1,
        width = 8,
        height = 6,
        dpi = 300
      )

# Marginal response curve
p2 <- plotResponse(
        model = final_model,
        var = var_name,
        type = "cloglog",
        only_presence = TRUE,
        marginal = TRUE,
        fun = mean,
        rug = TRUE,
        color = "red"
      )

ggsave(
        filename = file.path(response_dir, paste0("marginal_response_", var_name, ".png")),
        plot = p2,
        width = 8,
        height = 6,
        dpi = 300
      )

    }, error = function(e) {
      warning("Failed to create response curve for ", var_name, ": ", e$message)
    })
  }

  message("Response curves saved to: ", response_dir)


  # ================================================================================
  # CHAPTER 4.16: SUMMARY AND REPORTING
  # ================================================================================

  # ---------- 4.16.1 Comprehensive Model Summary ----------
  # Compiles all metrics, parameters, and results into single Excel file

# ---------- Model Report (Disabled) ----------
# Note: modelReport() disabled due to long processing time
# tryCatch({
#    modelReport(final_model, type = "cloglog",
#               folder = file.path(sdm_folder, paste0("report_", species_key)),
#              test = test_final, response_curves = TRUE, only_presence = TRUE,
#              jk = TRUE, env = final_predictors)
#
# }, error = function(e) warning("Model report failed: ", e$message))

  # Calculate data source breakdown
  source_summary <- occ_df %>%
    group_by(source) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(desc(count))

  source_breakdown <- paste(
    paste0(source_summary$source, ": ", source_summary$count),
    collapse = "; "
  )

  # Calculate processing time
  species_end_time <- Sys.time()
  species_elapsed <- as.numeric(difftime(species_end_time, species_start_time, units = "mins"))

  model_summary <- data.frame(
    species = species_name,
    code = species_key,
    species_key = species_key,
    date = as.character(Sys.Date()),
    start_time = format(species_start_time, "%Y-%m-%d %H:%M:%S"),
    end_time = format(species_end_time, "%Y-%m-%d %H:%M:%S"),
    processing_time_mins = round(species_elapsed, 2),
    raw_records = nrow(occ_df),
    raw_online = nrow(occ_df) - gabi_records,
    raw_gabi_antmaps = gabi_records,
    cleaned_records = nrow(occ_df_clean),
    final_records = nrow(p_coords),
    background_records = nrow(bg_coords),
    data_sources = source_breakdown,
    num_variables = length(names(final_predictors)),
    variables = paste(names(final_predictors), collapse = "; "),
    test_auc = round(auc_val, 3),
    test_tss = round(tss_val, 3),
    cv_auc = round(cv_auc_test, 3),
    cv_tss = round(cv_tss_test, 3),
    threshold = round(optimal_threshold, 3),
    sensitivity = round(sensitivity, 3),
    specificity = round(specificity, 3),
    true_positives = cm_stats$tp,
    false_positives = cm_stats$fp,
    true_negatives = cm_stats$tn,
    false_negatives = cm_stats$fn,
    fc = final_model@model@fc,
    reg = final_model@model@reg,
    cv_folds = CV_FOLDS,
    climate_mean_change = round(mean_change, 4),
    climate_gain_pct = round(pos_change_area, 1),
    climate_loss_pct = round(neg_change_area, 1)
  )

  rio::export(model_summary, file.path(tables_folder, paste0("summary_", species_key, ".xlsx")))

  # ---------- 4.16.2 Species Completion Message ----------
  message("\n", paste(rep("=", 60), collapse = ""))
  message("COMPLETED: ", species_name)
  message(paste(rep("=", 60), collapse = ""))
  message("Results: ", sdm_folder)
  message("Test - AUC: ", round(auc_val, 3), " | TSS: ", round(tss_val, 3))
  message("CV   - AUC: ", round(cv_auc_test, 3), " | TSS: ", round(cv_tss_test, 3))
  message("Finished: ", format(species_end_time, "%H:%M:%S"))
  message("Processing time: ", round(species_elapsed, 2), " minutes")

  # Clean up memory
  gc()
  flush.console()
  cat("\014")

}  # End of species loop


# ================================================================================
# CHAPTER 5: PIPELINE COMPLETION
# ================================================================================

# ---------- 5.1 Overall Processing Summary ----------
overall_end_time <- Sys.time()
overall_elapsed <- as.numeric(difftime(overall_end_time, overall_start_time, units = "mins"))
overall_elapsed_hours <- overall_elapsed / 60

message("\n", paste(rep("=", 80), collapse = ""))
message("ALL SPECIES PROCESSING COMPLETED!")
message(paste(rep("=", 80), collapse = ""))
message("Total species processed: ", nrow(test_species))
message("Results directory: ", base_dir)
message("")
message("TIMING SUMMARY:")
message("  Start time: ", format(overall_start_time, "%Y-%m-%d %H:%M:%S"))
message("  End time:   ", format(overall_end_time, "%Y-%m-%d %H:%M:%S"))
message("  Total time: ", round(overall_elapsed, 2), " minutes (", round(overall_elapsed_hours, 2), " hours)")
message("  Average per species: ", round(overall_elapsed / nrow(test_species), 2), " minutes")
message(paste(rep("=", 80), collapse = ""))

gc()

