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
# Script: 10_updated_sdmtune_loop_merged.R
# Purpose: Production SDM pipeline with full features and comprehensive outputs
# Author: SDM Pipeline Project
# Last Updated: 2026-03-31
#
# Description:
# Complete species distribution modeling pipeline using MaxEnt (via SDMtune).
# Predicts current and future distributions under climate change (SSP585 2021-2040).
# Includes robust error handling, automatic fallbacks, and comprehensive outputs.
#
# Key Features:
# - Circle-based background selection with automatic fallback (EPSG:4087)
# - Two-stage variable selection (importance + correlation filtering)
# - Conservative hyperparameter optimization (H, Q, HQ features; reg 2-10)
# - Cross-validation with uncertainty quantification
# - Clamped and unclamped predictions
# - MESS analysis for extrapolation detection
# - modEvA extended model evaluation (Boyce, calibration, HLfit)
# - Norwegian climate histograms (species niche vs Norwegian background)
# - Three-format model summary (JSON, Markdown, plain text)
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
# >> MERGED FROM Script 2: additional libraries for extended evaluation and output
library(modEvA)             # Extended model evaluation (AUC, calibration, thresholds)
library(jsonlite)           # JSON output for model_summary.json
library(scales)             # Scale functions (squish for ggplot2 colour limits)

# ---------- 1.2 Configuration Parameters ----------
MIN_RECORDS <- 50           # Minimum presence records required after cleaning
DOWNLOAD_LIMIT <- 100000    # Maximum records per data source (spocc)
CV_FOLDS <- 5               # Cross-validation folds (increased from 2 for robustness)
RECENT_YEARS <- 1960        # Temporal filter: only records from 1960 onwards

# Background sampling buffer distance (meters, projected EPSG:4087)
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
  pull(species)

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

    available_vars <- names(current_predictors)
    missing_vars <- setdiff(selected_vars, available_vars)

    if (length(missing_vars) > 0) {
      warning("  ", length(missing_vars), " selected variables not found, using available ones")
      selected_vars <- intersect(selected_vars, available_vars)
    }

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

# >> MERGED FROM Script 2: high-resolution land polygon for background land-filtering
# world_full: ne_countries(scale=50) unioned — full southern coverage (ymin ~-90)
# Used ONLY for st_filter (land/ocean discrimination) inside background selection.
# Do NOT use world_full for plotting — Antarctica extent causes ggplot banding artifacts.
world_full <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf") |>
  sf::st_union() |>
  sf::st_make_valid()

# >> MERGED FROM Script 2: bioclim variable label lookup
# Maps raw bioclim layer names to human-readable axis labels.
# Defined here (before the species loop) so it is available in the Norwegian
# climate histogram chapter and the named-predictor-stack export.
bioclim_labels <- list(
  bio1  = "Annual Mean Temperature (deg C x 10)",
  bio2  = "Mean Diurnal Range (deg C x 10)",
  bio3  = "Isothermality (%)",
  bio4  = "Temperature Seasonality (SD x 100)",
  bio5  = "Max Temperature of Warmest Month (deg C x 10)",
  bio6  = "Min Temperature of Coldest Month (deg C x 10)",
  bio7  = "Temperature Annual Range (deg C x 10)",
  bio8  = "Mean Temperature of Wettest Quarter (deg C x 10)",
  bio9  = "Mean Temperature of Driest Quarter (deg C x 10)",
  bio10 = "Mean Temperature of Warmest Quarter (deg C x 10)",
  bio11 = "Mean Temperature of Coldest Quarter (deg C x 10)",
  bio12 = "Annual Precipitation (mm)",
  bio13 = "Precipitation of Wettest Month (mm)",
  bio14 = "Precipitation of Driest Month (mm)",
  bio15 = "Precipitation Seasonality (CV)",
  bio16 = "Precipitation of Wettest Quarter (mm)",
  bio17 = "Precipitation of Driest Quarter (mm)",
  bio18 = "Precipitation of Warmest Quarter (mm)",
  bio19 = "Precipitation of Coldest Quarter (mm)"
)


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
  # Strip parenthetical author citations e.g. "(Wheeler, 1903) Wheeler, 1903",
  # sanitize remaining special characters, and cap at 40 chars to prevent
  # Windows MAX_PATH (260 char) truncation that silently drops ".png" endings.
  safe_species_name <- gsub("\\([^)]*\\)", "", species_name)   # remove (author, year)
  safe_species_name <- trimws(safe_species_name)
  safe_species_name <- gsub("[^a-zA-Z0-9]", "_", safe_species_name)  # replace special chars
  safe_species_name <- gsub("_+", "_", safe_species_name)       # collapse multiple underscores
  safe_species_name <- sub("_$", "", safe_species_name)          # trim trailing underscore
  safe_species_name <- substr(safe_species_name, 1, 40)          # cap length
  species_folder_name <- paste0(species_key, "_", safe_species_name)
  species_folder <- file.path(base_dir, species_folder_name)
  sdm_folder <- file.path(species_folder, "updated_SDM_maxnet")

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
  # >> MERGED FROM Script 2: create summary_folder for JSON/MD/TXT outputs
  summary_folder <- file.path(sdm_folder, "model_summary")
  if (!dir.exists(summary_folder)) dir.create(summary_folder, recursive = TRUE)


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

    species_name_clean <- gsub("\\s*\\([^)]+\\)\\s*$", "", species_name)
    species_name_clean <- trimws(species_name_clean)

    message("Matching species name: '", species_name_clean, "'")

    gabi_species <- gabi_data %>%
      filter(tolower(trimws(valid_species_name)) == tolower(species_name_clean))

    if (nrow(gabi_species) > 0) {
      message("Found ", nrow(gabi_species), " records in GABI database")

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

      gabi_occ <- gabi_occ %>%
        filter(!is.na(decimalLongitude) & !is.na(decimalLatitude))

      message("GABI records with valid coordinates: ", nrow(gabi_occ))

      char_cols <- c("species", "decimalLongitude", "decimalLatitude")
      occ_df <- occ_df %>% mutate(across(all_of(char_cols), as.character))
      gabi_occ <- gabi_occ %>% mutate(across(all_of(char_cols), as.character))

      occ_df <- bind_rows(occ_df, gabi_occ)
      gabi_records <- nrow(gabi_occ)

      message("Total records after adding GABI data: ", nrow(occ_df))
    } else {
      message("No matching records found in GABI database")
    }
  }

  message("Raw occurrence records (total): ", nrow(occ_df),
          " (online: ", nrow(occ_df) - gabi_records, ", GABI: ", gabi_records, ")")

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
      tryCatch({
        temporal_bias_lat <- lm(decimalLatitude ~ year, data = recent_with_year)
        temporal_bias_lon <- lm(decimalLongitude ~ year, data = recent_with_year)

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
      cc_sea(verbose = TRUE) %>%
      cc_zero(verbose = TRUE) %>%
      cc_dupl(verbose = TRUE) %>%
      cc_outl(method = "quantile", mltpl = 3, verbose = TRUE, value = "clean")
  }, error = function(e) {
    message("Error during coordinate cleaning: ", e$message)
    message("Attempting simplified cleaning without cc_sea...")
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
    coords_clean <- as.data.frame(occ_df_clean)

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
  # Fallback: dismo::randomPoints if circle method fails
  #
  # CONFLICT: Script 1 used EPSG:3857 (Web Mercator) + sf::st_simplify + world for land mask.
  # Script 2 uses EPSG:4087 (World Equidistant Cylindrical, distance-accurate globally) +
  # st_wrap_dateline + world_full for land mask. Applying Script 2's implementation.
  # Rationale: EPSG:3857 distorts heavily at high latitudes; antimeridian-crossing rings
  # corrupt polygons for species with circumpolar distributions.
  # Reference: r-spatial/sf Issue #1710
  message("Generating background points...")

  p_coords <- occ_df_clean[, c("decimalLongitude", "decimalLatitude")]

  bg_coords <- tryCatch({
    message("Attempting circle-based background selection (EPSG:4087)...")

    # Disable S2 for all polygon operations (avoids spherical validity errors)
    sf::sf_use_s2(FALSE)

    # Project to EPSG:4087 (World Equidistant Cylindrical) — distance-accurate globally
    presence_sf        <- sf::st_as_sf(p_coords, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
    presence_projected <- sf::st_transform(presence_sf, crs = 4087)
    circles            <- sf::st_buffer(presence_projected, dist = BUFFER_DISTANCE)
    merged_polygon     <- sf::st_union(circles)

    # Transform back to 4326 and repair antimeridian-crossing rings
    merged_polygon_geo <- sf::st_transform(merged_polygon, crs = 4326)
    merged_polygon_geo <- sf::st_wrap_dateline(merged_polygon_geo,
                                               options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
    merged_polygon_geo <- sf::st_make_valid(merged_polygon_geo)

    # Sample and filter to land (still S2 off)
    bg_points      <- sf::st_sample(merged_polygon_geo, size = 25000, type = "random")
    bg_coords_land <- sf::st_filter(sf::st_sf(geometry = bg_points), world_full)

    sf::sf_use_s2(TRUE)   # restore only after all operations done

    bg_coords_temp <- sf::st_coordinates(bg_coords_land)
    bg_coords_temp <- as.data.frame(bg_coords_temp[, c("X", "Y")])
    colnames(bg_coords_temp) <- c("x", "y")

    if (nrow(bg_coords_temp) < 100) stop("circle method returned < 100 land points")

    message("Circle-based background selection successful: ", nrow(bg_coords_temp), " points")
    bg_coords_temp

  }, error = function(e) {
    sf::sf_use_s2(TRUE)   # always restore S2 on error
    message("Circle-based background selection FAILED (", e$message, ") — falling back to dismo::randomPoints")

    bg_out_g <- as.data.frame(dismo::randomPoints(raster::raster(current_predictors[[1]]), n = 25000))
    colnames(bg_out_g) <- c("x", "y")
    message("dismo::randomPoints fallback successful: ", nrow(bg_out_g), " points")
    bg_out_g
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
  distribution_plot <- ggplot(world) +
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

  data_swd <- addSamplesToBg(data_swd, all = FALSE)
  message("Background locations after adding unique presence: ", sum(data_swd@pa == 0))

  c(train_data, test_data) %<-% trainValTest(data_swd, test = 0.2, only_presence = TRUE, seed = 25)


  # ================================================================================
  # CHAPTER 4.8: VARIABLE SELECTION (TWO-STAGE)
  # ================================================================================

  # ---------- 4.8.1 Stage 1: Permutation Importance (Top 15 Variables) ----------
  message("Variable selection...")

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

  vi <- SDMtune::varImp(initial_model, permut = 5)

  rio::export(vi, file.path(tables_folder, paste0("initial_variable_importance_", species_key, ".xlsx")))

  vi_plot <- plotVarImp(vi[, c("Variable", "Permutation_importance")], color = "steelblue")

  ggsave(file.path(sdm_folder, paste0("initial_variable_importance_", species_key, ".png")),
         vi_plot, width = 10, height = max(6, nrow(vi) * 0.3), dpi = 300)

  top_vars <- vi %>%
    arrange(desc(Permutation_importance)) %>%
    slice_head(n = min(15, nrow(vi))) %>%
    pull(Variable)

  current_predictors_subset <- current_predictors[[top_vars]]

  # ---------- 4.8.2 Stage 2: Correlation-Based Removal ----------
  background <- prepareSWD(species = species_name, a = bg_coords, env = current_predictors_subset)
  background <- addSamplesToBg(background, all = TRUE)

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
    tryCatch({
      varSel(
        model = train(method = "Maxnet", data = train_subset),
        metric = "auc",
        test = test_subset,
        bg4cor = background,
        method = "spearman",
        cor_th = 0.85,
        env = current_predictors_subset,
        use_pc = FALSE,
        progress = TRUE,
        permut = 5
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

  selected_vars_df <- data.frame(
    variable = selected_vars,
    rank = seq_along(selected_vars)
  )
  rio::export(selected_vars_df, file.path(tables_folder, paste0("selected_variables_", species_key, ".xlsx")))

  final_vi <- SDMtune::varImp(selected_variables_model, permut = 10)
  rio::export(final_vi, file.path(tables_folder, paste0("final_variable_importance_", species_key, ".xlsx")))

  final_vi_plot <- plotVarImp(final_vi[, c("Variable", "Permutation_importance")], color = "darkgreen")

  ggsave(file.path(sdm_folder, paste0("final_variable_importance_", species_key, ".png")),
         final_vi_plot, width = 10, height = max(6, nrow(final_vi) * 0.3), dpi = 300)

  final_swd <- prepareSWD(species = species_name, p = p_coords, a = bg_coords, env = final_predictors)
  c(train_final, test_final) %<-% trainValTest(final_swd, test = 0.2, only_presence = TRUE, seed = 25)


  # ================================================================================
  # CHAPTER 4.9: HYPERPARAMETER OPTIMIZATION
  # ================================================================================

  message("Optimizing hyperparameters...")

  h <- list(reg = seq(2, 10, 2),
            fc = c("h", "q", "hq"))

  optimized_model <- tryCatch({
    message("PRIMARY OPTIMIZATION: Testing FC=h,q,hq with reg=2,4,6,8,10...")
    result <- optimizeModel(
      model = train(method = "Maxnet", data = train_final),
      hypers = h,
      metric = "auc",
      test = test_final,
      pop = 10,
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

    h_retry <- list(
      reg = seq(0.5, 5, 0.5),
      fc = c(
        "l", "q", "h",
        "lq", "lh", "qh",
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
  message("Calculating evaluation metrics...")

  auc_val <- SDMtune::auc(final_model, test = test_final)
  tss_val <- SDMtune::tss(final_model, test = test_final)

  ths <- thresholds(final_model, type = "cloglog")
  optimal_threshold <- ths[3, 2]

  message("AUC: ", round(auc_val, 3), " | TSS: ", round(tss_val, 3),
          " | Threshold: ", round(optimal_threshold, 3))

  cm <- confMatrix(final_model, test = test_final, th = optimal_threshold, type = "cloglog")
  cm_stats <- cm[1, ]

  sensitivity <- cm_stats$tp / (cm_stats$tp + cm_stats$fn)
  specificity <- cm_stats$tn / (cm_stats$tn + cm_stats$fp)

  message("Confusion Matrix - TP: ", cm_stats$tp, " | FP: ", cm_stats$fp,
          " | TN: ", cm_stats$tn, " | FN: ", cm_stats$fn)
  message("Sensitivity: ", round(sensitivity, 3), " | Specificity: ", round(specificity, 3))

  # ---------- 4.10.2 Cross-Validation ----------
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
  message("Performing Jackknife test...")
  jk_results <- tryCatch({
    doJk(final_model, metric = "auc", test = test_final, with_only = TRUE)
  }, error = function(e) {
    message("Jackknife test failed: ", e$message)
    NULL
  })

  if (!is.null(jk_results)) {
    rio::export(jk_results, file.path(tables_folder, paste0("jackknife_results_", species_key, ".xlsx")))

    jk_plot <- plotJk(jk_results, type = "test", ref = auc_val)

    ggsave(file.path(sdm_folder, paste0("jackknife_plot_", species_key, ".png")),
           jk_plot, width = 10, height = max(6, nrow(jk_results) * 0.3), dpi = 300)
  }


  # ================================================================================
  # CHAPTER 4.11: MODEL PREDICTIONS
  # ================================================================================

  message("Generating predictions (clamped and unclamped)...")

  pred_current_clamped <- predict(cv_model, data = final_predictors, type = "cloglog",
                                  fun = "mean", progress = TRUE, clamp = TRUE)

  pred_current_unclamped <- predict(cv_model, data = final_predictors, type = "cloglog",
                                    fun = "mean", progress = TRUE, clamp = FALSE)

  pred_uncertainty <- predict(cv_model, data = final_predictors, type = "cloglog",
                              fun = c("max", "sd"), progress = TRUE, clamp = TRUE)

  pred_max <- pred_uncertainty$max
  pred_sd <- pred_uncertainty$sd

  clamp_diff_current <- pred_current_unclamped - pred_current_clamped


  # ================================================================================
  # CHAPTER 4.12: MESS ANALYSIS (EXTRAPOLATION DETECTION)
  # ================================================================================

  message("Computing MESS analysis...")

  # Initialize MESS stats to NA in case MESS fails (used later in summaries)
  pct_extrapolation <- NA
  mean_mess <- NA
  mess_raster <- NULL

  tryCatch({
    train_env <- final_model@data@data[final_model@data@pa == 1, ]

    if (!require("raster", quietly = TRUE)) {
      message("Installing raster package for MESS analysis...")
      install.packages("raster")
      library(raster)
    }

    final_predictors_raster <- raster::stack(final_predictors)

    mess_result <- dismo::mess(final_predictors_raster, train_env, full = FALSE)

    mess_raster <- terra::rast(mess_result)
    names(mess_raster) <- "MESS"

    terra::writeRaster(mess_raster,
                       file.path(raster_folder, paste0("mess_", species_key, ".tif")),
                       overwrite = TRUE)

    mess_vals <- terra::values(mess_raster, na.rm = TRUE)
    pct_extrapolation <- sum(mess_vals < 0, na.rm = TRUE) / length(mess_vals) * 100
    mean_mess <- mean(mess_vals, na.rm = TRUE)

    message("MESS - Mean: ", round(mean_mess, 2), " | Extrapolation: ", round(pct_extrapolation, 1), "%")

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

  message("Future projections (clamped and unclamped)...")

  future_predictors_subset <- future_predictors[[names(final_predictors)]]

  future_pred_clamped <- predict(cv_model, data = future_predictors_subset, type = "cloglog",
                                 fun = "mean", progress = TRUE, clamp = TRUE)

  future_pred_unclamped <- predict(cv_model, data = future_predictors_subset, type = "cloglog",
                                   fun = "mean", progress = TRUE, clamp = FALSE)

  if (!terra::compareGeom(pred_current_clamped, future_pred_clamped, stopOnError = FALSE)) {
    future_pred_clamped <- terra::resample(future_pred_clamped, pred_current_clamped, method = "bilinear")
    future_pred_unclamped <- terra::resample(future_pred_unclamped, pred_current_clamped, method = "bilinear")
  }

  clamp_diff_future <- future_pred_unclamped - future_pred_clamped

  climate_change_clamped <- future_pred_clamped - pred_current_clamped
  climate_change_unclamped <- future_pred_unclamped - pred_current_unclamped

  clamp_stats_current <- terra::values(clamp_diff_current, na.rm = TRUE)
  clamp_stats_future <- terra::values(clamp_diff_future, na.rm = TRUE)

  mean_clamp_current <- mean(abs(clamp_stats_current), na.rm = TRUE)
  mean_clamp_future <- mean(abs(clamp_stats_future), na.rm = TRUE)
  pct_clamped_current <- sum(abs(clamp_stats_current) > 0.01, na.rm = TRUE) / length(clamp_stats_current) * 100
  pct_clamped_future <- sum(abs(clamp_stats_future) > 0.01, na.rm = TRUE) / length(clamp_stats_future) * 100

  message("Clamping effect - Current: ", round(mean_clamp_current, 4), " (", round(pct_clamped_current, 1), "% affected)")
  message("Clamping effect - Future: ", round(mean_clamp_future, 4), " (", round(pct_clamped_future, 1), "% affected)")

  diff_values <- terra::values(climate_change_clamped, na.rm = TRUE)
  mean_change <- mean(diff_values, na.rm = TRUE)
  pos_change_area <- sum(diff_values > 0, na.rm = TRUE) / length(diff_values) * 100
  neg_change_area <- sum(diff_values < 0, na.rm = TRUE) / length(diff_values) * 100

  message("Climate change (clamped) - Mean: ", round(mean_change, 4), " | Gain: ", round(pos_change_area, 1),
          "% | Loss: ", round(neg_change_area, 1), "%")


  # ================================================================================
  # CHAPTER 4.14: OUTPUT RASTERS
  # ================================================================================

  message("Saving rasters (clamped and unclamped versions)...")

  terra::writeRaster(pred_current_clamped,
                     file.path(raster_folder, paste0("current_clamped_", species_key, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(pred_current_unclamped,
                     file.path(raster_folder, paste0("current_unclamped_", species_key, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(future_pred_clamped,
                     file.path(raster_folder, paste0("future_clamped_", species_key, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(future_pred_unclamped,
                     file.path(raster_folder, paste0("future_unclamped_", species_key, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(climate_change_clamped,
                     file.path(raster_folder, paste0("change_clamped_", species_key, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(climate_change_unclamped,
                     file.path(raster_folder, paste0("change_unclamped_", species_key, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(clamp_diff_current,
                     file.path(raster_folder, paste0("clamp_effect_current_", species_key, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(clamp_diff_future,
                     file.path(raster_folder, paste0("clamp_effect_future_", species_key, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(pred_sd,
                     file.path(raster_folder, paste0("uncertainty_", species_key, ".tif")),
                     overwrite = TRUE)

  # >> MERGED FROM Script 2: save named predictor stacks with human-readable layer names
  # rename_bioclim_layers() operates on a deep copy so downstream name lookups are unaffected
  tryCatch({
    rename_bioclim_layers <- function(rast_obj) {
      new_names <- sapply(names(rast_obj), function(v) {
        lbl <- bioclim_labels[[tolower(v)]]
        if (!is.null(lbl)) lbl else v
      }, USE.NAMES = FALSE)
      names(rast_obj) <- new_names
      rast_obj
    }

    final_pred_named  <- rename_bioclim_layers(terra::deepcopy(final_predictors))
    future_pred_named <- rename_bioclim_layers(terra::deepcopy(future_predictors_subset))

    terra::writeRaster(final_pred_named,
                       file.path(raster_folder, paste0("final_predictors_named_", species_key, ".tif")),
                       overwrite = TRUE)
    terra::writeRaster(future_pred_named,
                       file.path(raster_folder, paste0("future_predictors_named_", species_key, ".tif")),
                       overwrite = TRUE)

    message("  Saved named predictor stacks (",
            terra::nlyr(final_pred_named), " layers: current | ",
            terra::nlyr(future_pred_named), " layers: future)")
  }, error = function(e) {
    message("Warning: Could not save named predictor stacks: ", e$message)
  })


  # ================================================================================
  # CHAPTER 4.15: VISUALIZATIONS
  # ================================================================================

  message("Creating visualizations...")

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
    geom_sf(data = world, fill = NA, color = "gray50", linewidth = 0.15) +
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
    geom_sf(data = world, fill = NA, color = "gray50", linewidth = 0.15) +
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
    geom_sf(data = world, fill = NA, color = "gray50", linewidth = 0.15) +
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
  if (!is.null(mess_raster)) {
    mess_plot <- ggplot() +
      tidyterra::geom_spatraster(data = mess_raster, maxcell = 5e+07) +
      scale_fill_gradient2(
        low = "#d73027", mid = "#ffffbf", high = "#1a9850",
        midpoint = 0, na.value = "transparent", name = "MESS",
        limits = c(-100, 100), oob = scales::squish
      ) +
      labs(title = paste("MESS Analysis (Extrapolation Risk) -", species_name),
           subtitle = "Red = novel climate (extrapolation, less reliable); Green = within training range",
           x = "Longitude", y = "Latitude") +
      map_theme()

    ggsave(file.path(sdm_folder, paste0("mess_analysis_", species_key, ".png")),
           mess_plot, width = 14, height = 10, dpi = 300, bg = "white")
  } else {
    message("  Skipping MESS plot (MESS analysis failed earlier)")
  }

  message("Visualization plots complete")

  # ---------- 4.15.2 Clamping Comparison Plots (3-Panel) ----------
  message("Creating clamping comparison plots...")

  current_clamped_plt <- ggplot() +
    tidyterra::geom_spatraster(data = pred_current_clamped, maxcell = 5e+07) +
    scale_fill_viridis_c(name = "Suitability", option = "plasma", limits = c(0, 1), na.value = "white") +
    labs(title = "Current (Clamped)", subtitle = paste(species_name)) +
    map_theme() + theme(legend.position = "right")

  current_unclamped_plt <- ggplot() +
    tidyterra::geom_spatraster(data = pred_current_unclamped, maxcell = 5e+07) +
    scale_fill_viridis_c(name = "Suitability", option = "plasma", limits = c(0, 1), na.value = "white") +
    labs(title = "Current (Unclamped)") +
    map_theme() + theme(legend.position = "right")

  current_diff_plt <- ggplot() +
    tidyterra::geom_spatraster(data = clamp_diff_current, maxcell = 5e+07) +
    scale_fill_gradient2(name = "Difference", low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                         midpoint = 0, limits = c(-0.5, 0.5), oob = scales::squish) +
    labs(title = "Difference") +
    map_theme() + theme(legend.position = "right")

  future_clamped_plt <- ggplot() +
    tidyterra::geom_spatraster(data = future_pred_clamped, maxcell = 5e+07) +
    scale_fill_viridis_c(name = "Suitability", option = "plasma", limits = c(0, 1), na.value = "white") +
    labs(title = "Future (Clamped)", subtitle = paste(species_name, "- SSP585 2021-2040")) +
    map_theme() + theme(legend.position = "right")

  future_unclamped_plt <- ggplot() +
    tidyterra::geom_spatraster(data = future_pred_unclamped, maxcell = 5e+07) +
    scale_fill_viridis_c(name = "Suitability", option = "plasma", limits = c(0, 1), na.value = "white") +
    labs(title = "Future (Unclamped)") +
    map_theme() + theme(legend.position = "right")

  future_diff_plt <- ggplot() +
    tidyterra::geom_spatraster(data = clamp_diff_future, maxcell = 5e+07) +
    scale_fill_gradient2(name = "Difference", low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                         midpoint = 0, limits = c(-0.5, 0.5), oob = scales::squish) +
    labs(title = "Difference") +
    map_theme() + theme(legend.position = "right")

  # ---------- 4.15.3 Binary Threshold Maps ----------
  # CONFLICT: Script 1 used SDMtune::plotPA() which is capped at 50,000 cells.
  # Script 2 uses make_binary_map() with tidyterra::geom_spatraster() (maxcell=5e+07).
  # Applying Script 2's implementation for full-resolution binary maps.
  message("Creating threshold maps...")

  ths <- thresholds(final_model, type = "cloglog")

  message("Thresholds:")
  message("  Minimum training presence: ", round(ths[1, 2], 4))
  message("  Equal sens/spec: ", round(ths[2, 2], 4))
  message("  Max TSS: ", round(ths[3, 2], 4))

  # >> MERGED FROM Script 2: make_binary_map() — full-resolution binary maps
  make_binary_map <- function(suitability_rast, threshold, title_label) {
    pa <- terra::ifel(suitability_rast >= threshold, 1L, 0L)
    levels(pa) <- data.frame(value = c(0, 1), label = c("Absent", "Suitable"))
    ggplot() +
      tidyterra::geom_spatraster(data = pa, maxcell = 5e+07) +
      geom_sf(data = world, fill = NA, color = "gray50", linewidth = 0.15) +
      scale_fill_manual(
        values   = c("Absent" = "#5ec962", "Suitable" = "#fde725"),
        na.value = "transparent",
        name     = NULL
      ) +
      coord_sf(expand = FALSE) +
      labs(
        title    = paste(title_label, "-", species_name),
        subtitle = paste("Threshold (maxTSS):", round(threshold, 3),
                         "| AUC:", round(auc_val, 3),
                         "| TSS:", round(tss_val, 3)),
        x = "Longitude", y = "Latitude"
      ) +
      map_theme()
  }

  current_binary_maxTSS <- make_binary_map(pred_current_clamped, ths[3, 2],
                                           "Current Suitable Habitat (Clamped)")
  ggsave(file.path(sdm_folder, "th_current_binary_maxTSS_clamped.png"),
         current_binary_maxTSS, width = 40, height = 30, units = "cm", dpi = 300)

  future_binary_maxTSS <- make_binary_map(future_pred_clamped, ths[3, 2],
                                          "Future Suitable Habitat SSP585 (Clamped)")
  ggsave(file.path(sdm_folder, "th_future_binary_maxTSS_clamped.png"),
         future_binary_maxTSS, width = 40, height = 30, units = "cm", dpi = 300)

  # ---------- 4.15.4 Response Curves ----------
  message("Creating response curves...")

  response_dir <- file.path(sdm_folder, "response_curves")
  if (!dir.exists(response_dir)) dir.create(response_dir, recursive = TRUE)

  for (var_name in names(final_predictors)) {
    tryCatch({
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
             plot = p1, width = 8, height = 6, dpi = 300)

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

      ggsave(filename = file.path(response_dir, paste0("marginal_response_", var_name, ".png")),
             plot = p2, width = 8, height = 6, dpi = 300)

    }, error = function(e) {
      warning("Failed to create response curve for ", var_name, ": ", e$message)
    })
  }

  message("Response curves saved to: ", response_dir)


  # ================================================================================
  # CHAPTER 4.16: MODEVA MODEL EVALUATION
  # >> MERGED FROM Script 2: extended evaluation using modEvA package
  # Reference: https://www.r-bloggers.com/2022/05/model-evaluation-with-presence-points-and-raster-predictions/
  # ================================================================================

  message("Running modEvA model evaluation...")

  modeva_metrics <- NULL   # initialized here; used in Chapter 4.18 summary

  tryCatch({
    modeva_folder <- file.path(sdm_folder, "modeva")
    if (!dir.exists(modeva_folder)) dir.create(modeva_folder, recursive = TRUE)

    pres_pred <- terra::extract(pred_current_clamped, as.matrix(p_coords))
    pres_pred <- pres_pred[, 1]
    pres_pred <- pres_pred[!is.na(pres_pred)]

    bg_pred <- terra::extract(pred_current_clamped, as.matrix(bg_coords))
    bg_pred <- bg_pred[, 1]
    bg_pred <- bg_pred[!is.na(bg_pred)]

    obs_pa  <- c(rep(1, length(pres_pred)), rep(0, length(bg_pred)))
    pred_pa <- c(pres_pred, bg_pred)

    message("  Computing modEvA AUC...")
    modeva_auc <- tryCatch({
      modEvA::AUC(obs = obs_pa, pred = pred_pa, plot = FALSE)
    }, error = function(e) { message("  modEvA::AUC failed: ", e$message); NULL })

    message("  Computing threshold-based measures...")
    modeva_thresh <- tryCatch({
      modEvA::threshMeasures(obs = obs_pa, pred = pred_pa, thresh = "maxTSS",
                             measures = c("TSS", "kappa", "F1score", "Sensitivity", "Specificity"),
                             standardize = FALSE, plot = FALSE)
    }, error = function(e) { message("  modEvA::threshMeasures failed: ", e$message); NULL })

    message("  Computing Miller calibration...")
    modeva_calib <- tryCatch({
      modEvA::MillerCalib(obs = obs_pa, pred = pred_pa, plot = FALSE)
    }, error = function(e) { message("  modEvA::MillerCalib failed: ", e$message); NULL })

    message("  Computing Hosmer-Lemeshow fit...")
    modeva_hl <- tryCatch({
      modEvA::HLfit(obs = obs_pa, pred = pred_pa, bin.method = "quantiles", plot = FALSE)
    }, error = function(e) { message("  modEvA::HLfit failed: ", e$message); NULL })

    message("  Computing Boyce index...")
    modeva_boyce <- tryCatch({
      modEvA::Boyce(obs = obs_pa, pred = pred_pa, pbg = TRUE, plot = FALSE)
    }, error = function(e) { message("  modEvA::Boyce failed: ", e$message); NULL })

    message("  Computing prevalence and evenness...")
    modeva_prev <- tryCatch(modEvA::prevalence(obs = obs_pa), error = function(e) NULL)
    modeva_even <- tryCatch(modEvA::evenness(obs = obs_pa),   error = function(e) NULL)

    metrics_list <- list()
    if (!is.null(modeva_auc))    metrics_list[["AUC_modeva"]] <- modeva_auc$AUC
    if (!is.null(modeva_thresh) && !is.null(modeva_thresh$ThreshMeasures)) {
      thresh_df <- as.data.frame(modeva_thresh$ThreshMeasures)
      for (m in rownames(thresh_df)) metrics_list[[m]] <- thresh_df[m, "Value"]
    }
    if (!is.null(modeva_calib)) {
      metrics_list[["MillerCalib_slope"]]     <- modeva_calib$slope
      metrics_list[["MillerCalib_intercept"]] <- modeva_calib$intercept
    }
    if (!is.null(modeva_boyce) && !is.null(modeva_boyce$Boyce))
      metrics_list[["Boyce_index"]] <- modeva_boyce$Boyce
    if (!is.null(modeva_hl)) {
      if (!is.null(modeva_hl$chi.sq))  metrics_list[["HLfit_chi_sq"]]  <- modeva_hl$chi.sq
      if (!is.null(modeva_hl$df))      metrics_list[["HLfit_df"]]      <- modeva_hl$df
      if (!is.null(modeva_hl$p.value)) metrics_list[["HLfit_p_value"]] <- modeva_hl$p.value
    }
    if (!is.null(modeva_prev)) metrics_list[["prevalence"]] <- modeva_prev
    if (!is.null(modeva_even)) metrics_list[["evenness"]]   <- modeva_even

    if (length(metrics_list) > 0) {
      modeva_metrics <- data.frame(
        metric = names(metrics_list),
        value  = unlist(metrics_list),
        row.names = NULL
      )
      write.csv(modeva_metrics,
                file.path(tables_folder, paste0(species_key, "_modeva_metrics.csv")),
                row.names = FALSE)
      message("  modEvA metrics saved: ", nrow(modeva_metrics), " metrics")
    } else {
      message("  WARNING: all modEvA sub-calls returned NULL; no metrics collected")
    }

    modeva_plots <- list(
      list(file = paste0(species_key, "_modeva_01_roc.png"),
           fn   = function() modEvA::AUC(obs = obs_pa, pred = pred_pa, plot = TRUE,
                                         main = paste0("ROC Curve - ", species_name))),
      list(file = paste0(species_key, "_modeva_02_calibration.png"),
           fn   = function() modEvA::MillerCalib(obs = obs_pa, pred = pred_pa, plot = TRUE,
                                                 main = paste0("Miller Calibration - ", species_name))),
      list(file = paste0(species_key, "_modeva_03_hlfit.png"),
           fn   = function() modEvA::HLfit(obs = obs_pa, pred = pred_pa,
                                           bin.method = "quantiles", plot = TRUE,
                                           main = paste0("Hosmer-Lemeshow Fit - ", species_name))),
      list(file = paste0(species_key, "_modeva_04_threshmeasures.png"),
           fn   = function() modEvA::threshMeasures(obs = obs_pa, pred = pred_pa, thresh = "maxTSS",
                                                    measures = c("TSS", "kappa", "F1score",
                                                                 "Sensitivity", "Specificity"),
                                                    standardize = FALSE, plot = TRUE,
                                                    main = paste0("Threshold Measures - ", species_name))),
      list(file = paste0(species_key, "_modeva_05_preddensity.png"),
           fn   = function() modEvA::predDensity(pred = pred_pa, obs = obs_pa, pbg = TRUE,
                                                 main = paste0("Prediction Density - ", species_name))),
      list(file = paste0(species_key, "_modeva_06_boyce.png"),
           fn   = function() modEvA::Boyce(obs = obs_pa, pred = pred_pa, pbg = TRUE,
                                           bin.width = 0.1, plot = TRUE,
                                           main = paste0("Boyce Index - ", species_name))),
      list(file = paste0(species_key, "_modeva_07_predplot.png"),
           fn   = function() modEvA::predPlot(obs = obs_pa, pred = pred_pa, pbg = TRUE,
                                              thresh = "maxTSS",
                                              main = paste0("Predicted Values - ", species_name)))
    )

    saved_plots <- 0
    for (plt in modeva_plots) {
      tryCatch({
        png(file.path(modeva_folder, plt$file), width = 900, height = 750, res = 150)
        plt$fn()
        dev.off()
        saved_plots <- saved_plots + 1
      }, error = function(e) {
        try(dev.off(), silent = TRUE)
        message("  Plot failed (", plt$file, "): ", e$message)
      })
    }
    message("  modEvA plots saved: ", saved_plots, "/", length(modeva_plots))

  }, error = function(e) {
    message("WARNING: modEvA evaluation failed for ", species_name, ": ", e$message)
    message("Continuing without modEvA metrics...")
    modeva_metrics <<- NULL
  })


  # ================================================================================
  # CHAPTER 4.17: NORWEGIAN CLIMATE HISTOGRAMS
  # >> MERGED FROM Script 2: species climate niche vs Norwegian background
  # ================================================================================

  message("Generating Norwegian climate histograms...")

  tryCatch({
    hist_folder <- file.path(sdm_folder, "climate_histograms")
    if (!dir.exists(hist_folder)) dir.create(hist_folder, recursive = TRUE)

    norway_sf <- world[world$sovereignt == "Norway" & world$type == "Country", ]
    if (nrow(norway_sf) == 0) norway_sf <- world[world$sovereignt == "Norway", ]
    if (nrow(norway_sf) == 0) stop("Could not find Norway in rnaturalearth world polygon")

    norway_rast <- tryCatch({
      terra::mask(terra::crop(final_predictors, norway_sf), terra::vect(norway_sf))
    }, error = function(e) {
      message("  Could not mask to Norway polygon, using crop only: ", e$message)
      terra::crop(final_predictors, norway_sf)
    })

    set.seed(42)
    norway_sample <- tryCatch({
      terra::spatSample(norway_rast, size = 10000, method = "random",
                        na.rm = TRUE, as.df = TRUE)
    }, error = function(e) {
      message("  spatSample failed, using values(): ", e$message)
      as.data.frame(terra::values(norway_rast, na.rm = TRUE))
    })

    pres_env <- as.data.frame(terra::extract(final_predictors, as.matrix(p_coords)))
    pres_env <- pres_env[, names(final_predictors), drop = FALSE]

    message("  Norway sample: ", nrow(norway_sample), " points | Presence: ", nrow(pres_env), " points")

    hist_plots <- list()

    for (var in names(final_predictors)) {
      tryCatch({
        pres_vals <- pres_env[[var]]
        nor_vals  <- norway_sample[[var]]
        pres_vals <- pres_vals[!is.na(pres_vals)]
        nor_vals  <- nor_vals[!is.na(nor_vals)]

        if (length(pres_vals) < 5 || length(nor_vals) < 5) {
          message("  Skipping ", var, ": insufficient data"); next
        }

        x_label   <- if (!is.null(bioclim_labels[[tolower(var)]])) bioclim_labels[[tolower(var)]] else var
        train_min <- min(pres_vals, na.rm = TRUE)
        train_max <- max(pres_vals, na.rm = TRUE)

        plot_df <- rbind(
          data.frame(value = nor_vals,  group = "Norway background"),
          data.frame(value = pres_vals, group = "Species presence")
        )

        hist_p <- ggplot(plot_df, aes(x = value, fill = group, color = group)) +
          geom_histogram(aes(y = after_stat(density)), bins = 40,
                         alpha = 0.5, position = "identity") +
          geom_density(alpha = 0, linewidth = 0.8) +
          geom_vline(xintercept = train_min, linetype = "dashed",
                     color = "steelblue", linewidth = 0.6) +
          geom_vline(xintercept = train_max, linetype = "dashed",
                     color = "steelblue", linewidth = 0.6) +
          scale_fill_manual(values  = c("Norway background" = "grey60",
                                        "Species presence"  = "steelblue")) +
          scale_color_manual(values = c("Norway background" = "grey40",
                                        "Species presence"  = "steelblue4")) +
          labs(title    = paste(species_name, "-", var),
               subtitle = paste("Presence range:", round(train_min, 1), "-", round(train_max, 1)),
               x = x_label, y = "Density", fill = NULL, color = NULL) +
          theme_minimal() +
          theme(legend.position  = "bottom",
                plot.title       = element_text(size = 10, face = "bold"),
                plot.subtitle    = element_text(size = 8))

        ggsave(file.path(hist_folder, paste0(species_key, "_", var, ".png")),
               hist_p, width = 8, height = 6, dpi = 200)

        hist_plots[[var]] <- hist_p

      }, error = function(e) {
        message("  Histogram failed for ", var, ": ", e$message)
      })
    }

    if (length(hist_plots) > 0) {
      n_cols         <- min(3, length(hist_plots))
      composite_plot <- patchwork::wrap_plots(hist_plots, ncol = n_cols) +
        patchwork::plot_annotation(
          title    = paste("Climate Niche vs. Norwegian Background -", species_name),
          subtitle = "Grey = Norway background | Blue = Species presence | Dashed = training range",
          theme    = theme(
            plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 10, hjust = 0.5)
          )
        )

      ggsave(file.path(hist_folder, paste0(species_key, "_all_variables.png")),
             composite_plot,
             width     = n_cols * 8,
             height    = ceiling(length(hist_plots) / n_cols) * 6,
             dpi       = 200,
             limitsize = FALSE)

      message("  Norwegian climate histograms saved: ", length(hist_plots), " variables")
    }

  }, error = function(e) {
    message("WARNING: Norwegian climate histograms failed for ", species_name, ": ", e$message)
    message("Continuing without climate histograms...")
  })


  # ================================================================================
  # CHAPTER 4.18: COMPREHENSIVE MODEL SUMMARY (MD + JSON + TXT)
  # >> MERGED FROM Script 2: three-format machine/human-readable outputs
  # ================================================================================

  message("Writing comprehensive model summary...")

  tryCatch({
    # Climate ranges for each selected variable (from presence data)
    pres_env_for_summary <- tryCatch({
      as.data.frame(terra::extract(final_predictors, as.matrix(p_coords)))
    }, error = function(e) NULL)

    climate_ranges <- list()
    if (!is.null(pres_env_for_summary)) {
      for (var in names(final_predictors)) {
        vals <- pres_env_for_summary[[var]]
        vals <- vals[!is.na(vals)]
        if (length(vals) > 0)
          climate_ranges[[var]] <- list(min = round(min(vals), 2), max = round(max(vals), 2))
      }
    }

    vi_named <- setNames(round(final_vi$Permutation_importance, 4), final_vi$Variable)

    mess_pct_extrap <- tryCatch(round(pct_extrapolation, 1), error = function(e) NA)
    mess_mean_val   <- tryCatch(round(mean_mess, 2),         error = function(e) NA)

    get_modeva <- function(metric_name) {
      if (!is.null(modeva_metrics) && metric_name %in% modeva_metrics$metric)
        modeva_metrics$value[modeva_metrics$metric == metric_name]
      else NA
    }
    modeva_auc_val  <- get_modeva("AUC_modeva")
    modeva_tss_val  <- get_modeva("TSS")
    modeva_kappa    <- get_modeva("kappa")
    modeva_f1       <- get_modeva("F1score")
    miller_slope    <- get_modeva("MillerCalib_slope")
    miller_int      <- get_modeva("MillerCalib_intercept")
    boyce_index_val <- get_modeva("Boyce_index")
    prevalence_val  <- get_modeva("prevalence")
    evenness_val    <- get_modeva("evenness")
    hlfit_p_val     <- get_modeva("HLfit_p_value")

    # --- JSON ---
    json_data <- list(
      species           = species_name,
      species_key       = as.character(species_key),
      date_run          = as.character(Sys.Date()),
      pipeline_version  = "10_updated_sdmtune_loop_merged",
      n_presence        = nrow(p_coords),
      n_background      = nrow(bg_coords),
      variables         = names(final_predictors),
      n_variables       = length(names(final_predictors)),
      variable_importance     = vi_named,
      feature_class     = final_model@model@fc,
      regularization    = final_model@model@reg,
      auc_test          = round(auc_val, 4),
      tss_test          = round(tss_val, 4),
      auc_cv            = round(cv_auc_test, 4),
      tss_cv            = round(cv_tss_test, 4),
      optimal_threshold = round(optimal_threshold, 4),
      sensitivity       = round(sensitivity, 4),
      specificity       = round(specificity, 4),
      modeva_auc              = if (!is.na(modeva_auc_val))  round(modeva_auc_val, 4)  else NULL,
      modeva_tss              = if (!is.na(modeva_tss_val))  round(modeva_tss_val, 4)  else NULL,
      modeva_kappa            = if (!is.na(modeva_kappa))    round(modeva_kappa, 4)    else NULL,
      modeva_f1               = if (!is.na(modeva_f1))       round(modeva_f1, 4)       else NULL,
      miller_calib_slope      = if (!is.na(miller_slope))    round(miller_slope, 4)    else NULL,
      miller_calib_intercept  = if (!is.na(miller_int))      round(miller_int, 4)      else NULL,
      boyce_index             = if (!is.na(boyce_index_val)) round(boyce_index_val, 4) else NULL,
      prevalence              = if (!is.na(prevalence_val))  round(prevalence_val, 4)  else NULL,
      evenness                = if (!is.na(evenness_val))    round(evenness_val, 4)    else NULL,
      hlfit_p_value           = if (!is.na(hlfit_p_val))     round(hlfit_p_val, 4)     else NULL,
      climate_ranges    = climate_ranges,
      mess_pct_extrapolation  = mess_pct_extrap,
      mess_mean         = mess_mean_val,
      ssp585_mean_change      = round(mean_change, 4),
      ssp585_gain_pct         = round(pos_change_area, 1),
      ssp585_loss_pct         = round(neg_change_area, 1)
    )

    jsonlite::write_json(json_data,
                         file.path(summary_folder, "model_summary.json"),
                         pretty = TRUE, auto_unbox = TRUE)

    # --- Markdown ---
    var_table_rows <- paste0(
      sapply(names(vi_named), function(v) {
        imp <- vi_named[v]
        rng <- if (!is.null(climate_ranges[[v]]))
                 paste0(climate_ranges[[v]]$min, " - ", climate_ranges[[v]]$max)
               else "N/A"
        paste0("| ", v, " | ", round(imp, 3), " | ", rng, " |")
      }),
      collapse = "\n"
    )

    mess_text <- if (!is.na(mess_pct_extrap))
      paste0("- **Extrapolation area (MESS < 0):** ", mess_pct_extrap, "%\n",
             "- **Mean MESS score:** ", mess_mean_val)
    else
      "- MESS analysis not available for this run"

    if (!is.null(modeva_metrics) && nrow(modeva_metrics) > 0) {
      modeva_rows <- paste(
        apply(modeva_metrics, 1, function(r)
          paste0("| ", r["metric"], " | ", round(as.numeric(r["value"]), 4), " |")),
        collapse = "\n"
      )
    } else {
      modeva_rows <- "| (modEvA evaluation not available) | - |"
    }

    md_content <- paste0(
      "---\n",
      "species: \"", species_name, "\"\n",
      "date_run: \"", as.character(Sys.Date()), "\"\n",
      "pipeline_version: \"10_updated_sdmtune_loop_merged\"\n",
      "n_presence: ", nrow(p_coords), "\n",
      "n_background: ", nrow(bg_coords), "\n",
      "auc_test: ", round(auc_val, 4), "\n",
      "tss_test: ", round(tss_val, 4), "\n",
      "auc_cv: ", round(cv_auc_test, 4), "\n",
      "variables: [\"", paste(names(final_predictors), collapse = "\", \""), "\"]\n",
      "ssp585_mean_change: ", round(mean_change, 4), "\n",
      "mess_extrapolation_pct: ", ifelse(is.na(mess_pct_extrap), "null", mess_pct_extrap), "\n",
      "boyce_index: ", ifelse(is.na(boyce_index_val), "null", round(boyce_index_val, 4)), "\n",
      "prevalence: ", ifelse(is.na(prevalence_val), "null", round(prevalence_val, 4)), "\n",
      "evenness: ", ifelse(is.na(evenness_val), "null", round(evenness_val, 4)), "\n",
      "---\n\n",
      "# SDM Model Summary: ", species_name, "\n\n",
      "## 1. Species & Data\n\n",
      "- **Species:** *", species_name, "* (key: ", species_key, ")\n",
      "- **Run date:** ", as.character(Sys.Date()), "\n",
      "- **Pipeline version:** 10_updated_sdmtune_loop_merged\n",
      "- **Presence records used:** ", nrow(p_coords), " (after spatial thinning)\n",
      "- **Background points:** ", nrow(bg_coords), " (circle method, 3000 km radius, EPSG:4087)\n",
      "- **Occurrence source:** GBIF / multi-source via spocc + GABI AntMaps\n\n",
      "## 2. Variable Selection\n\n",
      "Two-stage selection: permutation importance (top 15) then Spearman correlation filtering (rho < 0.7).\n\n",
      "| Variable | Importance | Training Range |\n",
      "|----------|-----------|----------------|\n",
      var_table_rows, "\n\n",
      "- **Feature class (FC):** ", final_model@model@fc, "\n",
      "- **Regularization (reg):** ", final_model@model@reg, "\n\n",
      "## 3. Model Performance\n\n",
      "| Metric       | Value |\n",
      "|-------------|-------|\n",
      "| AUC (test)  | ", round(auc_val, 3), " |\n",
      "| TSS (test)  | ", round(tss_val, 3), " |\n",
      "| AUC (CV)    | ", round(cv_auc_test, 3), " |\n",
      "| TSS (CV)    | ", round(cv_tss_test, 3), " |\n",
      "| Sensitivity | ", round(sensitivity, 3), " |\n",
      "| Specificity | ", round(specificity, 3), " |\n",
      "| Threshold   | ", round(optimal_threshold, 3), " |\n\n",
      "### modEvA Extended Metrics\n\n",
      "| Metric  | Value |\n",
      "|---------|-------|\n",
      modeva_rows, "\n\n",
      "## 4. Extrapolation Risk (MESS)\n\n",
      mess_text, "\n\n",
      "## 5. Future Projections (SSP585 2021-2040)\n\n",
      "- **Mean suitability change:** ", round(mean_change, 4), "\n",
      "- **Area with increasing suitability:** ", round(pos_change_area, 1), "%\n",
      "- **Area with decreasing suitability:** ", round(neg_change_area, 1), "%\n\n",
      "## 6. Interpretation for Risk Assessment\n\n",
      "The model for *", species_name, "* achieved a test AUC of ", round(auc_val, 3),
      " and cross-validated AUC of ", round(cv_auc_test, 3), ", indicating ",
      ifelse(cv_auc_test >= 0.8, "good", ifelse(cv_auc_test >= 0.7, "acceptable", "limited")),
      " discriminatory ability. ",
      ifelse(!is.na(mess_pct_extrap) && mess_pct_extrap > 20,
        paste0("Note: ", mess_pct_extrap,
               "% of the prediction area falls outside the training climate envelope (MESS < 0), ",
               "indicating elevated extrapolation risk. Interpret predictions in these zones with caution. "),
        "The majority of the prediction area falls within the training climate envelope (MESS >= 0). "),
      "Under SSP585 2021-2040, the model projects a mean suitability change of ", round(mean_change, 4),
      ", with ", round(pos_change_area, 1), "% of the area showing increasing and ",
      round(neg_change_area, 1), "% showing decreasing suitability.\n\n",
      "## 7. Caveats & Limitations\n\n",
      "- Predictions are based on correlative climate-niche modeling (MaxEnt/Maxnet)\n",
      "- Model trained on ", nrow(p_coords), " presence records",
      ifelse(nrow(p_coords) < 100, " -- low sample size may reduce reliability\n", "\n"),
      "- Future projections use SSP585 (high-emission scenario) only\n",
      "- Biotic interactions, dispersal barriers, and land-use change are not modeled\n",
      "- MESS analysis identifies geographic areas with novel climates not in training data\n"
    )

    writeLines(md_content, file(file.path(summary_folder, "model_summary.md"), encoding = "UTF-8"))

    # --- Plain text ---
    model_quality <- if (cv_auc_test >= 0.8) "good" else if (cv_auc_test >= 0.7) "acceptable" else "limited"

    mess_sentence <- if (!is.na(mess_pct_extrap) && mess_pct_extrap > 20)
      paste0("However, ", mess_pct_extrap,
             "% of the prediction area falls outside the climate conditions the model was ",
             "trained on (extrapolation zones), so predictions in these areas should be treated with caution. ")
    else
      "The model's predictions are largely based on climate conditions represented in the training data. "

    human_text <- paste0(
      "Model Summary for ", species_name, " (", as.character(Sys.Date()), ")\n",
      paste(rep("=", nchar(paste0("Model Summary for ", species_name))), collapse = ""), "\n\n",
      "OVERVIEW\n",
      "This document summarises the results of a species distribution model (SDM) for ",
      species_name, ". The model was run using the SDMtune pipeline ",
      "(10_updated_sdmtune_loop_merged) on ", as.character(Sys.Date()), ".\n\n",
      "DATA\n",
      "The model was built using ", nrow(p_coords), " occurrence records downloaded from GBIF ",
      "and other biodiversity databases (including GABI AntMaps). Background points (n = ",
      nrow(bg_coords), ") were sampled from a 3000 km radius around occurrence points ",
      "(EPSG:4087, distance-accurate projection).\n\n",
      "VARIABLES\n",
      "The model used ", length(names(final_predictors)), " bioclimatic variables: ",
      paste(names(final_predictors), collapse = ", "), ". ",
      "Variables were selected in two stages: first by permutation importance (top 15), then by ",
      "removing correlated variables (Spearman rho > 0.7).\n\n",
      "MODEL PERFORMANCE\n",
      "The model showed ", model_quality, " performance. The test AUC was ", round(auc_val, 3),
      " and the cross-validated AUC was ", round(cv_auc_test, 3), " (5-fold cross-validation). ",
      "The True Skill Statistic (TSS) at the optimal threshold was ", round(tss_val, 3), ". ",
      "Sensitivity was ", round(sensitivity, 3), " and specificity was ", round(specificity, 3), ". ",
      if (!is.na(boyce_index_val)) paste0(
        "The continuous Boyce index was ", round(boyce_index_val, 3),
        ifelse(boyce_index_val > 0.5,
               " (positive: model predicts presence better than random).",
               ifelse(boyce_index_val < 0,
                      " (negative: model may perform worse than random in some areas).",
                      " (near zero: limited discriminatory power along suitability gradient)."))
      ) else "",
      "\n\n",
      "EXTRAPOLATION RISK\n",
      mess_sentence, "\n\n",
      "FUTURE PROJECTIONS\n",
      "Under the high-emission climate scenario (SSP585, 2021-2040), the model predicts a mean ",
      "suitability change of ", round(mean_change, 4), ". Suitability is projected to increase in ",
      round(pos_change_area, 1), "% of the prediction area and decrease in ",
      round(neg_change_area, 1), "% of the area.\n\n",
      "LIMITATIONS\n",
      "This model is correlative and does not account for biotic interactions, dispersal ability, or ",
      "land-use change. The model is based on climate variables only. Results should be interpreted ",
      "in conjunction with expert knowledge and other lines of evidence.\n"
    )

    writeLines(human_text, file(file.path(summary_folder, "model_summary_human.txt"), encoding = "UTF-8"))

    message("  Model summary files written to: ", summary_folder)
    message("  - model_summary.md  (YAML frontmatter + structured markdown)")
    message("  - model_summary.json (machine-readable JSON)")
    message("  - model_summary_human.txt (plain-language narrative)")

  }, error = function(e) {
    message("WARNING: Could not write comprehensive model summary for ", species_name,
            ": ", e$message)
    message("Continuing -- Excel summary will still be written...")
  })


  # ================================================================================
  # CHAPTER 4.19: SUMMARY AND REPORTING (EXCEL)
  # ================================================================================

  # ---------- 4.19.1 Comprehensive Model Summary ----------

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

  # >> MERGED FROM Script 2: Boyce index added to Excel summary
  boyce_for_excel <- tryCatch({
    if (!is.null(modeva_metrics) && "Boyce_index" %in% modeva_metrics$metric)
      round(modeva_metrics$value[modeva_metrics$metric == "Boyce_index"], 4)
    else NA_real_
  }, error = function(e) NA_real_)

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
    boyce_index = boyce_for_excel,   # >> MERGED FROM Script 2
    fc = final_model@model@fc,
    reg = final_model@model@reg,
    cv_folds = CV_FOLDS,
    climate_mean_change = round(mean_change, 4),
    climate_gain_pct = round(pos_change_area, 1),
    climate_loss_pct = round(neg_change_area, 1)
  )

  rio::export(model_summary, file.path(tables_folder, paste0("summary_", species_key, ".xlsx")))

  # ---------- 4.19.2 Species Completion Message ----------
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
