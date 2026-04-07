# Efficient raster resampling script
# Processes each input file and saves to output folder

library(terra)
library(parallel)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Target resolution (degrees)
target_resolution <- 0.04166667  # 2.5 arc-minutes

# Target CRS
target_crs <- "EPSG:4326"  # WGS84

# Target extent: terra::ext(xmin, xmax, ymin, ymax)
target_extent <- terra::ext(-180, 180, -60, 90)

# Input folder containing rasters to process
input_folder <- "./1_raw_data/bioclim/future_raw"

# Output folder
output_folder <- "./1_raw_data/bioclim/future"

                                     
# ============================================================================
# SETUP
# ============================================================================

# Validation
if (!dir.exists(input_folder)) {
  stop("Input folder not found: ", input_folder)
}

# Create output directory
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# List all .tif files
input_files <- list.files(path = input_folder, pattern = "\\.tif$",
                          full.names = TRUE, recursive = FALSE)

if (length(input_files) == 0) {
  stop("No .tif files found in: ", input_folder)
}

# Setup parallel processing (CPU cores - 1)
n_cores <- max(1, detectCores() - 1)

message("\n=== CONFIGURATION ===")
message("Input folder: ", input_folder)
message("Output folder: ", output_folder)
message("Files to process: ", length(input_files))
message("Resolution: ", target_resolution, " degrees")
message("CRS: ", target_crs)
message("Parallel cores: ", n_cores)
message("\n")

# ============================================================================
# PROCESSING
# ============================================================================

# Process single file
process_file <- function(file_path, extent_vec, target_resolution, target_crs, output_folder) {
  # Recreate extent and template
  target_extent <- ext(extent_vec[1], extent_vec[2], extent_vec[3], extent_vec[4])
  template <- rast(target_extent, resolution = target_resolution, crs = target_crs)

  # Load and process
  input_rast <- rast(file_path)
  resampled <- resample(input_rast, template, method = "bilinear")
  cropped <- crop(resampled, target_extent, snap = "out")

  # Save to output folder with same filename
  output_file <- file.path(output_folder, basename(file_path))
  writeRaster(cropped, output_file,
              filetype = "GTiff",
              datatype = "FLT4S",
              overwrite = TRUE,
              gdal = c("COMPRESS=LZW"))

  return(basename(file_path))
}

# Convert extent to numeric vector
extent_vec <- as.vector(target_extent)

# Process in parallel
message("Processing ", length(input_files), " files in parallel...")
start_time <- Sys.time()

if (n_cores > 1 && length(input_files) > 1) {
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, library(terra))
  clusterExport(cl, c("extent_vec", "target_resolution", "target_crs", "output_folder", "process_file"),
                envir = environment())

  results <- parLapply(cl, input_files, function(f) {
    process_file(f, extent_vec, target_resolution, target_crs, output_folder)
  })

  stopCluster(cl)
} else {
  results <- lapply(input_files, function(f) {
    message("  Processing: ", basename(f))
    process_file(f, extent_vec, target_resolution, target_crs, output_folder)
  })
}

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

message("\n=== COMPLETE ===")
message("Processed files: ", length(results))
message("Output folder: ", output_folder)
message("Processing time: ", round(elapsed, 1), " seconds")
message("Done!")
