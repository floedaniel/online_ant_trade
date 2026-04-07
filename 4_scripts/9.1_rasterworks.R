library(terra)

# ---------- Configuration ----------
reference_path <- "./1_raw_data/bioclim/current/current.tif"
env_folder <- "./1_raw_data/bioclim/unprocessed_bioclim"
output_folder <- "./1_raw_data/bioclim/future"

if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# ---------- Load Reference ----------
reference <- rast(reference_path)
land_mask <- !is.na(reference[[1]])
message("Reference loaded (", nlyr(reference), " layers)")

# ---------- List Files ----------
all_files <- list.files(env_folder, pattern = "\\.tif$", full.names = TRUE)
message("Found ", length(all_files), " .tif files")

# Rename files with spaces
for (f in all_files) {
  if (grepl(" ", basename(f))) {
    new_name <- gsub(" ", "_", basename(f))
    new_path <- file.path(dirname(f), new_name)
    file.rename(f, new_path)
    message("Renamed: ", basename(f), " → ", new_name)
  }
}

# Re-list after renaming
all_files <- list.files(env_folder, pattern = "\\.tif$", full.names = TRUE)

# Exclude stack files
process_files <- all_files[!grepl("stack", basename(all_files), ignore.case = TRUE)]
message("Processing ", length(process_files), " files")

# ---------- Process Each File ----------
success_count <- 0
fail_count <- 0

for (i in seq_along(process_files)) {
  file_path <- process_files[i]
  file_name <- basename(file_path)
  output_name <- tools::file_path_sans_ext(file_name)
  
  cat("\n[", i, "/", length(process_files), "] ", file_name, "... ", sep = "")
  
  tryCatch({
    # Load
    r <- rast(file_path)
    
    # Take first layer if multi-layer
    if (nlyr(r) > 1) {
      cat("(using layer 1) ")
      r <- r[[1]]
    }
    
    # Set CRS if missing
    if (is.na(crs(r))) crs(r) <- "EPSG:4326"
    
    # Resample, crop, clean, mask
    r <- resample(r, reference, method = "bilinear")
    r <- crop(r, reference)
    r[r == 999] <- NA
    r[r == 9999] <- NA
    r[r == -9999] <- NA
    r <- mask(r, land_mask, maskvalues = 0)
    
    # Save
    writeRaster(r, file.path(output_folder, paste0(output_name, ".tif")), overwrite = TRUE)
    
    cat("✓")
    success_count <- success_count + 1
    
  }, error = function(e) {
    cat("✗ ", e$message)
    fail_count <- fail_count + 1
  })
}

# ---------- Summary ----------
message("\n\n", paste(rep("=", 60), collapse = ""))
message("SUCCESS: ", success_count, " | FAILED: ", fail_count)
message(paste(rep("=", 60), collapse = ""))
message("Output folder: ", output_folder)


# -------------------------------------------------------------------------

# =============================================================================
# Visual Inspection of All Rasters with Proper Color Scaling
# =============================================================================

library(terra)
library(tidyverse)
library(tidyterra)
library(viridis)

# Load all predictors
message("Loading environmental data...")
current_predictors <- rast("./1_raw_data/bioclim/current/current.tif")

env_dir <- "./1_raw_data/bioclim/env"
env_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)
if (length(env_files) > 0) {
  env_stack <- rast(env_files)
  current_predictors <- c(current_predictors, env_stack)
}

message("Total layers: ", nlyr(current_predictors))

# Create output directory
output_dir <- "./output/raster_visual_check"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------- Plot Each Layer ----------
message("\nCreating visual plots for all ", nlyr(current_predictors), " layers...")

for (i in 1:nlyr(current_predictors)) {
  layer_name <- names(current_predictors)[i]
  layer <- current_predictors[[i]]
  
  message("Plotting ", i, "/", nlyr(current_predictors), ": ", layer_name)
  
  # Get statistics
  vals <- values(layer, mat = FALSE, na.rm = TRUE)
  min_val <- min(vals, na.rm = TRUE)
  max_val <- max(vals, na.rm = TRUE)
  mean_val <- mean(vals, na.rm = TRUE)
  n_na <- sum(is.na(values(layer, mat = FALSE)))
  
  # Create plot
  p <- ggplot() +
    tidyterra::geom_spatraster(data = layer) +
    scale_fill_viridis_c(
      name = "Value",
      option = "plasma",
      na.value = "red"  # Show NAs in red
    ) +
    labs(
      title = paste0("Layer ", i, ": ", layer_name),
      subtitle = sprintf("Min: %.2f | Max: %.2f | Mean: %.2f | NAs: %d", 
                         min_val, max_val, mean_val, n_na)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "right"
    )
  
  # Save plot
  filename <- sprintf("%03d_%s.png", i, gsub("[^A-Za-z0-9]", "_", layer_name))
  ggsave(
    file.path(output_dir, filename),
    plot = p,
    width = 12,
    height = 8,
    dpi = 150,
    bg = "white"
  )
}

message("\n✓ All plots saved to: ", output_dir)
message("\nNAs are shown in RED - check for suspicious red patterns!")
