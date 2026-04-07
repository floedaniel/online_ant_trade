# Fix internal layer names in all env folder TIF files
# Replace hyphens (-) with underscores (_) in layer metadata

library(terra)

message("\n=== FIXING LAYER NAMES IN ENV FOLDER ===\n")

env_dir <- "./1_raw_data/bioclim/env"
env_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)

message("Found ", length(env_files), " TIF files in env folder\n")

# Create backup directory
backup_dir <- "./1_raw_data/bioclim/env_backup"
if (!dir.exists(backup_dir)) {
  dir.create(backup_dir, recursive = TRUE)
  message("Created backup directory: ", backup_dir, "\n")
}

files_fixed <- 0
files_with_hyphens <- character()

for (i in seq_along(env_files)) {
  file_path <- env_files[i]
  file_name <- basename(file_path)

  # Load raster
  r <- terra::rast(file_path)

  # Get original layer names
  original_names <- names(r)

  # Fix names (replace - with _)
  fixed_names <- gsub("-", "_", original_names)

  # Check if any changes needed
  if (!identical(original_names, fixed_names)) {
    message(sprintf("[%d/%d] Fixing: %s", i, length(env_files), file_name))
    message("  Original: ", paste(original_names, collapse = ", "))
    message("  Fixed:    ", paste(fixed_names, collapse = ", "))

    # Create backup
    backup_path <- file.path(backup_dir, file_name)
    if (!file.exists(backup_path)) {
      file.copy(file_path, backup_path)
    }

    # Apply fixed names
    names(r) <- fixed_names

    # Save with fixed names (overwrite original)
    terra::writeRaster(r, file_path, overwrite = TRUE, gdal = c("COMPRESS=LZW"))

    files_fixed <- files_fixed + 1
    files_with_hyphens <- c(files_with_hyphens, file_name)
    message("  ✓ Fixed and saved\n")
  }
}

message("\n=== SUMMARY ===")
message("Total files checked: ", length(env_files))
message("Files with hyphens fixed: ", files_fixed)

if (files_fixed > 0) {
  message("\nFiles that were fixed:")
  for (f in files_with_hyphens) {
    message("  - ", f)
  }
  message("\nBackups saved to: ", backup_dir)
} else {
  message("\nNo files needed fixing - all layer names are already correct!")
}

message("\n=== COMPLETE ===\n")
