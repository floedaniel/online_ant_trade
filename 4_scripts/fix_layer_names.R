# Fix layer names in current.tif and future.tif
# Replace hyphens (-) with underscores (_) permanently

library(terra)

message("\n=== FIXING LAYER NAMES IN BIOCLIM RASTERS ===\n")

# File paths
current_file <- "./1_raw_data/bioclim/current/current.tif"
future_file <- "./1_raw_data/bioclim/future/future.tif"

# Create backups
current_backup <- "./1_raw_data/bioclim/current/current_backup.tif"
future_backup <- "./1_raw_data/bioclim/future/future_backup.tif"

message("Creating backups...")
file.copy(current_file, current_backup, overwrite = FALSE)
file.copy(future_file, future_backup, overwrite = FALSE)
message("Backups created: current_backup.tif, future_backup.tif\n")

# Process CURRENT file
message("Processing CURRENT file...")
current <- terra::rast(current_file)
message("  Original layer names:")
print(names(current))

# Fix names
current_names_fixed <- gsub("-", "_", names(current))
names(current) <- current_names_fixed

message("\n  Fixed layer names:")
print(names(current))

# Save with fixed names
message("\n  Saving updated current.tif...")
terra::writeRaster(current, current_file, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
message("  ✓ Current file updated")

# Process FUTURE file
message("\nProcessing FUTURE file...")
future <- terra::rast(future_file)
message("  Original layer names:")
print(names(future))

# Fix names
future_names_fixed <- gsub("-", "_", names(future))
names(future) <- future_names_fixed

message("\n  Fixed layer names:")
print(names(future))

# Save with fixed names
message("\n  Saving updated future.tif...")
terra::writeRaster(future, future_file, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
message("  ✓ Future file updated")

message("\n=== LAYER NAME FIX COMPLETE ===")
message("\nBackups saved as:")
message("  - current_backup.tif")
message("  - future_backup.tif")
message("\nTo restore from backup if needed:")
message('  file.copy("current_backup.tif", "current.tif", overwrite = TRUE)')
message('  file.copy("future_backup.tif", "future.tif", overwrite = TRUE)')
message("\n")
