# Inspect current and future bioclim rasters
library(terra)
library(dplyr)

message("\n=== INSPECTING CURRENT AND FUTURE BIOCLIM RASTERS ===\n")

# Load rasters
current_file <- "./1_raw_data/bioclim/current/current.tif"
future_file <- "./1_raw_data/bioclim/future/future.tif"

current <- terra::rast(current_file)
future <- terra::rast(future_file)

message("=== CURRENT BIOCLIM ===")
message("File: ", current_file)
message("Number of layers: ", terra::nlyr(current))
message("Layer names:")
print(names(current))

message("\nExtent:")
print(terra::ext(current))

message("\nResolution (degrees):")
print(terra::res(current))

message("\nCRS:")
print(terra::crs(current))

message("\nDimensions (rows x cols x layers): ", nrow(current), " x ", ncol(current), " x ", nlyr(current))
message("Total cells per layer: ", format(ncell(current), big.mark = ","))

# Sample statistics for first layer
message("\nSample statistics (first layer: ", names(current)[1], "):")
vals <- terra::spatSample(current[[1]], size = 10000, method = "regular", na.rm = TRUE, as.df = TRUE)[,1]
message("  Min: ", round(min(vals, na.rm = TRUE), 2))
message("  Max: ", round(max(vals, na.rm = TRUE), 2))
message("  Mean: ", round(mean(vals, na.rm = TRUE), 2))

message("\n", paste(rep("=", 70), collapse = ""))
message("\n=== FUTURE BIOCLIM ===")
message("File: ", future_file)
message("Number of layers: ", terra::nlyr(future))
message("Layer names:")
print(names(future))

message("\nExtent:")
print(terra::ext(future))

message("\nResolution (degrees):")
print(terra::res(future))

message("\nCRS:")
print(terra::crs(future))

message("\nDimensions (rows x cols x layers): ", nrow(future), " x ", ncol(future), " x ", nlyr(future))
message("Total cells per layer: ", format(ncell(future), big.mark = ","))

# Sample statistics for first layer
message("\nSample statistics (first layer: ", names(future)[1], "):")
vals <- terra::spatSample(future[[1]], size = 10000, method = "regular", na.rm = TRUE, as.df = TRUE)[,1]
message("  Min: ", round(min(vals, na.rm = TRUE), 2))
message("  Max: ", round(max(vals, na.rm = TRUE), 2))
message("  Mean: ", round(mean(vals, na.rm = TRUE), 2))

message("\n", paste(rep("=", 70), collapse = ""))
message("\n=== COMPATIBILITY CHECK ===")

# Check if current and future have same structure
message("\nLayer names match: ", identical(names(current), names(future)))
message("Extent match: ", terra::compareGeom(current, future, stopOnError = FALSE))
message("Resolution match: ", all(terra::res(current) == terra::res(future)))
message("Number of layers match: ", nlyr(current) == nlyr(future))

# Create detailed layer report
layer_report <- data.frame(
  layer_number = 1:nlyr(current),
  layer_name_current = names(current),
  layer_name_future = names(future),
  match = names(current) == names(future)
)

message("\nLayer-by-layer comparison:")
print(layer_report)

# Save report
rio::export(layer_report, "./5_outputs/bioclim_layer_comparison.xlsx")
message("\nReport saved to: ./5_outputs/bioclim_layer_comparison.xlsx")

# Compare to env folder
message("\n", paste(rep("=", 70), collapse = ""))
message("\n=== COMPARISON WITH ENV FOLDER ===")

env_dir <- "./1_raw_data/bioclim/env"
env_files <- list.files(env_dir, pattern = "\\.tif$", full.names = FALSE)

message("\nCurrent/Future layers: ", nlyr(current))
message("Env folder layers: ", length(env_files))
message("Total environmental variables: ", nlyr(current) + length(env_files))

message("\n=== INSPECTION COMPLETE ===\n")
