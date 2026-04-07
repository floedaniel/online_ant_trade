# Quick validation script to check if layer name fixes work
# Run this before running the main SDM loop

library(terra)

message("=== Validating Environmental Layer Names ===\n")

# Load environmental layers (same as main script)
current_dir <- "./1_raw_data/bioclim/current"
future_dir <- "./1_raw_data/bioclim/future"
env_dir <- "./1_raw_data/bioclim/env"

current_files <- list.files(current_dir, pattern = "\\.tif$", full.names = TRUE)
future_files <- list.files(future_dir, pattern = "\\.tif$", full.names = TRUE)
env_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)

message("Loading rasters...")
current_predictors <- terra::rast(current_files)
future_predictors <- terra::rast(future_files)
env_stack <- terra::rast(env_files)

# Combine layers
current_predictors <- c(current_predictors, env_stack)
future_predictors <- c(future_predictors, env_stack)

message("Initial layer counts:")
message("  Current: ", terra::nlyr(current_predictors))
message("  Future: ", terra::nlyr(future_predictors))
message()

# Fix layer names (from script)
current_names <- names(current_predictors)
current_names_fixed <- gsub("-", "_", current_names)
names(current_predictors) <- current_names_fixed

future_names <- names(future_predictors)
future_names_fixed <- gsub("-", "_", future_names)
names(future_predictors) <- future_names_fixed

message("Fixing hyphens in names:")
message("  Current: ", sum(current_names != current_names_fixed), " names fixed")
message("  Future: ", sum(future_names != future_names_fixed), " names fixed")
message()

# Standardize future bioclim layer names
future_names_std <- names(future_predictors)
renamed_count <- 0
for (i in seq_along(future_names_std)) {
  if (grepl("layer\\d+", future_names_std[i])) {
    layer_num <- as.numeric(sub(".*layer(\\d+).*", "\\1", future_names_std[i]))
    future_names_std[i] <- paste0("wc2.1_2.5m_bio_", layer_num)
    renamed_count <- renamed_count + 1
  }
}
names(future_predictors) <- future_names_std

message("Standardizing future bioclim names:")
message("  Renamed ", renamed_count, " layers from 'layer#' to 'wc2.1_2.5m_bio_#'")
message()

# Verify alignment
current_only <- setdiff(names(current_predictors), names(future_predictors))
future_only <- setdiff(names(future_predictors), names(current_predictors))

message("=== Name Alignment Check ===")
if (length(current_only) > 0) {
  message("WARNING: ", length(current_only), " layers in CURRENT but not in FUTURE:")
  for (var in current_only) {
    message("  - ", var)
  }
} else {
  message("✓ All current layers exist in future")
}
message()

if (length(future_only) > 0) {
  message("WARNING: ", length(future_only), " layers in FUTURE but not in CURRENT:")
  for (var in future_only) {
    message("  - ", var)
  }
} else {
  message("✓ All future layers exist in current")
}
message()

# Ensure both stacks have identical layer names
common_layers <- intersect(names(current_predictors), names(future_predictors))
message("Common layers: ", length(common_layers))
current_predictors <- current_predictors[[common_layers]]
future_predictors <- future_predictors[[common_layers]]

message()
message("=== Final Validation ===")
message("Current predictors: ", terra::nlyr(current_predictors), " layers")
message("Future predictors: ", terra::nlyr(future_predictors), " layers")

if (identical(names(current_predictors), names(future_predictors))) {
  message("✓✓✓ SUCCESS: Layer names match perfectly!")
  message()
  message("Sample of layer names (first 10):")
  for (i in 1:min(10, length(names(current_predictors)))) {
    message("  ", i, ". ", names(current_predictors)[i])
  }
} else {
  message("✗✗✗ ERROR: Layer names still don't match!")
  message("This should not happen. Check the name standardization logic.")
}

message()
message("=== Validation Complete ===")
message("If you see 'SUCCESS' above, the layer name fixes are working correctly.")
message("You can now run the main SDM loop with confidence.")
