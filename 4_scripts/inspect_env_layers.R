# Inspect environmental layers
library(terra)
library(dplyr)

# Load all environmental layers
env_dir <- "./1_raw_data/bioclim/env"
env_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)

message("Found ", length(env_files), " environmental layers")
message("\n=== LOADING AND INSPECTING RASTER PROPERTIES ===\n")

# Create summary dataframe
layer_info <- data.frame(
  file = character(),
  layer_name = character(),
  extent_xmin = numeric(),
  extent_xmax = numeric(),
  extent_ymin = numeric(),
  extent_ymax = numeric(),
  ncols = numeric(),
  nrows = numeric(),
  resolution_x = numeric(),
  resolution_y = numeric(),
  crs = character(),
  min_value = numeric(),
  max_value = numeric(),
  mean_value = numeric(),
  na_count = numeric(),
  total_cells = numeric(),
  stringsAsFactors = FALSE
)

# Loop through files and extract properties
for (i in seq_along(env_files)) {
  file_path <- env_files[i]
  file_name <- basename(file_path)

  message("Processing ", i, "/", length(env_files), ": ", file_name)

  tryCatch({
    r <- terra::rast(file_path)

    layer_info[i, "file"] <- file_name
    layer_info[i, "layer_name"] <- names(r)[1]

    ext <- terra::ext(r)
    layer_info[i, "extent_xmin"] <- ext[1]
    layer_info[i, "extent_xmax"] <- ext[2]
    layer_info[i, "extent_ymin"] <- ext[3]
    layer_info[i, "extent_ymax"] <- ext[4]

    layer_info[i, "ncols"] <- ncol(r)
    layer_info[i, "nrows"] <- nrow(r)

    res <- terra::res(r)
    layer_info[i, "resolution_x"] <- res[1]
    layer_info[i, "resolution_y"] <- res[2]

    layer_info[i, "crs"] <- as.character(terra::crs(r))

    # Calculate statistics (sample if too large)
    if (ncell(r) > 1e7) {
      vals <- terra::spatSample(r, size = 1e6, method = "regular", na.rm = FALSE, as.df = TRUE)[,1]
    } else {
      vals <- terra::values(r, mat = FALSE)
    }

    layer_info[i, "min_value"] <- min(vals, na.rm = TRUE)
    layer_info[i, "max_value"] <- max(vals, na.rm = TRUE)
    layer_info[i, "mean_value"] <- mean(vals, na.rm = TRUE)
    layer_info[i, "na_count"] <- sum(is.na(vals))
    layer_info[i, "total_cells"] <- ncell(r)

  }, error = function(e) {
    message("Error processing ", file_name, ": ", e$message)
  })
}

# Save full report
rio::export(layer_info, "./5_outputs/environmental_layers_report.xlsx")
message("\n=== FULL REPORT SAVED ===")
message("File: ./5_outputs/environmental_layers_report.xlsx\n")

# Categorize layers
message("\n=== LAYER CATEGORIES ===\n")

soil_bioclim <- layer_info %>% filter(grepl("^SBIO", file))
soil_monthly_temp <- layer_info %>% filter(grepl("^soilT_", file))
soil_properties <- layer_info %>% filter(grepl("(bdod|clay|nitrogen|ocd|phh2o|sand|silt|soc)", file))
land_cover <- layer_info %>% filter(grepl("(Barren|Cultivated|Deciduous|Evergreen|Herbaceous|Flooded|Shrubs|Other_Trees)", file))
other_vars <- layer_info %>% filter(grepl("(Aridity|Evapo)", file))

message("Soil Bioclim Variables (SBIO): ", nrow(soil_bioclim))
message("  - These are bioclim-style temperature metrics AT SOIL DEPTH")
print(soil_bioclim %>% select(file, layer_name))

message("\nMonthly Soil Temperature (soilT): ", nrow(soil_monthly_temp))
message("  - Monthly soil temperature values (12 months x 2 depths = 24 layers)")
print(soil_monthly_temp %>% select(file, layer_name))

message("\nSoil Properties: ", nrow(soil_properties))
message("  - Bulk density, clay, nitrogen, pH, sand, silt, SOC")
print(soil_properties %>% select(file, layer_name))

message("\nLand Cover: ", nrow(land_cover))
print(land_cover %>% select(file, layer_name))

message("\nOther Variables: ", nrow(other_vars))
print(other_vars %>% select(file, layer_name))

# Check for resolution consistency
message("\n=== RESOLUTION CONSISTENCY CHECK ===\n")
unique_resolutions <- layer_info %>%
  select(resolution_x, resolution_y) %>%
  distinct() %>%
  arrange(resolution_x)

message("Found ", nrow(unique_resolutions), " unique resolution(s):")
print(unique_resolutions)

if (nrow(unique_resolutions) > 1) {
  message("\nWARNING: Multiple resolutions detected!")
  message("Layers by resolution:")
  for (i in 1:nrow(unique_resolutions)) {
    res_x <- unique_resolutions$resolution_x[i]
    res_y <- unique_resolutions$resolution_y[i]
    matching_layers <- layer_info %>%
      filter(resolution_x == res_x, resolution_y == res_y) %>%
      pull(file)
    message("\n  Resolution ", res_x, " x ", res_y, " (", length(matching_layers), " layers)")
    message("  ", paste(head(matching_layers, 5), collapse = ", "))
    if (length(matching_layers) > 5) message("  ... and ", length(matching_layers) - 5, " more")
  }
}

# Check extent consistency
message("\n=== EXTENT CONSISTENCY CHECK ===\n")
unique_extents <- layer_info %>%
  select(extent_xmin, extent_xmax, extent_ymin, extent_ymax) %>%
  distinct()

message("Found ", nrow(unique_extents), " unique extent(s)")
if (nrow(unique_extents) > 1) {
  message("WARNING: Multiple extents detected - layers may not align!")
}

message("\n=== ANALYSIS COMPLETE ===\n")
