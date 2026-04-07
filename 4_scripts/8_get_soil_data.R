
library(geodata)
library(terra)
library(jsonlite)

# ---------- Configuration ----------
output_dir <- "./1_raw_data/bioclim/unprocessed"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Soil chemical/physical properties
soil_vars <- c(
  "sand",      # Nesting substrate preference
  "clay",      # Nest construction difficulty
  "soc",       # Soil fertility / prey availability
  "phh2o",     # pH - physiological constraint
  "nitrogen",  # Soil fertility
  "bdod"       # Soil compaction / excavation difficulty
)

# Depth: 5 cm (0-5 cm) - primary ant nesting layer
depth <- 5
stat <- "mean"

# ---------- Download Soil Properties ----------
message("=== PART 1: Downloading Soil Properties ===")
message("Variables: ", paste(soil_vars, collapse = ", "))
message("Depth: 0-", depth, " cm")

success_count <- 0
fail_count <- 0

for (var in soil_vars) {
  message("\nDownloading: ", var)
  
  tryCatch({
    soil_raster <- soil_world(
      var = var,
      depth = depth,
      stat = stat,
      path = output_dir
    )
    
    # Save individual file
    output_file <- file.path(output_dir, paste0("soil_", var, "_", depth, "cm.tif"))
    terra::writeRaster(soil_raster, output_file, overwrite = TRUE)
    
    message("  ✓ Saved: ", basename(output_file))
    success_count <- success_count + 1
    
  }, error = function(e) {
    warning("  ✗ Failed to download ", var, ": ", e$message)
    fail_count <- fail_count + 1
  })
  
  Sys.sleep(2)
}

message("\nSoil properties: ", success_count, " successful, ", fail_count, " failed")
