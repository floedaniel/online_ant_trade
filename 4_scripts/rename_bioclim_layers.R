# Fix Internal Raster Names for All Bioclim Layers
library(terra)

bio_names <- c(
  "BIO1_Annual_Mean_Temperature",
  "BIO2_Mean_Diurnal_Range",
  "BIO3_Isothermality",
  "BIO4_Temperature_Seasonality",
  "BIO5_Max_Temperature_Warmest_Month",
  "BIO6_Min_Temperature_Coldest_Month",
  "BIO7_Temperature_Annual_Range",
  "BIO8_Mean_Temperature_Wettest_Quarter",
  "BIO9_Mean_Temperature_Driest_Quarter",
  "BIO10_Mean_Temperature_Warmest_Quarter",
  "BIO11_Mean_Temperature_Coldest_Quarter",
  "BIO12_Annual_Precipitation",
  "BIO13_Precipitation_Wettest_Month",
  "BIO14_Precipitation_Driest_Month",
  "BIO15_Precipitation_Seasonality",
  "BIO16_Precipitation_Wettest_Quarter",
  "BIO17_Precipitation_Driest_Quarter",
  "BIO18_Precipitation_Warmest_Quarter",
  "BIO19_Precipitation_Coldest_Quarter"
)

# CURRENT
message("Fixing CURRENT bioclim internal names...")
current_dir <- "./1_raw_data/bioclim/current"
for (i in 1:19) {
  file <- file.path(current_dir, paste0(bio_names[i], ".tif"))
  if (file.exists(file)) {
    r <- rast(file)
    names(r) <- bio_names[i]
    temp <- sub("\\.tif$", "_temp.tif", file)
    writeRaster(r, temp, overwrite = TRUE)
    file.remove(file)
    file.rename(temp, file)
    message("  ", bio_names[i])
  }
}

# FUTURE
message("\nFixing FUTURE bioclim internal names...")
future_dir <- "./1_raw_data/bioclim/future"
for (i in 1:19) {
  file <- file.path(future_dir, paste0(bio_names[i], ".tif"))
  if (file.exists(file)) {
    r <- rast(file)
    names(r) <- bio_names[i]
    temp <- sub("\\.tif$", "_temp.tif", file)
    writeRaster(r, temp, overwrite = TRUE)
    file.remove(file)
    file.rename(temp, file)
    message("  ", bio_names[i])
  }
}

# ENV
message("\nFixing ENV layer internal names...")
env_dir <- "./1_raw_data/bioclim/env"
env_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)
for (file in env_files) {
  filename <- tools::file_path_sans_ext(basename(file))
  r <- rast(file)
  names(r) <- filename
  temp <- sub("\\.tif$", "_temp.tif", file)
  writeRaster(r, temp, overwrite = TRUE)
  file.remove(file)
  file.rename(temp, file)
  message("  ", filename)
}

# SBIO
message("\nFixing SBIO layer internal names...")
sbio_dir <- "./1_raw_data/bioclim/sbio"
sbio_files <- list.files(sbio_dir, pattern = "\\.tif$", full.names = TRUE)
for (file in sbio_files) {
  filename <- tools::file_path_sans_ext(basename(file))
  r <- rast(file)
  names(r) <- filename
  temp <- sub("\\.tif$", "_temp.tif", file)
  writeRaster(r, temp, overwrite = TRUE)
  file.remove(file)
  file.rename(temp, file)
  message("  ", filename)
}

message("\nDONE - All internal raster names fixed")
