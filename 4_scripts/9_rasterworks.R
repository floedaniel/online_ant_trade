
library(doParallel)
library(foreach)
library(raster)
library(terra)

# Define the reference raster (using one bioclim file as reference)
reference_path <- "./1_raw_data/bioclim/current/current.tif"
reference <- raster::raster(reference_path)

# (Optional) Load a mask for the study region (e.g., Eurasia) if needed
# e.g., eurasia <- shapefile("path/to/eurasia.shp")

# Set up parallel processing -----------------------------------------------
UseCores <- parallel::detectCores() - 13
cl <- makeCluster(UseCores)
registerDoParallel(cl)

# Define the folder with the environmental rasters
env_folder <- "./1_raw_data/bioclim/unprocces_bioclim"
# List all .tif files in that folder (non-recursive)
stack_list <- list.files(path = env_folder, full.names = TRUE, pattern = "\\.tif$", recursive = FALSE)

# Process each raster in parallel -----------------------------------------
foreach(i = 1:length(stack_list), .packages = c("raster", "terra")) %dopar% {
  # Load the raster using the raster package
  rast <- raster::raster(stack_list[i])
  filename <- basename(stack_list[i])
  
  # Compute min and max values and set projection
  rast <- setMinMax(rast)
  raster::projection(rast) <- '+proj=longlat +datum=WGS84 +no_defs'
  
  # Resample the raster to match the reference (bilinear interpolation)
  rast_resampled <- raster::resample(rast, reference, method = "bilinear")
  
  # Convert to terra for further processing
  rast_terra <- terra::rast(rast_resampled)
  
  # Crop to the extent of the reference raster
  rast_cropped <- terra::crop(rast_terra, terra::ext(reference), snap = "out")

  rast_final <- rast_cropped
  
  # Save the processed raster to the output folder
  output_filename <- file.path("./1_raw_data/bioclim/env/", paste0(filename, "_", i, ".tif"))
  terra::writeRaster(rast_final, output_filename,
                     filetype = "GTiff",
                     datatype = "INT2S",
                     overwrite = TRUE)
}

# Stop the cluster
stopCluster(cl)


# soil fix ----------------------------------------------------------------


library(doParallel)
library(foreach)
library(raster)
library(terra)

# Define the reference raster (using one bioclim file as reference)
reference_path <- "C:/Users/dafl/OneDrive - Folkehelseinstituttet/VKM Data/27.02.2025_maur_forprosjekt_biologisk_mangfold/data/1_raw_data/bioclim/current/current.tif"
reference <- raster::raster(reference_path)

# (Optional) Load a mask for the study region (e.g., Eurasia) if needed
# e.g., eurasia <- shapefile("path/to/eurasia.shp")

# Set up parallel processing -----------------------------------------------
UseCores <- parallel::detectCores() - 13
cl <- makeCluster(UseCores)
registerDoParallel(cl)

# Define the folder with the environmental rasters
env_folder <- "C:/Users/dafl/OneDrive - Folkehelseinstituttet/VKM Data/27.02.2025_maur_forprosjekt_biologisk_mangfold/data/1_raw_data/bioclim/Soil"
# List all .tif files in that folder (non-recursive)
stack_list <- list.files(path = env_folder, full.names = TRUE, pattern = "\\.tif$", recursive = FALSE)

# Process each raster in parallel -----------------------------------------
foreach(i = 1:length(stack_list), .packages = c("raster", "terra")) %dopar% {
  # Load the raster using the raster package
  rast <- raster::raster(stack_list[i])
  filename <- basename(stack_list[i])
  
  # Compute min and max values and set projection
  rast <- setMinMax(rast)
  raster::projection(rast) <- '+proj=longlat +datum=WGS84 +no_defs'
  
  # Resample the raster to match the reference (bilinear interpolation)
  rast_resampled <- raster::resample(rast, reference, method = "bilinear")
  
  # Convert to terra for further processing
  rast_terra <- terra::rast(rast_resampled)
  
  # Crop to the extent of the reference raster
  rast_cropped <- terra::crop(rast_terra, terra::ext(reference), snap = "out")
  
  rast_final <- rast_cropped
  
  # Save the processed raster to the output folder
  output_filename <- file.path("./1_raw_data/bioclim/env/", paste0(filename, "_", i, ".tif"))
  terra::writeRaster(rast_final, output_filename,
                     filetype = "GTiff",
                     datatype = "INT2S",
                     overwrite = TRUE)
}

# Stop the cluster
stopCluster(cl)
