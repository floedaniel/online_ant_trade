library(fs)
library(readxl)
library(rio)
library(dplyr)

# --- 1. Source empty Endnote library ---
source_enl <- "./species/empty.enl"
source_data <- "./species/empty.Data"

# Check that the source items exist
if (!file_exists(source_enl)) {
  stop("Source .enl file does not exist: ", source_enl)
}
if (!dir_exists(source_data)) {
  stop("Source .Data folder does not exist: ", source_data)
}

# --- 2. Find all multisource.xlsx files recursively ---
base_dir <- "./Species"
multisource_files <- list.files(
  path = base_dir,
  pattern = "_multisource\\.xlsx$",
  recursive = TRUE,
  full.names = TRUE
)

message("Found ", length(multisource_files), " multisource.xlsx files")

if (length(multisource_files) == 0) {
  stop("No multisource.xlsx files found in ", base_dir)
}

# --- 3. Process each species file separately ---
for (species_file in multisource_files) {
  
  message("\n==============================================")
  message("Processing: ", species_file)
  
  # Extract species directory and key
  species_dir <- dirname(species_file)
  species_key <- basename(dirname(species_dir))  # Extract acceptedKey from path
  
  message("Species key: ", species_key)
  message("Species directory: ", species_dir)
  
  # Define target paths (in same folder as multisource.xlsx)
  target_enl_file <- file.path(species_dir, paste0(species_key, ".enl"))
  target_data_folder <- file.path(species_dir, paste0(species_key, ".Data"))
  
  # Copy and rename the .enl file
  if (file_exists(target_enl_file)) {
    message("EndNote library already exists, overwriting: ", target_enl_file)
  }
  file_copy(source_enl, target_enl_file, overwrite = TRUE)
  message("Copied .enl file to: ", target_enl_file)
  
  # Copy and rename the .Data folder along with all its subfolders and files
  if (dir_exists(target_data_folder)) {
    dir_delete(target_data_folder)
    message("Deleted existing .Data folder: ", target_data_folder)
  }
  dir_copy(source_data, target_data_folder, overwrite = TRUE)
  message("Copied .Data folder to: ", target_data_folder)
  
  message("EndNote library distributed for species: ", species_key)
}

message("\n==============================================")
message("EndNote libraries distributed successfully!")
message("==============================================")
