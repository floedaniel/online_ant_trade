# Load required libraries
library(rio)
library(dplyr)

# Define the path to the master data file
master_file <- "./2_processed_data/complete_ant_data.xlsx"

# Import the master data
master_data <- rio::import(master_file)

# rio::export(master_data, "./2_processed_data/master_data.xlsx")

# Define the base folder for species directories
base_dir <- "./species"

# Create the base folder if it doesn't exist
if (!dir.exists(base_dir)) {
  dir.create(base_dir, recursive = TRUE)
  message("Created base directory: ", base_dir)
}

# Iterate over each species key (using the 'key' column)
for (acceptedKey in master_data$acceptedKey) {
  # Build the path for the new folder using the acceptedKey
  folder_path <- file.path(base_dir, as.character(acceptedKey))
  
  # Check if the folder already exists; if not, create it
  if (!dir.exists(folder_path)) {
    dir.create(folder_path)
    message("Created folder for species: ", acceptedKey)
  } else {
    message("Folder for species ", acceptedKey, " already exists. Skipping.")
  }
}
