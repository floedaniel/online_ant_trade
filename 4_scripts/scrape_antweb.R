# Load required libraries
library(rvest)
library(rio)
library(stringr)

# Define the path to the master data file
master_file <- "./2_processed_data/master_data.xlsx"

# Import the master data from Excel
master_data <- rio::import(master_file)

# Define the base directory where species data will be stored
base_dir <- "C:/Users/dafl/OneDrive - Folkehelseinstituttet/VKM Data/27.02.2025_maur_forprosjekt_biologisk_mangfold/data/species"

# Loop over each species in the master data
for(i in 1:nrow(master_data)) {
  # Extract speciesKey and species name from the master data
  speciesKey <- as.character(master_data$speciesKey[i])
  species_name_full <- master_data$species[i]
  
  # Construct the URL for the antwiki page by replacing spaces with underscores
  species_name_url <- gsub(" ", "_", species_name_full)
  url <- paste0("https://www.antwiki.org/wiki/", species_name_url)
  
  # Try to read the HTML content from the URL
  webpage <- tryCatch({
    read_html(url)
  }, error = function(e) {
    message("Error reading URL for ", species_name_full, ": ", e)
    return(NULL)
  })
  
  # If the webpage was successfully read, extract the text from the element with class 'mw-page-container'
  if(!is.null(webpage)) {
    page_data <- webpage %>% html_nodes("main") %>% html_text(trim = TRUE)
  } else {
    page_data <- "No data retrieved."
  }
  
  # Define the directory for this species based on its speciesKey
  species_dir <- file.path(base_dir, speciesKey)
  if(!dir.exists(species_dir)) {
    dir.create(species_dir, recursive = TRUE)
  }
  
  # Define the output file name (spaces replaced with underscores)
  file_name <- file.path(species_dir, paste0(gsub(" ", "_", species_name_full), ".txt"))
  
  # Write the scraped data to the file
  writeLines(page_data, con = file_name)
  
  # Optional: pause for a second between requests to be polite
  #Sys.sleep(1)
}
