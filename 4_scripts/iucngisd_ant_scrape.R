library(rvest)
library(httr)
library(dplyr)
library(stringr)
library(rio)

# Define the base directory where species data is stored (species folders already exist)
base_dir <- "C:/Users/dafl/OneDrive - Folkehelseinstituttet/VKM Data/27.02.2025_maur_forprosjekt_biologisk_mangfold/data/species"

# Import the master species file
master_file <- "./2_processed_data/master_data.xlsx"
master_data <- rio::import(master_file)

# Function to clean species names for filenames
clean_filename <- function(name) {
  name %>%
    str_replace_all("[^a-zA-Z0-9_]", "_") %>%  # Replace spaces & special chars with "_"
    str_trim() %>%  # Remove leading/trailing spaces
    tolower()       # Convert to lowercase
}

# Base URL for species pages using species name
base_url <- "https://www.iucngisd.org/gisd/speciesname/"

# Loop over each species in the master data
for(i in 1:nrow(master_data)) {
  # Extract speciesKey and species name from the master data
  speciesKey <- as.character(master_data$speciesKey[i])
  species_name <- master_data$species[i]
  
  # Determine the target folder for this species
  species_dir <- file.path(base_dir, speciesKey)
  if (!dir.exists(species_dir)) {
    message("Species folder does not exist for ", species_name, " (", speciesKey, "). Skipping PDF download.")
    next
  }
  
  # Construct the URL for the species page by replacing spaces with "+"
  species_name_url <- gsub(" ", "+", species_name)
  species_url <- paste0(base_url, species_name_url)
  
  message("Processing: ", species_name, " - URL: ", species_url)
  
  # Read the webpage for the species
  page <- tryCatch({
    read_html(species_url)
  }, error = function(e) {
    message("Failed to read page: ", species_url)
    return(NULL)
  })
  
  if (is.null(page)) next  # Skip if the page couldn't be loaded
  
  # Optionally extract the species title from the page (fallback to original name)
  extracted_species_name <- page %>% 
    html_nodes("#spe-title") %>% 
    html_text(trim = TRUE)
  if(length(extracted_species_name) == 0 || extracted_species_name == "") {
    extracted_species_name <- species_name
  }
  
  # Clean the species name for use as a filename
  species_name_clean <- clean_filename(extracted_species_name)
  
  # Attempt to extract the PDF link from an element with id "fa"
  pdf_link <- page %>% 
    html_nodes("#fa") %>% 
    html_attr("href")
  
  # If no PDF link found using the primary selector, try to find any link ending with ".pdf"
  if(length(pdf_link) == 0 || is.na(pdf_link) || pdf_link == "") {
    pdf_link <- page %>% 
      html_nodes("a[href$='.pdf']") %>% 
      html_attr("href")
  }
  
  # Ensure we have a valid PDF link
  if (length(pdf_link) > 0 && !is.na(pdf_link) && pdf_link != "") {
    # Convert to absolute URL if necessary
    if (!grepl("^http", pdf_link)) {
      pdf_link <- paste0("https://www.iucngisd.org/gisd/", pdf_link)
    }
    
    # Define the PDF filename using the cleaned species name within its species folder
    pdf_filename <- file.path(species_dir, paste0(species_name_clean, ".pdf"))
    
    # Check if the PDF is downloadable by sending a HEAD request
    pdf_response <- HEAD(pdf_link)
    content_type <- http_type(pdf_response)
    
    if (!is.null(content_type) && content_type == "application/pdf") {
      # Download the PDF and save it to disk
      download_status <- tryCatch({
        GET(pdf_link, write_disk(pdf_filename, overwrite = TRUE))
        TRUE
      }, error = function(e) FALSE)
      
      if (download_status && file.exists(pdf_filename) && file.size(pdf_filename) > 1000) {
        message("Downloaded: ", pdf_filename)
      } else {
        message("Failed to download a valid PDF for species: ", extracted_species_name)
        if (file.exists(pdf_filename)) file.remove(pdf_filename)  # Remove empty files
      }
    } else {
      message("Invalid PDF URL or content type for species: ", extracted_species_name)
    }
  } else {
    message("No valid PDF found for species: ", extracted_species_name)
  }
}
