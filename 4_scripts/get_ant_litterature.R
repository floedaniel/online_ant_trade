library(europepmc)
library(rentrez)
library(rcrossref)
library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(stringr)
library(readxl)
library(rio)

# Define paths
target_folder <- "C:/Users/dafl/OneDrive - Folkehelseinstituttet/VKM Data/27.02.2025_maur_forprosjekt_biologisk_mangfold/data/species"
master_file   <- "C:/Users/dafl/OneDrive - Folkehelseinstituttet/VKM Data/27.02.2025_maur_forprosjekt_biologisk_mangfold/data/2_processed_data/master_data.xlsx"

# Read master species file (assumes columns 'speciesKey' and 'species')
master_data <- rio::import(master_file) %>% as_tibble()

# Initialize list to store results (if needed)
all_species_results <- list()

# Loop through each species in the master file
for (i in seq_len(nrow(master_data))) {
  species_key  <- as.character(master_data$speciesKey[i])
  species_name <- master_data$species[i]
  
  cat("\nFetching data for:", species_name, " (Key:", species_key, ")\n")
  
  # Fetch data from EuropePMC
  epmc_doi <- tryCatch({
    epmc_results <- epmc_search(query = species_name, limit = 10000, synonym = TRUE)
    epmc_results %>%
      filter(!is.na(doi) & str_detect(title, regex(species_name, ignore_case = TRUE))) %>%
      select(title, doi) %>%
      mutate(source = "EuropePMC", species = species_name)
  }, error = function(e) return(NULL))
  
  # Fetch data from PubMed
  pubmed_doi <- tryCatch({
    search <- entrez_search(db = "pubmed", term = species_name, retmax = 10000)
    articles <- entrez_summary(db = "pubmed", id = search$ids)
    df <- data.frame(
      title = sapply(articles, function(x) if (!is.null(x$title)) x$title else NA),
      doi   = sapply(articles, function(x) if (!is.null(x$elocationid)) x$elocationid else NA)
    ) %>%
      filter(!is.na(doi)) %>%
      filter(str_detect(title, regex(species_name, ignore_case = TRUE))) %>%
      mutate(source = "PubMed", species = species_name)
    if(nrow(df) == 0) return(NULL) else return(df)
  }, error = function(e) return(NULL))
  
  # Fetch data from CrossRef
  crossref_doi <- tryCatch({
    query_string <- paste0('"', species_name, '"')
    crossref_results <- cr_works(query = query_string, limit = 1000, .progress = "text")
    df <- crossref_results$data %>%
      as_tibble() %>%
      filter(!is.na(doi)) %>%
      mutate(title = map_chr(title, ~ if (!is.null(.x)) paste(.x, collapse = " ") else NA)) %>%
      filter(str_detect(title, regex(species_name, ignore_case = TRUE))) %>%
      select(title, doi) %>%
      mutate(source = "CrossRef", species = species_name)
    if(nrow(df) == 0) return(NULL) else return(df)
  }, error = function(e) return(NULL))
  
  # Combine results from all sources
  species_results <- bind_rows(epmc_doi, pubmed_doi, crossref_doi) %>%
    distinct(doi, .keep_all = TRUE)
  
  if (!is.null(species_results) && nrow(species_results) > 0) {
    all_species_results[[species_name]] <- species_results
    
    # Generate RIS content
    ris_content <- species_results %>%
      rowwise() %>%
      mutate(
        ris_entry = paste0(
          "TY  - JOUR\n",
          "TI  - ", title, "\n",
          "DO  - ", doi, "\n",
          "ER  - \n\n"
        )
      ) %>%
      pull(ris_entry) %>%
      paste(collapse = "")
    
    # Define the Endnote folder for this species:
    # target_folder/<speciesKey>/<speciesKey>_Endnote
    endnote_folder <- file.path(target_folder, species_key, paste0(species_key, "_Endnote"))
    
    # Create the Endnote folder if it doesn't exist
    if (!dir.exists(endnote_folder)) {
      dir.create(endnote_folder, recursive = TRUE)
      message("Created folder: ", endnote_folder)
    }
    
    # Use species name for the RIS file name (sanitize spaces and slashes)
    sanitized_species_name <- gsub("[ /]", "_", species_name)
    species_ris_filename <- file.path(endnote_folder, paste0(sanitized_species_name, ".ris"))
    
    # Write the RIS content to the file
    writeLines(ris_content, con = species_ris_filename)
    message("RIS file created for species: ", species_name, " at ", species_ris_filename)
  } else {
    message("No results found for species: ", species_name)
  }
}

# Combine all species results (if needed)
final_results <- bind_rows(all_species_results)
as_tibble(final_results)
