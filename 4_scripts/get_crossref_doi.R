library(httr)
library(jsonlite)
library(dplyr)
library(readxl)
library(rio)

# Set working directories and file paths ------------------------------------
base_dir    <- "./Species"
master_file <- "./2_processed_data/complete_ant_data.xlsx"

# Helper function to clean species names ------------------------------------
clean_species_name <- function(name) {
  # Remove author names and dates in parentheses
  name_clean <- gsub("\\s*\\([^)]+\\)", "", name)
  # Remove trailing author names and variety designations
  name_clean <- gsub("\\s+(Latreille|Forel|Christ|Santschi|var\\.|subsp\\.).*$", "", name_clean)
  # Trim whitespace
  name_clean <- trimws(name_clean)
  # Get just genus + species (first two words)
  name_parts <- strsplit(name_clean, "\\s+")[[1]]
  if (length(name_parts) >= 2) {
    return(paste(name_parts[1:2], collapse = " "))
  }
  return(name_clean)
}

# Read and group species by acceptedKey -------------------------------------
master_species <- rio::import(master_file) %>%
  as_tibble() %>%
  select(acceptedKey, name) %>%
  filter(!is.na(acceptedKey) & !is.na(name)) %>%
  group_by(acceptedKey) %>%
  summarise(
    names = list(name),
    n_synonyms = n(),
    .groups = "drop"
  )

# Define the base URL for the CrossRef API -------------------------------------
base_url <- "https://api.crossref.org/works"

# Process one acceptedKey at a time --------------------------------------------
for (i in seq_len(nrow(master_species))) {
  species_key <- master_species$acceptedKey[i]
  species_names <- master_species$names[[i]]
  n_synonyms <- master_species$n_synonyms[i]
  
  message("\n==============================================")
  message("Processing acceptedKey: ", species_key)
  message("Number of synonyms: ", n_synonyms)
  message("Synonyms: ", paste(species_names, collapse = " | "))
  
  # Create folder structure ----------------------------------------------------
  species_folder  <- file.path(base_dir, species_key)
  endnote_folder  <- file.path(species_folder, paste0(species_key, "_Endnote"))
  
  if (!dir.exists(endnote_folder)) {
    dir.create(endnote_folder, recursive = TRUE)
  }
  
  # Clean species names and remove duplicates after cleaning ------------------
  species_names_cleaned <- unique(sapply(species_names, clean_species_name))
  message("Unique search terms after cleaning: ", paste(species_names_cleaned, collapse = " | "))
  
  # Collect results from all synonyms -----------------------------------------
  all_articles <- list()
  
  for (j in seq_along(species_names_cleaned)) {
    species_name <- species_names_cleaned[j]
    message("  Searching for: ", species_name)
    
    # Set up query parameters for the CrossRef API -----------------------------
    query_params <- list(
      "query.title" = paste0('"', species_name, '"'),
      rows = 1000,
      filter = paste(
        "type:journal-article",
        "from-pub-date:2000-01-01",
        sep = ","
      ),
      sort = "relevance"
    )
    
    # Make the API request -----------------------------------------------------
    response <- GET(base_url, query = query_params)
    
    if (status_code(response) == 200) {
      content_text <- content(response, "text", encoding = "UTF-8")
      json_data    <- fromJSON(content_text, flatten = TRUE)
      
      if (!is.null(json_data$message$items)) {
        results <- json_data$message$items
        results_df <- as_tibble(results)
        
        # Extract and format publication date
        results_df <- results_df %>%
          mutate(
            published_date = sapply(`issued.date-parts`, function(x) {
              if (!is.null(x)) paste(x[[1]], collapse = "-") else NA
            })
          ) %>%
          mutate(title = sapply(title, function(x) {
            if (!is.null(x)) paste(x, collapse = " ") else NA
          }))
        
        # Dynamically select columns
        columns_to_select <- c("title", "doi" = "DOI", "publisher", "URL", 
                               "volume", "reference-count", "published_date")
        if ("abstract" %in% colnames(results_df)) {
          columns_to_select <- c(columns_to_select, "abstract")
        }
        
        # Select and filter results
        articles <- results_df %>%
          dplyr::select(any_of(columns_to_select)) %>%
          filter(!is.na(published_date), published_date >= "2000") %>%
          arrange(desc(`reference-count`))
        
        # Filter for species name in title or abstract
        if ("abstract" %in% colnames(articles)) {
          articles <- articles %>%
            filter(
              (!is.na(title) & grepl(paste0("\\b", species_name, "\\b"), title, ignore.case = TRUE)) |
                (!is.na(abstract) & grepl(paste0("\\b", species_name, "\\b"), abstract, ignore.case = TRUE))
            )
        } else {
          articles <- articles %>%
            filter(!is.na(title) & grepl(paste0("\\b", species_name, "\\b"), title, ignore.case = TRUE))
        }
        
        if (nrow(articles) > 0) {
          message("    Found ", nrow(articles), " articles")
          all_articles[[j]] <- articles
        } else {
          message("    No matching articles found")
        }
      } else {
        message("    No items returned from API")
      }
    } else {
      message("    API request failed. Status code: ", status_code(response))
    }
    
    Sys.sleep(0.5)  # Be nice to the API
  }
  
  # Combine all articles and remove duplicates based on DOI --------------------
  if (length(all_articles) > 0) {
    combined_articles <- bind_rows(all_articles) %>%
      distinct(doi, .keep_all = TRUE) %>%
      arrange(desc(`reference-count`))
    
    message("\nTotal unique articles found: ", nrow(combined_articles))
    
    # Generate .ris content ----------------------------------------------------
    ris_content <- combined_articles %>%
      rowwise() %>%
      mutate(
        ris_entry = paste0(
          "TY  - JOUR\n",
          "TI  - ", title, "\n",
          "JO  - ", ifelse(is.na(publisher), "", publisher), "\n",
          "VL  - ", ifelse(is.na(volume), "", volume), "\n",
          "UR  - ", ifelse(is.na(URL), "", URL), "\n",
          "DO  - ", doi, "\n",
          "PY  - ", published_date, "\n",
          "ER  - \n\n"
        )
      ) %>%
      pull(ris_entry) %>%
      paste(collapse = "")
    
    # Write the RIS file -------------------------------------------------------
    output_file <- file.path(endnote_folder, paste0(species_key, ".ris"))
    writeLines(ris_content, con = output_file)
    rio::export(combined_articles, file.path(endnote_folder, paste0(species_key, ".xlsx")))
    
    message("RIS file created at: ", output_file)
  } else {
    message("No articles found for acceptedKey: ", species_key)
  }
}

message("\n==============================================")
message("Processing complete!")