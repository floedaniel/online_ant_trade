library(europepmc)
library(rentrez)
library(rcrossref)
library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(stringr)
library(openalexR)
library(rplos)
library(rio)

# Set working directories and file paths ------------------------------------
base_dir    <- "./Species"
master_file <- "./2_processed_data/master_data.xlsx"

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

# -------------------- SEARCH FUNCTION ------------------------
query_sources <- function(term) {
  results <- list()
  
  # EuropePMC: Ensure DOI is character
  epmc_doi <- tryCatch({
    res <- epmc_search(query = term, limit = 10000, synonym = FALSE, sort = 'cited')
    res %>%
      filter(!is.na(doi)) %>%
      select(title, doi) %>%
      mutate(doi = as.character(doi), source = "EuropePMC")
  }, error = function(e) NULL)
  
  # PubMed: Ensure DOI is character
  pubmed_doi <- tryCatch({
    search_res <- entrez_search(db = "pubmed", term = term, retmax = 10000)
    if (length(search_res$ids) == 0) return(NULL)
    
    articles <- entrez_summary(db = "pubmed", id = search_res$ids)
    
    # Extract DOI from elocationid and filter for presence
    df_pubmed <- data.frame(
      title = sapply(articles, function(x) if (!is.null(x$title)) x$title else NA),
      doi   = sapply(articles, function(x) {
        loc <- x$elocationid
        if (!is.null(loc) && length(grep("doi:", loc, ignore.case = TRUE)) > 0) {
          return(str_extract(loc, "(?i)doi:[^\\s]*"))
        } else if (length(grep("doi", loc, ignore.case = TRUE)) > 0) {
          return(loc) 
        } else {
          return(NA)
        }
      }),
      stringsAsFactors = FALSE
    ) 
    
    df_pubmed %>%
      filter(!is.na(doi) & doi != "") %>%
      mutate(doi = as.character(doi), source = "PubMed") # Explicit character conversion
  }, error = function(e) NULL)
  
  # CrossRef: Ensure DOI is character
  crossref_doi <- tryCatch({
    cr_res <- cr_works(query = paste0('"', term, '"'), limit = 10000)
    cr_res$data %>%
      as_tibble() %>%
      filter(!is.na(doi)) %>%
      mutate(title = map_chr(title, ~ if (!is.null(.x)) paste(.x, collapse = " ") else NA)) %>%
      select(title, doi) %>%
      mutate(doi = as.character(doi), source = "CrossRef")
  }, error = function(e) NULL)
  
  # OpenAlex: Ensure DOI is character
  openalex_doi <- tryCatch({
    oa_res <- oa_fetch(
      entity = "works",
      output = "list",  # Get as list, not tibble
      search = term,
      verbose = FALSE
    )
    
    if (is.null(oa_res) || length(oa_res) == 0) return(NULL)
    
    # Extract only what we need from the list
    data.frame(
      title = sapply(oa_res, function(x) x$display_name %||% NA),
      doi = sapply(oa_res, function(x) x$doi %||% NA),
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(doi) & doi != "") %>%
      mutate(doi = as.character(doi), source = "OpenAlex")
    
  }, error = function(e) {
    message("    OpenAlex error: ", e$message)
    return(NULL)
  })
  
  # PLOS: Ensure DOI is character
  rplos_doi <- tryCatch({
    plos_res <- searchplos(q = term, fl = c("id", "title", "abstract"), limit = 10000)
    df <- plos_res$data
    if (nrow(df) == 0) return(NULL)
    df <- df %>%
      rename(doi = id) %>%
      filter(!is.na(doi)) %>%
      select(title, doi) %>%
      mutate(doi = as.character(doi), source = "PLOS") # Explicit character conversion
    df
  }, error = function(e) NULL)
  
  # Combine all results. `bind_rows` should now succeed.
  bind_rows(epmc_doi, pubmed_doi, crossref_doi, openalex_doi, rplos_doi)
}

# Read and group species by acceptedKey -------------------------------------
master_species <- rio::import(master_file) %>% 
  as_tibble(,1:20) %>%
  select(acceptedKey, name) %>%
  filter(!is.na(acceptedKey) & !is.na(name)) %>%
  group_by(acceptedKey) %>%
  summarise(
    names = list(name),
    n_synonyms = n(),
    .groups = "drop"
  )

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
  all_results <- list()
  
  for (j in seq_along(species_names_cleaned)) {
    species_name <- species_names_cleaned[j]
    message("  Searching for: ", species_name)
    
    res <- query_sources(species_name)
    
    if (!is.null(res) && nrow(res) > 0) {
      message("    Found ", nrow(res), " articles across all sources")
      all_results[[j]] <- res
    } else {
      message("    No articles found")
    }
    
    Sys.sleep(.3)  # Be nice to the APIs
  }
  
  # Combine all articles and remove duplicates based on DOI --------------------
  if (length(all_results) > 0) {
    combined_articles <- bind_rows(all_results) %>%
      distinct(doi, .keep_all = TRUE)
    
    message("\nTotal unique articles found: ", nrow(combined_articles))
    
    # Generate .ris content ----------------------------------------------------
    ris_content <- combined_articles %>%
      rowwise() %>%
      mutate(
        # Attempt to clean up DOIs by removing 'doi:' prefix if present
        cleaned_doi = gsub("^doi:", "", doi, ignore.case = TRUE),
        ris_entry = paste0(
          "TY  - JOUR\n",
          "TI  - ", title, "\n",
          "DO  - ", cleaned_doi, "\n",
          "N1  - Source: ", source, "\n",
          "ER  - \n\n"
        )
      ) %>%
      pull(ris_entry) %>%
      paste(collapse = "")
    
    # Write the RIS file -------------------------------------------------------
    output_file <- file.path(endnote_folder, paste0(species_key, "_multisource.ris"))
    writeLines(ris_content, con = output_file)
    
    # Write the Excel file -----------------------------------------------------
    rio::export(combined_articles, file.path(endnote_folder, paste0(species_key, "_multisource.xlsx")))
    
    message("RIS file created at: ", output_file)
  } else {
    message("No articles found for acceptedKey: ", species_key)
  }
}

message("\n==============================================")
message("Processing complete!")
