# Data from 
# https://antcheck.info/
# API https://antcheck.info/api/v2/docs#/
# AntCheck API Data Retrieval Script
# This script retrieves ant sales data from the AntCheck.info API
# AntCheck API - Simple data retrieval
library(httr)
library(jsonlite)
library(dplyr)
library(stringr)

# Read API key
api_key <- readLines("C:/Users/dafl/Desktop/API keys/antcheck.txt", warn = FALSE)[1]
base_url <- "https://antcheck.info/api/v2"

# Function to get all data (limit=-1 returns everything)
get_data <- function(endpoint) {
  url <- paste0(base_url, endpoint, "?limit=-1")
  response <- GET(url, add_headers("X-Api-Key" = api_key))
  data <- fromJSON(content(response, "text", encoding = "UTF-8"), flatten = TRUE)
  return(as.data.frame(data))
}

# Fetch all data
species_df <- get_data("/ants/species")
products_df <- get_data("/ecommerce/products")
variants_df <- get_data("/ecommerce/variants")
shops_df <- get_data("/ecommerce/shops")

# -------------------------------------------------------------------------
# Clean titles - remove special chars except parentheses
clean_title <- function(x) {
  x <- gsub("<[^>]+>", "", x)  # Remove HTML
  x <- gsub("[^A-Za-z0-9()\\s'-]", " ", x)  # Keep letters, numbers, (), -, '
  x <- gsub("\\s+", " ", x)  # Multiple spaces to single
  x <- trimws(x)  # Trim
  ifelse(x == "", NA, x)
}

# Clean descriptions - remove HTML
clean_description <- function(x) {
  x <- gsub("<[^>]+>", "", x)  # Remove HTML tags
  x <- gsub("&nbsp;", " ", x)  # HTML entities
  x <- gsub("&amp;", "&", x)
  x <- gsub("&lt;", "<", x)
  x <- gsub("&gt;", ">", x)
  x <- gsub("&quot;", '"', x)
  x <- gsub("&#39;", "'", x)
  x <- gsub("\\s+", " ", x)  # Multiple spaces to single
  x <- trimws(x)  # Trim
  ifelse(x == "", NA, x)
}

# Apply cleaning to products
products_df <- products_df %>%
  mutate(
    title = clean_title(title),
    description = clean_description(description),
    comment = clean_title(comment)
  ) %>%
  mutate(across(where(is.character), ~ifelse(. == "" | trimws(.) == "", NA, .)))

# -------------------------------------------------------------------------
# Filter for ants only
ants <- products_df %>% 
  filter(product_type == "ants") %>% 
  as_tibble()

ants

# -------------------------------------------------------------------------

# Extract species names using GBIF backbone
library(readr)

gbif_path <- "C:/Users/dafl/OneDrive - Folkehelseinstituttet/VKM Data/04.09.2025_Betula_Nordic_Network/data/1_raw_data/Gbif backbone/backbone/Taxon.tsv"
gbif <- read_tsv(gbif_path, show_col_types = FALSE)

# Function to extract species name (GBIF validated)
# Function to extract species name (GBIF validated)
extract_species_gbif <- function(title, lookup_lower) {
  if (is.na(title)) return(NA_character_)
  
  # Try text in parentheses
  in_parens <- str_extract(title, "\\(([A-Za-z]+ [a-z]+)\\)")
  if (!is.na(in_parens)) {
    candidate <- str_remove_all(in_parens, "[\\(\\)]")
    if (tolower(candidate) %in% lookup_lower) return(candidate)
  }
  
  # Try at start
  at_start <- str_extract(title, "^([A-Za-z]+ [a-z]+)")
  if (!is.na(at_start) && tolower(at_start) %in% lookup_lower) return(at_start)
  
  # Try any pattern
  patterns <- str_extract_all(title, "[A-Za-z]+ [a-z]+")[[1]]
  for (pattern in patterns) {
    if (tolower(pattern) %in% lookup_lower) return(pattern)
  }
  
  return(NA_character_)
}

# Function to extract fuzzy match
extract_species_fuzzy <- function(title) {
  if (is.na(title)) return(NA_character_)
  
  # Handle: Genus (Subgenus) species
  parts <- str_match(title, "^([A-Za-z]+) \\([A-Za-z]+\\) ([a-z]+)")
  if (!is.na(parts[1])) {
    return(paste(parts[2], parts[3]))
  }
  
  # Extract 2-3 word pattern in parentheses
  in_parens <- str_extract(title, "\\(([A-Za-z]+ [a-z]+(?: [a-z]+)?)\\)")
  if (!is.na(in_parens)) {
    return(str_remove_all(in_parens, "[\\(\\)]"))
  }
  
  # Try at start
  at_start <- str_extract(title, "^([A-Za-z]+ [a-z]+(?: [a-z]+)?)")
  if (!is.na(at_start)) return(at_start)
  
  # Try any 2-word pattern
  patterns <- str_extract_all(title, "[A-Za-z]+ [a-z]+")[[1]]
  if (length(patterns) > 0) return(patterns[1])
  
  return(NA_character_)
}

# Extract species names
ants <- ants %>%
  mutate(
    gbif_species = sapply(title, extract_species_gbif, 
                          lookup_lower = tolower(ants_gbif$genus_species)),
    fuzzy_match = sapply(title, extract_species_fuzzy)
  )

ants

ants <- ants %>%
  mutate(species_name = coalesce(gbif_species, fuzzy_match))

rio::export(ants, "./2_processed_data/AntCheck_data.xlsx")

ants <- ants %>% select(-created_at, -updated_at, -changed_at)

ants %>%  distinct(species_name, .keep_all = T)

ants %>% filter(in_stock=="TRUE") %>% distinct(species_name, .keep_all = T)

ants %>% filter(is_active=="TRUE") %>% distinct(species_name, .keep_all = T)

ants %>% filter(is_active=="TRUE" & in_stock=="TRUE") %>% distinct(species_name, .keep_all = T)


ants <- ants %>% 
  filter(str_detect(description, regex("invasive", ignore_case = TRUE)))

ants <- ants[1:10,]
