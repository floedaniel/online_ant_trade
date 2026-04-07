library(dplyr)
library(purrr)
library(rvest)
library(httr)

# loop over pages 0-20 and extract .title elements
scrape_antstore <- function(page_num) {
  url <- paste0("https://www.antstore.net/shop/de/Ameisen/?page=", page_num)
  
  tryCatch({
    response <- GET(url, timeout(10))
    
    if (status_code(response) != 200) return(NULL)
    
    page <- read_html(response)
    
    # extract all .title elements
    titles <- page %>%
      html_elements(".title") %>%
      html_text() %>%
      trimws()
    
    tibble(
      page = page_num,
      species_name = titles
    )
  }, error = function(e) {
    NULL
  })
}

# scrape pages 0-20
antstore_species <- map(
  0:20,
  ~{
    Sys.sleep(0.5) # rate limit
    scrape_antstore(.x)
  },
  .progress = TRUE
) %>%
  compact() %>%
  bind_rows()

# remove page column and keep only unique species
antstore_species <- antstore_species %>%
  select(species_name) %>%
  distinct()

antstore_species

# -------------------------------------------------------------------------

library(dplyr)
library(purrr)
library(stringr)
library(httr)
library(jsonlite)

# clean species names - remove common names in parentheses and subgenus
clean_species_names <- antstore_species %>%
  mutate(
    # remove text in parentheses (common names)
    species_clean = str_remove(species_name, "\\s*\\([^)]+\\)"),
    # trim whitespace
    species_clean = str_trim(species_clean)
  )

# gbif species match function
gbif_match_species <- function(species_name) {
  base_url <- "https://api.gbif.org/v1/species/match"
  
  tryCatch({
    response <- GET(
      base_url,
      query = list(name = species_name, verbose = TRUE)
    )
    
    if (status_code(response) == 200) {
      result <- content(response, as = "parsed")
      
      tibble(
        query_name = species_name,
        matched_name = result$canonicalName %||% NA_character_,
        scientific_name = result$scientificName %||% NA_character_,
        status = result$status %||% NA_character_,
        confidence = result$confidence %||% NA_integer_,
        match_type = result$matchType %||% NA_character_,
        kingdom = result$kingdom %||% NA_character_,
        phylum = result$phylum %||% NA_character_,
        class = result$class %||% NA_character_,
        order = result$order %||% NA_character_,
        family = result$family %||% NA_character_,
        genus = result$genus %||% NA_character_,
        species = result$species %||% NA_character_,
        rank = result$rank %||% NA_character_,
        gbif_id = result$usageKey %||% NA_integer_
      )
    } else {
      tibble(
        query_name = species_name,
        matched_name = NA_character_,
        scientific_name = NA_character_,
        status = "NO_MATCH",
        confidence = NA_integer_,
        match_type = NA_character_,
        kingdom = NA_character_,
        phylum = NA_character_,
        class = NA_character_,
        order = NA_character_,
        family = NA_character_,
        genus = NA_character_,
        species = NA_character_,
        rank = NA_character_,
        gbif_id = NA_integer_
      )
    }
  }, error = function(e) {
    tibble(
      query_name = species_name,
      matched_name = NA_character_,
      scientific_name = NA_character_,
      status = "ERROR",
      confidence = NA_integer_,
      match_type = NA_character_,
      kingdom = NA_character_,
      phylum = NA_character_,
      class = NA_character_,
      order = NA_character_,
      family = NA_character_,
      genus = NA_character_,
      species = NA_character_,
      rank = NA_character_,
      gbif_id = NA_integer_
    )
  })
}

# match all species to gbif
antstore_validated <- clean_species_names %>%
  pull(species_clean) %>%
  unique() %>%
  map(
    ~{
      Sys.sleep(0.3) # rate limit for gbif api
      gbif_match_species(.x)
    },
    .progress = TRUE
  ) %>%
  bind_rows()

# join back with original names
antstore_final <- antstore_species %>%
  mutate(species_clean = str_remove(species_name, "\\s*\\([^)]+\\)")) %>%
  mutate(species_clean = str_trim(species_clean)) %>%
  left_join(antstore_validated, by = c("species_clean" = "query_name"))

# save results
#write.csv(antstore_final, "antstore_species_validated.csv", row.names = FALSE)

# summary statistics
cat("\n=== VALIDATION SUMMARY ===\n")
cat("Total species queried:", nrow(antstore_final), "\n")
cat("Exact matches:", sum(antstore_final$match_type == "EXACT", na.rm = TRUE), "\n")
cat("Fuzzy matches:", sum(antstore_final$match_type == "FUZZY", na.rm = TRUE), "\n")
cat("No matches:", sum(is.na(antstore_final$matched_name)), "\n")
cat("Accepted names:", sum(antstore_final$status == "ACCEPTED", na.rm = TRUE), "\n")
cat("Synonyms:", sum(antstore_final$status == "SYNONYM", na.rm = TRUE), "\n\n")

# show problematic matches
low_confidence <- antstore_final %>%
  filter(confidence < 90 | is.na(matched_name)) %>%
  select(species_name, matched_name, confidence, match_type, status)

if (nrow(low_confidence) > 0) {
  cat("Low confidence or failed matches:\n")
  print(low_confidence, n = 20)
}

antstore_final

# save to csv
write.csv(antstore_final, "2_processed_data/antstore_species.csv", row.names = FALSE)


# -------------------------------------------------------------------------

antstore_final <- rio::import("2_processed_data/antstore_species.csv", row.names = FALSE)

antstore_final %>% filter(matched_name %in% ant_species$species )

antstore_final %>% filter(!matched_name %in% ant_species$species )



