library(jsonlite)
library(tidyverse)

# https://antapi.org/endpoints/ants


# Function to get all data from antapi.org
lookup_antapi_full <- function(species_name) {
  parts <- str_split(species_name, " ")[[1]]
  
  if (length(parts) < 2) {
    return(tibble(query = species_name, found = FALSE))
  }
  
  genus <- parts[1]
  species <- parts[2]
  
  url <- paste0("https://api.antapi.org/ants/species/", genus, "/", species)
  
  cat("Checking:", species_name, "\n")
  
  result <- tryCatch({
    data <- fromJSON(url, flatten = TRUE)
    
    # Convert the entire result to a tibble
    df <- as_tibble(data)
    df <- df %>% mutate(query = species_name, .before = 1)
    
    return(df)
  }, error = function(e) {
    tibble(query = species_name, found = FALSE)
  })
  
  Sys.sleep(0.5)
  return(result)
}


# Apply to all missing species
df_antapi_full <- map_dfr(ant_data$species, lookup_antapi_full)

df_antapi_full_na <- df_antapi_full %>%  filter(is.na(family))  

df_antapi_full <- df_antapi_full %>%  filter(!is.na(family))  

# Export
write_xlsx(df_antapi_full,     "./2_processed_data/df_antapi_full.xlsx")
write_xlsx(df_antapi_full_na,     "./2_processed_data/df_antapi_full_na.xlsx")

