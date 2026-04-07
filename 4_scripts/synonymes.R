library(rgbif)
library(rio)
library(dplyr)
library(purrr)

# Import the master data file
master_file <- "./2_processed_data/master_data.xlsx"
master_data <- rio::import(master_file)

# Helper function to retrieve synonyms for a given species key
get_synonyms <- function(species_key) {
  res <- tryCatch({
    name_usage(key = species_key, data = "synonyms")$data
  }, error = function(e) {
    message("Error retrieving synonyms for key ", species_key, ": ", e$message)
    return(NULL)
  })
  
  # If results exist and have at least one row, add the species key as a column.
  if (!is.null(res) && nrow(res) > 0) {
    res <- res %>% mutate(speciesKey = species_key) %>% select(speciesKey, everything())
    return(res)
  } else {
    # If no synonyms found, return a row with NA for synonym-related columns.
    return(data.frame(speciesKey = species_key, synonym = NA, stringsAsFactors = FALSE))
  }
}

# Loop over each species key from master_data and combine the results
all_synonyms <- master_data$speciesKey %>% 
  map_df(~ get_synonyms(.x))

# View the combined synonyms data frame
print(all_synonyms)
