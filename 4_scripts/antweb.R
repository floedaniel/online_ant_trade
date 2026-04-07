# main dir
# https://www.antweb.org/web/
# https://www.antweb.org/web/data/specimenData/
# https://www.antweb.org/web/workingdir/

# read ant web data 
library(readr)

library(stringr)
species_list <- read_delim("https://www.antweb.org/web/workingdir/worldants_speciesList.txt", delim = "\t")

# ============================================
# STEP 1: Create species name column
# ============================================

species_list$species_name <- NA

# If only genus (no species), use genus
species_list$species_name[is.na(species_list$species)] <- species_list$genus[is.na(species_list$species)]

# If has species but no subspecies, use genus + species
has_species <- !is.na(species_list$species) & is.na(species_list$subspecies)

species_list$species_name[has_species] <- paste(species_list$genus[has_species], 
                                                species_list$species[has_species])

# If has subspecies, use genus + species + subspecies
has_subspecies <- !is.na(species_list$species) & !is.na(species_list$subspecies)
species_list$species_name[has_subspecies] <- paste(species_list$genus[has_subspecies], 
                                                   species_list$species[has_subspecies],
                                                   species_list$subspecies[has_subspecies])

# ============================================
# STEP 2: Clean author_date_html column
# ============================================

# Remove HTML tags
species_list$author_clean <- str_replace_all(species_list$`author date html`, "<[^>]+>", "")

# Fix Unicode characters
species_list$author_clean <- str_replace_all(species_list$author_clean, "<U\\+00E9>", "é")
species_list$author_clean <- str_replace_all(species_list$author_clean, "<U\\+00E0>", "à")
species_list$author_clean <- str_replace_all(species_list$author_clean, "<U\\+00E7>", "ç")
species_list$author_clean <- str_replace_all(species_list$author_clean, "<U\\+00F4>", "ô")
species_list$author_clean <- str_replace_all(species_list$author_clean, "<U\\+00EA>", "ê")
species_list$author_clean <- str_replace_all(species_list$author_clean, "<U\\+2020>", "†")

# Remove extra whitespace
species_list$author_clean <- str_trim(species_list$author_clean)
species_list$author_clean <- str_squish(species_list$author_clean)

# ============================================
# STEP 3: Clean taxonomic_history_html column (if exists)
# ============================================

if("taxonomic history html" %in% names(species_list)) {
  species_list$history_clean <- str_replace_all(species_list$`taxonomic history html`, "<[^>]+>", "")
  species_list$history_clean <- str_replace_all(species_list$history_clean, "<U\\+00E9>", "é")
  species_list$history_clean <- str_replace_all(species_list$history_clean, "<U\\+00E0>", "à")
  species_list$history_clean <- str_replace_all(species_list$history_clean, "<U\\+00E7>", "ç")
  species_list$history_clean <- str_replace_all(species_list$history_clean, "<U\\+2020>", "†")
  species_list$history_clean <- str_trim(species_list$history_clean)
  species_list$history_clean <- str_squish(species_list$history_clean)
}

# -------------------------------------------------------------------------
# Clean taxonomic_history_html column

# Remove HTML tags
species_list$history_clean <- str_replace_all(species_list$`taxonomic history html`, "<[^>]+>", "")

# Fix HTML entities
species_list$history_clean <- str_replace_all(species_list$history_clean, "&quot;", '"')
species_list$history_clean <- str_replace_all(species_list$history_clean, "&amp;", "&")
species_list$history_clean <- str_replace_all(species_list$history_clean, "&lt;", "<")
species_list$history_clean <- str_replace_all(species_list$history_clean, "&gt;", ">")

# Fix Unicode characters
species_list$history_clean <- str_replace_all(species_list$history_clean, "<U\\+00E9>", "é")
species_list$history_clean <- str_replace_all(species_list$history_clean, "<U\\+00E0>", "à")
species_list$history_clean <- str_replace_all(species_list$history_clean, "<U\\+00E7>", "ç")
species_list$history_clean <- str_replace_all(species_list$history_clean, "<U\\+00F4>", "ô")
species_list$history_clean <- str_replace_all(species_list$history_clean, "<U\\+00EA>", "ê")
species_list$history_clean <- str_replace_all(species_list$history_clean, "<U\\+2020>", "†")

# Remove extra whitespace
species_list$history_clean <- str_squish(species_list$history_clean)
species_list$history_clean <- str_trim(species_list$history_clean)

# ============================================
# STEP 4: Clean country column
# ============================================

# Remove parenthetical info like "(Silhouette)"
species_list$country_clean <- str_replace_all(species_list$country, "\\s*\\([^)]+\\)", "")
species_list$country_clean <- str_trim(species_list$country_clean)

# ============================================
# STEP 5: Move species_name to first column
# ============================================

# Get all column names except species_name
other_cols <- setdiff(names(species_list), "species_name")

# Reorder: species_name first, then everything else
species_list <- species_list[, c("species_name", other_cols)]

# ============================================
# DONE! Check your data
# ============================================

# View first few rows
head(species_list[, c("species_name", "author_clean", "country_clean")])

species_list <-species_list %>% filter(!fossil=="TRUE")

species_list$`taxonomic history html`

issues <- rio::import("./2_processed_data/ant_issues.xlsx")

t <- left_join(issues, species_list, by=c("species"="species_name"))

t %>% select("original_name", "final_name", "species","current valid rank", "current valid name",  "history_clean", "original combination" )

# ditribution data ? -------------------------------------------------------------------------

library(readr)

dir_path <- "C:/Users/dafl/Desktop/antweb"
txt_files <- list.files(dir_path, pattern = "\\.txt$", full.names = TRUE)

# Read with error handling
data_list <- lapply(txt_files, function(file) {
  tryCatch({
    read_delim(file, delim = "\t", show_col_types = FALSE)
  }, error = function(e) {
    message(paste("Error reading", basename(file), ":", e$message))
    return(NULL)
  })
})

names(data_list) <- basename(txt_files)

# Remove any NULL entries (failed reads)
data_list <- data_list[!sapply(data_list, is.null)]

library(dplyr)
library(purrr)

data_list <- map(data_list, ~mutate(.x, across(everything(), as.character)))
combined_data <- bind_rows(data_list, .id = "file_name")


# -------------------------------------------------------------------------


          