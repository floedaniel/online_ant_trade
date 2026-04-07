library(rvest)
library(tidyverse)
library(writexl)

# Define URLs to scrape
urls <- c(
  "https://www.fourmishome.fr/Africa-bbeaaaaab.asp",
  "https://www.fourmishome.fr/America-bbeaaaaac.asp",
  "https://www.fourmishome.fr/Europe-bbeaaaaad.asp",
  "https://www.fourmishome.fr/Asia-bbeaaaaae.asp",
  "https://www.fourmishome.fr/Australia-bbeaaaaaf.asp"
)

# Initialize empty list to store results
results <- list()

# Loop through URLs
for (i in seq_along(urls)) {
  
  url <- urls[i]
  
  cat("Scraping", url, "...\n")
  
  # Add error handling in case a page doesn't exist
  page <- tryCatch({
    read_html(url)
  }, error = function(e) {
    cat("Error on", url, "\n")
    return(NULL)
  })
  
  if (!is.null(page)) {
    # Extract species titles
    titles <- page %>%
      html_nodes(".x3_nom") %>%
      html_text()
    
    # Store in list
    if (length(titles) > 0) {
      results[[i]] <- tibble(
        source = "fourmishome",
        species = titles
      )
    }
  }
  
  # Be polite - add a small delay between requests
  Sys.sleep(1)
}

# Combine all results into one data frame
df <- bind_rows(results)

df

# View results
print(df)

# Export to Excel
write_xlsx(df, "./1_raw_data/fourmishome_species.xlsx")

cat("\nTotal species scraped:", nrow(df), "\n")
