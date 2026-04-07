library(rvest)
library(tidyverse)
library(writexl)

# Initialize empty list to store results
results <- list()

# Loop through pages 1-20
for (pagenumber in 1:20) {
  
  url <- paste0("https://antsatan.com/shop/page/", pagenumber, "/")
  
  cat("Scraping page", pagenumber, "...\n")
  
  # Add error handling in case a page doesn't exist
  page <- tryCatch({
    read_html(url)
  }, error = function(e) {
    cat("Error on page", pagenumber, "\n")
    return(NULL)
  })
  
  if (!is.null(page)) {
    # Extract species titles
    titles <- page %>%
      html_nodes(".woocommerce-loop-product__title") %>%
      html_text()
    
    # Store in list
    if (length(titles) > 0) {
      results[[pagenumber]] <- tibble(
        source = "antsatan.com",
        species = titles
      )
    }
  }
  
  # Be polite - add a small delay between requests
  Sys.sleep(1)
}

# Combine all results into one data frame
df <- bind_rows(results)

df <- df %>%
  mutate(species = str_trim(species),
         species = word(species, 1, 2)) %>% 
  distinct(species, .keep_all = T)

# View results
print(df)
dir()
# Export to Excel
write_xlsx(df, "./1_raw_data/scraped/antsatan_species.xlsx")

cat("\nTotal species scraped:", nrow(df), "\n")
