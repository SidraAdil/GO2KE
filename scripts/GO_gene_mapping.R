############################################################
# GO_gene_mapping.R
#
# Purpose:
# Assign Gene Ontology biological process gene sets to
# Key Events (KEs) using the GO API and convert HGNC
# symbols to Ensembl gene identifiers.
#
# Input:
#   v4_brain.tsv
#   v3_kidney.tsv
#   v4_liver.tsv
#   v3_lungs.tsv
#
# Output:
#   g_brain.tsv
#   g_kidney.tsv
#   g_liver.tsv
#   g_lungs.tsv
############################################################
############################################################
# 1. Load input files
############################################################
 
# calling libraries
library(readr)
library(tidyverse)

# setting working dir

setwd('')
getwd()

# this data has same naming conventions, with no duplication (has distinct keids) and missing values handled (NA added)

list.files(path = '',
           pattern = ".tsv$")

path = ''

brain <- read.delim(file = file.path(path, 'v4_brain.tsv'),
                    sep = '\t')
kidney <- read.delim(file = file.path(path, 'v3_kidney.tsv'),
                     sep = '\t')
liver <- read.delim(file = file.path(path, 'v4_liver.tsv'),
                    sep = '\t')
lungs <- read.delim(file = file.path(path, 'v3_lungs.tsv'),
                    sep = '\t')


############################################################
# 2. Fetch ontology IDs
############################################################

# Keep only BP columns, and rows where ontology_id starts with "GO:"
filter_go_data <- function(df) {
  df %>%
    dplyr::select(1:8) %>%
    filter(BP_source == 'GO')
  
}

# Apply to all organs
data_go_brain <- filter_go_data(brain)

data_go_kidney <- filter_go_data(kidney)

data_go_liver <- filter_go_data(liver)

data_go_lungs <- filter_go_data(lungs)

# for simplicity add $ontology_id
make_ontology_id <- function(df, uri_col = "BP_uri") {
  # Extract the last segment of the URI (e.g. "GO_0038166" or "GO:0038166")
  id <- sub(".*/", "", df[[uri_col]])
  
  # Convert the first underscore to a colon, so "GO_0038166" becomes "GO:0038166"
  df$ontology_id <- sub("_", ":", id)
  
  df$BP_id <- NULL
  df$BP_source <- NULL
  return(df)
}

data_go_brain <- make_ontology_id(data_go_brain)
data_go_kidney <- make_ontology_id(data_go_kidney)
data_go_liver <- make_ontology_id(data_go_liver)
data_go_lungs <- make_ontology_id(data_go_lungs)

# save the files containing only GO ids
getwd()
write_tsv(data_go_brain, 'go_ids_brain.tsv')
write_tsv(data_go_kidney, 'go_ids_kidney.tsv')
write_tsv(data_go_liver, 'go_ids_liver.tsv')
write_tsv(data_go_lungs, 'go_ids_lungs.tsv')


############################################################
# 3. Retrieve gene sets from GO API
############################################################


# Loading libraries
library(httr)
library(jsonlite)
library(dplyr)
library(purrr)

############ Functions to get direct + indirect annotations ############ 

# Function to get human genes from GO API with retry logic
get_human_genes <- function(go_id, max_retries = 3, base_wait = 2) {
  url <- paste0("https://api.geneontology.org/api/bioentity/function/", go_id, "/genes")
  print(paste("Fetching:", url))
  
  params <- list(
    taxon = "NCBITaxon:9606",  # Homo sapiens
    rows = 10000,
    evidence = "IDA,IMP,IGI,IPI,IEP,EXP"  # Experimental evidence only
  )
  
  for (attempt in 1:max_retries) {
    tryCatch({
      cat(paste("  Attempt", attempt, "of", max_retries, "\n"))
      
      response <- GET(url, query = params, timeout(30))  # 30 second timeout
      status <- status_code(response)
      
      if (status == 200) {
        data <- fromJSON(content(response, "text"))
        if (!is.null(data$associations) && length(data$associations) > 0) {
          genes <- unique(data$associations$subject$label)
          cat(paste("  ✓ Success: Found", length(genes), "genes\n"))
          return(genes)
        } else {
          cat("  ✓ Success: No genes found for this GO term\n")
          return(character(0))
        }
      } else {
        # Handle different HTTP status codes
        if (status == 429) {
          cat(paste("  ⚠ Rate limit exceeded (HTTP 429). Waiting longer...\n"))
          wait_time <- base_wait * (2^attempt)  # Exponential backoff
        } else if (status >= 500) {
          cat(paste("  ⚠ Server error (HTTP", status, "). Retrying...\n"))
          wait_time <- base_wait * attempt
        } else if (status == 404) {
          cat(paste("  ⚠ GO term not found (HTTP 404). This may be expected.\n"))
          return(character(0))  # Don't retry for 404s
        } else {
          cat(paste("  ⚠ HTTP error", status, ". Retrying...\n"))
          wait_time <- base_wait * attempt
        }
        
        if (attempt < max_retries) {
          cat(paste("  Waiting", wait_time, "seconds before retry...\n"))
          Sys.sleep(wait_time)
        }
      }
      
    }, error = function(e) {
      if (grepl("Timeout", e$message) || grepl("timeout", e$message)) {
        cat(paste("  ⚠ Request timeout:", e$message, "\n"))
      } else if (grepl("resolve host", e$message) || grepl("network", e$message)) {
        cat(paste("  ⚠ Network error:", e$message, "\n"))
      } else {
        cat(paste("  ⚠ Unexpected error:", e$message, "\n"))
      }
      
      if (attempt < max_retries) {
        wait_time <- base_wait * attempt
        cat(paste("  Waiting", wait_time, "seconds before retry...\n"))
        Sys.sleep(wait_time)
      }
    })
  }
  
  # If we get here, all retries failed
  cat(paste("  ✗ FAILED: All", max_retries, "attempts failed for GO term:", go_id, "\n"))
  warning(paste("Failed to fetch genes for GO term:", go_id, "after", max_retries, "attempts"))
  return(character(0))
}

# Function to add gene sets to KE data with improved error handling
add_gene_sets <- function(ke_data, max_retries = 3, base_wait = 2) {
  
  ke_data$GOAPI_gene_set <- ""
  ke_data$GOAPI_total_genes <- 0
  failed_terms <- c()
  
  cat("Starting gene set enrichment for", nrow(ke_data), "entries...\n")
  cat("=================================================\n")
  
  for (i in 1:nrow(ke_data)) {
    cat(paste("\n[", i, "/", nrow(ke_data), "] Processing KEID:", ke_data$keid[i], "\n"))
    
    # Handle multiple GO terms (separated by semicolon or comma)
    go_terms <- trimws(ke_data$ontology_id[i])
    
    all_genes <- c()
    
    for (go_term in go_terms) {
      if (grepl("^GO:", go_term)) {
        cat(paste("Processing GO term:", go_term, "\n"))
        
        genes <- get_human_genes(go_term, max_retries, base_wait)
        
        if (length(genes) == 0) {
          # Check if this was a real failure or just no genes found
          # We'll track terms that completely failed
          failed_terms <- c(failed_terms, paste(ke_data$keid[i], go_term, sep = ":"))
        }
        
        all_genes <- c(all_genes, genes)
        
        # Rate limiting between requests
        Sys.sleep(1)
      } else {
        cat(paste("  Skipping non-GO term:", go_term, "\n"))
      }
    }
    
    # Remove duplicates
    unique_genes <- unique(all_genes[!is.na(all_genes) & all_genes != ""])
    
    ke_data$GOAPI_gene_set[i] <- paste(unique_genes, collapse = ",")
    ke_data$GOAPI_total_genes[i] <- length(unique_genes)
    
    cat(paste("  → Final result:", length(unique_genes), "unique genes\n"))
    
    # Progress indicator
    if (i %% 10 == 0) {
      cat(paste("\n*** Progress: Completed", i, "of", nrow(ke_data), "entries ***\n"))
    }
  }
  
  cat("\n=================================================\n")
  cat("Gene set enrichment completed!\n")
  
  # Report any failures
  if (length(failed_terms) > 0) {
    cat("\n⚠ WARNING: The following KEID:GO_term combinations had issues:\n")
    for (failed_term in failed_terms) {
      cat(paste("  -", failed_term, "\n"))
    }
    cat("You may want to re-run these manually or check if the GO terms are valid.\n")
  }
  
  # Summary statistics
  total_genes_found <- sum(ke_data$GOAPI_total_genes)
  entries_with_genes <- sum(ke_data$GOAPI_total_genes > 0)
  cat(paste("\nSUMMARY:\n"))
  cat(paste("  Total unique genes found:", total_genes_found, "\n"))
  cat(paste("  Entries with genes:", entries_with_genes, "out of", nrow(ke_data), "\n"))
  cat(paste("  Success rate:", round(entries_with_genes/nrow(ke_data)*100, 1), "%\n"))
  
  return(ke_data)
}


# We can adjust max_retries and base_wait as needed
all_kidney <- add_gene_sets(data_go_kidney, max_retries = 3, base_wait = 2)
all_brain <- add_gene_sets(data_go_brain, max_retries = 3, base_wait = 2)
all_liver <- add_gene_sets(data_go_liver, max_retries = 3, base_wait = 2)
all_lungs <- add_gene_sets(data_go_lungs, max_retries = 3, base_wait = 2)

# save files 
write_tsv(all_brain, 'all_brain.tsv')
write_tsv(all_lungs, 'all_lungs.tsv')
write_tsv(all_liver, 'all_liver.tsv')
write_tsv(all_kidney, 'all_kidney.tsv')

############################################################
# 4. Convert HGNC symbols to Ensembl IDs
############################################################
library(biomaRt)
library(stringr)

# Initialize Ensembl connection
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# FUNCTION to convert HGNC symbols to Ensembl IDs + add additional columns
get_ensembl_ids <- function(gene_symbols, mart) {
  if (is.na(gene_symbols) || gene_symbols == "") return(data.frame(
    ensembl_ids = NA, ensembl_API_count = 0, unmapped_ids = NA, num_unmapped_ids = 0
  ))
  
  symbols <- unlist(str_split(gene_symbols, ",\\s*"))
  results <- getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id", "gene_biotype"),
    filters = "hgnc_symbol",
    values = symbols,
    mart = mart
  )
  
  # keep only protein coding
  results <- results[results$gene_biotype == "protein_coding", ]
  
  unique_ids <- unique(results$ensembl_gene_id)
  mapped_symbols <- unique(results$hgnc_symbol)
  
  unmapped <- setdiff(symbols, mapped_symbols)
  
  data.frame(
    ensembl_ids = ifelse(length(unique_ids) > 0, paste(unique_ids, collapse = ", "), NA),
    ensembl_API_count = length(unique_ids),
    unmapped_ids = ifelse(length(unmapped) > 0, paste(unmapped, collapse = ", "), NA),
    num_unmapped_ids = length(unmapped)
  )
}

# Wrapper to apply Ensembl mapping to any organ dataframe
add_ensembl_columns <- function(df, mart, gene_col = "GOAPI_gene_set") {
  
  # Apply mapping first
  mapped <- lapply(df[[gene_col]], function(x) get_ensembl_ids(x, mart = mart))
  
  # Convert list of named vectors to a dataframe
  mapped_df <- do.call(rbind, lapply(mapped, as.data.frame))
  rownames(mapped_df) <- NULL
  
  # Bind back to original df
  df <- bind_cols(df, mapped_df)
  
  return(df)
}

# make a named list of organ dataframes
organ_list <- list(
  brain  = all_brain,
  kidney = all_kidney,
  liver  = all_liver,
  lungs  = all_lungs
)

# apply the add_ensembl_columns function to each organ df
organ_list <- lapply(organ_list, add_ensembl_columns, mart = mart)

############################################################
# 5. Export results
############################################################

# Assign back to individual objects with "g_" prefix
list2env(setNames(organ_list, paste0("g_", names(organ_list))), 
         envir = .GlobalEnv)

# save for each organ
write_tsv(g_kidney, 'KE_gene_kidney.tsv')
write_tsv(g_brain, 'KE_gene_brain.tsv')
write_tsv(g_liver, 'KE_gene_liver.tsv')
write_tsv(g_lungs, 'KE_gene_lungs.tsv')

