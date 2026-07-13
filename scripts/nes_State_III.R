############################################################
# KE perturbation analysis using ssGSEA
#
# This script computes KE perturbation scores for the enhanced
# liver AOP network (State III) using differential expression
# results from the Open TG-GATEs primary human hepatocyte dataset.
#
# Gene ranking:
#   log2FC × -log10(p-value)
#
# Gene weighting:
#   |log2FC|^0.25
#
# Output:
#   - Compound-specific KE perturbation tables
#   - Combined KE perturbation results
############################################################

############################################################
# Load libraries and input data
############################################################

# Set working directory
setwd("")

library(dplyr)
library(tidyr)
library(readr)

list.files(pattern = '.tsv$')
path = './state_III_nodes.tsv'
# Read State III node table and retain only KEs with associated gene sets
nodes <- read_tsv(path)

nodes <- nodes %>%
  filter(!is.na(KE_geneset))

###################################################
# KE perturbation scores
###################################################

# Directory containing differential expression results

data_path <- "./corrected_data"

# Pesticide compounds included in the case study

all_compounds <- c("Hexachlorobenzene", "Rotenone", "Cycloheximide", "2,4-Dinitrophenol")

############################################################
# Function: Calculate enrichment score (ES)
############################################################

# Computes the enrichment score (ES) for a single KE gene set.
# Genes are ranked using:
#     log2FC × -log10(p-value)
# while enrichment weighting is based on:
#     |log2FC|^alpha

calculate_enrichment_score_pvalue_ranked <- function(ranking_scores, log2fc_vector, 
                                                     gene_set_genes, alpha = 0.25) {
  N <- length(ranking_scores)
  N_h <- length(gene_set_genes)
  P_miss_constant <- 1 / (N - N_h)
  N_R <- sum(abs(log2fc_vector[gene_set_genes])^alpha)
  
  # Sort genes by ranking score (log2FC × -log10(p-value))
  sorted_genes <- names(sort(ranking_scores, decreasing = TRUE))
  gene_set_lookup <- sorted_genes %in% gene_set_genes
  
  running_sum <- 0
  max_deviation <- 0
  genes_found <- 0
  n_genes <- length(sorted_genes)
  
  for (i in seq_len(n_genes)) {
    if (gene_set_lookup[i]) {
      gene <- sorted_genes[i]
      # Still use log2FC for weighting (not ranking score)
      P_hit <- (abs(log2fc_vector[gene])^alpha) / N_R
      running_sum <- running_sum + P_hit
      genes_found <- genes_found + 1
    } else {
      running_sum <- running_sum - P_miss_constant
    }
    
    if (abs(running_sum) > abs(max_deviation)) {
      max_deviation <- running_sum
    }
    
    if (genes_found == N_h) {
      remaining <- n_genes - i
      potential_final <- running_sum - (remaining * P_miss_constant)
      if (abs(potential_final) <= abs(max_deviation)) break
    }
  }
  
  return(max_deviation)
}

############################################################
# Function: Calculate normalized enrichment score (NES)
############################################################

# Calculates KE perturbation scores using single-sample GSEA (ssGSEA)
# and estimates statistical significance through permutation testing.

calculate_ssGSEA_NES_pvalue_ranked <- function(de_results_list, gene_sets, alpha = 0.25, 
                                               nperm = 1000, seed = 123) {
  set.seed(seed)
  
  n_conditions <- length(de_results_list)
  n_gene_sets <- length(gene_sets)
  
  nes_matrix <- matrix(NA, nrow = n_gene_sets, ncol = n_conditions)
  pval_matrix <- matrix(NA, nrow = n_gene_sets, ncol = n_conditions)
  rownames(nes_matrix) <- rownames(pval_matrix) <- names(gene_sets)
  colnames(nes_matrix) <- colnames(pval_matrix) <- names(de_results_list)
  
  cat("Calculating ssGSEA NES with p-value ranked genes\n")
  cat("Gene sets:", n_gene_sets, "| Conditions:", n_conditions, "\n")
  cat("Seed:", seed, "| Permutations:", nperm, "\n")
  
  all_genes <- de_results_list[[1]]$gene_id[!is.na(de_results_list[[1]]$log2FC)]
  gene_overlaps <- lapply(gene_sets, function(gs) intersect(gs, all_genes))
  
  for (cond_idx in 1:n_conditions) {
    condition_name <- names(de_results_list)[cond_idx]
    cat("\nCondition:", condition_name, "\n")
    
    de_data <- de_results_list[[cond_idx]]
    
    # Extract log2FC and p-values
    log2fc_vector <- setNames(de_data$log2FC, de_data$gene_id)
    p_values <- setNames(de_data$p_value, de_data$gene_id)
    
    # Remove NAs
    valid_genes <- !is.na(log2fc_vector) & !is.na(p_values)
    log2fc_vector <- log2fc_vector[valid_genes]
    p_values <- p_values[valid_genes]
    
    # Calculate ranking scores: log2FC × -log10(p_value)
    # Add small value to avoid log10(0)
    p_values_safe <- pmax(p_values, 1e-300)
    ranking_scores <- log2fc_vector * (-log10(p_values_safe))
    
    cat("  Genes with valid data:", length(ranking_scores), "\n")
    
    # Generate null distributions for each geneset size
    gene_set_sizes <- sapply(gene_overlaps, length)
    unique_sizes <- unique(gene_set_sizes)
    null_distributions <- list()
    
    for (size in unique_sizes) {
      null_es <- replicate(nperm, {
        random_genes <- sample(names(ranking_scores), size)
        calculate_enrichment_score_pvalue_ranked(ranking_scores, log2fc_vector, 
                                                 random_genes, alpha)
      })
      null_distributions[[as.character(size)]] <- null_es
    }
    
    # Calculate ES for each gene set
    pb <- txtProgressBar(min = 0, max = n_gene_sets, style = 3)
    for (gs_idx in 1:n_gene_sets) {
      common_genes <- gene_overlaps[[gs_idx]]
      
      # Calculate observed ES
      es_observed <- calculate_enrichment_score_pvalue_ranked(ranking_scores, log2fc_vector,
                                                              common_genes, alpha)
      null_es <- null_distributions[[as.character(length(common_genes))]]
      
      # Calculate NES and p-value
      if (es_observed >= 0) {
        pos_null <- null_es[null_es >= 0]
        nes <- if (length(pos_null) > 0) es_observed / mean(pos_null) else 0
        pval <- sum(null_es >= es_observed) / nperm
      } else {
        neg_null <- null_es[null_es < 0]
        nes <- if (length(neg_null) > 0) es_observed / abs(mean(neg_null)) else 0
        pval <- sum(null_es <= es_observed) / nperm
      }
      
      if (!is.finite(nes)) nes <- 0
      if (pval == 0) pval <- 1 / nperm
      
      nes_matrix[gs_idx, cond_idx] <- nes
      pval_matrix[gs_idx, cond_idx] <- pval
      setTxtProgressBar(pb, gs_idx)
    }
    close(pb)
  }
  
  cat("\n\nComplete!\n")
  return(list(nes = nes_matrix, pval = pval_matrix, gene_overlaps = gene_overlaps))
}

############################################################
# Function: Run KE perturbation analysis
############################################################

# Performs KE perturbation analysis for one compound across all
# exposure conditions.
#
# Workflow:
#   1. Build KE gene sets
#   2. Filter gene sets (10–800 genes)
#   3. Compute ES, NES and permutation p-values
#   4. Determine significantly perturbed KEs
#   5. Export KE-level results

run_KE_analysis <- function(nodes, de_data_list, compound_name = "Compound", seed = 123) {
  cat("=== KE ssGSEA Analysis (P-value Ranked) ===\n")
  cat("Compound:", compound_name, "\n")
  cat("Number of conditions:", length(de_data_list), "\n\n")
  
  # Count total nodes
  cat("Total nodes:", nrow(nodes), "\n")
  
  # Build gene sets from nodes with genesets
  gene_sets <- list()
  for (i in 1:nrow(nodes)) {
    if (is.na(nodes$KE_geneset[i])) {
      next
    }
    
    genes <- trimws(unlist(strsplit(nodes$KE_geneset[i], ",")))
    genes <- genes[genes != ""]
    if (length(genes) > 0) {
      gene_sets[[nodes$id[i]]] <- genes
    }
  }
  
  cat("Nodes with genesets:", length(gene_sets), "\n\n")
  
  # Calculate overlaps with DE data
  all_genes <- de_data_list[[1]]$gene_id[!is.na(de_data_list[[1]]$log2FC)]
  gene_overlaps <- lapply(gene_sets, function(gs) intersect(gs, all_genes))
  
  # Calculate overlap sizes
  overlap_sizes <- sapply(gene_overlaps, length)
  
  # Apply filtering (10-800 genes)
  too_small <- overlap_sizes < 10
  too_large <- overlap_sizes > 800
  valid_ke <- overlap_sizes >= 10 & overlap_sizes <= 800
  
  cat("Geneset filtering (10-800 genes):\n")
  cat("Included:", sum(valid_ke), "\n")
  cat("Excluded (too small < 10):", sum(too_small), "\n")
  cat("Excluded (too large > 800):", sum(too_large), "\n\n")
  
  # Filter gene sets and overlaps
  gene_sets <- gene_sets[valid_ke]
  gene_overlaps <- gene_overlaps[valid_ke]
  
  # Run ssGSEA with p-value ranking
  results <- calculate_ssGSEA_NES_pvalue_ranked(de_data_list, gene_sets, seed = seed)
  
  # Build overlap info
  overlap_info <- data.frame(
    KE_ID = names(results$gene_overlaps),
    KE_genes_in_DE_count = sapply(results$gene_overlaps, length),
    KE_genes_in_DE = sapply(results$gene_overlaps, function(x) paste(x, collapse = ",")),
    stringsAsFactors = FALSE
  )
  
  overlap_info$KE_genes_NOT_in_DE_count <- sapply(names(results$gene_overlaps), function(ke_id) {
    all_genes_in_ke <- gene_sets[[ke_id]]
    genes_in_de <- results$gene_overlaps[[ke_id]]
    length(setdiff(all_genes_in_ke, genes_in_de))
  })
  
  overlap_info$KE_genes_NOT_in_DE <- sapply(names(results$gene_overlaps), function(ke_id) {
    all_genes_in_ke <- gene_sets[[ke_id]]
    genes_in_de <- results$gene_overlaps[[ke_id]]
    genes_not_in_de <- setdiff(all_genes_in_ke, genes_in_de)
    paste(genes_not_in_de, collapse = ",")
  })
  
  # Reshape results
  nes_long <- as.data.frame(results$nes) %>%
    tibble::rownames_to_column("KE_ID") %>%
    pivot_longer(cols = -KE_ID, names_to = "Condition", values_to = "NES")
  
  pval_long <- as.data.frame(results$pval) %>%
    tibble::rownames_to_column("KE_ID") %>%
    pivot_longer(cols = -KE_ID, names_to = "Condition", values_to = "p_value")
  
  combined <- nes_long %>%
    left_join(pval_long, by = c("KE_ID", "Condition"))
  
  nes_wide <- combined %>%
    select(KE_ID, Condition, NES) %>%
    pivot_wider(names_from = Condition, values_from = NES, names_prefix = "NES_")
  
  pval_wide <- combined %>%
    select(KE_ID, Condition, p_value) %>%
    pivot_wider(names_from = Condition, values_from = p_value, names_prefix = "pval_")
  
  node_cols <- intersect(c("id", "name", "type", "KE_geneset_count", "KE_geneset", "BP_title"), 
                         colnames(nodes))
  
  final_df <- nes_wide %>%
    left_join(pval_wide, by = "KE_ID") %>%
    left_join(nodes %>% select(all_of(node_cols)), by = c("KE_ID" = "id")) %>%
    left_join(overlap_info, by = "KE_ID") %>%
    select(KE_ID, name, type, KE_geneset_count, any_of("BP_title"),
           any_of("KE_geneset"), 
           KE_genes_in_DE_count, KE_genes_in_DE, 
           KE_genes_NOT_in_DE_count, KE_genes_NOT_in_DE, 
           everything())
  
  # Add significance columns
  for (cond in names(de_data_list)) {
    nes_col <- paste0("NES_", cond)
    pval_col <- paste0("pval_", cond)
    sig_col <- paste0("significant_", cond)
    
    final_df[[sig_col]] <- ifelse(
      !is.na(final_df[[pval_col]]) & final_df[[pval_col]] < 0.05 & abs(final_df[[nes_col]]) > 1.5,
      ifelse(final_df[[nes_col]] > 0, "Active+", "Active-"),
      "Not significant"
    )
  }
  
  return(final_df)
}

############################################################
# Process all pesticide compounds
############################################################

# Create results directory
nes_results_path <- file.path(getwd(), "nes_results")

if (!dir.exists(nes_results_path)) {
  dir.create(nes_results_path)
}

all_results <- list()

for (compound in all_compounds) {
  cat("\n\n========================================\n")
  cat("Processing compound:", compound, "\n")
  cat("========================================\n\n")
  
  files <- list.files(data_path, pattern = paste0("^", compound, "_.*\\.tsv$"), full.names = TRUE)
  
  if (length(files) == 0) {
    cat("No files found for", compound, "\n")
    next
  }
  
  cat("Found", length(files), "files\n")
  
  compound_data <- list()
  for (file in files) {
    condition <- gsub(paste0(data_path, "/", compound, "_"), "", file)
    condition <- gsub(".tsv", "", condition)
    condition_name <- paste0(compound, "_", condition)
    
    compound_data[[condition_name]] <- read_tsv(file, show_col_types = FALSE)
    cat("  Loaded:", condition_name, "\n")
  }
  
  results <- run_KE_analysis(
    nodes = nodes,
    de_data_list = compound_data,
    compound_name = compound,
    seed = 123
  )
  
  # Add compound column
  results$compound <- compound
  all_results[[compound]] <- results
  
  # Save individual file
  output_file <- file.path(nes_results_path, paste0("STATE3_", compound, "_NES_pvalue_ranked.tsv"))
  write_tsv(results, output_file)
  cat("\nSaved:", output_file, "\n")
  
  # Summary
  cat("\n=== Summary of Significant KEs ===\n")
  for (cond in names(compound_data)) {
    sig_col <- paste0("significant_", cond)
    n_active_plus <- sum(results[[sig_col]] == "Active+", na.rm = TRUE)
    n_active_minus <- sum(results[[sig_col]] == "Active-", na.rm = TRUE)
    cat(sprintf("%s: %d Active+, %d Active-\n", cond, n_active_plus, n_active_minus))
  }
}

# export combined results
combined_results <- bind_rows(all_results, .id = "compound_id")
combined_file <- file.path(nes_results_path, "STATE3_ALL_compounds_NES_pvalue_ranked.tsv")
write_tsv(combined_results, combined_file)
cat("\n\nSaved combined results:", combined_file, "\n")






