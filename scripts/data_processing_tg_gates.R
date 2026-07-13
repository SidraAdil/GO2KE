############################################################
# 1. Load libraries and Open TG-GATEs transcriptomics data
############################################################
### rds object from Open TG-GATES (PHH) ###

setwd('')
library(ToxicoGx)
library(dplyr)
library(tidyr)

# --- Load data ---
tg_data <- readRDS("ToxicoSet_TGGATES_human.rds")

############################################################
# 2. Select and process pesticide compounds for the liver case study
############################################################

# --- Define compounds of interest ---
target_compounds <- c("Hexachlorobenzene", "Rotenone", "Cycloheximide", "2,4-Dinitrophenol")

# --- Extract RNA expression data ---
expr <- molecularProfiles(tg_data, mDataType = "rna")
pheno_df <- as.data.frame(phenoInfo(tg_data, mDataType = "rna"))


# (Open TG-GATEs data are already log-normalized)
range(expr)


# Subset expression data for selected pesticides
samples_of_interest <- pheno_df[pheno_df$drugid %in% target_compounds, ]

# Keep expression and pheno data SEPARATE
expr_subset <- expr[, rownames(samples_of_interest)]
expr_df <- as.data.frame(t(expr_subset))  # Samples as rows, genes as columns

# Verify structure
cat("Expression data dimensions:", dim(expr_df), "\n")
cat("Pheno data dimensions:", dim(samples_of_interest), "\n")
cat("Dose levels present:", paste(unique(samples_of_interest$dose_level), collapse = ", "), "\n\n")

# --- Differential expression function ---
library(limma)
library(dplyr)
library(readr)

# --- Differential expression function using limma ---
# Define differential expression analysis function

perform_de_analysis_limma <- function(expr_data, pheno_data, compound_name) {
  compound_samples <- pheno_data[pheno_data$drugid == compound_name, ]
  available_samples <- intersect(compound_samples$samplename, rownames(expr_data))
  
  if (length(available_samples) == 0) {
    warning(paste("No samples found for", compound_name))
    return(NULL)
  }
  
  compound_samples <- compound_samples[compound_samples$samplename %in% available_samples, ]
  pheno_cols <- colnames(pheno_data)
  gene_cols <- setdiff(colnames(expr_data), pheno_cols)
  
  results_list <- list()
  
  time_points <- unique(compound_samples$duration)
  dose_levels <- c("Low", "Middle", "High")
  
  for (time_point in time_points) {
    for (dose_level in dose_levels) {
      control_samples <- compound_samples[compound_samples$duration == time_point & 
                                            compound_samples$dose_level == "Control", ]
      treated_samples <- compound_samples[compound_samples$duration == time_point & 
                                            compound_samples$dose_level == dose_level, ]
      
      if (nrow(control_samples) == 0 | nrow(treated_samples) == 0) next
      
      # Combine samples
      all_samples <- rbind(control_samples, treated_samples)
      
      # Extract expression matrix (genes in rows, samples in columns)
      expr_matrix <- t(expr_data[all_samples$samplename, gene_cols, drop = FALSE])
      
      # Create design matrix
      group <- factor(c(rep("Control", nrow(control_samples)), 
                        rep("Treated", nrow(treated_samples))))
      design <- model.matrix(~0 + group)
      colnames(design) <- c("Control", "Treated")
      
      # Fit linear model
      fit <- lmFit(expr_matrix, design)
      
      # Create contrast matrix (Treated vs Control)
      contrast.matrix <- makeContrasts(Treated - Control, levels = design)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      
      # Extract results
      results <- topTable(fit2, number = Inf, sort.by = "none")
      
      # Format output
      gene_results <- data.frame(
        gene_id = rownames(results),
        time_point = time_point,
        dose_level = dose_level,
        concentration = unique(treated_samples$concentration),
        log2FC = results$logFC,
        p_value = results$P.Value,
        adj_p_value = results$adj.P.Val,
        t_statistic = results$t,
        B_statistic = results$B,
        stringsAsFactors = FALSE
      )
      
      results_list[[paste(compound_name, time_point, dose_level, sep = "_")]] <- gene_results
    }
  }
  
  all_results <- do.call(rbind, results_list)
  rownames(all_results) <- NULL
  return(all_results)
}

# --- Perform analysis for all compounds ---
all_compound_results <- list()
for (compound in target_compounds) {
  cat("Processing:", compound, "\n")
  compound_pheno <- samples_of_interest[samples_of_interest$drugid == compound, ]
  results <- perform_de_analysis_limma(expr_df, compound_pheno, compound)
  all_compound_results[[compound]] <- results
}

# --- Save results ---
if (!dir.exists("corrected_data")) {
  dir.create("corrected_data")
}

for (compound in target_compounds) {
  results <- all_compound_results[[compound]]
  if (is.null(results)) next
  
  safe_name <- gsub(" ", "_", compound)
  
  for (dose in c("Low", "Middle", "High")) {
    for (time in unique(results$time_point)) {
      filtered <- results %>%
        filter(dose_level == dose, time_point == time) %>%
        mutate(gene_id = sub("_at$", "", gene_id))
      
      if (nrow(filtered) > 0) {
        filename <- file.path("corrected_data", paste0(safe_name, "_", tolower(dose), "_vs_cont_", time, ".tsv"))
        write_tsv(filtered, filename)
      }
    }
  }
}

cat("Analysis complete for all compounds!\n")

# --- Check results ---
# Example: Check one compound
example_result <- all_compound_results[[1]]
cat("\nNumber of genes with adj_p < 0.05:", sum(example_result$adj_p_value < 0.05, na.rm = TRUE), "\n")
cat("Number of genes with p < 0.05:", sum(example_result$p_value < 0.05, na.rm = TRUE), "\n")
summary(example_result$log2FC)

##################################

# see the dose and time info of compounds
pheno_info <- pheno_df %>%
  select(7:9, 11, 13, 14) %>%
  filter(drugid %in% target_compounds) 

# Export phenotype information
write_tsv(pheno_info, 'pheno_info.tsv')
















