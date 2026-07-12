# Adding-gene-sets-to-Key-Events

This repository contains the R scripts used to systematically assign Gene Ontology (GO) biological process gene sets to Key Events (KEs) in the Adverse Outcome Pathway (AOP)-Wiki. The workflow supports transcriptomics-based quantification of AOPs by generating KE-associated gene sets. We demonstrate its application through a liver toxicity case study using the Open TG-GATEs primary human hepatocyte (PHH) transcriptomics dataset.

## Repository Structure

### 1. GO-Based Gene Set Assignment

**Script**

- `GO_gene_mapping.R`

**Description**

Retrieves genes associated with GO biological process terms and assigns gene sets to both pre-annotated and newly annotated Key Events. Gene identifiers are standardized to Ensembl IDs to facilitate downstream transcriptomics analyses.

---

### 2. Liver AOP Network Construction

**Scripts**

- `AOP_selection.R`
- `AOP_network_generation.R`

**Description**

Retrieves liver-associated AOPs from AOP-Wiki using SPARQL queries, performs manual curation of AOPs, constructs the liver AOP network, and visualizes the resulting network in Cytoscape.

---

### 3. Case Study: Pesticide-Induced Liver Toxicity

**Scripts**

- `data_processing_tg_gates.R`
- `nes_State_III.R`

**Description**

Processes the Open TG-GATEs primary human hepatocyte transcriptomics dataset, computes KE perturbation scores using single-sample Gene Set Enrichment Analysis (ssGSEA), and evaluates pesticide-induced perturbations in the enhanced liver AOP network.

---

## Input Data

The manually curated Gene Ontology biological process annotations for previously unannotated Key Events are available on Zenodo:

**DOI:** https://doi.org/10.5281/zenodo.21264953

These annotations include organ-specific mappings for:

- Brain
- Kidney
- Liver
- Lung

---

## Requirements

- R (version 4.4 or later)
- Cytoscape
- Required R packages are specified within the individual scripts.


