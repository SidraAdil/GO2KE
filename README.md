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

### `AOP_selection.R`
Retrieves liver-associated Adverse Outcome Pathways (AOPs) from the AOP-Wiki SPARQL endpoint using predefined liver-related keywords. Following manual curation of the retrieved AOPs, the script extracts the corresponding Key Events (KEs) and Key Event Relationships (KERs) required for network construction.

**Outputs**
- `selected_liver_aops.tsv`
- `kes.tsv`
- `ker.tsv`

---

### `AOP_network_preparation.R`
Prepares the liver AOP network by generating node and edge tables and integrating GO-based KE gene sets. The script produces network files for both **State I** (pre-annotated KEs) and **State III** (improved pre-annotated + newly annotated KEs), including annotation status and KE gene set information.

**Outputs**
- `nodes.tsv`
- `edges.tsv`
- `state_I_nodes.tsv`
- `state_III_nodes.tsv`

---

### `AOP_network_visualization.R`
Generates and visualizes the liver AOP network in Cytoscape using the prepared node and edge tables. Nodes are styled according to annotation status, node size is scaled by KE gene set size, and network statistics are calculated to identify highly connected (hub) KEs.

**Outputs**
- Cytoscape session (`.cys`)
- `top_10_hub_nodes.tsv`

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


