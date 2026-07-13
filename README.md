# Adding-gene-sets-to-Key-Events

This repository contains the R scripts used to systematically assign Gene Ontology (GO) biological process gene sets to Key Events (KEs) in the Adverse Outcome Pathway (AOP)-Wiki. The workflow supports transcriptomics-based quantification of AOPs by generating KE-associated gene sets. We demonstrate its application through a liver toxicity case study using the Open TG-GATEs primary human hepatocyte (PHH) transcriptomics dataset.

---

## Workflow

### 1. Ontology Mapping and GO-Based Gene Set Assignment

**Script**

**`scripts/GO_gene_mapping.R`**

**Description**

Retrieves genes associated with Gene Ontology (GO) biological process terms and assigns KE-associated gene sets to both pre-annotated and newly annotated Key Events. Gene identifiers are standardized to Ensembl IDs to facilitate downstream transcriptomics analyses.

---

### 2. Liver Case Study: AOP Network Construction

#### Script

**`scripts/AOP_selection.R`**

Retrieves liver-associated Adverse Outcome Pathways (AOPs) from the AOP-Wiki SPARQL endpoint using predefined liver-related keywords. Following manual curation of the retrieved AOPs, the script extracts the corresponding Key Events (KEs) and Key Event Relationships (KERs) required for network construction.

**Outputs**

- `selected_liver_aops.tsv`
- `kes.tsv`
- `ker.tsv`

---

#### Script

**`scripts/AOP_network_preparation.R`**

Constructs the liver AOP network by generating node and edge tables and integrating GO-based KE gene sets. The script prepares network files for both:

- **State I:** pre-annotated KEs
- **State III:** improved pre-annotated and newly annotated KEs

Annotation status and KE-associated gene sets are added to each network node.

**Outputs**

- `nodes.tsv`
- `edges.tsv`
- `state_I_nodes.tsv`
- `state_III_nodes.tsv`

---

#### Script

**`scripts/AOP_network_visualization.R`**

Generates the liver AOP network in Cytoscape using the prepared node and edge tables. Nodes are colored according to annotation status, node size is scaled by KE gene set size, and network statistics are computed to identify highly connected (hub) Key Events.

**Outputs**

- Cytoscape session (`.cys`)
- `top_10_hub_nodes.tsv`

---

### 3. Transcriptomics Data Processing

#### Script

**`scripts/data_processing_tg_gates.R`**

Processes transcriptomics data from the Open TG-GATEs primary human hepatocyte (PHH) dataset. The script:

- extracts the selected pesticide exposure experiments,
- performs differential gene expression analysis using **limma**,
- calculates log₂ fold changes and statistical significance,
- exports one expression file for each dose–time combination.

**Outputs**

- Differential expression tables (`.tsv`)
- `pheno_info.tsv`

---

### 4. KE Perturbation Analysis

#### Script

**`scripts/nes_State_III.R`**

Computes Key Event perturbation scores using single-sample Gene Set Enrichment Analysis (ssGSEA) and evaluates perturbations in the enhanced liver AOP network (State III).

---

## Input Data

The manually curated Gene Ontology biological process annotations for previously unannotated Key Events are available on Zenodo:

**DOI:** https://doi.org/10.5281/zenodo.21264953

The dataset includes GO biological process mappings for:

- Brain
- Kidney
- Liver
- Lung

---

## Requirements

- R (version 4.4 or later)
- Cytoscape
- Required R packages are listed within the individual scripts.
