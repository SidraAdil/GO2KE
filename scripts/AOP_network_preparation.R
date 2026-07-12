############################################################
# AOP Network Preparation
#
# This script prepares the liver AOP network for downstream
# transcriptomics-based analysis. It:
#   (1) constructs the network nodes and edges,
#   (2) assigns KE-associated gene sets,
#   (3) generates State I and State III node tables, and
#   (4) exports the processed network files.
############################################################

############################################################
# 1. Set Working Directory and Load Required Libraries
############################################################
setwd('')

# call libraries
library(readr)
library(tidyverse)
library(dplyr)


############################################################
# 2. Read Input Files
############################################################

list.files(pattern = '.tsv$')

aops <- read_tsv("selected_liver_aops.tsv")

kes <- read_tsv('kes.tsv')

ker <- read_tsv('ker.tsv')

############################################################
# 3. Identify Multi-role Key Events
#
# Some KEs appear in multiple biological roles (e.g., both
# Molecular Initiating Event and Key Event). These are
# identified and standardized to ensure each KE has a single
# node representation in the network.
############################################################

# Identify duplicated KE entries resulting from multiple roles

# see duplicates
table(kes$keLabel)[table(kes$keLabel) > 1]

# see full rows of duplicates
duplicated_kes <- kes %>%
  group_by(keLabel) %>%
  filter(n() > 1) %>%
  arrange(keLabel)

##### Identify multi-role categories #####
multi_role_types <- kes %>%
  group_by(keLabel) %>%
  summarise(types = paste(sort(unique(type)), collapse = "/")) %>%
  filter(grepl("/", types))   # only those with >1 type


############################################################
# 4. Prepare Network Nodes
############################################################

# preparing nodes 
nodes <- data.frame(
  id = kes$keLabel,
  name = kes$keTitle,
  type = kes$type,
  stringsAsFactors = FALSE
)

# Standardize multi-role KEs to a single node type
for(i in 1:nrow(multi_role_types)) {
  if(multi_role_types$types[i] == "KE/MIE") {
    nodes$type[nodes$id == multi_role_types$keLabel[i]] <- "MIE"
  } else if(multi_role_types$types[i] == "AO/KE") {
    nodes$type[nodes$id == multi_role_types$keLabel[i]] <- "AO"
  }
}


# Remove duplicate node entries
nodes <- unique(nodes)

############################################################
# 5. Prepare Network Edges
############################################################

edges <- data.frame(
  kerLabel = ker$kerLabel,
  source = ker$upstreamLabel,
  target = ker$downstreamLabel,
  interaction = "leads_to",
  stringsAsFactors = FALSE
)

# Remove duplicate edge entries
edges <- unique(edges)

############################################################
# 6. Export Base Network Files
############################################################

# save nodes and edges 
write_tsv(nodes, 'nodes.tsv')
write_tsv(edges, 'edges.tsv')


###########################################################

############################################################
# 7. Generate State I Node Table
#    (Pre-annotated Key Events)
############################################################

# prepare state I nodes (pre-annotated data only)
# Integrate pre-existing KE gene sets into the node table

pre_anno_data <- read_tsv('already_anno_liver_geneset.tsv')
unique(pre_anno_data$keid)

# Inspect input columns
names(pre_anno_data)

library(dplyr)

state_I_nodes <- nodes %>%
  left_join(
    pre_anno_data %>%
      dplyr::select(keid, ensembl_ids, ensembl_API_count, BP_title, ontology_id),
    by = c("id" = "keid") 
  ) %>%
  mutate(
    anno_status = ifelse(!is.na(ensembl_ids), 'pre_annotated', NA)
  )

# Rename columns for consistency across network states
names(state_I_nodes)
state_I_nodes <- state_I_nodes %>%
  rename(KE_geneset = ensembl_ids,
         KE_geneset_count = ensembl_API_count) 

# Export State I node table
write_tsv(state_I_nodes, 'state_I_nodes.tsv')

############################################################
# 8. Generate State III Node Table
#    (Improved Pre-annotated + Newly Annotated Key Events)
############################################################

# Read newly annotated KE gene sets
newly_anno <- read_tsv('KE_gene_liver.tsv')

# Read improved pre-annotated KE gene sets
improved_liver <- read_tsv('improved_anno_liver_geneset.tsv')

# Inspect input columns
names(newly_anno)
names(improved_liver)

# Merge newly annotated KE gene sets into the node table
# add $ensembl_ids, $ensembl_API_count based on $keid in nodes and name in state III
library(dplyr)
state_III_nodes <- nodes %>%
  left_join(
    newly_anno %>%
      dplyr::select(keid, ensembl_ids, ensembl_API_count, ontology_id, BP_title),
    by = c("id" = "keid")
  ) %>%
  mutate(
    anno_status = ifelse(!is.na(ensembl_ids), 'newly_annotated', NA)
  )


# Verify successful integration of newly annotated KEs
length(newly_anno$keid %in% state_III_nodes$id)

# Merge improved pre-annotated KE gene sets
names(state_III_nodes)

state_III_nodes <- state_III_nodes %>%
  left_join(
    improved_liver %>%
      dplyr::select(keid, ensembl_ids, ensembl_API_count, ontology_id, BP_title),
    by = c("id" = "keid")
  ) %>%
  mutate(
    ensembl_ids = ifelse(is.na(ensembl_ids.x), ensembl_ids.y, ensembl_ids.x),
    ensembl_API_count = ifelse(is.na(ensembl_API_count.x), ensembl_API_count.y, ensembl_API_count.x),
    ontology_id = ifelse(is.na(ontology_id.x), ontology_id.y, ontology_id.x),
    BP_title = ifelse(is.na(BP_title.x), BP_title.y, BP_title.x),
    #anno_status = ifelse(is.na(anno_status), 'improved_pre_annotated', anno_status)
    
  ) %>%
  dplyr::select(-ends_with(".x"), -ends_with(".y"))   # drop helper columns


# Verify successful integration of improved annotations
length(improved_liver$keid %in% state_III_nodes$id)

# Assign annotation status to each Key Event
# Annotation categories:
#   pre_annotated
#   improved_pre_annotated
#   newly_annotated
#   NA (no associated gene set)

state_III_nodes <- state_III_nodes %>%
  left_join(
    pre_anno_data %>%
      dplyr::select(keid, ontology_id),
    by = c("id" = "keid"),
    suffix = c("", "_pre")
  ) %>%
  mutate(
    anno_status = case_when(
      !is.na(anno_status) ~ anno_status,  # keep 'newly_annotated'
      is.na(ontology_id) ~ NA_character_,  # no annotation at all
      !is.na(ontology_id_pre) & ontology_id == ontology_id_pre ~ "pre_annotated",  # exists in pre_annotated data AND matches
      !is.na(ontology_id_pre) & ontology_id != ontology_id_pre ~ "improved_pre_annotated",  # exists in pre_annotated data BUT differs
      is.na(ontology_id_pre) ~ "improved_pre_annotated",  # NOT in pre_annotated data = new improvement
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(-ontology_id_pre)

# Summarize annotation status
table(state_III_nodes$anno_status, useNA = "ifany")

length(unique(state_III_nodes$id))

############################################################
# 9. Summarize Annotation Coverage
############################################################

# Calculate the proportion of KEs associated with gene sets
percent_annotated <- sum(!is.na(state_III_nodes$ensembl_ids)) / nrow(state_III_nodes) * 100
percent_annotated

# rename columns
names(state_III_nodes)
state_III_nodes <- state_III_nodes %>%
  rename(KE_geneset = ensembl_ids,
         KE_geneset_count = ensembl_API_count) 

############################################################
# 10. Resolve Multiple GO Annotations per Key Event
############################################################

# see multi kec (key event component) in both files
names(state_I_nodes)
state_I_nodes %>%
  group_by(id) %>%
  filter(n()>1)

# No multiple GO annotations detected in State I


# Resolve multiple GO annotations in State III

# Identify duplicated annotation records
state_III_nodes %>%
  group_by(across(everything())) %>%
  filter(n() > 1) %>%
  arrange(id)

# Remove duplicated annotation records
refined_state_III_nodes <- state_III_nodes %>%
  distinct()

multi_kec_state_III <- refined_state_III_nodes %>%
  group_by(id) %>%
  filter(n() > 1)

# remove too generic # KE 295 catalytic activity, KE 457 signaling
refined_state_III_nodes <- refined_state_III_nodes %>%
  filter(!(id == "KE 295" & BP_title == "catalytic activity")) %>%
  filter(!(id == "KE 457" & BP_title == "signaling"))

# Recheck remaining multi-annotation KEs
multi_kec_state_III <- refined_state_III_nodes %>%
  group_by(id) %>%
  filter(n() > 1)

tab <- multi_kec_state_III %>%
  select(1,2,4,6:8)

names(refined_state_III_nodes)
# take union of $KE_geneset to handle multi KEC (multiple GO annotations per KE)
refined <- refined_state_III_nodes %>%
  group_by(id, name, type) %>%
  summarise(
    KE_geneset = paste(unique(unlist(strsplit(KE_geneset, ", "))), collapse = ", "),
    KE_geneset_count = sum(KE_geneset_count),
    BP_title = paste(unique(BP_title), collapse = ", "),
    ontology_id = paste(unique(ontology_id) , collapse = ', '),
    anno_status = paste(unique(anno_status), collapse = ','),
    .groups = "drop"
  ) 

names(refined)

# Convert character "NA" values to missing values
refined$KE_geneset[refined$KE_geneset == "NA"] <- NA
refined$BP_title[refined$BP_title == "NA"] <- NA
refined$ontology_id[refined$ontology_id == "NA"] <- NA

############################################################
# 11. Export Final State III Node Table
############################################################

write_tsv(refined, 'state_III_nodes.tsv')

############################################################
# 13. Compare Annotation Coverage Across Network States
############################################################

# Calculate annotation coverage for State I 
names(state_III_nodes)
percent_annotated <- sum(!is.na(state_III_nodes$KE_geneset)) / nrow(state_III_nodes) * 100
percent_annotated

# Calculate annotation coverage for State I
percent_annotated <- sum(!is.na(state_I_nodes$KE_geneset)) / nrow(state_I_nodes) * 100
percent_annotated



