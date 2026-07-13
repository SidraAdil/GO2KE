############################################################
# AOP network visualization in Cytoscape
############################################################

# Set node size proportional to the associated KE gene set size
# (State III: improved pre-annotated + newly annotated KEs)

############################################################
# 1. Read network files
############################################################
nodes <- read_tsv('state_III_nodes.tsv')
edges <- read_tsv('edges.tsv')

unique(nodes$anno_status)

############################################################
# 2. Prepare node attributes
############################################################
# Replace missing gene set counts with zero for node size mapping

nodes <- nodes %>%
  mutate(KE_geneset_count = ifelse(is.na(KE_geneset_count), 0, KE_geneset_count))

############################################################
# 3. Create the Cytoscape network
############################################################
# Create a new network from node and edge tables
# Remove an existing network with the same name if present
network_name <- "Network_geneset_coverage_StateIII"

# Delete network if it already exists
if(network_name %in% getNetworkList()) {
  deleteNetwork(network_name)
}

# Create new network
createNetworkFromDataFrames(nodes, edges, title = network_name)

# IMPORTANT: Import the node type data into Cytoscape
loadTableData(nodes, data.key.column = "id", 
              table = "node", 
              table.key.column = "id")

# Step 5: Apply layout
layoutNetwork("force-directed")

# Then manually adjust spacing by scaling positions
scaleLayout(scaleFactor = 2.0, axis = 'x')

# Step 6: Create and apply visual style
style_name <- "AOP_Style_StateIII"
createVisualStyle(style_name)

# SET NODE SHAPE FIRST
setNodeShapeDefault("ELLIPSE", style.name = style_name)

# Lock node dimensions
lockNodeDimensions(TRUE, style.name = style_name)

# Node COLOR mapping by anno_status
setNodeColorDefault("#D3D3D3", style.name = style_name)  # grey default for NA

setNodeColorMapping(
  "anno_status",
  c("pre_annotated", "newly_annotated", "improved_pre_annotated", "pre_annotated,improved_pre_annotated"),
  c("#CCECFF", "#BCEE68", "#DDA0DD", "#DDA0DD"),  # blue, pink, light purple, light purple
  mapping.type = "d",
  style.name = style_name
)

# Node SIZE mapping by geneset count
setNodeSizeDefault(40, style.name = style_name)

setNodeSizeMapping(
  table.column = "KE_geneset_count",
  table.column.values = c(0, max(nodes$KE_geneset_count, na.rm = TRUE)),
  sizes = c(40, 120),
  mapping.type = "c",
  style.name = style_name
)

# Node BORDER COLOR mapping by type (only MIE=green, AO=red, no border for KE)
setNodeBorderColorDefault("#FFFFFF", style.name = style_name)  # white (invisible) default

setNodeBorderColorMapping(
  "type",
  c("MIE", "AO"),
  c("#228B22", "#DC143C"),  # green, red
  mapping.type = "d",
  style.name = style_name
)

# Set border width
setNodeBorderWidthDefault(5, style.name = style_name)

# Edge arrow
setEdgeTargetArrowShapeDefault("ARROW", style.name = style_name)

# Node labels
setNodeLabelMapping(
  table.column = "id",
  style.name = style_name
)

# Apply the style
setVisualStyle(style_name, network_name)

# Save the full Cytoscape session
saveSession("state_III_coverage.cys")

###### identifying hub KEs and KERs #########
# In R/RCy3
library(RCy3)

# Calculate network statistics
analyzeNetwork()

# Get node degree (number of connections)
node_stats <- getTableColumns('node', c('name', 'id', 'type', 'Degree', 'BetweennessCentrality'))

# Find hub nodes (high degree)
hubs <- node_stats[order(-node_stats$Degree), ]
head(hubs, 10)  # Top 10 hubs

# save top 10 hubs
library(dplyr)
top_10_hubs <- node_stats %>%
  arrange(desc(BetweennessCentrality)) %>%
  head(10)
getwd()
write_tsv(top_10_hubs, "top_10_hub_nodes.tsv")
