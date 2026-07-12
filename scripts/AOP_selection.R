############################################################
# AOP_selection.R
#
# Purpose:
# Retrieve liver-associated AOPs from AOP-Wiki,
# manually curate the selected AOPs,
# and extract Key Events (KEs) and
# Key Event Relationships (KERs)
# for liver AOP network construction.
############################################################

############################################################
# 1. Retrieve liver-associated AOPs
############################################################

# setting working directory
setwd("")

# call libraries
library(httr)
library(tidyverse)
library(jsonlite)

# SPARQL endpoint (API)
# https://aopwiki.rdf.bigcat-bioinformatics.org/sparql


###### get all AOPs associated with Liver based on defined key words ###### 

# define query
query <- '
SELECT DISTINCT ?aop AS ?link 
       (CONCAT("AOP ", ?aopNumber) AS ?aopNumber) 
       ?title 
       ?aotitle 
       (GROUP_CONCAT(DISTINCT ?keid; SEPARATOR=", ") AS ?KEids) 
       (COUNT(DISTINCT ?keid) AS ?totalKEs)
WHERE {
    ?aop a aopo:AdverseOutcomePathway ;
         dc:title ?title ;
         aopo:has_adverse_outcome ?ao ;
         aopo:has_key_event ?ke .
    
    ?ao dc:title ?aotitle .
    ?ke rdfs:label ?keid .

    # Merge title and AO title
    BIND(CONCAT(LCASE(?title), " ", LCASE(?aotitle)) AS ?searchText)

    # AOP number
    BIND(STRAFTER(STR(?aop), "aop/") AS ?aopNumber)

    # Use regex for multiple keywords at once
    FILTER regex(?searchText, "liver|hepatic|hepato|steatosis|hepatitis|cholestasis|cirrhosis|kupffer cells|hepatic sinusoidal|bile canaliculi|fatty liver|lipidosis", "i")
}
GROUP BY ?aop ?aopNumber ?title ?aotitle
ORDER BY ?aopNumber

'

# Get response
response <- httr::POST("https://aopwiki.rdf.bigcat-bioinformatics.org/sparql", 
                       body = list(query = query, format = "json"), 
                       encode = "form")
# read the results
data <- jsonlite::fromJSON(httr::content(response, "text"))
results <- data$results$bindings

# Inspect returned results
head(results,1)
str(results)
names(results)

# Extract values from the nested JSON response
results$title$value

results <- data.frame(
  link = results$link$value,
  aopNumber = results$aopNumber$value,
  title = results$title$value,
  aotitle = results$aotitle$value,
  KEids = results$KEids$value,
  totalKEs = results$totalKEs$value
  
)

all_aops <- results

# save liver aops
write.table(all_aops, file = "all_liver_aops.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)



############################################################
# 2. Manual curations
############################################################

###### manually select the AOPs ###### 
###### Exclude non-liver, multi-organ and non-mammalian AOPs 

excluded_aops <- c("AOP 8", "AOP 537", "AOP 458", "AOP 194")

selected_liver_aops <-
  all_aops %>%
  filter(!(aopNumber %in% excluded_aops))

# save results
write.table(selected_liver_aops, file = 'selected_liver_aops.tsv',
            sep = '\t',
            row.names = F,
            quote = F)

######################################

######################################

# see the total num of KEs from selected_liver_aops
sum(as.integer(selected_liver_aops$totalKEs)) 

# there are 2 AOPs [273, 278] are appeared twice bc these have 2 AOs

# Generate a VALUES list for the downstream SPARQL queries
aop_list <- paste0('"', selected_liver_aops$aopNumber, '"', collapse = " ")
cat(aop_list)


############################################################
#3. Extract Key Events
############################################################

##### Get all the KEs (nodes) for liver network #####

query_kes <- paste0('

SELECT DISTINCT ?aopNumber ?ke ?keLabel ?keTitle ?type
WHERE {

    VALUES ?aopNumber { ', aop_list, ' }

    ?aop a aopo:AdverseOutcomePathway ;
         rdfs:label ?aopNumber .

    {
        ?aop aopo:has_molecular_initiating_event ?ke .
        ?ke rdfs:label ?keLabel ;
            dc:title ?keTitle .
        BIND("MIE" AS ?type)
    }

    UNION

    {
        ?aop aopo:has_adverse_outcome ?ke .
        ?ke rdfs:label ?keLabel ;
            dc:title ?keTitle .
        BIND("AO" AS ?type)
    }

    UNION

    {
        ?aop aopo:has_key_event ?ke .
        ?ke rdfs:label ?keLabel ;
            dc:title ?keTitle .

        FILTER NOT EXISTS {
            ?aop aopo:has_molecular_initiating_event ?ke
        }

        FILTER NOT EXISTS {
            ?aop aopo:has_adverse_outcome ?ke
        }

        BIND("KE" AS ?type)
    }
}

ORDER BY ?aopNumber ?type ?keLabel

')

# Get response
response_kes <- httr::POST("https://aopwiki.rdf.bigcat-bioinformatics.org/sparql", 
                           body = list(query = query_kes, format = "json"), 
                           encode = "form")
# read the results
data_kes <- jsonlite::fromJSON(httr::content(response_kes, "text"))
results_kes <- data_kes$results$bindings

# Inspect returned results
head(results_kes,1)
str(results_kes)
names(results_kes)

# Extract only the value columns from nested dataframes
results_kes <- data.frame(
  aopNumber = results_kes$aopNumber$value,
  ke = results_kes$ke$value,
  keLabel = results_kes$keLabel$value,
  keTitle = results_kes$keTitle$value,
  type = results_kes$type$value
)


write.table(results_kes, file = "kes.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# see the total number of KEs
length(results_kes$keLabel) 

# see the number of UNIQUE KEs/ nodes
length(unique(results_kes$keLabel)) 
# Total number of duplicated entries
sum(duplicated(results_kes$keLabel)) 

# see how many times each KE appears
table(results_kes$keLabel)

# see only those with more than 1 occurrence
table(results_kes$keLabel)[table(results_kes$keLabel) > 1]

# the reason for repetition: AOPs 273 and 278 are repeated twice as these have 2 AOs each. Moreover, some KEs are present in more than 1 AOPs

############################################################
#4. Extract Key Event Relationships
############################################################

##### Get the KERs (edges) for liver network #####

# define query
query_ker <- paste0('

SELECT DISTINCT
       ?ker
       ?kerLabel
       ?upstreamKE
       ?downstreamKE
       ?upstreamLabel
       ?downstreamLabel
WHERE {

    VALUES ?aopNumber { ', aop_list, ' }

    ?aop a aopo:AdverseOutcomePathway ;
         rdfs:label ?aopNumber ;
         aopo:has_key_event_relationship ?ker .

    ?ker a aopo:KeyEventRelationship ;
         aopo:has_upstream_key_event ?upstreamKE ;
         aopo:has_downstream_key_event ?downstreamKE ;
         rdfs:label ?kerLabel .

    ?upstreamKE rdfs:label ?upstreamLabel .
    ?downstreamKE rdfs:label ?downstreamLabel .

}

')

# Get response
response_ker <- httr::POST("https://aopwiki.rdf.bigcat-bioinformatics.org/sparql", 
                           body = list(query = query_ker, format = "json"), 
                           encode = "form")
# read the results
data_ker <- jsonlite::fromJSON(httr::content(response_ker, "text"))
results_ker <- data_ker$results$bindings

# Inspect returned results
head(results_ker,1)
str(results_ker)
names(results_ker)

# Extract only the value columns from nested dataframes
results_ker <- data.frame(
  ker = results_ker$ker$value,
  kerLabel = results_ker$kerLabel$value,
  upstreamKE = results_ker$upstreamKE$value,
  downstreamKE = results_ker$downstreamKE$value,
  upstreamLabel = results_ker$upstreamLabel$value,
  downstreamLabel = results_ker$downstreamLabel$value
)

# save results
write.table(results_ker, file = "ker.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)


# save session 


