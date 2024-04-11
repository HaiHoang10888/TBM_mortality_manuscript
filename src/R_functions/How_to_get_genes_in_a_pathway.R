# using ReactomeContentService4R  to get Human reactome pathway database
library(ReactomeContentService4R)
pathways <- getSchemaClass(class = "Pathway", species = "human", all = TRUE) # there are totally 2546 pathways in this database.

#name of 2546 pathways
pathways$displayName

#ID of 2546 pathways
pathways$stId

# get genes involved in a target pathway (3 steps)
## 1. get pathway ID (for example first pathway)
target_pathway_ID <- pathways$stId[1]

## 2. extraction information of this pathway 
t <- event2Ids(event.id = target_pathway_ID)

## 3. get all genes involved in first pathway
Pathway_genes <- t$geneSymbol
