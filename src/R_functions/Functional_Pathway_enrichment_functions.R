# 10. Functional Pathway

## get genelist in pathways

get_genelist_pathway <- function(name=NULL, go_id=NULL){
  if(!is.null(name)) go_id = GOID( GOTERM[ Term(GOTERM) == name])
  allegs = get(go_id, org.Hs.egGO2ALLEGS)
  gene_list = unlist(mget(allegs,org.Hs.egSYMBOL))
  gene_list <- unname(gene_list)  
  gene_list <- unique(gene_list)
  return(gene_list)
}

get_genelist_pathway_kegg <- function(hsa_ID){
  library("KEGGREST")
  names <- keggGet(hsa_ID)[[1]]$GENE
  namesodd <-  names[seq(0,length(names),2)]
  namestrue <- gsub("\\;.*","",namesodd)
}


## reactome pathway
GSEREACTOME_function <- function(res,p.cutoff){
  DEGs=subset(res,res$padj<p.cutoff)
  DEGs$ENTREZID <-  mapIds(x = org.Hs.eg.db,
                           keys=row.names(DEGs),
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first",
                           ifnotfound = NA)
  DEGs <- DEGs[!is.na(DEGs$ENTREZID),]
  foldchanges  <-DEGs$log2FoldChange
  names(foldchanges) <- DEGs$ENTREZID
  foldchanges <- sort(foldchanges,decreasing = T)
  GSEREACTOME <- gsePathway(foldchanges, 
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH", 
                            verbose = FALSE)
  return(GSEREACTOME)
}

View_REACTOME_function <- function(res,p.cutoff,name,layout){
  DEGs=subset(res,res$padj<p.cutoff)
  DEGs$ENTREZID <-  mapIds(x = org.Hs.eg.db,
                           keys=row.names(DEGs),
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first",
                           ifnotfound = NA)
  DEGs <- DEGs[!is.na(DEGs$ENTREZID),]
  foldchanges  <-DEGs$log2FoldChange
  names(foldchanges) <- DEGs$ENTREZID
  foldchanges <- sort(foldchanges,decreasing = T)
  viewPathway(name, 
              readable = TRUE, 
              foldChange = foldchanges, layout = layout)
}


View_REACTOME_function_pvalue <- function(res,p.cutoff,name,layout){
  DEGs=subset(res,res$pvalue<p.cutoff)
  DEGs$ENTREZID <-  mapIds(x = org.Hs.eg.db,
                           keys=row.names(DEGs),
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first",
                           ifnotfound = NA)
  DEGs <- DEGs[!is.na(DEGs$ENTREZID),]
  foldchanges  <-DEGs$log2FoldChange
  names(foldchanges) <- DEGs$ENTREZID
  foldchanges <- sort(foldchanges,decreasing = T)
  viewPathway(name, 
              readable = TRUE, 
              foldChange = foldchanges, layout = layout)
}

## 10.1 GOEnrich_function
##----------function 10.1
GOEnrich_function <- function(res,p.cutoff,universe=NULL){
  DEGs=subset(res,res$padj<p.cutoff)
  DEGs$ENTREZID <-  mapIds(x = org.Hs.eg.db,
                           keys=row.names(DEGs),
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first",
                           ifnotfound = NA)
  DEGs <- DEGs[!is.na(DEGs$ENTREZID),]
  if(!is.null(universe)){
    universe$ENTREZID <-  mapIds(x = org.Hs.eg.db,
                                 keys=row.names(universe),
                                 column="ENTREZID",
                                 keytype="SYMBOL",
                                 multiVals="first",
                                 ifnotfound = NA)
    universe <- universe[!is.na(universe$ENTREZID),]
    universe <- universe$ENTREZID
  }
  
  
  GOenrich <- enrichGO(gene = DEGs$ENTREZID,
                       #universe = row.names(res),
                       OrgDb = org.Hs.eg.db,
                       universe= universe,
                       keyType = "ENTREZID",
                       ont="BP",
                       minGSSize    = 10,
                       maxGSSize    = 500,
                       readable      = TRUE)
  return(GOenrich)
}



## 10.2
##----------function 10.2
KEGGEnrich_function <- function(res,p.cutoff,universe=NULL){
  DEGs=subset(res,res$padj<p.cutoff)
  DEGs$ENTREZID <-  mapIds(x = org.Hs.eg.db,
                           keys=row.names(DEGs),
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first",
                           ifnotfound = NA)
  DEGs <- DEGs[!is.na(DEGs$ENTREZID),]
  if(!is.null(universe)){
    universe$ENTREZID <-  mapIds(x = org.Hs.eg.db,
                                 keys=row.names(universe),
                                 column="ENTREZID",
                                 keytype="SYMBOL",
                                 multiVals="first",
                                 ifnotfound = NA)
    universe <- universe[!is.na(universe$ENTREZID),]
    universe <- universe$ENTREZID
  }
  KEGGenrich <- enrichKEGG(gene = DEGs$ENTREZID,
                           organism     = "hsa",
                           keyType = "kegg",
                           universe=universe,
                           minGSSize    = 10,
                           maxGSSize    = 800)
  return(KEGGenrich)
}


## 10.3 GOGSE_function
##----------function 10.3
GOGSE_function <- function(res,p.cutoff){
  DEGs=subset(res,res$padj<p.cutoff)
  DEGs$ENTREZID <-  mapIds(x = org.Hs.eg.db,
                           keys=row.names(DEGs),
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first",
                           ifnotfound = NA)
  DEGs <- DEGs[!is.na(DEGs$ENTREZID),]
  foldchanges  <-DEGs$log2FoldChange
  names(foldchanges) <- rownames(DEGs)
  foldchanges <- sort(foldchanges,decreasing = T)
  GOGSE <- gseGO(geneList     = foldchanges,
                 OrgDb          = org.Hs.eg.db,
                 ont = "BP",
                 keyType     = "SYMBOL",
                 by = "fgsea",
                 minGSSize    = 10,
                 pvalueCutoff = 0.05,
                 verbose      = TRUE)
  return(GOGSE)
}

## 10.4 KEGGGSE_function
##----------function 10.4
KEGGGSE_function <- function(res,p.cutoff){
  DEGs=subset(res,res$padj<p.cutoff)
  DEGs$ENTREZID <-  mapIds(x = org.Hs.eg.db,
                           keys=row.names(DEGs),
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first",
                           ifnotfound = NA)
  DEGs <- DEGs[!is.na(DEGs$ENTREZID),]
  foldchanges  <-DEGs$log2FoldChange
  names(foldchanges) <- DEGs$ENTREZID
  foldchanges <- sort(foldchanges,decreasing = T)
  KEGGGSE <- gseKEGG(geneList     = foldchanges,
                     organism     = "hsa",
                     by = "fgsea",
                     minGSSize    = 10,
                     maxGSSize    = 800,
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     keyType       = "kegg")
  return(KEGGGSE)
}

## 10.5 Final pathway function
##----------function Final 10 Enrich function
Enrich_function <- function(res,p.cutoff,type,universe=NULL,export=F,name=NULL){
  if (type=="GOEnrich") {
    result <- GOEnrich_function(res,p.cutoff)
  }
  else if (type=="KEGGEnrich") {
    result <- KEGGEnrich_function(res,p.cutoff)
  }
  else if (type=="GOGSE") {
    result <- GOGSE_function(res,p.cutoff)
  }
  else if (type=="KEGGGSE") {
    result <- KEGGGSE_function(res,p.cutoff)
  }
  if (export==T) {
    result_EX <- data.frame(result)
    result_EX <- result_EX[order(result_EX$p.adjust),]
    write.csv(result_EX,paste(name,".csv",sep = ""))
  }
  return(result)
}



## 10.6 GOGSE_function using all genes
##----------function 10.3
GOGSE_function_All <- function(res,p.cutoff){
  DEGs=subset(res,res$padj<p.cutoff)
  DEGs$ENTREZID <-  mapIds(x = org.Hs.eg.db,
                           keys=row.names(DEGs),
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first",
                           ifnotfound = NA)
  DEGs <- DEGs[!is.na(DEGs$ENTREZID),]
  DEGs$foldchanges <- -log10(DEGs$pvalue)*sign(DEGs$log2FoldChange)
  foldchanges  <-DEGs$foldchanges
  names(foldchanges) <- rownames(DEGs)
  foldchanges <- sort(foldchanges,decreasing = T)
  GOGSE <- gseGO(geneList     = foldchanges,
                 OrgDb          = org.Hs.eg.db,
                 ont = "BP",
                 keyType     = "SYMBOL",
                 by = "fgsea",
                 minGSSize    = 20,
                 pvalueCutoff = 0.05,
                 verbose      = TRUE)
  return(GOGSE)
}

## 10.7 KEGGGSE_function using all genes
##----------function 10.4
KEGGGSE_function_All <- function(res,p.cutoff){
  DEGs=subset(res,res$padj<p.cutoff)
  DEGs$ENTREZID <-  mapIds(x = org.Hs.eg.db,
                           keys=row.names(DEGs),
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first",
                           ifnotfound = NA)
  DEGs <- DEGs[!is.na(DEGs$ENTREZID),]
  DEGs$foldchanges <- -log10(DEGs$pvalue)*sign(DEGs$log2FoldChange)
  foldchanges  <-DEGs$foldchanges
  names(foldchanges) <- DEGs$ENTREZID
  foldchanges <- sort(foldchanges,decreasing = T)
  KEGGGSE <- gseKEGG(geneList     = foldchanges,
                     organism     = "hsa",
                     by = "fgsea",
                     minGSSize    = 10,
                     maxGSSize    = 500,
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     keyType       = "kegg")
  return(KEGGGSE)
}

## 10.8 Final pathway function using all genes
##----------function Final 10 Enrich function
Enrich_function_All <- function(res,p.cutoff,type,export,name=NULL){
  if (type=="GOEnrich") {
    result <- GOEnrich_function(res,p.cutoff)
  }
  else if (type=="KEGGEnrich") {
    result <- KEGGEnrich_function(res,p.cutoff)
  }
  else if (type=="GOGSE") {
    result <- GOGSE_function_All(res,p.cutoff)
  }
  else if (type=="KEGGGSE") {
    result <- KEGGGSE_function_All(res,p.cutoff)
  }
  if (export==T) {
    result_EX <- data.frame(result)
    result_EX <- result_EX[order(result_EX$p.adjust),]
    write.csv(result_EX,paste(name,".csv",sep = ""))
  }
  return(result)
}
