# Performs Over Representation Analysis for a list of genes and a GMT
#
# @keywords internal
#
# @param topgenes a vector of genes
# @param gmt.list a gmt from prepare.gmt function
# @param allgenes a vector containing all genes to be considered as universe
#
# @return a data.frame containing the results
#
#
# res_list <- lapply(names(mods), ora_GO, mods=mods)
# 
# mods <- split(cem@module[, "genes"], cem@module[, "modules"])
# mod_name <- "black" 
# result <- ora_GO(mod_name,mods=mods)

# load("WGCNA_Dex_consensus/cem_WGCNA_Dex_D14_4041.RData")
# cem <- mod_ora_GO(cem)
# save(cem,file = "WGCNA_Dex_consensus/cem_ORA_GO_WGCNA_Dex_D14_4041.RData")

library(clusterProfiler)
library(org.Hs.eg.db)

pval2symbol <- function(matrix) {
  modtraitsymbol <- matrix
  modtraitsymbol[modtraitsymbol < 0.001] <- "***"
  modtraitsymbol[modtraitsymbol >= 0.001 & modtraitsymbol < 0.01] <- "**"
  modtraitsymbol[modtraitsymbol >= 0.01 & modtraitsymbol < 0.05] <- "*"
  modtraitsymbol[modtraitsymbol >= 0.05] <- ""
  return(modtraitsymbol)
}


ora_KEGG_single <- function(genelist, universe=NULL, OrgDb = org.Hs.eg.db){
  genelist <- bitr(genelist, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  genelist <- genelist$ENTREZID
  enriched <- clusterProfiler::enrichKEGG(gene = genelist,
                                        organism     = "hsa",
                                        keyType = "kegg",
                                        universe=universe,
                                        minGSSize    = 20,
                                        maxGSSize    = 800)
  enriched <- clusterProfiler::setReadable(enriched, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  return(enriched)
}

ora_GO_single <- function(genelist, universe=NULL, OrgDb = org.Hs.eg.db){
  genelist <- bitr(genelist, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  genelist <- genelist$ENTREZID
  enriched <- clusterProfiler::enrichGO(gene = genelist,
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.05,
                                        OrgDb = OrgDb,
                                        keyType = "ENTREZID",
                                        ont="BP",
                                        minGSSize    = 20,
                                        maxGSSize    = 800,
                                        readable      = T)
  return(enriched)
}





ora_GO <- function(mod_name, universe=NULL, OrgDb = org.Hs.eg.db, mods){
  topgenes <- mods[[mod_name]]
  topgenes <- bitr(topgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  topgenes <- topgenes$ENTREZID
  enriched <- clusterProfiler::enrichGO(gene = topgenes,
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.05,
                                        OrgDb = OrgDb,
                                        keyType = "ENTREZID",
                                        ont="BP",
                                        minGSSize    = 20,
                                        maxGSSize    = 500,
                                        readable      = T)
  if (!is.null(enriched) && !is.logical(enriched)) {
    result <- enriched@result
    #result <- subset(result, p.adjust < 0.05)
  } else {
    if(mod_name != "Not.Correlated"){
      warning("Enrichment for module ", mod_name, " is NULL")
    }
    result <- data.frame(Module=character(), ID=character(),
                         Description=character(),
                         GeneRatio=numeric(), BgRatio=numeric(),
                         pvalue=numeric(), p.adjust=numeric(),
                         qvalue=numeric(), geneID=character(),
                         Count=numeric(), stringsAsFactors=FALSE)
  }
  return(result)
}



mod_ora_GO <- function(cem, verbose=FALSE){
  if (verbose) {
    message('Running ORA')
    message("Using all genes in org.Hs.eg.db file as universe.")
  }
  if(is.null(module_genes(cem))){
    warning("No modules in CEMiTool object! Did you run find_modules()?")
    return(cem)
  }
  mods <- split(cem@module[, "genes"], cem@module[, "modules"])
  res_list <- lapply(names(mods),ora_GO, mods=mods)
  if (all(lapply(res_list, nrow) == 0)){
    warning("Enrichment is NULL. Either your gene list is inadequate or your modules really aren't enriched for any of the pathways in GO.")
    return(cem)
  }
  names(res_list) <- names(mods)
  res <- lapply(names(res_list), function(x){
    if(nrow(res_list[[x]]) > 0){
      as.data.frame(cbind(x, res_list[[x]]))
    }
  })
  res <- do.call(rbind, res)
  names(res)[names(res) == "x"] <- "Module"
  
  rownames(res) <- NULL
  cem@ora <- res
  return(cem)
}




#' ORA visualization for one module
#'
#' @keywords internal
#'
#' @param es a data.frame from ora function containing only one module
#' @param ordr_by column to order the data.frame
#' @param max_length max length of a gene set name
#' @param pv_cut p-value cuttoff
#' @param graph_color color of bars
#' @param title title of the graph
#'
#' @return a list with ggplot2 object and the number of significant gene sets

plot_ora_single <- function(es, ordr_by='p.adjust', max_length=70, pv_cut=0.05,
                            graph_color="#4169E1", title="Over Representation Analysis"){
  #x=test[dupes]
  comsub <- function(x){
    #split the first and last element by character
    d_x <- strsplit(x[c(1, length(x))], "")
    #search for the first not common element, and so, get the last matching one
    der_com <- match(FALSE, do.call("==", d_x))-1
    return(substr(x, 1, 60))
  }
  es[, "GeneSet"] <- es[, "Description"]
  
  # limits name length
  ovf_rows <- which(nchar(es[, "GeneSet"]) > max_length) # overflow
  ovf_data <- es[ovf_rows, "GeneSet"]
  test <- strtrim(ovf_data, max_length)
  dupes <- duplicated(test) | duplicated(test, fromLast=TRUE)
  if(sum(dupes) > 0){
    test[dupes] <- ovf_data[dupes]
    test[dupes] <- comsub(test[dupes])
    max_length <- max(nchar(test))
  }
  
  es[ovf_rows, "GeneSet"] <-  paste0(strtrim(test, max_length), "...")
  es[, "GeneSet"] <- stringr::str_wrap(es[, "GeneSet"], width = 20)
  
  # order bars
  lvls <- es[order(es[, ordr_by], decreasing=TRUE), "GeneSet"]
  es[, "GeneSet"] <- factor(es[, "GeneSet"], levels=unique(lvls))
  es[, "alpha"] <- 1
  es[es[, ordr_by] > pv_cut, "alpha"] <- 0
  
  # Avoid 0's
  es[es[, ordr_by] > 0.8, ordr_by] <- 0.8
  my_squish <- function(...){
    return(scales::squish(..., only.finite=FALSE))
  }
  
  # plot
  y_axis <- paste('-log10(', ordr_by, ')')
  pl <- ggplot(es, aes_string(x="GeneSet", y=y_axis, alpha="alpha", fill=y_axis)) +
    geom_bar(stat="identity") +
    theme(axis.text=element_text(size=8), legend.title=element_blank()) +
    coord_flip() +
    scale_alpha(range=c(0.4, 1), guide="none") +
    labs(y="-log10(adjusted p-value)", title=title, x="") +
    geom_hline(yintercept=-log10(pv_cut), colour="grey", linetype="longdash") +
    scale_fill_gradient(low="gray", high=graph_color, limits=c(2, 5), oob=my_squish)
  #res <- list('pl'=pl, numsig=sum(es[, ordr_by] < pv_cut, na.rm=TRUE))
  return(pl)
}
        
