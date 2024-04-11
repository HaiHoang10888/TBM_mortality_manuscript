

## How to use

#GOCompare <- get(load(paste0(pw_out,"TBM_PTB_compareCluster_validation100times.RData")))
# PPI <- read.csv("Pathways_output/PPI_Subnetwork1_Minimum_node_table.csv", header = T)
# res_anno_DE <- Load_res(data_name = res_file,filter = T)
# res_anno_DE <- as.data.frame(res_anno_DE)
# PPI <- subset(PPI,Expression!=0)
# PPI <- reannotate_genesymbol_PPI(PPI=PPI, anno=anno)
# res_export <- merge_res_PPI_pathway(res = res_anno_DE,pathway = GOCompare,PPI = PPI)
# write.csv(res_export,"output_folder/1237DEGenes_Merged_TBM_vs_PTB_D0.csv")


reannotate_genesymbol_PPI <- function(PPI,anno){
  index <- which(PPI$Label %nin% anno$gene_symbol)
  if (length(index)>0) {
    for (i in 1:length(index)){
      temp_index <- index[i]
      temp_entrez <- PPI$Id[temp_index]
      eg = bitr(temp_entrez, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
      gene <- eg$SYMBOL
      PPI$Label[temp_index] <- gene
    }
  }
  return(PPI)
  
}


merge_res_PPI_pathway <- function(res,pathway,PPI){
  for (i in 1:nrow(res)) {
    if (is(pathway, "gseaResult")|is(pathway, "enrichResult")) {
      Pw_data <- pathway@result
    }
    if (is(pathway, "compareClusterResult")) {
      Pw_data <- pathway@compareClusterResult
    }
    Pw_data <- Pw_data[order(Pw_data$p.adjust),]
    colnames(Pw_data)[which(colnames(Pw_data)=="core_enrichment")] <- "geneID"
    temp_genes <- rownames(res)[i]
    pwt <- NULL
    Pval <- NULL
    pwt_oder <- NULL
    for (pt in 1:nrow(Pw_data)){
      Tem_pw_gene <-  Pw_data$geneID[pt]
      Tem_pw <- Pw_data$Description[pt]
      tem_p <- Pw_data$p.adjust[pt]
      Tem_pw_gene <- scan(text = Tem_pw_gene, sep = "/", what = "")
      if (temp_genes %in% Tem_pw_gene) {
        pwt <- paste(pwt,Tem_pw,sep = "/")
        Pval <- paste(Pval,tem_p,sep = "/")
        pwt_oder <- paste(pwt_oder,pt,sep = "/")
      }
      
      res$pawthay[i] <- NA
      res$Pvalue_pawthay[i] <- NA
      res$Order_pawthay[i] <- NA
      if (!is.null(pwt)) {
        res$pawthay[i] <- pwt
        res$Pvalue_pawthay[i] <- Pval
        res$Order_pawthay[i] <- pwt_oder
      }
      
    }
  }
  res1 <- res[PPI$Label,]
  res1$PPI_degree <- PPI$Degree
  res1$PPI_betweenness <- PPI$Betweenness
  res2 <- subset(res, rownames(res) %nin% PPI$Label)
  res2$PPI_degree <- NA
  res2$PPI_betweenness <- NA
  res <- rbind(res1,res2)
  return(res)
}


