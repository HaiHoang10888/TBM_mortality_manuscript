

# source("C:/Users/haiht/Dropbox/Project2020/WholeGenomeSeqHuman27TB/Hai_function/Mixomic_function.R")


subset_GOGSE<-function(x, list_G0)
{
  if (is(x, "gseaResult")) {
    final_result <- x
    ele_result <- final_result@result
    ele_result <- subset(ele_result,ID %in% list_G0)
    gene_Set <- final_result@geneSets[list_G0]
    core_gene <- unique(unlist(gene_Set))
    gene_list <- final_result@geneList[core_gene]
    final_result@result <- ele_result
    final_result@geneSets <- gene_Set
    final_result@geneList <- gene_list
    return(final_result)
  }
  if (is(x, "enrichResult")) {
    final_result <- x
    ele_result <- final_result@result
    ele_result <- subset(ele_result,ID %in% list_G0)
    gene_Set <- final_result@geneSets[list_G0]
    core_gene <- unique(unlist(gene_Set))
    final_result@result <- ele_result
    final_result@geneSets <- gene_Set
    
  }
  if (is(x, "compareClusterResult")) {
    final_result <- x
    ele_result <- final_result@compareClusterResult
    ele_result <- subset(ele_result,ID %in% list_G0)
    final_result@compareClusterResult <- ele_result
  }
  return(final_result)
}



# p <- draw_pheatmap_pathway_All(dds_anno=vst_data,class_goup="TB_type",DE_gene=DE_gene
#                                ,GOGSE=GOCompare@compareClusterResult$ID,
#                                gene_set=go.bp.gs,cluster_cols=F)
# data_tree <- Select_genes_for_GO_pathway(heatmap = p$pathway_heatmap,k=2,GOGSE = GOCompare,vst_data = vst_data ,
#                                          res_anno = res_anno_R,type = "spca",group ="TB_type" )


## create cluster for heatmap tree
##input is heatmap, and number of clusters
create_heat_tree_cluster <- function(heatmap,k){
  cluster = as.data.frame(cutree(heatmap$tree_row,k = k))
  colnames(cluster) <- "tree"
  cluster$GO <- rownames(cluster)
  cluster$ID <- substr(cluster$GO,1,10)
  data_tree <- list()
  for (i in 1:k){
    str <- paste0("cluster", i, " <-", " subset(cluster, tree==", i, ")")
    GO_cluster <-  eval(parse( text=str))
    str <- paste0("GO_cluster", i, " <-", "GO_cluster$GO")
    GO_cluster <-  eval(parse( text=str))
    GO_cluster <- as.vector(GO_cluster)
    str <- paste0("data_tree$tree$","GO_cluster",i, " <- GO_cluster")
    eval(parse( text=str))
  }
  rownames(cluster) <- NULL
  cluster <- cluster[order(cluster$tree),]
  data_tree$heatmap_tree <- cluster
  return(data_tree)
} 

## create list gene for 
create_list_genes_GO_cluster <- function(GOGSE, data_tree){
  for (i in names(data_tree$tree)){
    temp_list <- unname(unlist(data_tree$tree[i]))
    temp_list <- substr(temp_list,1,10)
    temp_GOGSE <- subset_GOGSE(x=GOGSE,list_G0=temp_list)
    if (is(temp_GOGSE, "gseaResult")) {
      temp_genes <- temp_GOGSE@result$core_enrichment
    }
    if (is(temp_GOGSE, "enrichResult")) {
      temp_genes <- temp_GOGSE@result$geneID
    }
    if (is(temp_GOGSE, "compareClusterResult")) {
      temp_genes <- temp_GOGSE@compareClusterResult$geneID
    }
    
    temp_genes <- scan(text = temp_genes, sep = "/", what = "")
    temp_genes <- unique(temp_genes)
    str <- paste0("data_tree$genes$",i, " <- temp_genes")
    eval(parse( text=str))
    #print(i)
    # print(temp_GOGSE@result$core_enrichment)
  }
  return(data_tree)
}


create_spca_for_heatmap_GO_cluster <- function(vst_data,data_tree, group="trial_arm",type="spca", plot=T){
  for (i in names(data_tree$genes)){
    temp_genes <- unname(unlist(data_tree$genes[i]))
    spca <- create_plsda_plot(data = vst_data,
                                genelist =temp_genes,group = group,
                                title = NULL,type=type, plot=plot)
    str <- paste0("data_tree$spca$",i," <-  spca")
    eval(parse( text=str))
  }
  return(data_tree)
}


create_Loading_genes_for_cluster_GO <- function(data_tree, res_anno){
  All <- data.frame(value.var=NA,comp=NA,GO_cluster=NA,genes=NA,log2FoldChange=NA,padj=NA)
  res_anno <- as.data.frame(res_anno)
  for(i in names(data_tree$spca)){
    #print(i)
    str <- paste0("temp_pca  <- ", "data_tree$spca$",i )
    eval(parse( text=str))
    type_pca <- class(temp_pca)
    temp_data1 <- selectVar(temp_pca, comp = 1)$value
    temp_data1$comp <- 1
    temp_data1$genes <- rownames(temp_data1)
    temp_data1_res <- subset(res_anno, rownames(res_anno) %in% rownames(temp_data1))
    temp_data1_res <-  temp_data1_res[,c("log2FoldChange","padj")]
    temp_data1_res <- temp_data1_res[rownames(temp_data1),]
    temp_data1 <- cbind(temp_data1,temp_data1_res)
    temp_data2 <- selectVar(temp_pca, comp = 2)$value
    temp_data2$comp <- 2
    temp_data2$genes <- rownames(temp_data2)
    temp_data2_res <- subset(res_anno, rownames(res_anno) %in% rownames(temp_data2))
    temp_data2_res <-  temp_data2_res[,c("log2FoldChange","padj")]
    temp_data2_res <- temp_data2_res[rownames(temp_data2),]
    temp_data2 <- cbind(temp_data2,temp_data2_res)
    if("mixo_plsda" %in% type_pca){
      temp_data1 <- temp_data1[1:5,]
      temp_data2 <- temp_data2[1:5,]
    }
    temp_data <- rbind(temp_data1, temp_data2)
    temp_data[!duplicated(temp_data$genes), ]
    temp_data$GO_cluster <- i
    str <- paste0("data_tree$spca_loading$",i," <-  temp_data")
    eval(parse( text=str))
    All <- rbind(All, temp_data)
    All <- All[!is.na(All$value.var),]
  }
  data_tree$spca_loading$All <- All
  return(data_tree)
}

Select_genes_for_GO_pathway <- function(heatmap,k,GOGSE,GO_target=NULL,vst_data,res_anno=res_anno, group="trial_arm",type="spca", plot=T){
  data_tree1 <- create_heat_tree_cluster(heatmap,k)
  data_tree2 <- create_list_genes_GO_cluster(GOGSE,data_tree1)
  data_tree3 <- create_spca_for_heatmap_GO_cluster(vst_data,data_tree2, 
                                                  group,type,plot)
  data_tree <- create_Loading_genes_for_cluster_GO(data_tree3, res_anno = res_anno)
  if(!is.null(GO_target)){
    GOGSE <- subset(GOGSE, GOGSE@result$ID %in% GO_target)
  }
  data_tree$GOGSE <- as.data.frame(GOGSE)
  return(data_tree)
  
}


# data_tree <- Select_genes_for_GO_pathway(heatmap = p$pathway_heatmap,k=3,
#                                         GOGSE = GOGSE_2,vst_data = vst_data_ratio ,res_anno = res_anno,
#                                         type = "plsda")




# ------other plot for PCA
# identical(rownames(vst_data), rownames(data_tree$spca$GO_cluster2$x))
# 
# col.sideColors <- palette()[vst_data$trial_arm]
# col.sideColors
# cim(data_tree$spca$GO_cluster3,transpose =T, row.sideColors = col.sideColors, cluster = "column")
# 
# plotVar(data_tree$spca$GO_cluster2,cex=3) 
# 
# 
# selectVar(data_tree$spca$GO_cluster2, comp = 1)$value
# 
# plotLoadings(data_tree$spca$GO_cluster3, comp=1)
# 
# biplot(data_tree$spca$GO_cluster2, ind.names = F, group=vst_data$trial_arm)


## create an excell file with each GO target in one sheet and list of all genes in each GO target
create_excel_output_Pathway_old <- function(GO_target, GO_data_base=go.bp.gs, res_anno,
                                        file){
  library(stringr)
  library(org.Hs.eg.db)
  library(xlsx)
  if (isClass(GO_target)) {
    GO_target <-  GO_target@result$ID
  }
  PW_gene <-  subset(GO_data_base, substr(names(GO_data_base),1,10) %in% GO_target)
  n =1 
  for (i in GO_target){
    temp_gene_ID <- unname(unlist(PW_gene[i]))
    temp_gene <-  mapIds(x = org.Hs.eg.db,
                         keys=temp_gene_ID,
                         column="SYMBOL",
                         keytype="ENTREZID",
                         multiVals="first",
                         ifnotfound = NA)
    temp_gene <- unname(temp_gene)
    res_anno_temp <- subset(res_anno, rownames(res_anno) %in% temp_gene)
    res_anno_temp <- res_anno_temp[order(res_anno_temp$padj),]
    i_name <- substr(i,12, nchar(i))
    if (nchar(i_name)>31) {
      i_name <- str_replace_all(i_name, ":", "")
      i_name <- str_replace_all(i_name, " ", "")
      i_name <- str_replace_all(i_name, "/", "or")
    }
    if (nchar(i_name)<31) {
      i_name <- str_replace_all(i_name, ":", "_")
      i_name <- str_replace_all(i_name, " ", "_")
      #i_name <- str_replace_all(i_name, "|", "or")
      i_name <- str_replace_all(i_name, "/", "or")
    }
    print(i_name)
    if (n<2) {
      write.xlsx(res_anno_temp, file = paste0(file,".xlsx"), 
                 sheetName=i_name, append=F)
    }
    else write.xlsx(res_anno_temp, paste0(file,".xlsx"), 
                    sheetName=i_name, append=TRUE)
    n = n+1
    
    #print(i_name)
    
  }
  
}


create_excel_output_Pathway <- function(GO_target, GO_data_base=go.bp.gs, res_anno,
                                        file){
  #options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
  library(stringr)
  library(org.Hs.eg.db)
  library(xlsx)
  library(dplyr)
  if (isClass(GO_target)) {
    GO_target <-  GO_target@result$ID
  }
  PW_gene <-  subset(GO_data_base, substr(names(GO_data_base),1,10) %in% GO_target)
  names(PW_gene) <- substr(names(PW_gene),1,10)
  t =1 
  wb = createWorkbook()
  for (i in GO_target){
    temp_gene_ID <- unname(unlist(PW_gene[i]))
    temp_gene <-  mapIds(x = org.Hs.eg.db,
                         keys=temp_gene_ID,
                         column="SYMBOL",
                         keytype="ENTREZID",
                         multiVals="first",
                         ifnotfound = NA)
    temp_gene <- unname(temp_gene)
    res_anno$gene_symbol <- rownames(res_anno)
    order_col <- move_column(names(res_anno), "gene_symbol first")
    res_anno <- res_anno[,order_col]
    res_anno_temp <- subset(res_anno, rownames(res_anno) %in% temp_gene)
    res_anno_temp <- res_anno_temp[order(res_anno_temp$padj),]
    i_name <- i
    i_name <- str_replace_all(i_name, ":", "_")
    sheet = createSheet(wb, i_name)
    addDataFrame(res_anno_temp, sheet=sheet, startColumn=1, row.names=FALSE)
    #print(i)
  }
  saveWorkbook(wb, paste0(file,".xlsx"))
}

## create res_ano for extract gene list pathway excell file
create_res_anno_ptrend <- function(res_anno, result_trend){
  result_trend <- result_trend[ncol(result_trend)]
  result_trend$ID <- rownames(result_trend)
  
  res_anno$p_and_FC_cutoff <- with(res_anno,ifelse(padj<p.cutoff & abs(log2FoldChange)>=LFC.cutoff, "Pass","Fail"))
  res_anno_M <- res_anno
  res_anno_M$ID <- rownames(res_anno_M)
  res_anno_M <- as.data.frame(res_anno_M)
  res_anno_M <- merge(res_anno_M,result_trend,by = "ID", all.x=TRUE)
  res_anno <- res_anno[res_anno_M$ID,]
  if(all(rownames(res_anno) == res_anno_M$ID)){
    res_anno$p.trend.adj <- res_anno_M$p.trend.adj
  }
  return(res_anno)
}

## create list of genes to validate based on pathway heatmap cluster.
write_csv_selected_genes_pathway_heatmap <- function(data_tree,file="selected_gene_to_validate_pathway_heatmap_D14"){
  library(xlsx)
  write.xlsx(data_tree$GOGSE, file = paste0(file,".xlsx"), 
             sheetName="GOGSE", append=F)
  write.xlsx(data_tree$heatmap_tree, file = paste0(file,".xlsx"), 
             sheetName="heatmap_tree", append=T)
  write.xlsx(data_tree$spca_loading$All, file = paste0(file,".xlsx"), 
             sheetName="spca_loading", append=T)
}
  

