
# I. General functions for pre-process count data and clinical data

## 1. Process count expression data
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

create_count_data_all <- function(workdir,Filter_rRNA=F,Filter_Sex=F){
  count_data <- paste0(workdir,"count_expression/raw_data/Run1/count_old_Run1vsReRun1_separate.csv")
  anno <- paste0(workdir,"count_expression/raw_data/GeneAnnotationFile_EnsembltoGeneIds.txt")
  #----------------load Count_data
  count_data <-  read.csv(count_data)
  #----------------Load anno file
  anno <- read.table(anno,header = T)
  anno$class_new <- with(anno, ifelse(class=="protein_coding", "protein_coding", ifelse(class=="rRNA_pseudogene"|class=="rRNA", "rRNA","non protein_coding")))
  count_data_all_X <- merge(anno, count_data, by.x = "gene_id", by.y = "Geneid")
  rownames(count_data_all_X) <- count_data_all_X$gene_id
  
  # subset for rRNA and hemoglobin  (n = 557 and 11)
  anno_rRNA <- anno[which(anno$class %in% c("rRNA_pseudogene", "rRNA")),]
  globin_genes <- c("HBA1","HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ", "MB")
  anno_globin <- anno[which(anno$gene_symbol %in% globin_genes),]
  
  # subset for for gene that not encode protein (n = 40708)
  anno_not_protein <- anno[which(anno$class !="protein_coding"),]
  anno_protein <- anno[which(anno$class =="protein_coding"),]
  
  # subset for genes in chromosome X or Y (n = 2946)
  anno_sex_chr_gene <- anno[which(anno$chr =="Y"|anno$chr =="X"),]
  
  #remove reads mapping to globin and rRNA genes (n=568) and sex chromosome
  Exlcuded_rRNA <- c(as.character(anno_globin$gene_id),
                     as.character(anno_rRNA$gene_id))
  Exlcuded_sex <- as.character(anno_sex_chr_gene$gene_id)
  Exlcuded_rRNA_sex <- c(Exlcuded_rRNA,Exlcuded_sex)
  
  if (Filter_rRNA==T & Filter_Sex==T) {
    count_data_all_X <- subset(count_data_all_X,
                               rownames(count_data_all_X) %nin% Exlcuded_rRNA_sex)
  } else if (Filter_rRNA==T) {
    count_data_all_X <- subset(count_data_all_X,
                               rownames(count_data_all_X) %nin% Exlcuded_rRNA)
  } else if (Filter_Sex==T) {
    count_data_all_X <- subset(count_data_all_X,
                               rownames(count_data_all_X) %nin% Exlcuded_sex)
  }
  
  return(count_data_all_X)
  
}

## ---- can create count data check to check for merge steps
# rownames(count_data) <- count_data$Geneid
# rownames(anno) <- anno$gene_id
# count_data <- count_data[rownames(anno),]
# count_data_all_check <- cbind(anno,count_data)
# count_data_all_check <- count_data_all_check[rownames(count_data_all),]
# count_data_all_check <- count_data_all_check[,colnames(count_data_all)]
# identical(count_data_all,count_data_all_check)
## check the result of function create countdata
# nrow(count_data_all_X)
# nrow(count_data)
# nrow(anno)
# rownames(count_data) <- count_data$Geneid
# identical(rownames(count_data), anno$gene_id)
# count_data_all_X <- count_data_all_X[anno$gene_id,]
# identical(rownames(count_data_all_X), anno$gene_id)
# identical(count_data_all_X$chr_position, anno$chr_position)
# 


## 2. create annotation data for gene
create_anno_data <- function(workdir){
  anno <- paste0(workdir,"count_expression/raw_data/GeneAnnotationFile_EnsembltoGeneIds.txt")
  #----------------Load anno file
  anno <- read.table(anno,header = T)
  anno$class_new <- with(anno, ifelse(class=="protein_coding", "protein_coding", ifelse(class=="rRNA_pseudogene"|class=="rRNA", "rRNA","non protein_coding")))
  rownames(anno) <- anno$gene_id
  return(anno)
}


## 3. subset for housekeeping genes
#------------------ Load list of houseskepping genes list
create_HK_genes_11_list <- function(anno){
  HK_genes_11 <- c("C1orf43", "CHMP2A", "EMC7", "GPI", "PSMB2", "PSMB4", "RAB7A", "REEP5", "SNRPD3", "VCP", "VPS29")
  HK_genes_11_anno <- anno[which(anno$gene_symbol %in% HK_genes_11),]
  HK_genes_11_list <- as.character(HK_genes_11_anno$gene_id)
  return(HK_genes_11_list)
}

create_HK_genes_HUPO <- function(anno){
  HK_genes <- c("RPLPO")
  HK_genes_anno <- anno[which(anno$gene_symbol %in% HK_genes),]
  HK_genes_11_list <- as.character(HK_genes_anno$gene_id)
  return(HK_genes_11_list)
}


## 4. general function to read in TBM and PTB data and exclude low quality samples
create_data <- function(workdir, data){
  ## list of sample with RIn <4 and Read depth <10
  Low_quality_sample <- c("HOA7818A211", "HOA7818A1", "HOA7818A127", "HOA7818A92")
  ## Sample pair with low quality samples
  paired_sample <- c("HOA7818A128", "HOA7812A286", "HOA7814A282", "HOA7818A210", "HOA7814A185")
  ## samples need to exclude
  Excluded_sample <- c(Low_quality_sample, paired_sample)
  data <- paste0(workdir,data)
  data <- get(load(data))
  rownames(data) <- data$LIMS_ID
  ## filter for low quality samples
  data <- subset(data,LIMS_ID %nin% Excluded_sample)
  data$LTA4H <- factor(data$LTA4H, levels = c("CC","CT","TT"))
  return(data)
}


## 5. subset for clinical data
create_TB_TBM_data <- function(workdir, data_PTB, data_TBM, output){
  ## add information about day stop using study drugs
  STDR_stop <-  paste0(workdir,"clinical_data/process_data/STDRSTOP_27TB_27MAY2020.csv") 
  STDR_stop <- read.csv(STDR_stop)
  STDR_stop_D14 <- subset(STDR_stop, Exclude_D14 =="Y"| Exclude_D14 =="UK")
  STDR_stop_D60 <- subset(STDR_stop, Exclude_D60 =="Y")
  Excluded_sample_D14 <- c(STDR_stop_D14$RNA.seq_base, 
                           STDR_stop_D14$RNA.seq_D14)
  Excluded_sample_D60 <- c(STDR_stop_D60$RNA.seq_base,
                           STDR_stop_D60$RNA.seq_D14,STDR_stop_D60$RNA.seq_D60)
  Antidrug_use_D14 <- subset(STDR_stop, Anti_D14 =="EX")
  Antidrug_use_D60 <- subset(STDR_stop, Anti_D60 =="EX")
  
  # create list of samples to exclude b/c stop using study drugs
  Antidrug_exclude_D14 <- c(Antidrug_use_D14$RNA.seq_base, 
                            Antidrug_use_D14$RNA.seq_D14)
  Antidrug_exclude_D60 <- c(Antidrug_use_D60$RNA.seq_base, 
                            Antidrug_use_D60$RNA.seq_D14)
  
  # read in PTB and TBM data
  data_PTB <- create_data(workdir=workdir,data = data_PTB)
  data_TBM <- create_data(workdir=workdir,data = data_TBM)
  
  ## subset condition for data and check
  cond_out <- c("TBM_D0", "PTB_D0", "TB_D0", 
                "All_TBM_D0", "All_TB_D0",
                "TBM_D14", "TBM_D0_D14", "TBM_D60",
                "TBM_D0_D60", "TBM_D0_D14_D60", "TBM_D14_D60", "All_27TB")
  if(output %nin% cond_out) 
    stop("Error: output should be in TBM_D0, PTB_D0, TB_D0, All_TBM_D0, All_PTB_D0, All_TB_D0, TBM_D14, TBM_D0_D14, TBM_D60, TBM_D0_D60, TBM_D0_D14_D60, TBM_D14_D60")
  if (output == "TBM_D0") {
    data <- subset(data_TBM,trial_number=="27TB" & Timepoint=="D0")
  }
  else if (output == "All_27TB") {
    data <- subset(data_TBM,trial_number=="27TB")
  }
  else if (output == "PTB_D0") {
    data <- subset(data_PTB,trial_number=="MDR_PTB" | trial_number=="DS_PTB")
  }
  else if (output == "TB_D0") {
    data_PTB$TB_type <- "PTB"
    data_TBM$TB_type <- "TBM"
    select_var <- colnames(data_PTB)[which(colnames(data_PTB) %in% colnames(data_TBM))]
    data_PTB_S <- subset(data_PTB,trial_number=="MDR_PTB" | trial_number=="DS_PTB")
    data_TBM_S <- subset(data_TBM,trial_number=="27TB" & Timepoint=="D0")
    select_PTB <- data_PTB_S[,select_var]
    select_TBM <- data_TBM_S[,select_var]
    data <- rbind(select_PTB,select_TBM)
  }
  else if (output == "All_TBM_D0") {
    data <- subset(data_TBM, Timepoint=="D0")
  }
  else if (output == "All_TB_D0") {
    data_PTB$TB_type <- "PTB"
    data_TBM$TB_type <- "TBM"
    select_var <- colnames(data_PTB)[which(colnames(data_PTB) %in% colnames(data_TBM))]
    select_PTB <- data_PTB[,select_var]
    select_TBM <- data_TBM[,select_var]
    data <- rbind(select_PTB,select_TBM)
  }
  #subset for data and exclude samples that stop using study drugs
  else if (output == "TBM_D14") {
    data <- subset(data_TBM,trial_number=="27TB" & Timepoint=="D14")
    data <- subset(data, LIMS_ID %nin% Excluded_sample_D14)
    data <- subset(data, LIMS_ID %nin% Antidrug_exclude_D14)
  }
  else if (output == "TBM_D0_D14") {
    data_D14 <- subset(data_TBM,trial_number=="27TB" & Timepoint=="D14")
    data_D14 <- subset(data_D14, LIMS_ID %nin% Excluded_sample_D14)
    data_D14 <- subset(data_D14, LIMS_ID %nin% Antidrug_exclude_D14)
    data <- subset(data_TBM, Patient_ID %in% data_D14$Patient_ID)
    data <- subset(data, Timepoint=="D0"|Timepoint=="D14")
  }
  else if (output == "TBM_D60") {
    data <- subset(data_TBM,trial_number=="27TB" & Timepoint=="D60")
    data <- subset(data, LIMS_ID %nin% Excluded_sample_D60)
    data <- subset(data, LIMS_ID %nin% Antidrug_exclude_D60)
  }
  else if (output == "TBM_D0_D60") {
    data_D60 <- subset(data_TBM,trial_number=="27TB" & Timepoint=="D60")
    data_D60 <- subset(data_D60, LIMS_ID %nin% Excluded_sample_D60)
    data_D60 <- subset(data_D60, LIMS_ID %nin% Antidrug_exclude_D60)
    data <- subset(data_TBM, Patient_ID %in% data_D60$Patient_ID)
    data <- subset(data, Timepoint != "D14")
  }
  # subset for data and exclude samples that stop using study drugs
  else if (output == "TBM_D0_D14_D60") {
    data_D60 <- subset(data_TBM,trial_number=="27TB" & Timepoint=="D60")
    data_D60 <- subset(data_D60, LIMS_ID %nin% Excluded_sample_D60)
    data_D60 <- subset(data_D60, LIMS_ID %nin% Antidrug_exclude_D60)
    data <- subset(data_TBM, Patient_ID %in% data_D60$Patient_ID)
  }
  else if (output == "TBM_D14_D60") {
    data_D60 <- subset(data_TBM,trial_number=="27TB" & Timepoint=="D60")
    data_D60 <- subset(data_D60, LIMS_ID %nin% Excluded_sample_D60)
    data_D60 <- subset(data_D60, LIMS_ID %nin% Antidrug_exclude_D60)
    data <- subset(data_TBM, Patient_ID %in% data_D60$Patient_ID)
    data <- subset(data, Timepoint != "D0")
  }
  else  data <- NULL
  
  return(data)
}

## 6. Function for subset expression gene data for subset of study or sample
create_exp_data <- function(count_data, clinical_data){
  rownames(count_data) <- count_data$gene_id
  # set rowname for clinical data using ID of sequecing
  rownames(clinical_data) <- clinical_data$LIMS_ID
  exp_data <- count_data[,which(colnames(count_data) %in% rownames(clinical_data))]
  exp_data <- exp_data[,rownames(clinical_data)]
  return(exp_data)
}


# II. Process data for Deseq2 packages
## 1. Function to select variable for colData in Deseq2
create_coldata <- function(data,Select_Var){
  coldata <- data[, Select_Var]
  # renames timepoint and treatment if data contains these vars
  if (!is.null(coldata$Timepoint)){
    coldata$Time <- with(coldata, ifelse(Timepoint == "D0",0,ifelse(Timepoint=="D14", 14, 60)))
    coldata$Time <- as.factor(coldata$Time)
  }
  if (!is.null(coldata$Timepoint) & !is.null(coldata$trial_arm)) {
    coldata$Dex_time <- as.factor(paste(coldata$trial_arm,coldata$Time,sep = "_"))
  }
  return(coldata)
}


create_dds <- function(count,coldata,group,number,filter,genelist,design1,design2,fitType){
  cts <- count
  str <- paste("min(table(coldata$",group,"))","/number",sep = "")
  filter_min_group <-  eval(parse( text=str ))
  str2 <- paste0("coldata$",group)
  str2 <- eval(parse( text=str2))
  if(is.numeric(str2)) filter_min_group <- length(str2)/5
  ## filter=5, exclude gene with count less than 5 in min of group samples size,
  # if treatment is binary, if treatment is continuous, exclude genes with count < 5 in 20% sample 
  filter <- apply(cts, 1, function(x) length(x[x>filter]) >= filter_min_group)
  cts <- cts[filter,]
  cts <- as.matrix(cts)
  ## create coldata
  set <- newSeqExpressionSet(as.matrix(cts),
                             phenoData = coldata)
  ## normalized data
  set_N <- betweenLaneNormalization(set, which="upper",offset=F)
  set_RUV <- RUVg(set_N, genelist, k=1)
  str1 <- paste("DESeqDataSetFromMatrix(countData = counts(set_RUV),colData= pData(set_RUV), design = ~", design1, ")",sep = "")
  dds <- eval(parse( text=str1))
  dds <- DESeq(dds,minReplicatesForReplace=7, parallel = T,fitType = c("parametric"))
  dds <- replaceOutliersWithTrimmedMean(dds)
  str2 <- paste("DESeqDataSetFromMatrix(countData = counts(dds),colData= pData(set_RUV), design = ~", design2,")",sep = "")
  dds <- eval(parse( text=str2))
  dds <- DESeq(dds,minReplicatesForReplace=7, parallel = T,fitType = c("parametric"))
  return(dds)
}



create_dds_1times <- function(count,coldata,group,number,filter,gene_ref,design,fitType="parametric",
                              parallel=T){
  cts <- count
  # filter for low expressed genes
  str <- paste("min(table(coldata$",group,"))","/number",sep = "")
  filter_min_group <-  eval(parse( text=str ))
  str2 <- paste0("coldata$",group)
  str2 <- eval(parse( text=str2))
  ## filter=5, exclude gene with count less than 5 in min of group samples size,
  # if treatment is binary, if treatment is continuous, exclude genes with count < 5 in 20% sample
  if(is.numeric(str2)) filter_min_group <- length(str2)/5
  filter <- apply(cts, 1, function(x) length(x[x>filter]) >= filter_min_group)
  cts <- cts[filter,]
  cts <- as.matrix(cts)
  ## create coldata
  set <- newSeqExpressionSet(as.matrix(cts),
                             phenoData = coldata)
  ## normalized data
  set_N <- betweenLaneNormalization(set, which="upper",offset=F)
  set_RUV <- RUVg(set_N, gene_ref, k=1)
  str1 <- paste("DESeqDataSetFromMatrix(countData = counts(set_RUV),colData= pData(set_RUV), design = ~", design, ")",sep = "")
  dds <- eval(parse( text=str1))
  dds <- DESeq(dds,minReplicatesForReplace=7, parallel = parallel,fitType = fitType)
  return(dds)
}


# 2. Deseq2 function
##---- Function 2
create_dds_test <- function(data_Ruv,design1,design2){
  str1 <- paste("DESeqDataSetFromMatrix(countData = counts(data_Ruv),colData= pData(data_Ruv), design = ~", design1, ")",sep = "")
  dds <- eval(parse( text=str1))
  dds <- DESeq(dds,minReplicatesForReplace=7, parallel = T)
  dds <- replaceOutliersWithTrimmedMean(dds)
  str2 <- paste("DESeqDataSetFromMatrix(countData = counts(dds),colData= pData(data_Ruv), design = ~", design2,")",sep = "")
  dds <- eval(parse( text=str2))
  dds <- DESeq(dds,minReplicatesForReplace=7, parallel = T)
  return(dds)
}



anno_gene_symbol <- function(data,data_anno,export, name=NULL){
  anno_DESeq2 <- data_anno[,c("gene_id", "gene_symbol")]
  anno_DESeq2 <- subset(anno_DESeq2, gene_id %in% rownames(data))
  rownames(anno_DESeq2) <- anno_DESeq2$gene_id
  anno_DESeq2 <- anno_DESeq2[rownames(data),]
  rownames(data) <- anno_DESeq2$gene_symbol
  if (export==T){
    data_ex <- as.data.frame(data)
    write.csv(data_ex,paste(name,".csv",sep = ""))
  }
  return(data)
}



# 3. create_vsd_t data 
##---- Function 3
create_vsd_t <- function(data_dds,data_coldata, blind=F){
  if (blind==F) {
    vsd <- varianceStabilizingTransformation(data_dds, blind = F)
  } else  vsd <- varianceStabilizingTransformation(data_dds, blind = T)
  vsd_t <- t(assay(vsd))
  col_name <- colnames(vsd_t)
  vsd_t <- as.data.frame(vsd_t)
  colnames(vsd_t) <- col_name
  vsd_t <- cbind(data_coldata,vsd_t)
  col_name <- colnames(vsd_t)
  vsd_t <- as.data.frame(vsd_t)
  colnames(vsd_t) <- col_name
  return(vsd_t)
}

### create vst_ratio data 
create_vst_ratio <- function(genelist,vst_data, time1=0,time,group=NULL,subtract="no",log="yes"){
  vst_data_D0 <- subset(vst_data, Time==time1)
  if(is.null(group)){
    group= colnames(vst_data_D0)[1:which(colnames(vst_data_D0)=="sizeFactor")]
    vst_data_D0_clinical <- vst_data_D0[,c(group)]
  }
  if(!is.null(group)){
    vst_data_D0_clinical <- vst_data_D0[,c("Patient_ID",group,"Time")]
  }
  vst_data_D0 <- vst_data_D0[,genelist]
  vst_data_D <- subset(vst_data, Time==time)
  vst_data_D <- vst_data_D[,genelist]
  vst_data_ratio <- as.matrix(vst_data_D)/as.matrix(vst_data_D0)
  if(subtract=="yes"){
    vst_data_ratio <- as.matrix(vst_data_D) - as.matrix(vst_data_D0)
  }
  if (log=="no"){
    vst_data_ratio <- 2^vst_data_ratio
  } 
  vst_data_ratio <- data.frame(vst_data_ratio)
  colnames(vst_data_ratio) <- colnames(vst_data_D0)
  vst_data_ratio <- cbind(vst_data_D0_clinical,vst_data_ratio)
  return(vst_data_ratio)
}



# 4. create vsd_t_plot data for group of genes
create_vsd_t_plot <- function(X,data_vsd_t){
  plot_var <- c("Patient_ID", "Time", "TBMDEX", X)
  vsd_t_plot <- data_vsd_t[,plot_var]
  return(vsd_t_plot)
} 


# 5. featurePlot for multiple genes
##---- Function 5
Function_featurePlot <- function(X,Y,data){
  library(caret)
  featurePlot(x = data[,X], 
              y = Y,
              plot = "box",
              scales = list(x = list(relation="free"), 
                            y = list(relation="free")), 
              pch = "|")
}



# 6. Plot_change of count of interest genes
##---- Function 6
Plot_change <- function(X,data){
  for (i in X) {
    gene <- i
    data_plot <- data[,c("Patient_ID", "Time", "TBMDEX",gene)]
    data_plot$TBMDEX <- revalue(data_plot$TBMDEX, c("Placebo"="P","Dexamethasone"="Dex"))
    colnames(data_plot)[ncol(data_plot)] <- "count"
    P <- ggplot(data_plot, aes(x=as.numeric(Time), y=count,colour=TBMDEX,group=Patient_ID)) +
      geom_line(lwd=1) +
      labs(title = gene, x="TimePoint", y="Log2 count")
    print(P)
  }
}



# 7. create_ttest_table
##---- Function 7
create_ttest_table <- function(data_vsd_t, group){
  data_vsd_t$group <- droplevels(group)
  res <- compareGroups(group ~ ., data = data_vsd_t, method = 4)
  #plot(res[c(10, 4)], type = "png")
  tab <- createTable(res, show.n=T)
  pvals <- getResults(tab, "p.overall")
  tab_data <- as.data.frame(tab$descr)
  tab_data$padj_BH <- round(p.adjust(pvals, method = "BH"),5)
  tab_data <- tab_data[order(tab_data$padj_BH),]
  #write.csv(tab_data,"T.test_Dex_genes_D14.csv")
  export2md(tab)
}

# 8. draw_pheatmap
##---- Function 8
draw_pheatmap_matrix_T <- function(data,group, subgroup=NULL,genelist,cluster_cols=T,
                                   cluster_rows=T,scale="row", kmeans_k=NA,show_rownames=F,
                                   cutree=NA,fontsize_row=4,export=F,filename=NULL)
  {
  library(gplots)
  library(pheatmap)
  str <- paste0("data$",group)
  group <-  eval(parse( text=str))
  df_col <- as.data.frame(group, row.names=rownames(data))
  if (!is.null(subgroup)) {
    str <- paste0("data$",subgroup)
    subgroup <-  eval(parse( text=str))
    
    
  }
  data <- data[,genelist]
  data <- t(data)
  mycol <- colorpanel(1000,"blue","white","red")
  heat_map <- pheatmap(data, scale=scale,cluster_cols =cluster_cols,cluster_rows = cluster_rows,
                       show_rownames = show_rownames, kmeans_k=kmeans_k,color= mycol,
                       show_colnames = F,annotation_col = df_col,clustering_method="complete",
                       height=15,fontsize = 15,fontsize_row=fontsize_row,width=15,cutree_rows = cutree,border_color =NA)
  print(heat_map)
  if (export) save_pheatmap_pdf(heat_map, filename)
  return(heat_map)
}


save_pheatmap_pdf <- function(x, filename, width=700, height=700) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  tiff(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

draw_pheatmap <- function (res, dds,n_gene,p,L2FC,group,scale="row",
                           cluster_cols=T,cluster_rows=T,kmeans_k=NA,cutree=NA,
                           show_rownames=F,border_color=NA
                           ,export=F,filename=NULL){
  library(gplots)
  library(pheatmap)
  resOrdered <- subset(res, padj <p & abs(log2FoldChange)>=L2FC)
  resOrdered <- resOrdered[order(resOrdered$padj),]
  select_genes <- rownames(resOrdered)[1:n_gene]
  vst_dds <- varianceStabilizingTransformation(dds, blind=F)
  mat <-  assay(vst_dds)[select_genes, ]
  mycol <- colorpanel(1000,"blue","white","red")
  df <- colData(dds)
  df_col <- as.data.frame(df[,group], row.names=rownames(df))
  colnames(df_col)[1] <- "group"
  df_col <- df_col[]
  df_col <- df_col[order(df_col[, 1]), , drop = FALSE]
  mat <- mat[,rownames(df_col)]
  heat_map <- pheatmap(mat, scale=scale,cluster_cols =cluster_cols,cluster_rows=cluster_rows,
                       show_rownames = show_rownames,
                       kmeans_k=kmeans_k,color= mycol,show_colnames = F,
                       annotation_col = df_col,clustering_method="complete",
                       height=14,fontsize = 12,width=14, cutree_rows = cutree,
                       border_color=border_color)
  #print(heat_map)
  if (export) save_pheatmap_pdf(heat_map, filename)
  return(heat_map)
}



Draw_pheatmap.type <- function (Data, annCol, type = colnames(annCol)[1], 
                                conditions = "Auto",cluster=T,...) 
{
  library(gplots)
  library(pheatmap)
  mycol <- colorpanel(1000,"blue","white","red")
  Data <- t(Data)
  res <- list()
  annCol <- annCol[, type, drop = FALSE]
  if (is.null(rownames(annCol))) 
    stop("annCol must have row names!")
  if (any(!rownames(annCol) %in% colnames(Data))) 
    stop("annCol has col that are not present in col of Data!")
  annCol <- annCol[order(annCol[, 1]), , drop = FALSE]
  samplesOriginalOrder <- colnames(Data)
  if (conditions[1] == "Auto") 
    conditions <- unique(as.character(annCol[, 1]))
  if (any(!conditions %in% unique(as.character(annCol[, 1])))) {
    warning("Some of the conditions are not in annCol.")
    conditions <- intersect(conditions, unique(as.character(annCol[, 
                                                                   1])))
  }
  pheatmapS <- list()
  dataPlot <- c()
  ann1 <- c()
  #cond = "Dexamethasone"
  for (cond in conditions) {
    condSamples <- rownames(annCol)[which(annCol == cond)]
    if (length(condSamples) > 1) {
      pa <- pheatmap(Data[, condSamples, drop = FALSE], 
                     cluster_rows = F, scale = "column",silent = T)
      pheatmapS[[as.character(cond)]] <- pa
      o2 <- pa$tree_col$order
      #print(o2)
    }
    else {
      o2 <- 1
    }
    dataPlot <- cbind(dataPlot, Data[, condSamples[o2], drop = FALSE])
    ann1 <- rbind(ann1, annCol[condSamples[o2], , drop = FALSE])
  }
  pAll <- pheatmap(dataPlot, annotation_col = ann1, cluster_cols = F,
                   show_rownames=FALSE, show_colnames=F,scale = "row",color= mycol,
                   ...)
  res[["pheatmapS"]] <- pheatmapS
  res[["pheat"]] <- pAll
  res[["ordering"]] <- match(rownames(dataPlot), samplesOriginalOrder)
  res[["annColAll"]] <- ann1
  invisible(res)
}

## how to use
#vst_data_matrix <- vst_data[,DE_gene_pvsFC]
#vst_group <- vst_data[,1:17]

#Draw_pheatmap.type(Data = vst_data_matrix,annCol=vst_group,type="trial_number")


## return cluster by k
# mtcars.clust <- cbind(mtcars, 
#                       cluster = cutree(res$tree_row, 
#                                        k = 10))

draw_pheatmap_pathway <- function(dds_anno,class_goup,DE_gene,GOGSE,gene_set,cluster_cols=T, method="gsva"){
  library(GSVA)
  if (is.data.frame(dds_anno)) {
    group <- dds_anno %>%
      dplyr::select(class_goup)
    vst_data <- dds_anno[,DE_gene]
    vst_data_gage <- vst_data %>% 
      t() %>% as.data.frame()
    } 
  else if(isClass(dds_anno)){
    group <- as.data.frame(colData(dds_anno))
    group <- group %>%
      dplyr::select(class_goup)
    vst_data <- assay(varianceStabilizingTransformation(dds_anno,blind = F))
    vst_data_gage <- vst_data[which(rownames(vst_data) %in% DE_gene),]
    vst_data_gage <- as.data.frame(vst_data_gage)
  }
  
  vst_data_gage$ENTREZID <-  mapIds(x = org.Hs.eg.db,
                                    keys=row.names(vst_data_gage),
                                    column="ENTREZID",
                                    keytype="SYMBOL",
                                    multiVals="first",
                                    ifnotfound = NA)
  vst_data_gage <- vst_data_gage[!is.na(vst_data_gage$ENTREZID),]
  rownames(vst_data_gage) <- vst_data_gage$ENTREZID
  vst_data_gage <- vst_data_gage %>%
    dplyr::select(-ENTREZID) %>%
    as.matrix()
  
  gsva_results <- gsva(
    vst_data_gage,
    gene_set,
    parallel.sz=1L,
    # Appropriate for our vst transformed data  
    kcdf = "Gaussian",
    # Minimum gene set size
    min.sz = 1,
    # Maximum gene set size
    max.sz = Inf,
    # Compute Gaussian-distributed scores
    mx.diff = TRUE,
    method=method,
    # Don't print out the progress bar
    verbose = F
  )
  gsva_results_d <- as.data.frame(gsva_results)
  gsva_results_d$ID <- substr(row.names(gsva_results_d),1,10)
  if (isClass(GOGSE)) {
    gsva_results_H <- subset(gsva_results_d, ID %in% GOGSE@result$ID)
  }
  if (is.vector(GOGSE)) {
    gsva_results_H <- subset(gsva_results_d, ID %in% GOGSE)
  }
  gsva_results_H <- gsva_results_H %>%
    dplyr::select(-ID) %>%
    as.matrix()
  pathway_heatmap <- pheatmap::pheatmap(gsva_results_H,
                                        annotation_col = group, # Add metadata labels!
                                        show_colnames = F, # Don't show sample labels
                                        fontsize_row = 6,
                                        cluster_cols = cluster_cols,
                                        scale="row"
                                        # Shrink the pathway labels a tad
  )
  gsva_data <- list()
  gsva_data$gsva_results <- gsva_results
  gsva_data$pathway_heatmap <- pathway_heatmap
  return(gsva_data)
}

# gene_set = go.bp.gs
# 
# names(gene_set)
# gene_set <-  subset(gene_set, substr(names(gene_set),1,10) %in% temp_GO_ID)
# names(gene_set)
# GOGSE@result$ID

### calculate score using all background genes
draw_pheatmap_pathway_All <- function(dds_anno,class_goup,DE_gene,GOGSE,gene_set,cluster_cols=T){
  library(GSVA)
  if (is.data.frame(dds_anno)) {
    group <- dds_anno %>%
      dplyr::select(class_goup)
    vst_data <- dds_anno[,DE_gene]
    vst_data_gage <- vst_data %>% 
      t() %>% as.data.frame()
  } 
  else if(isClass(dds_anno)){
    group <- as.data.frame(colData(dds_anno))
    group <- group %>%
      dplyr::select(class_goup)
    vst_data <- assay(varianceStabilizingTransformation(dds_anno,blind = F))
    vst_data_gage <- vst_data[which(rownames(vst_data) %in% DE_gene),]
    vst_data_gage <- as.data.frame(vst_data_gage)
  }
  
  vst_data_gage$ENTREZID <-  mapIds(x = org.Hs.eg.db,
                                    keys=row.names(vst_data_gage),
                                    column="ENTREZID",
                                    keytype="SYMBOL",
                                    multiVals="first",
                                    ifnotfound = NA)
  vst_data_gage <- vst_data_gage[!is.na(vst_data_gage$ENTREZID),]
  rownames(vst_data_gage) <- vst_data_gage$ENTREZID
  vst_data_gage <- vst_data_gage %>%
    dplyr::select(-ENTREZID) %>%
    as.matrix()
  
  if (isClass(GOGSE)) {
    temp_GO_ID <- GOGSE@result$ID
  }
  if (is.vector(GOGSE)) {
    temp_GO_ID <- GOGSE
  }
  gene_set <-  subset(gene_set, substr(names(gene_set),1,10) %in% temp_GO_ID)
  gsva_results <- gsva(
    vst_data_gage,
    gene_set,
    parallel.sz=1L,
    # Appropriate for our vst transformed data  
    kcdf = "Gaussian",
    # Minimum gene set size
    min.sz = 1,
    # Maximum gene set size
    max.sz = Inf,
    # Compute Gaussian-distributed scores
    mx.diff = TRUE,
    method="gsva",
    # Don't print out the progress bar
    verbose = F
  )
  gsva_results_d <- as.data.frame(gsva_results)
  gsva_results_d$ID <- substr(row.names(gsva_results_d),1,10)
  if (isClass(GOGSE)) {
    gsva_results_H <- subset(gsva_results_d, ID %in% GOGSE@result$ID)
  }
  if (is.vector(GOGSE)) {
    gsva_results_H <- subset(gsva_results_d, ID %in% GOGSE)
  }
  gsva_results_H <- gsva_results_H %>%
    dplyr::select(-ID) %>%
    as.matrix()
  pathway_heatmap <- pheatmap::pheatmap(gsva_results_H,
                                        annotation_col = group, # Add metadata labels!
                                        show_colnames = F, # Don't show sample labels
                                        fontsize_row = 6,
                                        cluster_cols = cluster_cols,
                                        scale="row"
                                        # Shrink the pathway labels a tad
  )
  gsva_data <- list()
  gsva_data$gsva_results <- gsva_results
  gsva_data$pathway_heatmap <- pathway_heatmap
  gsva_data$gene_set <- gene_set
  return(gsva_data)
}



## Venn diagram
create_venn_data <- function(data_x,
                             main,
                             category.names =NULL,
                             output=F, filename=NULL){
  x = data_x
  if (length(x)==2){
    col=c("#440154ff", '#21908dff')
    fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3))
    cat.pos = c(10, 180)
    cat.dist = c(0.055, 0.055)
  }
  else if (length(x)==3){
    col=c("#440154ff", '#21908dff', '#fde725ff')
    fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3))
    cat.pos = c(-27, 27, 135)
    cat.dist = c(0.055, 0.055, 0.085)
    }
  else if (length(x)==4){
    col=c("#440154ff", '#21908dff', '#fde725ff', 'F1B6DA')
    fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3),alpha('#F1B6DA',0.3))
    cat.pos = c(-27, 27, 135, 1)
    cat.dist = c(0.055, 0.055, 0.085, 0.5)
    }
  t <- venn.diagram(
    x = x,
    main=main,
    main.cex = 1.3,
    main.col="blue4",
    category.names = category.names,
    filename = filename,
    output = output,
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    lwd = 1,
    col = col,
    fill = fill,
    cex = 0.8,
    fontfamily = "sans",
    cat.cex = 1,
    cat.default.pos = "outer",
    cat.fontfamily = "sans",
    cat.pos = cat.pos,
    cat.dist = cat.dist,
    scaled = F,
  )
  grid.newpage()
  grid.draw(t)
}




# 9. glimma_plot
##----------function 9

glimma_plot <- function(res,dds,group,name,folder=NULL){
  library(Glimma)
  res$status <- ifelse(is.na(res$padj)|res$padj>0.05,0,ifelse(res$log2FoldChange > 0, 1, -1))
  anno_res<- as.data.frame(rownames(res))
  colnames(anno_res) <- "gene_symbol"
  data_group <- as.data.frame(colData(dds))
  group = data_group[,which(colnames(data_group)==group)]
  glMDPlot(res,
           anno=anno_res,
           groups = group,
           counts = counts(dds,normalized=TRUE),
           transform = T,
           status=res$status,
           side.main = "gene_symbol",
           folder=folder,
           html = name)
}



## PCA plot for matrix data
create_PCA_matrix_plot <- function(data, genelist, group,title){
  str <- paste0("data$",group)
  group <-  eval(parse( text=str ))
  data_matrix <- data[,genelist]
  PCA_T <- prcomp(data_matrix, scale. = T )
  PCA <- data.frame(PCA_T$x)
  PCA_1 <- PCA[,1:5]
  PCA_1$ID <- rownames(PCA_1)
  PCA_1$group <- group
  s <- summary(PCA_T)
  #str(s)
  PVE <-s$importance[2,] * 100
  i1 <- round(PVE[1],2)
  i2 <- round(PVE[2],2)
  i3 <- round(PVE[3],2)
  library(ggplot2)
  p <- ggplot(PCA_1, aes(x = PC1, y = PC2, colour = group)) + 
    geom_point(size=3, alpha=0.95) +
    theme (axis.text.x  = element_text(size=16), axis.text.y = element_text(size=16)) + 
    theme(axis.title=element_text(size=16,face="bold"), plot.title = element_text(color="blue", size=14, face="bold")) + xlab(paste0("PC1 (", i1, "%)")) + ylab(paste0("PC2 (", i2, "%)")) + ggtitle(title)
  return(p)
  
  
}



create_plsda_plot <- function(data,genelist,group,title,type, plot=T){
  library(mixOmics)
  legend_title <- group
  str <- paste0("data$",group)
  group <-  eval(parse( text=str ))
  data_matrix <- data[,genelist]
  data_matrix <- as.matrix(data_matrix)
  if (type=="pca") {
    miOxmic_DEX <- mixOmics::pca(X=data_matrix, ncomp=4,center = T, scale = TRUE)
    if(plot==T){
      plotIndiv(miOxmic_DEX, group = group, ind.names = FALSE, 
                legend = TRUE,ellipse = TRUE, 
                legend.title=legend_title,title = title,guide =T)
      }
  }
  if(type=="spca"){
    miOxmic_DEX <- mixOmics::spca(X=data_matrix,ncomp = 3, keepX = c(10,5,5))
    if(plot==T){
      plotIndiv(miOxmic_DEX, legend=T, title=title, 
                group = group,
                ellipse=T, 
                legend.title=legend_title, 
                ind.names = F,ellipse.level=0.95)
      }
  }
  if (type== "plsda"){
    miOxmic_DEX <- mixOmics::plsda(X=data_matrix, Y= group, ncomp=5)
    if(plot==T){
      plotIndiv(miOxmic_DEX, legend=T, title=title, 
                ellipse=T, 
                legend.title=legend_title, 
                ind.names = F, star = F
                , ellipse.level=0.95)
    }

  }
  return(miOxmic_DEX)
}

download_pathview <- function(res,pathway.id,out.suffix="test",
                              species = "hsa",gene_score=1,kegg.dir = "."){
  DEGs=res
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
  pathview(gene.data = foldchanges, 
           pathway.id = pathway.id,
           kegg.dir = kegg.dir,
           species = species,
           out.suffix = out.suffix,
           limit = list(gene=gene_score, cpd=1))
}


create_distribution_pvalue_res <- function(data,genelist1,genelist2){
  All_De <- subset(data, padj<0.05)
  All_De <- as.data.frame(All_De)
  All_De$group <- "All_DE_Genes"
  Ico <- subset(All_De, 
                rownames(All_De) %in% genelist1)
  Ico$group <- "Eicosanoids"
  IF <- subset(All_De, 
              rownames(All_De) %in% genelist2)
  IF$group <- "Inflammatory"
  data_res_plot <- do.call(rbind,list(All_De,
                                      Ico,
                                      IF))
  data_res_plot$Log_padj <- -log10(data_res_plot$padj)
  p <- ggplot(data_res_plot, aes(x=padj, fill=group)) +
    geom_histogram( alpha=.5, position="identity",aes(y = ..count..)) + 
    facet_grid(group ~ .,scales="free_y")+
    ggtitle("adjusted P distribution")+
    xlab("P adjusted") +
    ylab("Density count")
  print(p)
  return(p)
}

create_topGO <- function(res,p.cutoff=0.05, LFC=0.57,signode=10,prefix="test"){
  library(org.Hs.eg.db)
  DEGs=subset(res,res$padj<p.cutoff)
  geneList <- DEGs$log2FoldChange
  names(geneList) <- rownames(DEGs)
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                geneSelectionFun = function(x)x >=LFC,
                nodeSize = 10,
                ID = "SYMBOL",
                annot = annFUN.org , mapping = "org.Hs.eg.db")
  resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes = signode, useInfo = 'all')
  printGraph(GOdata, resultKS.elim, firstSigNodes = signode, fn.prefix = prefix, useInfo = "all", pdfSW = TRUE)
  
  return(resultKS.elim)
}

# 
# library(future)
# plan(multiprocess, gc=TRUE)

subset_gsea<-function(x, list_G0)
{
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



Load_dds <- function(data_name,p.cutoff = 0.05, LFC.cutoff= 0.57,data_anno = anno,export = F,name = NULL){
  library(dplyr)
  dds <- get(load(paste0("RData/",data_name)))
  dds_anno <- anno_gene_symbol(data = dds,data_anno = anno,export = export,name = name)
  colnames(colData(dds_anno))[which(colnames(colData(dds_anno))=="LTA4H")] <- "LTA4H.rs"
  if ("TBMDEX" %in% colnames(colData(dds_anno))) {
    colnames(colData(dds_anno))[which(colnames(colData(dds_anno))=="TBMDEX")] <- "trial_arm"
  }
  return(dds_anno)
}
Load_res <- function(data_name,p.cutoff = 0.05, LFC.cutoff= 0.57,data_anno = anno,filter=T,export = F,name = NULL){
  res <- get(load(paste0("RData/",data_name)))
  res_anno <- anno_gene_symbol(data = res,data_anno = anno,export = export,name = name)
  if (filter==T) {
    res_anno <- subset(res_anno, padj<p.cutoff & abs(log2FoldChange)>=LFC.cutoff)
    res_anno <- res_anno[order(res_anno$pvalue),]
  }
  
  return(res_anno)
}

