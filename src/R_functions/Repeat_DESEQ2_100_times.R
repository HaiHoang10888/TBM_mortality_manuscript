

## create clinical and expression data
clinical_data <- create_TB_TBM_data(workdir = workdir,data_PTB = data_PTB,data_TBM =data_TBM ,output = "PTB_D0")

## input batch data
Batch_data <- read.csv("RData/clinical_data_All_TB_D0_batch_selection.csv")
rownames(Batch_data) <- Batch_data$LIMS_ID
Batch_data <- Batch_data[rownames(clinical_data),]

if (all(rownames(Batch_data)==rownames(clinical_data))) {
  clinical_data$Batch_Patient_100 <- Batch_data$Batch_Patient_100
}

Xpert_value <- read.csv("RData/28TB_29TB_meta_data_Xpert_value.csv", row.names = "LMS_ID", na.strings = "")

Xpert_value <- Xpert_value[row.names(clinical_data),]
if (all(rownames(Xpert_value)==rownames(clinical_data))) {
  clinical_data$GenXpert <- Xpert_value$GenXpert
  clinical_data$CT_mean <- Xpert_value$CT_mean
}

clinical_data <- subset(clinical_data, !is.na(CT_mean))






Res_result <- data.frame(Gene=character(0),baseMean = numeric(0), log2FoldChange = numeric(0),
                       lfcSE = numeric(0),stat = numeric(0),pvalue = numeric(0),
                       padj = numeric(0))




clinical_data$CT_mean
data_import = clinical_data


for (i in 1:100){
  #
  p.cutoff = 0.1; LFC.cutoff= 0.57
  print(i)
  Temp_index <- sample(1:nrow(data_import),nrow(data_import)/2)
  data_Dis <- data_import[Temp_index,]
  data_Val <- data_import[-Temp_index,]

  exp_Dis <- create_exp_data(count_data = coun_data_all, clinical_data = data_Dis)
  exp_Val <- create_exp_data(count_data = coun_data_all, clinical_data = data_Val)
  
  design1 <- group <- "CT_mean"
  design2 <- paste0("W_1 + ", group)
  
  dds_Dis <- create_dds_1times(count = exp_Dis,coldata = data_Dis,
                               group = group,number = 1,
                               filter = 5,genelist = HK_genes_11_list,
                               design1 = design1,design2 = design2)
  
  res_Dis <- results(dds_Dis, name="CT_mean",alpha = 0.05,parallel = T)
  res_anno_Dis <- anno_gene_symbol(data = res_Dis,data_anno = anno,export=F,name=NULL)
  res_anno_DE_Dis <- subset(res_anno_Dis, padj<p.cutoff)
  
  dds_Val <- create_dds_1times(count = exp_Val,coldata = data_Val,
                               group = group,number = 1,
                               filter = 5,genelist = HK_genes_11_list,
                               design1 = design1,design2 = design2)
  
  res_Val <- results(dds_Val, name="CT_mean",alpha = 0.05,parallel = T)
  res_anno_Val <- anno_gene_symbol(data = res_Val,data_anno = anno,export=F,name=NULL)
  res_anno_DE_Val <- subset(res_anno_Val, padj<p.cutoff)
  
  res_anno_DE_Up_Dis <- subset(res_anno_DE_Dis, log2FoldChange>0)
  res_anno_DE_Down_Dis <- subset(res_anno_DE_Dis, log2FoldChange<0)
  
  res_anno_DE_Up_Val <- subset(res_anno_DE_Val, log2FoldChange>0)
  res_anno_DE_Down_Val <- subset(res_anno_DE_Val, log2FoldChange<0)
  
  Repeat_Gene_Up <- intersect(rownames(res_anno_DE_Up_Dis),rownames(res_anno_DE_Up_Val))
  Repeat_Gene_Down <- intersect(rownames(res_anno_DE_Down_Dis),rownames(res_anno_DE_Down_Val))
  Repeat_Gene <- c(Repeat_Gene_Up,Repeat_Gene_Down)
  
  res_Sub_Dis <- as.matrix(res_anno_DE_Dis[Repeat_Gene,])
  res_sub_Val <- as.matrix(res_anno_DE_Val[Repeat_Gene,])
  
  res_temp <- (res_Sub_Dis + res_sub_Val)/2
  res_temp <- as.data.frame(res_temp)
  res_temp$Gene <- rownames(res_temp)
  if(length(rownames(res_temp))>0){
    rownames(res_temp) <- 1:nrow(res_temp)
    rownames(res_temp) <- 1:nrow(res_temp)
    res_temp <- res_temp[,move_column(names(res_temp), "Gene first")]
    Res_result <- rbind(Res_result,res_temp)
  }

}

getwd()
save(Res_result, file = "RData/6_PTB_CT_Xpert_Validation_100times_1.RData")
# dds_Dis_1 <- create_dds_1times(count = exp_Dis,coldata = Clinical_data_Dis,
#                       group = group,number = 1,
#                       filter = 5,genelist = HK_genes_11_list,
#                       design1 = design1,design2 = design2)
# 
# res_Dis_1 <- results(dds_Dis_1, contrast=c("TB_type", "TBM", "PTB"),alpha = 0.05,parallel = T)
# res_Dis <- results(dds_Dis, contrast=c("TB_type", "TBM", "PTB"),alpha = 0.05,parallel = T)
# 
# res_DE_Dis <- subset(res_Dis, padj<p.cutoff & abs(log2FoldChange)>=LFC.cutoff)
# res_DE_Dis_1 <- subset(res_Dis_1, padj<p.cutoff & abs(log2FoldChange)>=LFC.cutoff)
# summary(res_DE_Dis)
# summary(res_DE_Dis_1)



t1 <- sample(1:nrow(TBM_data),nrow(TBM_data)/2)
t2 <- sample(1:nrow(TBM_data),nrow(TBM_data)/2)

all(t1 %in% t2)
set.seed(1234)

for (i in 1:100){
  t1 <- sample(1:nrow(TBM_data),nrow(TBM_data)/2)
  t2 <- sample(1:nrow(TBM_data),nrow(TBM_data)/2)
  print(i)
  print(all(t1 %in% t2))
}

for (i in 1:100){
  #
  p.cutoff = 0.05; LFC.cutoff= 0.57
  print(i)
  TBM_index <- sample(1:nrow(TBM_data),nrow(TBM_data)/2)
  PTB_index <- sample(1:nrow(PTB_data),round(nrow(PTB_data)/2))
  TBM_Dis <- TBM_data[TBM_index,]
  TBM_Val <- TBM_data[-TBM_index,]
  PTB_Dis <- PTB_data[PTB_index,]
  PTB_Val <- PTB_data[-PTB_index,]
  
  Clinical_data_Dis <- rbind(TBM_Dis,PTB_Dis)
  Clinical_data_Dis$GeneXpert_Base <- factor(Clinical_data_Dis$GeneXpert_Base)
  
  Clinical_data_Val <- rbind(TBM_Val,PTB_Val)
  Clinical_data_Val$GeneXpert_Base <- factor(Clinical_data_Val$GeneXpert_Base)
  
  exp_Dis <- create_exp_data(count_data = coun_data_all, clinical_data = Clinical_data_Dis)
  exp_Val <- create_exp_data(count_data = coun_data_all, clinical_data = Clinical_data_Val)
  
  design1 <- group <- "GeneXpert_Base"
  design2 <- paste0("W_1 + ", group)
  
  dds_Dis <- create_dds_1times(count = exp_Dis,coldata = Clinical_data_Dis,
                    group = group,number = 1,
                    filter = 5,genelist = HK_genes_11_list,
                    design1 = design1,design2 = design2)
  
  res_Dis <- results(dds_Dis, contrast=c("GeneXpert_Base", "High", "Low"),alpha = 0.05,parallel = T)
  res_anno_Dis <- anno_gene_symbol(data = res_Dis,data_anno = anno,export=F,name=NULL)
  res_anno_DE_Dis <- subset(res_anno_Dis, padj<p.cutoff & abs(log2FoldChange)>=LFC.cutoff)
  
  dds_Val <- create_dds_1times(count = exp_Val,coldata = Clinical_data_Val,
                        group = group,number = 1,
                        filter = 5,genelist = HK_genes_11_list,
                        design1 = design1,design2 = design2)
  
  res_Val <- results(dds_Val, contrast=c("GeneXpert_Base", "High", "Low"),alpha = 0.05,parallel = T)
  res_anno_Val <- anno_gene_symbol(data = res_Val,data_anno = anno,export=F,name=NULL)
  res_anno_DE_Val <- subset(res_anno_Val, padj<p.cutoff & abs(log2FoldChange)>=LFC.cutoff)
  
  res_anno_DE_Up_Dis <- subset(res_anno_DE_Dis, log2FoldChange>0)
  res_anno_DE_Down_Dis <- subset(res_anno_DE_Dis, log2FoldChange<0)
  
  res_anno_DE_Up_Val <- subset(res_anno_DE_Val, log2FoldChange>0)
  res_anno_DE_Down_Val <- subset(res_anno_DE_Val, log2FoldChange<0)
  
  Repeat_Gene_Up <- intersect(rownames(res_anno_DE_Up_Dis),rownames(res_anno_DE_Up_Val))
  Repeat_Gene_Down <- intersect(rownames(res_anno_DE_Down_Dis),rownames(res_anno_DE_Down_Val))
  Repeat_Gene <- c(Repeat_Gene_Up,Repeat_Gene_Down)
  
  res_Sub_Dis <- as.matrix(res_anno_DE_Dis[Repeat_Gene,])
  res_sub_Val <- as.matrix(res_anno_DE_Val[Repeat_Gene,])
  
  res_temp <- (res_Sub_Dis + res_sub_Val)/2
  res_temp <- as.data.frame(res_temp)
  res_temp$Gene <- rownames(res_temp)
  rownames(res_temp) <- 1:nrow(res_temp)
  res_temp <- res_temp[,move_column(names(res_temp), "Gene first")]
  
  Res_result <- rbind(Res_result,res_temp)
  
  
}



save(Res_result, file = "RData/6_PTB_HighvsLow_Xpert_Validation_100times_1.RData")


nrow(Res_result) 

min(Res_result$padj)

Res_result$Gene[which(Res_result$padj==min(Res_result$padj))]

Res_result$Gene

Res_result_gen1 <- subset(Res_result,Gene=="TRPV5")

nrow(Res_result_gen1)
head(Res_result_gen1)
a = c(1,23,5,6)

b = NA

intersect(a,b)
res_Sub_Dis <- as.matrix(res_anno_DE[1:10,])
res_sub_Val <- as.matrix(res_anno_DE[1:10,])
res_temp <- (res_Sub_Dis + res_sub_Val)/2

res_anno_1 <- as.matrix(res_anno_DE)




rownames(res_anno_1)
colnames(res_anno_1)

clinical_data$trial_number
