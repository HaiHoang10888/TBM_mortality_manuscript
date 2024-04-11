## get library
library(knitr)
library(reshape)
library(Hmisc)
library(ggplot2)
library(gridExtra)
library(limma)
library(edgeR)
library(reshape)
library(DESeq2)
library(tidyr)
library(kableExtra)
library(dplyr)
library(plyr)
library(sva)
library(compareGroups)
library(expss)
library(VennDiagram)
library(RUVSeq)
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(GO.db)
library(DOSE)
library(pathview)
library("BiocParallel")
library(mixOmics)
library(DaMiRseq)
library(EDASeq)
library(enrichplot)
library(GOSemSim)
library(ReactomePA)

## ---------------Run on server-----------------------
# workdir='/home/ubuntu/tb_volume/Hai_Cpath_Analysis/RNA_output/PAXGENE_27TB/Hai_function/'
# source("Function_data_process_and_analysis_8MAR2021.R")
rm(list=ls())
workdir='/home/ubuntu/tb_volume/Hai_Cpath_Analysis/RNA_output/PAXGENE_27TB/count/count_old_Run1vsRun2_separate/Bootstrap_dex'
setwd(workdir)
source("bootstrap_and_vst_data_ratio_function.R")

#==============================================================================
# 1.1 Bootstrap random 1 without replace
#==============================================================================

## Bootrap all sample and create RData
#dds <- get(load("RData/dds_All_TBM_Dex_Feb_2022.RData"))
set.seed(1)
#load vst data all sample
dds <- get(load("RData/16ddsTC_D0_D14_AllTBM_19Mar21.RData"))

gene_list = rownames(dds)
length(gene_list)
vst_data <- create_vsd_t(data_dds = dds,data_coldata = colData(dds))

## create data D14-D0 and take log transform
#num_clinicalVar = which(colnames(vst_data)=="sizeFactor")
vst_data_ratio <- create_vst_ratio(genelist = gene_list,vst_data = vst_data,time = 14,group = "TBMDEX", subtract = "yes",log="yes")

n= nrow(vst_data_ratio)
## split data in to 90 discovery and 89 validation
## first time with set.seed(1), second time with set.seed(2)...

vst_dex <- subset(vst_data_ratio,TBMDEX=="Dexamethasone")
vst_pla <- subset(vst_data_ratio,TBMDEX=="Placebo")

index_dex_D <- sample(1:nrow(vst_dex), 50, replace=F)
index_pla_D <- sample(1:nrow(vst_pla), 40, replace=F)

vst_dex_D <- vst_dex[index_dex_D,]
vst_pla_D <- vst_pla[index_pla_D,]


index_dex_V <-(1:nrow(vst_dex))[!(1:nrow(vst_dex) %in% index_dex_D)]
index_pla_V <-(1:nrow(vst_pla))[!(1:nrow(vst_pla) %in% index_pla_D)]

B=1000
bbetas <- NULL
branks <- NULL
n_D <- nrow(vst_dex_D)
n_P <- nrow(vst_pla_D)
rate <- 0.632
m_D<-ceiling(n_D*rate)
m_P<-ceiling(n_P*rate)
replace=FALSE

for(i in 1:B){
  #set.seed(12345+i)
  print(i)
  vst_data_1 <- vst_dex_D[sample(n_D,m_D,replace=replace),]
  vst_data_2 <- vst_pla_D[sample(n_P,m_P,replace=F),]
  data_combine <- rbind(vst_data_1,vst_data_2)
  temp_beta <- delta_treatmentvscontrol(data_combine,gene_list = gene_list)
  temp_rank <- as.matrix(order(sort(abs(temp_beta),decreasing = T,index.return=T)$ix))
  rownames(temp_rank) <- rownames(temp_beta)
  bbetas <- cbind(bbetas,temp_beta)
  branks <- cbind(branks,temp_rank)
}

save(bbetas, file="bootstrap/bbetas_Discovery_D14_random_split_15Mar22.RData")
save(branks, file="bootstrap/branks_Discovery_D14_random_split_15Mar22.RData")


#==============================================================================
# 1.1 Bootstrap random 1 without replace
#==============================================================================

## Bootrap all sample and create RData
#dds <- get(load("RData/dds_All_TBM_Dex_Feb_2022.RData"))
set.seed(1)
#load vst data all sample
dds <- get(load("RData/16ddsTC_D0_D14_AllTBM_19Mar21.RData"))

gene_list = rownames(dds)
length(gene_list)
vst_data <- create_vsd_t(data_dds = dds,data_coldata = colData(dds))

## create data D14-D0 and take log transform
#num_clinicalVar = which(colnames(vst_data)=="sizeFactor")
vst_data_ratio <- create_vst_ratio(genelist = gene_list,vst_data = vst_data,time = 14,group = "TBMDEX", subtract = "yes",log="yes")

n= nrow(vst_data_ratio)
## split data in to 90 discovery and 89 validation
## first time with set.seed(1), second time with set.seed(2)...
set.seed(1)
vst_dex <- subset(vst_data_ratio,TBMDEX=="Dexamethasone")
vst_pla <- subset(vst_data_ratio,TBMDEX=="Placebo")

index_dex_D <- sample(1:nrow(vst_dex), 50, replace=F)
index_pla_D <- sample(1:nrow(vst_pla), 40, replace=F)

vst_dex_D <- vst_dex[index_dex_D,]
vst_pla_D <- vst_pla[index_pla_D,]


index_dex_V <-(1:nrow(vst_dex))[!(1:nrow(vst_dex) %in% index_dex_D)]
index_pla_V <-(1:nrow(vst_pla))[!(1:nrow(vst_pla) %in% index_pla_D)]

B=1000
bbetas <- NULL
branks <- NULL
n_D <- nrow(vst_dex_D)
n_P <- nrow(vst_pla_D)
rate <- 0.632
m_D<-ceiling(n_D*rate)
m_P<-ceiling(n_P*rate)
replace=FALSE

for(i in 1:B){
  #set.seed(12345+i)
  print(i)
  vst_data_1 <- vst_dex_D[sample(n_D,m_D,replace=replace),]
  vst_data_2 <- vst_pla_D[sample(n_P,m_P,replace=F),]
  data_combine <- rbind(vst_data_1,vst_data_2)
  temp_beta <- delta_treatmentvscontrol(data_combine,gene_list = gene_list)
  temp_rank <- as.matrix(order(sort(abs(temp_beta),decreasing = T,index.return=T)$ix))
  rownames(temp_rank) <- rownames(temp_beta)
  bbetas <- cbind(bbetas,temp_beta)
  branks <- cbind(branks,temp_rank)
}

save(bbetas, file="bootstrap/bbetas_Discovery_D14_random_split_15Mar22.RData")
save(branks, file="bootstrap/branks_Discovery_D14_random_split_15Mar22.RData")



## Bootstrap Discovery
## split samples at baseline discovery n=100, validation n=107
dds <- get(load("RData/1ddsTC_D0_D14_DEX_NoLRT_discovery_14Apr21.RData"))
gene_list = rownames(dds)
vst_data <- create_vsd_t(data_dds = dds,data_coldata = colData(dds))

colnames(vst_data)[1] <- "TBMDEX"
num_clinicalVar = which(colnames(vst_data)=="sizeFactor")
vst_data_ratio <- create_vst_ratio(genelist = gene_list,vst_data = vst_data,time = 14,group = "TBMDEX", subtract = "yes",log="yes")
n= nrow(vst_data_ratio)

vst_data_dex <- subset(vst_data_ratio,TBMDEX=="Dexamethasone")
vst_data_pla <- subset(vst_data_ratio,TBMDEX=="Placebo")

B=1000
bbetas <- NULL
branks <- NULL
n_D <- nrow(vst_data_dex)
n_P <- nrow(vst_data_pla)
rate <- 1 #0.632
m_D<-ceiling(n_D*rate)
m_P<-ceiling(n_P*rate)
replace = TRUE

set.seed(1)
for(i in 1:B){
  #set.seed(12345+i)
  print(i)
  vst_data_1 <- vst_data_dex[sample(nrow(vst_data_dex),m_D,replace=replace),]
  vst_data_2 <- vst_data_pla[sample(nrow(vst_data_pla),m_P,replace=replace),]
  data_combine <- rbind(vst_data_1,vst_data_2)
  temp_beta <- delta_treatmentvscontrol(data_combine,gene_list = gene_list)
  temp_rank <- as.matrix(order(sort(abs(temp_beta),decreasing = T,index.return=T)$ix))
  rownames(temp_rank) <- rownames(temp_beta)
  bbetas <- cbind(bbetas,temp_beta)
  branks <- cbind(branks,temp_rank)
}
save(bbetas, file="bootstrap/bbetas_Discovery_D14_predefined_split_withreplace_15Mar22.RData")
save(branks, file="bootstrap/branks_Discovery_D14_predefined_split_withreplace_15Mar22.RData")



###-----------------change to private laptop---------------------------------------------------
rm(list=ls())
workdir='D:/RNA_output/PAXGENE 27TB/count/'
output_dir = "count_old_Run1vsRun2_separate/Summary_04Mar21_TBM_DEX_Discovery&Validate"
setwd(paste0(workdir))
source("C:/Users/haiht/Dropbox/Project2020/WholeGenomeSeqHuman27TB/Hai_function/bootstrap_and_vst_data_ratio_function.R")
source("C:/Users/haiht/Dropbox/Project2020/WholeGenomeSeqHuman27TB/Hai_function/Function_data_process_and_analysis_8MAR2021.R")
set.seed(12345)
# create annotation object
anno <- create_anno_data(workdir=workdir)
# create list 11 houseskeeping genes
HK_genes_11_list <- create_HK_genes_11_list(anno=anno)

Antidrug_used <- read.csv("Anti_Drug_before_Random_19Mar2021.csv")

## condition to create clinical data
cond_out <- c("TBM_D0", "PTB_D0", "TB_D0", 
              "All_TBM_D0", "All_TB_D0",
              "TBM_D14", "TBM_D0_D14", "TBM_D60",
              "TBM_D0_D60", "TBM_D0_D14_D60", "TBM_D14_D60")

## set file data for PTB and TBM
data_PTB <- 'Meta_data_PTB_8MAR2021.RData'
data_TBM <- 'TBM_ALL_mergedbyrbind.RData_13MAR2021.RData'
setwd(paste0(workdir,output_dir))

dds <- get(load("RData/1ddsTC_D0_D14_DEX_NoLRT_discovery_14Apr21.RData"))
res <- get(load("RData/1resTC_D0_D14_DEX_NoLRT_discovery_14Apr21.RData"))
res_anno <- anno_gene_symbol(data = res,data_anno = anno,export=F,name=NULL)
p.cutoff = 0.05; LFC.cutoff= 0.57
res_DE <- subset(res, padj<p.cutoff & abs(log2FoldChange)>=LFC.cutoff)


head(res_DE)

summary(res_DE)
common_DE <- rownames(res_DE)

## -------------------First approach----------------
res_dis <- get(load("RData/1resTC_D0_D14_DEX_NoLRT_discovery_14Apr21.RData"))
#res_anno_dis <- anno_gene_symbol(data = res_dis,data_anno = anno,export=F,name=NULL)
#subset for DE gene discovery
p.cutoff = 0.05; LFC.cutoff= 0.57
res_DE_dis <- subset(res_dis, padj<p.cutoff & abs(log2FoldChange)>=LFC.cutoff)
summary(res_DE_dis)


res_val <- get(load("RData/2resTC_D0_D14_DEX_NoLRT_validate_14Apr21.RData"))
#res_anno_val <- anno_gene_symbol(data = res_val,data_anno = anno,export=F,name=NULL)
res_DE_val_O <- subset(res_val, padj<p.cutoff & abs(log2FoldChange)>=LFC.cutoff)



## subset for common gene in validation res
common_gene <- intersect(rownames(res_DE_dis),rownames(res_val))
length(common_gene)
not_common <- common_gene[which(rownames(res_DE_dis) %nin% common_gene)]
res_DE_dis[not_common,]

## subset for validation DE
res_DE_val <- res_val[common_gene,]
res_DE_val$padjBY <- p.adjust(res_DE_val$pvalue,method ="BY")


## produce vst ratio data for validation and perform lm for calculating p value
dds_val <- get(load("RData/2ddsTC_D0_D14_DEX_NoLRT_validate_14Apr21.RData"))
vst_data_val <- create_vsd_t(data_dds = dds_val,data_coldata = colData(dds_val))
colnames(vst_data_val)[1]  <-"TBMDEX" 


## create ratio vst data for D14 - D0
vst_data_ratio <- create_vst_ratio(genelist = rownames(res_DE_val),vst_data = vst_data_val,time = 14,group = "TBMDEX", subtract = "yes",log="yes")
vst_data_lm <- vst_data_ratio

#
library(stringr)
gene_index <- colnames(vst_data_lm)[1:50]
#i=1
res_DE_val$lm_EF <- NA
res_DE_val$p_lm <- NA
for(i in 1:length(rownames(res_DE_val))){
  print(i)
  temp_gene <- rownames(res_DE_val)[i]
  selected_var <- c("TBMDEX",temp_gene)
  temp_vst_data <- vst_data_lm[,selected_var]
  colnames(temp_vst_data)[2] <- "x1"
  fit <- lm(x1~TBMDEX,data = temp_vst_data)
  t <- summary(fit)$coefficients
  res_DE_val$lm_EF[i] <- t[2,1]
  res_DE_val$p_lm[i] <- t[2,4]
}

## adjust p linear regression by BY
res_DE_val$padj_lm <- p.adjust(res_DE_val$p_lm,method ="BY")
res_DE_val_BY <- subset(res_DE_val,padjBY<p.cutoff & abs(log2FoldChange)>=LFC.cutoff)
summary(res_DE_val_BY)
nrow(res_DE_val_BY)

res_DE_val_lm <- subset(res_DE_val,padj_lm<p.cutoff & abs(log2FoldChange)>=LFC.cutoff)
summary(res_DE_val_lm)

nrow(res_DE_val_lm)

res_DE_val_lm

### 125 gene in old Analysis
all_genes <- list("Discovery" = rownames(res_DE_dis), "Validation" = rownames(res_DE_val_O))
jpeg("bootstrap/Figure/Venn_dis_vs_val_DEgene.jpeg",units = "in", width = 5, height = 5, res = 300)
create_venn_data(all_genes,main="Repeated DE Genes in Discovery and validation",category.names = names(all_genes))
dev.off()

### DE gene in New Analysis by Deseq2
all_genes <- list("Discovery" = rownames(res_DE_dis), "Validation" = rownames(res_DE_val_BY))
jpeg("bootstrap/Figure/Venn_dis_vs_val_DEgene_deseq2.jpeg",units = "in", width = 5, height = 5, res = 300)
create_venn_data(all_genes,main="Repeated DE Genes in Discovery and validation",category.names = names(all_genes))
dev.off()

### DE gene in New Analysis by lm
all_genes <- list("Discovery" = rownames(res_DE_dis), "Validation" = rownames(res_DE_val_lm))
jpeg("bootstrap/Figure/Venn_dis_vs_val_DEgene_lm.jpeg",units = "in", width = 5, height = 5, res = 300)
create_venn_data(all_genes,main="Repeated DE Genes in Discovery and validation",category.names = names(all_genes))
dev.off()

### DE gene in New Analysis between lm and deseq2
all_genes <- list("Deseq2" = rownames(res_DE_val_BY), "Linear regression" = rownames(res_DE_val_lm))
jpeg("bootstrap/Figure/Venn_dis_vs_val_DEgene_lmvsDeseq2.jpeg",units = "in", width = 5, height = 5, res = 300)
create_venn_data(all_genes,main="Repeated DE Genes in Discovery and validation",category.names = names(all_genes))
dev.off()


Common_125 <- intersect(rownames(res_DE_dis), rownames(res_DE_val_O))
### DE gene in New Analysis between lm and deseq2
all_genes <- list("New analysis" = rownames(res_DE_val_BY), "Old analysis" = Common_125)
jpeg("bootstrap/Figure/Venn_New_vs_Old_DEgene.jpeg",units = "in", width = 5, height = 5, res = 300)
create_venn_data(all_genes,main="Repeated DE Genes in New and Old analysis",category.names = names(all_genes))
dev.off()


### histogram basemean 136 DE gene
jpeg("bootstrap/Figure/hist_136_DE_basemean.jpeg",units = "in", width = 5, height = 5, res = 300)
hist(log2(res_DE_val_BY$baseMean), xlab="Log2 Basemean", main = " histogram basemean")
abline(v=4,col="red")
dev.off()



### ===================second approach Ranking boot strap==================
rm(list=ls()) ## then run function and data input again
p.cutoff = 0.05; LFC.cutoff= 0.57
dds <- get(load("RData/1ddsTC_D0_D14_DEX_NoLRT_discovery_14Apr21.RData"))
res <- get(load("RData/1resTC_D0_D14_DEX_NoLRT_discovery_14Apr21.RData"))
res_val <- get(load("RData/2resTC_D0_D14_DEX_NoLRT_validate_14Apr21.RData"))
dds_val <- get(load("RData/2ddsTC_D0_D14_DEX_NoLRT_validate_14Apr21.RData"))
vst_data_val <- create_vsd_t(data_dds = dds_val,data_coldata = colData(dds_val))
colnames(vst_data_val)[1]  <-"TBMDEX" 

## load bootstrap result
vst_data_delta <- get(load("bootstrap/vst_data_delta_Discovery_D14_21Feb22_N.RData"))
vst_data_delta_order <- get(load("bootstrap/vst_data_delta_order_Discovery_D14_21Feb22_N.RData"))

out_put_data <- data.frame(Gene_symbol=rownames(vst_data_delta), 
                           Mean_LFC=rowMeans(vst_data_delta),
                           L_CI_mean=rowQuantiles(vst_data_delta, probs = 0.25),
                           U_CI_mean=rowQuantiles(vst_data_delta, probs = 0.75),
                           Mean_order=rowMeans(vst_data_delta_order),
                           L_CI_order=rowQuantiles(vst_data_delta_order, probs = 0.25),
                           U_CI_order=rowQuantiles(vst_data_delta_order, probs = 0.75))
L_CI_cutoff <- nrow(out_put_data)-150
U_CI_cutoff <- nrow(out_put_data)-100

out_put_data$Top_rank <- with(out_put_data, 
                                 ifelse(L_CI_order>=L_CI_cutoff & U_CI_order>=L_CI_cutoff,"Yes","No"))
out_put_data$FC_direction <- with(out_put_data, 
                                 ifelse(Mean_LFC>0, "up-regulated","down-regulated"))
out_put_data$significant <- with(out_put_data, 
                                  ifelse(sign(L_CI_mean)==sign(U_CI_mean), "sig","non-sig"))

out_put_data_DE <- subset(out_put_data, significant=="sig" & abs(Mean_LFC) >0.57)
out_put_data_toprank <- subset(out_put_data, Top_rank=="Yes")
nrow(out_put_data_DE)
nrow(out_put_data_toprank)
##summary DE genes
table(out_put_data_DE$FC_direction)
table(out_put_data_toprank$FC_direction)



## merge basemean to data
comon_gene <- rownames(res)[which(rownames(res) %in% rownames(out_put_data))]
out_put_data_S <- out_put_data[comon_gene,]
res <- res[comon_gene,]
nrow(res)

if(all(rownames(res)==rownames(out_put_data_S))){
  out_put_data_S$baseMean <- res$baseMean
}
out_put_data_S$Log_baseMean <- log2(out_put_data_S$baseMean)
out_put_data_S$ABS_Mean_LFC <- abs(out_put_data_S$Mean_LFC)

## draw scater plot to see relationship between LFC and basemean
jpeg("bootstrap/Figure/LFCvsBasemean.jpeg",units = "in", width = 5, height = 5, res = 300)
ggplot(out_put_data_S, aes(x=Log_baseMean, y=ABS_Mean_LFC)) + 
  geom_point()+
  geom_smooth()+
  xlab("Log2 Mean normalized count")+
  ylab("Absolute LFC")

dev.off()

all(rownames(out_put_data_DE) %in% rownames(res_val))

### Subset for significant gene or toprank
res_DE_val <- res_val[rownames(out_put_data_DE),]
res_toprank_val <- res_val[rownames(out_put_data_toprank),]

## adjust p linear regression by BY
res_DE_val$padjBY <- p.adjust(res_DE_val$padj,method ="BY")
res_toprank_val$padjBY <- p.adjust(res_toprank_val$padj,method ="BY")

colnames(res_DE_val)

res_DE_val_BY <- subset(res_DE_val,padjBY<p.cutoff )

summary(res_DE_val)
nrow(res_DE_val)
table(res_DE_val$)

res_DE_val_FC_CI <- subset(res_DE_val,padjBY < 0.05)
summary(res_DE_val_FC_CI)


all_genes <- list("Up-regulated" = up_gene, "Down-regulated" = down_gene)
create_venn_data(all_genes,main="boostrap analysis",category.names = names(all_genes))


out_put_data_rank_DE <- subset(out_put_data,Top_rank=="Yes")
res_DE_val_rank_DE <- res_DE_val[rownames(out_put_data_rank_DE),]






















gene_3 <- c(rownames(anno)[which(anno$gene_symbol=="LCN2")],
            rownames(anno)[which(anno$gene_symbol=="OLFM4")],
            rownames(anno)[which(anno$gene_symbol=="MARCO")])

vst_data_ratio_X <- vst_data_ratio[,gene_3]

write.csv(vst_data_ratio_X, "3_gene_death_substract.csv")


data_plot <- read.csv("3_gene_death_ratio_addclinical.csv")
#data_plot <- read.csv("3_gene_death_substract_addclinical.csv")


colnames(data_plot)
wilcox.test(data_plot$DEATH_EVE_3M_o,data_plot$OLFM4)
wilcox.test(OLFM4 ~ DEATH_EVE_3M_o, data = data_plot,
            exact = FALSE)


data_plot$TBMDEX <- factor(data_plot$TBMDEX,
                              levels = c("Placebo1","Dex1",
                                         "Placebo2","Dex2",
                                         "Placebo3","Dex3"))

data_plot$DEATH_EVE_3M <- factor(data_plot$DEATH_EVE_3M,
                                 levels = c("No1","Yes1",
                                            "No2","Yes2",
                                            "No3","Yes3"))


colnames(data_plot)

jpeg(paste0("Figure_high_reso/Expression3genes.jpeg"), units="in", width=8, height=8, res=300)

boxplot(data_plot$Gene ~ data_plot$TBMDEX , xlab="",
        ylab = expression(Log^2 ~ (expression ~ D14/D0)))

# Add data points
stripchart(data_plot$Gene ~ data_plot$TBMDEX, vertical = TRUE, method = "jitter",
           pch = 19, add = TRUE, col =c("black","blue") )
dev.off()


for(i in c("LCN2","OLFM4","MARCO")){
  gene=i
  group <- "TBMDEX"
  var <- c(gene,group)
  data <- data_plot[,var] 
  colnames(data)  <- c("x","group")
  
  # Basic boxplot
  jpeg(paste0("Figure_high_reso/",gene,"_",group,".jpeg"), units="in", width=5, height=5, res=300)
  
  boxplot(data$x ~ data$group , xlab="",
          ylab = expression(Log^2 ~ (expression ~ D14/D0)))
  
  # Add data points
  stripchart(data$x ~ data$group, vertical = TRUE, method = "jitter",
             pch = 19, add = TRUE, col =c("black","blue") )
  dev.off()
  
}



gene="LCN2"
group <- "DEATH_EVE_3M"
var <- c(gene,group)
data <- data_plot[,var] 
colnames(data)  <- c("x","group")

# Basic boxplot
jpeg(paste0("Figure_high_reso/",gene,"_",group,".jpeg"), units="in", width=5, height=5, res=300)

boxplot(data$x ~ data$group , xlab="",
        ylab = expression(Log^2 ~ (expression ~ D14/D0)))

# Add data points
stripchart(data$x ~ data$group, vertical = TRUE, method = "jitter",
           pch = 19, add = TRUE, col =c("black","blue") )
dev.off()



min(test)

vst_data_1_mean[which(rownames(vst_data_1_mean)=="ENSG00000276566")]
vst_data_2_mean[which(rownames(vst_data_2_mean)=="ENSG00000276566")]


vst_data_ratio[which(colnames(vst_data_ratio)=="ENSG00000102837")]

vst_data$Time
vst_data_ratio$TBMDEX

mean(vst_data_ratio$ENSG00000102837[1:79])
mean(vst_data_ratio$ENSG00000102837[80:179])
library(Lock5Data)
data(CommuteAtlanta)
str(CommuteAtlanta)


age.mean = with(vst_data, mean(Age))


B = 1000
n = nrow(vst_data)
boot.samples = matrix(sample(vst_data$Age, size = B * n, replace = TRUE),
                      B, n)

dim(boot.samples)
boot.statistics = apply(boot.samples, 1, mean)



orde




sim <- function(y, prev , or) {
  n <- length(y)
  p <- length(prev)
  if(p != length(or)) stop('prev and or must have the same length')
# prev = Pr(x = 1|y = 0) ; let the odds for this be oprev = prev/( 1 - prev)
# or = odds(x = 1 | y = 1 )/oprev
# Pr(x = 1 | y = 1) = oprev /((1/or) + oprev)
  oprev <- prev/(1 - prev)
  p1 <- oprev/((1/or) + oprev)
  n0 <- sum(y == 0)
  n1 <- sum(y == 1)
  # For n0 observations sample x so that Pr( x = 0 | y = 0) = prev
  nxy0 <- rbinom(p, size =n0 , prob = prev)
  nxy1 <- rbinom(p, size =n1 , prob =p1)
  # Compute p sample odds ratios
  sor <- (n0 - nxy0) * nxy1/(nxy0 * (n1 - nxy1))
  g <- function(x) ifelse(x >= 1, x, 1/x)
  r1 <- rank(sor)[which.max(or)]/p
  r2 <- rank(or)[which.max(sor)]/p
  data.frame(prev , or , nx= nxy0/n0 , obsprev0 =(nxy0 + nxy1)/n,
               obsprev = nxy1/(nxy0 + nxy1 ), obsor =sor , n=n,
               N = paste('n', n, sep=':'),
               Features = paste('Features', p, sep=':'),
               mmoe = quantile(g(sor/or), 0.90 , na.rm = TRUE ),
               obsranktrue =r1 , truerankobs =r2 ,
               rho=cor(sor , or , method ='spearman', use='pair'))
}




U <- NULL
set.seed(1)

for(n in c(50 , 100 , 250 , 500 , 1000 , 2000)) {
  for(p in c(10 , 50, 500 , 1000 , 2000)) {
    for(yprev in c(.1 , .3)) {
      y <- rbinom(n, 1, yprev)
      prev <- runif(p, .05 , .5)
      or <- exp(rnorm(p, 0, .25))
      u <- cbind(sim(y, prev , or),
                 Yprev = paste ('Prevalence of Outcome ', yprev , sep=':'))
      U <- rbind(U, u)
    }
  }
}


require(ggplot2)
pl <- function(yprev) {
  br <- c(.01 , .1 , .5 , 1, 2.5 , 5, 25, 100)
  ggplot(subset(U, Yprev = yprev),
         aes(x=or , y= obsor)) + geom_point() + facet_grid( Features~N) +
    ggtitle(paste('Prevalence of Outcome', yprev , sep=':')) +
    xlab('True ORs') + ylab('Estimated ORs') +
    scale_x_log10(breaks =br) + scale_y_log10(breaks =br) +
    theme( axis.text.x = element_text(size = rel(0.8), angle =-45 ,
                                         hjust =0, vjust =1)) +
    geom_abline(col='red')
}
pl(0.1)
pl(0.3)


ggplot(U, aes(x=n, y=mmoe)) + geom_point() + facet_grid(Features ~ Yprev) +
  geom_hline(aes(yintercept =1.5 , col='red')) +
  ylim(1, 10) +
  ylab('0.9 Quantile of Multiplicative Margin of Error in OR Across Features ')


ggplot(U, aes(x=n, y=rho)) + geom_point() +
  facet_grid(Features ~ Yprev) +
  ylab(expression(paste('Spearman', rho , 'Rank Correlation Between',
                        OR , 'and', hat(OR), 'Across Features')))

#=================================================================================
#=================================================================================



# Function to simulate the raw data
# prev is the vector of prevalences of x when y=0 as before
# y prev is the overall prevalence of y
# n is the sample size
# or is the vector of true odds ratios
sim <- function(n, yprev , prev , or) {
  y <- rbinom(n, 1, yprev)
  p <- length(prev)
  if(p != length(or)) stop ('prev and or must have the same length ')
  # p r e v = P r ( x = 1 | y = 0 ) ; l e t t h e o d d s f o r t h i s b e o p r e v = p r e v / ( 1 - p r e v )
  # o r = o d d s ( x = 1 | y = 1 ) / o p r e v
  # P r ( x = 1 | y = 1 ) = o p r e v / ( ( 1 / o r ) + o p r e v )
  oprev <- prev/(1 - prev)
  p1 <- oprev/((1/or) + oprev)
  x <- matrix(NA , nrow =n, ncol =p)
  for(j in 1 : p)
    x[, j] <- ifelse(y== 1, rbinom(n, 1, prob = p1[j]),
                     rbinom(n, 1, prob = prev[j]))
  list (x=x, y=y)
}
# F u n c t i o n t o c o m p u t e t h e s a m p l e o d d s r a t i o s g i v e n x m a t r i x a n d y v e c t o r
ors <- function(x, y) {
  p <- ncol(x)
  or <- numeric(p)
  #j=1
  for(j in 1:p) {
    f <- table(x[, j], y)
    or[j] <- f[2, 2] * f[1, 1]/(f[1, 2] * f[2, 1])
  }
  or
}


# Generate sample of size 600 with 300 features
# Log odds ratios have a normal distribution with mean 0 SD 0.3
# x have a random prevalence uniform[0.05, 0.5]
# y has prevalence 0.3
set.seed(188)
n <- 600; p <- 300
prev <- runif(p, .05 , .5)
or <- exp(rnorm(p, 0, .3))
z <- sim(n, 0.3, prev , or)

z$x
head(z)
# Compute estimated ORs
x <- z$x; y <- z$y
sor <- ors(x, y)
# Show how estimates related to true ORs

ggplot(data.frame(or , sor), aes(x=or , y=sor)) + geom_point() +
  xlab('True OR') + ylab('Estimated OR')

# Print the largest estimated OR and its column number,
# and corresponding true OR , and similarly for the smallest.
largest <- max(sor)
imax <- which.max(sor)
true.imax <- or[imax]
mmoe.imax <- largest/true.imax
smallest <- min(sor)
imin <- which.min(sor)
true.imin <- or[imin]
mmoe.imin <- smallest/true.imin
cat('\nLargest observed OR\n')

cat('OR:', round(largest , 2) , ' Feature #', imax , ' True OR:',
    round(true.imax , 2) , ' MMOE :', round(mmoe.imax , 2) , '\n')

cat('Rank of winning feature among true ORs:', sum(or <= or[imax]) , '\n\n')

cat('Smallest observed OR\n')

cat('OR:', round(smallest , 2) , ' Feature #', imin , ' True OR:',
    round(true.imin , 2) , ' MMOE :', round(mmoe.imin , 2) , '\n')


## use bootstrap to estimate

set.seed (11)
B <- 1000
ranksS <- ranksL <- mmoeS <- mmoeL <- numeric(B)
k=1
for(k in 1:B) {
  # Draw a sample of size n with replacement
  i <- sample(1:n, n, replace = TRUE)
  # Compute sample ORs on the new sample
  bor <- ors(x[i,], y[i])
  blargest <-  max(bor)
  bmax <- which.max(bor)
  ranksL[k] <- sum(bor <= largest)
  mmoeL[k] <- blargest/sor[bmax]
  bsmallest <- min(bor)
  bmin <- which.min(bor)
  ranksS[k] <- sum(bor <= smallest)
  mmoeS[k] <- bsmallest/sor[bmin]
}


pr <- function(which, ranks, mmoe, mmoe.true, estor, or.true) {
  gm <- exp(mean(log(mmoe)))
  cat(which , 'OR\n')
  cat('CL for rank :', quantile(ranks, c(0.025 , 0.975)),
      'Median MMOE :', round(median(mmoe),2) ,
      'Geometric mean MMOE :', round(gm,2) ,
      '\nTrue MMOE :', round(mmoe.true , 2) , '\n')
  bmmoe <- if(which == 'Largest') gm else median(mmoe)
  cat('Bootstrap bias-corrected ', tolower(which), 'OR:',
      round(estor/bmmoe , 2) ,
      'Original OR:', round(estor , 2) ,
      'True OR:', round(or.true , 2) ,
      '\n\n')
}

which='Largest'; ranks=ranksL; mmoe=mmoeL; mmoe.true=mmoe.imax; estor=largest;
or.true=true.imax
pr('Largest', ranksL , mmoeL , mmoe.imax , largest , true.imax )

pr('Smallest ', ranksS , mmoeS , mmoe.imin , smallest , true.imin)
