

library(reshape)
library(Hmisc)
library(ggplot2)
library(gridExtra)
library(reshape)
library(tidyr)
library(kableExtra)
library(dplyr)
library(plyr)
library(compareGroups)
library(AnnotationDbi)
library(GO.db)
library(org.Hs.eg.db)


rm(list=ls()) # clear all previous object

source("C:/Users/haiht/Desktop/Dropbox/Project2020/WholeGenomeSeqHuman27TB/Hai_function/Function_data_process_and_analysis_8MAR2021_deskop.R")

source("C:/Users/haiht/Dropbox/Project2020/WholeGenomeSeqHuman27TB/Hai_function/Function_data_process_and_analysis_8MAR2021.R")
#--------------define working directory

output <- "count_old_Run1vsRun2_separate/Deconvolution_6May21"
workdir='C:/Users/haiht/Desktop/RNA_output/PAXGENE 27TB/count/'

##set the working directory
setwd(workdir)

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
data_TBM <- 'TBM_ALL_mergedbyrbind.RData_8MAR2021.RData'


## create countdata_Run1
# coun_data_all <- create_count_Run1(workdir = workdir,Filter_rRNA = F)

## create countdata_all
count_data_all <- create_count_data_all(workdir = workdir,Filter_rRNA = F)

sample_list <- colnames(count_data_all)[which(colnames(count_data_all)=="HOA7812A1"):ncol(count_data_all)]

exp_data <- count_data_all[,sample_list]
genes_length <- count_data_all$gene_length

# move working dir to output folder
setwd(paste0(workdir,output))


clinical_data <- create_TB_TBM_data(workdir = workdir,data_PTB = data_PTB,data_TBM =data_TBM ,output = "TBM_D0")
TBM_D0_P <- subset(clinical_data,trial_arm=="Placebo"& trial_number=="27TB"&Timepoint=="D0")
TBM_D0_D <- subset(clinical_data,trial_arm=="Dexamethasone"& trial_number=="27TB"&Timepoint=="D0")
list_TBM_D0 <- c(rownames(TBM_D0_P), rownames(TBM_D0_D))

clinical_data <- create_TB_TBM_data(workdir = workdir,data_PTB = data_PTB,data_TBM =data_TBM ,output = "TBM_D0_D14")
TBM_D14_P <- subset(clinical_data,trial_arm=="Placebo"& trial_number=="27TB"&Timepoint=="D14")
TBM_D14_D <- subset(clinical_data,trial_arm=="Dexamethasone"& trial_number=="27TB"&Timepoint=="D14")
list_TBM_D14 <- c(rownames(TBM_D14_P), rownames(TBM_D14_D))

clinical_data <- create_TB_TBM_data(workdir = workdir,data_PTB = data_PTB,data_TBM =data_TBM ,output = "TBM_D0_D60")
TBM_D60_P <- subset(clinical_data,trial_arm=="Placebo"& trial_number=="27TB"&Timepoint=="D60")
TBM_D60_D <- subset(clinical_data,trial_arm=="Dexamethasone"& trial_number=="27TB"&Timepoint=="D60")
list_TBM_D60 <- c(rownames(TBM_D60_P), rownames(TBM_D60_D))

clinical_data <- create_TB_TBM_data(workdir = workdir,data_PTB = data_PTB,data_TBM =data_TBM ,output = "All_TBM_D0")
TBM_26TB_D0 <- subset(clinical_data,trial_number=="26TB"&Timepoint=="D0")
TBM_27TB_D0 <- subset(clinical_data,trial_number=="27TB"&Timepoint=="D0")

clinical_data <- create_TB_TBM_data(workdir = workdir,data_PTB = data_PTB,data_TBM =data_TBM ,output = "All_TB_D0")
PTB_D0 <- subset(clinical_data,trial_arm=="PTB")
PTB_28TB <- subset(clinical_data,trial_number=="MDR_PTB")
PTB_29TB <- subset(clinical_data,trial_number=="DS_PTB")
PTB_AFBneg <- subset(clinical_data,trial_number=="AFBneg_PTB")


## script for deconvolution 
## create function for  convert raw count to rpkm or tpm
rpkm <- function(counts, lengths) {
  rate <- counts / lengths 
  rate / sum(counts) * 1e6
}

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}



## convert raw count to tpm

#rpkms <- apply(counts, 2, function(x) rpkm(x, genes$Length))

exp_tpms <- apply(exp_data, 2, function(x) tpm(x, genes_length))

if(identical(rownames(exp_tpms), count_data_all$gene_id)){
  rownames(exp_tpms) <- count_data_all$gene_symbol
}


gene_duplicate <- rownames(exp_tpms)[which(duplicated(rownames(exp_tpms)))] 
gene_duplicate <- paste0(gene_duplicate,"_R")

## renames rows which are duplicate
rownames(exp_tpms)[which(duplicated(rownames(exp_tpms)))]  <- gene_duplicate

exp_tpms_test <- as.data.frame(exp_tpms)

x= c()
for(i in colnames(exp_tpms_test)){
  t <- max(exp_tpms_test[which(colnames(exp_tpms_test)==i)])
  print(i)
  print(t)
  x = c(x,t)
}

write.csv(exp_tpms, "exp_tpms_Allsamples_6May21.csv")


# Put in your actual path where the text files are saved
temp_dir <- "C:/Users/haiht/Desktop/10Class_cells"
setwd(temp_dir)

## -1 to remove current dir, which we dont want to include
sub_dir <- list.dirs(path = temp_dir, full.names = TRUE, recursive = TRUE)[-1] 

GEP_cell_population <- list()

library(stringr)
for( i in sub_dir){
  temp_cells <- str_replace(i,paste0(temp_dir,"/"),"")
  #temp_file <- paste0(i,"")
  #print(temp_cells)
  txt_files = list.files(path=i, pattern="*.txt")
  pos <- 1 
  for (k in txt_files){
    file_link <- paste0(i,"/",k)
    df_name <- str_replace(k,paste0("_",str_replace_all(temp_cells,"_",""),".txt"),"")
    #print(temp_cells)
    #print(file_link)
    temp_df <- read.table(file_link, row.names = 1,header = T, na.strings="NA")
    temp_text <- paste0("GEP_cell_population$",temp_cells, "$",df_name, " <- temp_df")
    eval(parse( text=temp_text))
    #print(df_name)
    if(pos==1) temp_All_TBM <- temp_df
    else if (pos>1 & all(rownames(temp_All_TBM)==rownames(temp_df))) temp_All_TBM <- cbind(temp_All_TBM, temp_df) 
    else  stop("temp_All_TBM and temp_df must have same row names!")
    pos <- pos +1
  }
  temp_text2 <- paste0("GEP_cell_population$",temp_cells, "$All_TBM", " <- temp_All_TBM")
  eval(parse( text=temp_text2))
}


save(GEP_cell_population, file= "GEP_cell_population.RData")

## impute
Neutrophil <- GEP_cell_population$Neutrophils$All_TBM

D <- rownames(TBM_D0_D)
Neutrophil_D0 <- subset(Neutrophil, select=c(D))

# inpute NA with 0.5
Neutrophil_D0[is.na(Neutrophil_D0)] <- 1
Neutrophil_D0 <- as.matrix(Neutrophil_D0)




Neutro_degran <- get_genelist_pathway("neutrophil degranulation")
IFN_I <- get_genelist_pathway("type I interferon signaling pathway")

Neutrophil_D0 <- subset(Neutrophil_D0, rownames(Neutrophil_D0) %in% IFN_I)

idx <- which(apply(Neutrophil_D0, 1, function(x) all(x == 1)) == T)
Neutrophil_D0 <- Neutrophil_D0[-idx,]
View(Neutrophil_D0)

df_col <- c(rep("Dex", nrow(TBM_D0_D)), rep("Placebo", nrow(TBM_D0_P)))

df_col <- data.frame(Treatment=df_col)
df_col <- list_TBM_D0


library(pheatmap)
pheatmap(Neutrophil_D0, scale = "row", cluster_cols=F,annotation_col = df_col)

library(org.Hs.eg.db)
library(GO.db)
IFN_I <- get_genelist_pathway("type I interferon signaling pathway")

Neutrophil_D0 <- subset(Neutrophil_D0, rownames(Neutrophil_D0) %in% IFN_I)

df_col <- c(rep("Dex", nrow(TBM_D0_P)), rep("Placebo", nrow(TBM_D0_D)))
df_col <- data.frame(Treatment=df_col)
list_TBM_D0 <- c(rownames(TBM_D0_P), rownames(TBM_D0_D))

getwd()
Neutrophil <- read.table("RData/All_TBM_Monocytes.txt", header = T, row.names = 1)
Neutrophil_D0 <- subset(Neutrophil, select=c(list_TBM_D0))
Neutrophil_D0[is.na(Neutrophil_D0)] <- 1


Neutrophil_D0 <- subset(Neutrophil_D0, rownames(Neutrophil_D0) %in% IFN_I)

idx <- which(apply(Neutrophil_D0, 1, function(x) all(x == 1)) == T)
Neutrophil_D0 <- Neutrophil_D0[-idx,]
#View(Neutrophil_D0)

rownames(df_col) <- list_TBM_D0

Neutrophil_D0 <- as.matrix(Neutrophil_D0)


library(pheatmap)
pheatmap(Neutrophil_D0, scale = "row", cluster_cols=F,annotation_col = df_col)



#-------------------------------------------
list_TBM_D14 <- c(rownames(TBM_D14_P), rownames(TBM_D14_D))
df_col <- c(rep("Dex", nrow(TBM_D14_P)), rep("Placebo", nrow(TBM_D14_D)))
df_col <- data.frame(Treatment=df_col)
rownames(df_col) <- list_TBM_D14

Neutrophil_D14 <- subset(Neutrophil, select=c(list_TBM_D14))
Neutrophil_D14[is.na(Neutrophil_D14)] <- 1

IFN_I <- get_genelist_pathway("type I interferon signaling pathway")

Neutrophil_D14 <- subset(Neutrophil_D14, rownames(Neutrophil_D14) %in% IFN_I)

idx <- which(apply(Neutrophil_D14, 1, function(x) all(x == 1)) == T)
Neutrophil_D14 <- Neutrophil_D14[-idx,]
# View(Neutrophil_D14)

Neutrophil_D14 <- as.matrix(Neutrophil_D14)
library(pheatmap)
pheatmap(Neutrophil_D14, scale = "row", cluster_cols=F,annotation_col = df_col)



Neutrophil_D14 <- subset(Neutrophil, select=c(list_TBM_D14))
Neutrophil_D14[is.na(Neutrophil_D14)] <- 1
Neutrophil_D14 <- subset(Neutrophil_D14, rownames(Neutrophil_D14) %in% Neutro_degran)
idx <- which(apply(Neutrophil_D14, 1, function(x) all(x == 1)) == T)
Neutrophil_D14 <- Neutrophil_D14[-idx,]

library(pheatmap)
pheatmap(Neutrophil_D14, scale = "row", cluster_cols=F,annotation_col = df_col)








