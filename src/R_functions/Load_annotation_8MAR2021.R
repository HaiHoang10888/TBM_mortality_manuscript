
# I. Set up input data
#--------------define working directory
workdir='D:/RNA_output/PAXGENE 27TB/count/'
##set the working directory
setwd(workdir)

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
Select_var <- c("Patient_ID","LIMS_ID","Timepoint",
                "trial_arm","TBMGRADE_Base","LTA4H","RIN",
                "GeneXpert_Base","SEX","Age","DIAGNOSIS_SCORE", 
                "DEATH_EVE_3M","DEATH_EVE_2M")

## set file data for PTB and TBM
data_PTB <- 'Meta_data_PTB_8MAR2021.RData'
data_TBM <- 'TBM_ALL_mergedbyrbind.RData_8MAR2021.RData'