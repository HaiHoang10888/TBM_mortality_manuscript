
# source all functions
library(stringr)
gitdir <- paste0(strsplit(getwd(), split = "/github_repos", fixed=T)[[1]][1], "/github_repos/")

source(paste0(gitdir,"/R_functions/Load_library_8MAR2021.R"))
source(paste0(gitdir,"/R_functions/Function_data_process_8MAR2021.R"))
source(paste0(gitdir,"/R_functions/Functional_Pathway_enrichment_functions.R"))
source(paste0(gitdir,"/R_functions/Emapplot_cluster_Hai.R"))
source(paste0(gitdir,"/R_functions/EnhancedVocanoplot_hai.R"))
source(paste0(gitdir,"/R_functions/Select_genes_for_Pathway_heatmap.R"))
source(paste0(gitdir,"/R_functions/Deconvolution_function.R"))
source(paste0(gitdir,"/R_functions/move_column_Function.R"))
source(paste0(gitdir,"/R_functions/Function_extract_res_output.R"))
source(paste0(gitdir,"/R_functions/EnhancedVocanoplot_continuous_Outcome_new.R"))
source(paste0(gitdir,"/R_functions/WGCNA_CEMitool_Allfunctions_5May2022.R"))

