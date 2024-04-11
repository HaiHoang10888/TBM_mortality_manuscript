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


create_vsd_t <- function(data_dds,data_coldata){
  vsd <- varianceStabilizingTransformation(data_dds, blind = F)
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

delta_treatmentvscontrol <- function(vst_data_ratio,gene_list){
  vst_data_1 <- subset(vst_data_ratio,TBMDEX=="Dexamethasone")
  vst_data_1 <- vst_data_1[,gene_list]
  vst_data_1_mean <- colMeans(vst_data_1)
  vst_data_1_mean <- as.matrix(vst_data_1_mean)
  
  vst_data_2 <- subset(vst_data_ratio,TBMDEX=="Placebo")
  vst_data_2 <- vst_data_2[,gene_list]
  vst_data_2_mean <- colMeans(vst_data_2)
  vst_data_2_mean <- as.matrix(vst_data_2_mean)
  vst_data_delta <- vst_data_1_mean-vst_data_2_mean
  return(vst_data_delta)
}
