Create_GEP_for_each_celltypes <- function(workdir="D:/RNA_output/PAXGENE 27TB/count/",target,file){
  sub_dir <- paste0(workdir,target)
  Neutrophils <- read.table(paste0(sub_dir,"/Neutrophils.txt"), header = T, row.names = 1)
  Monocytes <- read.table(paste0(sub_dir,"/Monocytes.txt"), header = T, row.names = 1)
  Eosinophils <- read.table(paste0(sub_dir,"/Eosinophils.txt"), header = T, row.names = 1)
  Bcells <- read.table(paste0(sub_dir,"/Bcells.txt"), header = T, row.names = 1)
  TcellsCD4 <- read.table(paste0(sub_dir,"/TcellsCD4.txt"), header = T, row.names = 1)
  TcellsCD8 <- read.table(paste0(sub_dir,"/TcellsCD8.txt"), header = T, row.names = 1)
  NKcells <- read.table(paste0(sub_dir,"/NKcells.txt"), header = T, row.names = 1)
  Dendriticcells <- read.table(paste0(sub_dir,"/Dendriticcells.txt"), header = T, row.names = 1)
  Mastcells <- read.table(paste0(sub_dir,"/Mastcells.txt"), header = T, row.names = 1)
  Plasmacells <- read.table(paste0(sub_dir,"/Plasmacells.txt"), header = T, row.names = 1)
  save(Neutrophils,Monocytes,Eosinophils,Bcells,TcellsCD4,
       TcellsCD8,NKcells,Dendriticcells,Mastcells,Plasmacells,
       file = paste0(sub_dir,"/",file,".RData"))
}
cells <- c("Neutrophils","Monocytes","Eosinophils",
           "Bcells", "TcellsCD4", "TcellsCD8",
           "NKcells", "Dendriticcells", "Mastcells", "Plasmacells")
extract_GEP_specific_cells <- function(cells,Sample_ID,Pathway_gene){
  data_output <- data.frame(gene=Pathway_gene)
  for (data_name in cells) {
    #print(data_name)
    str1 <- paste0("Temp_data <- ", data_name)
    eval(parse( text=str1))
    Temp_data <- subset(Temp_data, select=c(Sample_ID))
    Temp_data <- Temp_data[Pathway_gene,]
    Temp_data[is.na(Temp_data)] <- 1
    str2 <- paste0("data_output$", data_name, "<-  rowMeans(Temp_data)")
    eval(parse( text=str2))
  }
  return(data_output)
}

get_genelist_RA <- function(pathway){
  library(ReactomeContentService4R)
  infor_data <- event2Ids(event.id = pathway)
  gene_list <- infor_data$geneSymbol
  return(gene_list)
}


draw_corplot_for_cibersort_old <- function(cells=c("Neutrophils","Monocytes",
                                               "Eosinophils","Bcells",
                                               "TcellsCD4","TcellsCD8",
                                               "NKcells","Dendriticcells",
                                               "Mastcells","Plasmacells"),
                                       Sample_ID,
                                       pathway,
                                       DE_gene,
                                       type="GO",
                                       win.asp=1){
  library(corrplot)
  if(type=="GO") Pathway_gene <- get_genelist_pathway(pathway)
  if(type=="RA") Pathway_gene <- get_genelist_RA(pathway)
  Pathway_gene <- Pathway_gene[which(Pathway_gene %in% DE_gene)]
  data_GEP <- extract_GEP_specific_cells(cells = cells,Sample_ID = Sample_ID,
                                    Pathway_gene = Pathway_gene)
  
  m <- as.matrix(data_GEP[-1])
  m <- t(apply(m, 1,function(x) x / sum(x,na.rm = T)*100))
  row.names(m) <- data_GEP$gene
  m <- log2(m)
  d <- dist(as.matrix(m))
  hc <- hclust(d)
  gene_order <- hc$order
  m <- m[gene_order,]
  m <- t(m)
  corrplot(m, method = "circle",is.corr = FALSE,mar = c(0, 0, 0, 0),
           win.asp=win.asp, tl.cex=0.6)
  
  
}


# load(paste0(paste0(workdir,target,"/",file,".RData")))
# Create_GEP_for_each_celltypes
# file="test_GEP"
# target <- "count_old_Run1vsRun2_separate/Deconvolution_6May21/RData/Allsample_All_IFN_Neutro_DE179TBM_D14"
# Create_GEP_for_each_celltypes(target=target,file=file)


draw_corplot_matrix<- function(data_GEP,win.asp=1){
  library(corrplot)
  m <- as.matrix(data_GEP[-1])
  m <- log2(m)
  row.names(m) <- data_GEP$gene
  corrplot(m, method = "circle",is.corr = FALSE,mar = c(0, 0, 0, 0), win.asp=win.asp, tl.cex=0.6)
}

draw_corplot_for_wholeblood <- function(vst_data,
                                        DE_gene,
                                        cells_blood_type=c("WBC","NEUTLE","LYMLE","MONOLE"),
                                        pathway,
                                        tl.cex=0.6){
  library(corrplot)
  target_gene <- get_genelist_pathway(pathway)
  target_gene <- target_gene[which(target_gene %in% DE_gene)]
  target_var <- c(cells_blood_type,target_gene)
  vst_data_p <- as.matrix(vst_data[,target_var])
  vst_data_p <- vst_data_p[complete.cases(vst_data_p),]
  M <- cor(vst_data_p)
  M <- M[cells_blood_type,target_gene]
  d <- dist(t(as.matrix(M)))
  hc <- hclust(d)
  gene_order <- hc$order
  M <- M[,gene_order]
  #ncol(vst_data_p)
  testRes = cor.mtest(vst_data_p, conf.level = 0.95)
  testRes$p <- testRes$p[1:4,5:length(target_var)]
  testRes$p <- testRes$p[,gene_order]
  corrplot(M, p.mat = testRes$p,insig='blank', tl.cex=tl.cex)
}


corr_simple <- function(data=df,sig=0.5){
  #convert data to numeric in order to run correlations
  #convert to factor first to keep the integrity of the data - each value will become a number rather than turn into NA
  df_cor <- data %>% mutate_if(is.character, as.factor)
  df_cor <- df_cor %>% mutate_if(is.factor, as.numeric)
  #run a correlation and drop the insignificant ones
  corr <- cor(df_cor)
  #prepare to drop duplicates and correlations of 1     
  corr[lower.tri(corr,diag=TRUE)] <- NA 
  #drop perfect correlations
  corr[corr == 1] <- NA 
  #turn into a 3-column table
  corr <- as.data.frame(as.table(corr))
  #remove the NA values from above 
  corr <- na.omit(corr) 
  #select significant values  
  corr <- subset(corr, abs(Freq) > sig) 
  #sort by highest correlation
  corr <- corr[order(-abs(corr$Freq)),] 
  #print table
  print(corr)
  #turn corr back into matrix in order to plot with corrplot
  mtx_corr <- reshape2::acast(corr, Var1~Var2, value.var="Freq")
  
  #plot correlations visually
  corrplot(mtx_corr, is.corr=FALSE, tl.col="black", na.label=" ")
}



draw_corplot_for_cibersort <- function(cells=c("Neutrophils","Monocytes",
                                               "Eosinophils","Bcells",
                                               "TcellsCD4","TcellsCD8",
                                               "NKcells","Dendriticcells",
                                               "Mastcells","Plasmacells"),
                                       Sample_ID,
                                       pathway=NULL,
                                       DE_gene,
                                       type="GO",
                                       win.asp=1){
  library(corrplot)
  if(type=="GO") Pathway_gene <- get_genelist_pathway(pathway)
  if(type=="RA") Pathway_gene <- get_genelist_RA(pathway)
  if(type=="random") Pathway_gene <- DE_gene
  Pathway_gene <- Pathway_gene[which(Pathway_gene %in% DE_gene)]
  data_GEP <- extract_GEP_specific_cells(cells = cells,Sample_ID = Sample_ID,
                                         Pathway_gene = Pathway_gene)
  
  m <- as.matrix(data_GEP[-1])
  row.names(m) <- data_GEP$gene
  m <- log2(m)
  if (length(row.names(m))>3) {
    d <- dist(as.matrix(m))
    hc <- hclust(d)
    gene_order <- hc$order
    m <- m[gene_order,]
  }
  m <- t(m)
  corrplot(m, method = "circle",is.corr = FALSE,mar = c(0, 0, 0, 0),
           win.asp=win.asp, tl.cex=0.6)
  
  
}

draw_corplot_for_cibersort_tmod <- function(cells=c("Neutrophils","Monocytes",
                                               "Eosinophils","Bcells",
                                               "TcellsCD4","TcellsCD8",
                                               "NKcells","Dendriticcells",
                                               "Mastcells","Plasmacells"),
                                       Sample_ID,
                                       module,
                                       DE_gene,
                                       mset,
                                       win.asp=1){
  library(corrplot)
  genelist <- paste0("mset$MODULES2GENES$",module)
  genelist <- eval(parse( text=genelist ))
  genelist <- unname(genelist)
  module_gene <- genelist[which(genelist %in% DE_gene)]
  data_GEP <- extract_GEP_specific_cells(cells = cells,Sample_ID = Sample_ID,
                                         Pathway_gene = module_gene)
  
  m <- as.matrix(data_GEP[-1])
  row.names(m) <- data_GEP$gene
  m <- log2(m)
  if (length(row.names(m))>3) {
    d <- dist(as.matrix(m))
    hc <- hclust(d)
    gene_order <- hc$order
    m <- m[gene_order,]
  }
  m <- t(m)
  corrplot(m, method = "circle",is.corr = FALSE,mar = c(0, 0, 0, 0),
           win.asp=win.asp, tl.cex=0.6)
  
  
}





draw_corplot_for_cibersort_pt <- function(cells=c("Neutrophils","Monocytes",
                                               "Eosinophils","Bcells",
                                               "TcellsCD4","TcellsCD8",
                                               "NKcells","Dendriticcells",
                                               "Mastcells","Plasmacells"),
                                       Sample_ID,
                                       pathway=NULL,
                                       DE_gene,
                                       type="GO",
                                       win.asp=1){
  library(corrplot)
  if(type=="GO") Pathway_gene <- get_genelist_pathway(pathway)
  if(type=="RA") Pathway_gene <- get_genelist_RA(pathway)
  if(type=="random") Pathway_gene <- DE_gene
  Pathway_gene <- Pathway_gene[which(Pathway_gene %in% DE_gene)]
  data_GEP <- extract_GEP_specific_cells(cells = cells,Sample_ID = Sample_ID,
                                         Pathway_gene = Pathway_gene)
  
  m <- as.matrix(data_GEP[-1])
  m <- t(apply(m, 1,function(x) x / sum(x,na.rm = T)*100))
  row.names(m) <- data_GEP$gene
  if (length(row.names(m))>3) {
    d <- dist(as.matrix(m))
    hc <- hclust(d)
    gene_order <- hc$order
    m <- m[gene_order,]
  }
  
  m <- t(m)
  corrplot(m, method = "circle",is.corr = FALSE,mar = c(0, 0, 0, 0),
           win.asp=win.asp, tl.cex=0.6)
  
  
}



draw_corplot_for_cibersort_pt_tmod <- function(cells=c("Neutrophils","Monocytes",
                                                  "Eosinophils","Bcells",
                                                  "TcellsCD4","TcellsCD8",
                                                  "NKcells","Dendriticcells",
                                                  "Mastcells","Plasmacells"),
                                          Sample_ID,
                                          module,
                                          DE_gene,
                                          mset,
                                          win.asp=1){
  library(corrplot)
  genelist <- paste0("mset$MODULES2GENES$",module)
  genelist <- eval(parse( text=genelist ))
  genelist <- unname(genelist)
  module_gene <- genelist[which(genelist %in% DE_gene)]
  data_GEP <- extract_GEP_specific_cells(cells = cells,Sample_ID = Sample_ID,
                                         Pathway_gene = module_gene)
  
  m <- as.matrix(data_GEP[-1])
  m <- t(apply(m, 1,function(x) x / sum(x,na.rm = T)*100))
  row.names(m) <- data_GEP$gene
  if (length(row.names(m))>3) {
    d <- dist(as.matrix(m))
    hc <- hclust(d)
    gene_order <- hc$order
    m <- m[gene_order,]
  }
  
  m <- t(m)
  corrplot(m, method = "circle",is.corr = FALSE,mar = c(0, 0, 0, 0),
           win.asp=win.asp, tl.cex=0.6)
  
  
}

# load(paste0(paste0(workdir,target,"/",file,".RData")))
# Create_GEP_for_each_celltypes
# file="test_GEP"
# target <- "count_old_Run1vsRun2_separate/Deconvolution_6May21/RData/Allsample_All_IFN_Neutro_DE179TBM_D14"
# Create_GEP_for_each_celltypes(target=target,file=file)


draw_corplot_matrix<- function(data_GEP,win.asp=1){
  library(corrplot)
  m <- as.matrix(data_GEP[-1])
  m <- log2(m)
  row.names(m) <- data_GEP$gene
  corrplot(m, method = "circle",is.corr = FALSE,mar = c(0, 0, 0, 0), win.asp=win.asp, tl.cex=0.6)
}



extract_GEP_specific_cells_raw <-
  function(cells,Sample_ID,Gene_list="LTA4H"){
    data_output <- data.frame(Sample_ID=Sample_ID)
    for (data_name in cells) {
      #print(data_name)
      str1 <- paste0("Temp_data <- ", data_name)
      eval(parse( text=str1))
      Temp_data <- subset(Temp_data, select=c(Sample_ID))
      Temp_data <- Temp_data[Gene_list,]
      Temp_data <- t(as.matrix(Temp_data))
      Temp_data[is.na(Temp_data)] <- 1
      Temp_data <- as.data.frame(Temp_data)
      str2 <- paste0("data_output$", Gene_list,".",data_name, "<-  Temp_data[,1]")
      eval(parse( text=str2))
    }
    return(data_output)
  }



.getmodules2 <- function(modules=NULL, mset="all", known.only=FALSE, skipcheck=FALSE) {
  
  # user provided mset
  if(is(mset, "list")) {
    mset <- .mset_sanity_check(mset, modules)
    if(is.null(modules)) modules <- mset[["MODULES"]]$ID
  } else {
    
    tmod <- .gettmod()
    
    mset <- match.arg( mset, c( "all", unique( tmod$MODULES$SourceID )) )
    if( mset != "all" ) { 
      if( is.null( modules ) ) modules <- tmod$MODULES$ID
      modules <- modules[ tmod$MODULES[modules,]$SourceID == mset ]
    }
    
    mset <- tmod
  }
  
  # filter the modules if hand-picked
  if(!is.null(modules)) {
    mset <- mset[modules,]
  }
  
  if(known.only) {
    mset <- mset[ ! is.na(mset$MODULES$Title) & ! mset$MODULES$Title %in% c( "TBA", "Undetermined", ""), ]
  }
  
  mset
}


## we need to use the tmod data set internally
.myDataEnv <- new.env(parent=emptyenv())

if(!exists("tmod", .myDataEnv)) {
  data("tmod", package="tmod", envir=.myDataEnv)
}

.gettmod <- function() {
  .myDataEnv[["tmod"]]
}


## check user provided mset for sanity
.mset_sanity_check <- function(mset, modules=NULL) {
  
  # sanity checks
  if(!all( c( "MODULES", "MODULES2GENES", "GENES") %in% names(mset)))
    stop("Required members missing from the list mset parameter")
  
  for( i in c( "MODULES", "GENES" )) {
    if(is.null(mset[[i]])) 
      stop(sprintf("Member %s of mset is NULL", i))
    if(!is(mset[[i]], "data.frame")) 
      stop(sprintf("Member %s of mset is not a data frame", i))
  }
  
  if(!is(mset[["MODULES2GENES"]], "list"))
    stop("Member MODULES2GENES of mset is not a list")
  
  if(!all(c("ID", "Title") %in% colnames(mset[["MODULES"]])))
    stop("Required columns missing from member MODULES")
  
  if(any(duplicated(mset[["MODULES"]]$ID))) 
    stop("Module IDs must not be duplicated")
  
  mset[["MODULES"]]$ID <- as.character(mset[["MODULES"]]$ID)
  rownames(mset[["MODULES"]]) <- mset[["MODULES"]]$ID
  
  if(!"ID" %in% colnames(mset[["GENES"]]))
    stop("Required column ID missing from member GENES")
  mset[["GENES"]]$ID <- as.character(mset[["GENES"]]$ID)
  
  missing <- !mset[["MODULES"]]$ID %in% names(mset[["MODULES2GENES"]])
  if(any(missing)) {
    stop(sprintf("Modules from MODULES member are missing in MODULES2GENES, first is %s", mset[["MODULES"]]$ID[missing][1] ))
  }
  
  if(!is.null(modules)) {
    
    if(!all(modules %in% mset[["MODULES"]]$ID )) {
      stop("Modules specified with the modules parameter are missing from definitions in the mset parameter")
    }
  }
  
  if(!is(mset, "tmod")) mset <- new("tmod", mset)
  mset
}
