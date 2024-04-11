

# source("C:/Users/haiht/Dropbox/Project2020/WholeGenomeSeqHuman27TB/Hai_function/Mixomic_function.R")


subset_GOGSE<-function(x, list_G0)
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
  return(data_tree)
} 


creat_list_genes_GO_cluster <- function(GOGSE, data_tree){
  for (i in names(data_tree$tree)){
    temp_list <- unname(unlist(data_tree$tree[i]))
    temp_list <- substr(temp_list,1,10)
    temp_GOGSE <- subset_GOGSE(x=GOGSE,list_G0=temp_list)
    temp_genes <- temp_GOGSE@result$core_enrichment
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


create_Loading_genes_for_cluster_GO <- function(data_tree){
  All <- data.frame(value.var=NA,comp=NA,GO_cluster=NA)
  for(i in names(data_tree$spca)){
    #print(i)
    str <- paste0("temp_pca  <- ", "data_tree$spca$",i )
    eval(parse( text=str))
    type_pca <- class(temp_pca)
    temp_data1 <- selectVar(temp_pca, comp = 1)$value
    temp_data1$comp <- 1
    temp_data2 <- selectVar(temp_pca, comp = 2)$value
    temp_data2$comp <- 2
    if("mixo_plsda" %in% type_pca){
      temp_data1 <- temp_data1[1:5,]
      temp_data2 <- temp_data2[1:5,]
    }
    temp_data <- rbind(temp_data1, temp_data2)
    temp_data[!duplicated(rownames(temp_data)), ]
    temp_data$GO_cluster <- i
    str <- paste0("data_tree$spca_loading$",i," <-  temp_data")
    eval(parse( text=str))
    All <- rbind(All, temp_data)
    All <- All[!is.na(All$value.var),]
  }
  data_tree$spca_loading$All <- All
  return(data_tree)
}


Select_genes_for_GO_pathway <- function(heatmap,k,GOGSE,vst_data,group="trial_arm",type="spca", plot=T){
  data_tree1 <- create_heat_tree_cluster(heatmap,k)
  data_tree2 <- creat_list_genes_GO_cluster(GOGSE,data_tree1)
  data_tree3 <- create_spca_for_heatmap_GO_cluster(vst_data,data_tree2, 
                                                  group,type,plot)
  data_tree <- create_Loading_genes_for_cluster_GO(data_tree3)
  return(data_tree)
  
}




data_tree <- Select_genes_for_GO_pathway(heatmap = p1$pathway_heatmap,k=3,
                                         GOGSE = GOGSE_2,vst_data = vst_data_ratio, type = "plsda")







test <- subset_gsea(x=GOGSE_2,cluster_1$GO_ID )

k=3
heatmap=p1$pathway_heatmap
gsea <- GOGSE_2
vst_data<- vst_data_ratio


plotIndiv(data_tree$spca$GO_cluster1, legend=T, title="title", 
          group = vst_data$trial_arm,
          ellipse=T, 
          legend.title=NULL, 
          ind.names = F,ellipse.level=0.95)


identical(rownames(vst_data), rownames(data_tree$spca$GO_cluster2$x))

col.sideColors <- palette()[vst_data$trial_arm]
col.sideColors
cim(data_tree$spca$GO_cluster3,transpose =T, row.sideColors = col.sideColors, cluster = "column")

plotVar(data_tree$spca$GO_cluster2,cex=3) 










class(data_tree$spca$GO_cluster2)
t <- selectVar(data_tree$spca$GO_cluster2, comp = 1)$value
selectVar(data_tree$spca$GO_cluster1, comp = 2)$value



plotLoadings(data_tree$spca$GO_cluster3, comp=1)

biplot(data_tree$spca$GO_cluster2, ind.names = F, group=vst_data$trial_arm)


.biplot <- 
  function(x,
           comp = c(1,2),
           block = NULL,
           ind.names = TRUE,
           group = NULL,
           cutoff = 0,
           col.per.group=NULL,
           col = NULL,
           ind.names.size = 3,
           ind.names.col = color.mixo(4),
           ind.names.repel = TRUE,
           pch = 19,
           pch.levels=NULL,
           pch.size = 2,
           var.names = TRUE,
           var.names.col = 'grey40',
           var.names.size = 4,
           var.names.angle = FALSE,
           var.arrow.col = 'grey40',
           var.arrow.size = 0.5,
           var.arrow.length = 0.2,
           ind.legend.title = NULL,
           vline = FALSE,
           hline = FALSE,
           legend = if (is.null(group)) FALSE else TRUE,
           legend.title = NULL,
           pch.legend.title = NULL,
           cex = 1.05,
           ...
  ){
    object <- x
    rm(x)
    ## for implicit support of non-pca objects - experimental
    block <- .change_if_null(block, 'X')
    comp <- .check_comp(comp, ncomp = object$ncomp)
    block <- match.arg(block, choices = names(object$variates))
    hide <- 'none'
    selection <- rowSums(object$loadings[[block]][, comp]) != 0 
    loadings <- object$loadings[[block]][selection, ]
    loadings <- data.frame(loadings)
    
    ## scale check
    if (isFALSE(object$call$scale))
      warning("The 'tune.spca' algorithm has used scale=FALSE. We recommend scaling the data",
              " to improve orthogonality in the sparse components.")
    ## cutoff correlation
    cutoff <- .check.cutoff(cutoff)
    cors <- cor(object[[block]][, selection], object$variates[[block]][, comp], use = 'pairwise' )
    # cors <- apply(cors, 1, function(x) (sqrt(sum(x^2))))
    # above.cutoff <- cors >= cutoff
    above.cutoff <- apply(cors, 1, function(x) any(abs(x) >= cutoff))
    loadings <- loadings[above.cutoff,]
    
    variates <- object$variates[[block]]
    variates <- data.frame(variates)
    ## scaler of var vs sample coordinates
    scaler <- max(variates, na.rm = TRUE)/max(loadings, na.rm = TRUE)
    
    PCs <- paste0('component_', comp)
    expl_vars <- round(object$prop_expl_var[[block]]*100)[comp]
    axes.titles <- sprintf("%s   (%s%%)", PCs, expl_vars)
    ind.names <- .get.character.vector(ind.names, vec = rownames(variates))
    
    variates$ind.names <- ind.names
    col.group <-
      .get.cols.and.group(
        col.per.group = col.per.group,
        group = group,
        col = col,
        object = object,
        n_ind = nrow(variates)
      )
    group <- col.group$group
    col.per.group <- col.group$col.per.group
    if (length(col.per.group) == 1)
    {
      legend <- FALSE
    }
    
    ## ------------- outline
    gg_biplot <- 
      ggplot() + 
      theme_classic() +  
      labs(x = axes.titles[1], 
           y = axes.titles[2])
    ## vline and hline
    if (vline)
    {
      gg_biplot <- gg_biplot + geom_vline(xintercept = 0, size = 0.3, col = 'grey75')
    }
    if (hline)
    {
      gg_biplot <- gg_biplot +  geom_hline(yintercept = 0, size = 0.3, col = 'grey75')
    }
    
    
    ## ------------- inds
    if (! 'ind' %in% hide) 
    {
      if (!is.null(pch)) 
      {
        ## ------------- advanced user args
        fill <- ifelse(is.null(list(...)$fill), 'black', list(...)$fill)
        alpha <- ifelse(is.null(list(...)$alpha), 1, list(...)$alpha)
        
        pch.res <- .get.pch(pch, pch.levels, n_ind = nrow(variates))
        pch <- pch.res$pch
        pch.levels <- pch.res$pch.levels
        pch.legend <- pch.res$pch.legend
        
        ## get 'pch' and 'group' arg used for legends so we can handle
        ## legends whether needed or not in a unified way (see scale_*_manual)
        pch.legend.title <- .change_if_null(pch.legend.title, as.character(as.list(match.call())['pch']))
        if (is.null(legend.title))
        {
          legend.title <- ifelse(is(object, 'DA'), yes = 'Y', no = as.character(as.list(match.call())['group']))
        }
        gg_biplot <- gg_biplot + 
          geom_point(aes(x = variates[, comp[1]], 
                         y = variates[, comp[2]],
                         col = group,
                         shape = pch),
                     fill = fill,
                     alpha = alpha,
                     size = pch.size)
        
        pch_legend <- NULL
        if (isTRUE(pch.legend)) {
          pch_legend <- guide_legend(title = pch.legend.title, override.aes = list(size = 5))
        }
        gg_biplot <- gg_biplot + 
          scale_shape_manual(values = pch.levels, guide = pch_legend)
      }
      else
      {
        ind.names.repel <- FALSE
      }
      if (!is.null(ind.names))
      {
        if (isTRUE(ind.names.repel)) {
          gg_biplot <- gg_biplot + 
            geom_text_repel(mapping = aes(x = variates[, comp[1]],
                                          y = variates[, comp[2]],
                                          label = ind.names,
                                          col = group
            ), 
            size = ind.names.size,
            show.legend = FALSE)
        } else {
          gg_biplot <- gg_biplot + 
            geom_text(mapping = aes(x = variates[, comp[1]],
                                    y = variates[, comp[2]],
                                    label = ind.names,
                                    col = group
            ), 
            size = ind.names.size,
            show.legend = FALSE)
        }
        
      }
      col_legend <- NULL
      if (isTRUE(legend)) {
        col_legend <- guide_legend(title = legend.title, override.aes = list(size = 5))
      }
      
      gg_biplot <- gg_biplot + 
        scale_color_manual(values = col.per.group, guide = col_legend)
      
      gg_biplot <-
        gg_biplot +
        theme(
          legend.text = element_text(size = rel(cex)),
          legend.title = element_text(size = rel(cex)),
          axis.title =  element_text(size = rel(cex)),
          axis.text =  element_text(size = rel(cex))
        )
      
    }
    
    ## ------------- vars
    
    if (! 'var' %in% hide) 
    {
      loadings <- loadings*scaler
      var.names.col <- .get.ind.colors(group = NULL, 
                                       col = var.names.col,
                                       col.per.group = NULL, 
                                       n_ind = nrow(loadings))
      if (!is.null(var.arrow.col))
      {
        var.arrow.col <- .get.ind.colors(group = NULL, 
                                         col = var.arrow.col,
                                         col.per.group = NULL, 
                                         n_ind = nrow(loadings))
        loadings$var.names.col <- var.names.col
        loadings$var.arrow.col <- var.arrow.col
        ## lines and arrows
        gg_biplot <-
          gg_biplot + geom_segment(
            aes(
              x = 0,
              y = 0,
              xend = loadings[,comp[1]],
              yend = loadings[,comp[2]],
            ),
            col = var.arrow.col,
            arrow = arrow(length = unit(var.arrow.length, "cm")),
            size = var.arrow.size,
            show.legend = FALSE
          )
      }
      
      ## labels
      var.labels <- .get.character.vector(arg = var.names, vec = rownames(loadings))
      ## label angles
      angle <- rep(0, nrow(loadings))
      
      if (!is.null(var.names)) 
      {
        angle <- rep(0, nrow(loadings))
        if (var.names.angle == TRUE)
        {
          angle <- atan(loadings[, comp[2]]/loadings[, comp[1]]) * 360/(2 * pi)
        }
        
        gg_biplot <-
          gg_biplot + geom_text_repel(
            aes(
              x = loadings[, comp[1]],
              y = loadings[, comp[2]],
              label = var.labels,
              angle = angle,
              hjust = ifelse(loadings[, comp[1]] > 0, 1, 0),
              vjust = ifelse(loadings[, comp[2]] > 0, 1, 0)
            ),
            col = var.names.col,
            size = var.names.size,
            box.padding = 0.1,
            point.padding = 0.1
          )
      } 
      
      ## second set of axes
      gg_biplot <- gg_biplot + scale_y_continuous(sec.axis = sec_axis(~.*1/scaler)) +
        scale_x_continuous(sec.axis = sec_axis(~.*1/scaler)) 
    }
    gg_biplot
  }








setClass("employee", slots=list(name="data.frame", id="numeric", contact="character"))


obj <- new("employee",name=cluster, id=1002, contact="West Avenue")
obj


dotplot_internal <- function(object, x = "geneRatio", color = "p.adjust",
                             showCategory=10, size=NULL, split = NULL,
                             font.size=12, title = "", orderBy="x",
                             decreasing=TRUE, label_format = 30) {
  
  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
    if (is.null(size))
      size <- "Count"
  } else if (x == "count" || x == "Count") {
    x <- "Count"
    if (is.null(size))
      size <- "GeneRatio"
  } else if (is(x, "formula")) {
    x <- as.character(x)[2]
    if (is.null(size))
      size <- "Count"
  } else {
    ## message("invalid x, setting to 'GeneRatio' by default")
    ## x <- "GeneRatio"
    ## size <- "Count"
    if (is.null(size))
      size  <- "Count"
  }
  
  df <- fortify(object, showCategory = showCategory, split=split)
  ## already parsed in fortify
  ## df$GeneRatio <- parse_ratio(df$GeneRatio)
  
  if (orderBy !=  'x' && !orderBy %in% colnames(df)) {
    message('wrong orderBy parameter; set to default `orderBy = "x"`')
    orderBy <- "x"
  }
  
  if (orderBy == "x") {
    df <- dplyr::mutate(df, x = eval(parse(text=x)))
  }
  
  label_func <- default_labeller(label_format)
  if(is.function(label_format)) {
    label_func <- label_format
  }
  
  idx <- order(df[[orderBy]], decreasing = decreasing)
  df$Description <- factor(df$Description,
                           levels=rev(unique(df$Description[idx])))
  
  print(df)
  ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
    geom_point() +
    scale_color_continuous(low="red", high="blue", name = color,
                           guide=guide_colorbar(reverse=TRUE)) +
    scale_y_discrete(labels = label_func) +
    ylab(NULL) + ggtitle(title) + theme_dose(font.size) +
    scale_size(range=c(3, 8))
  
}