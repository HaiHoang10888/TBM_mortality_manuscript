##' Get the similarity matrix
##'
##' @param y A data.frame of enrichment result
##' @param geneSets A list, the names of geneSets are term ids,
##' and every object is a vertor of genes.
##' @param method Method of calculating the similarity between nodes,
##' one of "Resnik", "Lin", "Rel", "Jiang" , "Wang"  and
##' "JC" (Jaccard similarity coefficient) methods
##' @param semData GOSemSimDATA object
##' @noRd
get_similarity_matrix <- function(y, geneSets, method, semData = NULL) {
  id <- y[, "ID"]
  geneSets <- geneSets[id]
  n <- nrow(y)
  y_id <- unlist(strsplit(y$ID[1], ":"))[1]
  ## Choose the method to calculate the similarity
  if (method == "JC") {
    w <- matrix(NA, nrow=n, ncol=n)
    colnames(w) <- rownames(w) <- y$Description
    for (i in seq_len(n-1)) {
      for (j in (i+1):n) {
        w[i,j] <- overlap_ratio(geneSets[id[i]], geneSets[id[j]])
      }
    }
    return(w)
  }
  
  if (y_id == "GO") {
    if(is.null(semData)) {
      stop("The semData parameter is missing,
                and it can be obtained through godata function in GOSemSim package.")
    }
    w <- GOSemSim::mgoSim(id, id, semData=semData, measure=method,
                          combine=NULL)
  }
  
  if (y_id == "DOID") w <- DOSE::doSim(id, id, measure=method)
  rownames(y) <- y$ID
  rownames(w) <- colnames(w) <- y[colnames(w), "Description"]
  return(w)
}


##' Check whether the similarity matrix exists
##'
##' @param x result of enrichment analysis
##'
##' @noRd
has_pairsim <- function(x) {
  if (length(x@termsim) == 0) {
    error_message <- paste("Term similarity matrix not available.",
                           "Please use pairwise_termsim function to",
                           "deal with the results of enrichment analysis.")
    stop(error_message)
  }
  
}


##' Get graph.data.frame() result
##'
##' @importFrom igraph graph.empty
##' @importFrom igraph graph.data.frame
##' @param y A data.frame of enrichment result.
##' @param geneSets A list gene sets with the names of enrichment IDs
##' @param color a string, the column name of y for nodes colours
##' @param cex_line Numeric, scale of line width
##' @param min_edge The minimum similarity threshold for whether 
##' two nodes are connected, should between 0 and 1, default value is 0.2.
##' @param pair_sim Semantic similarity matrix.
##' @param method Method of calculating the similarity between nodes,
##' one of "Resnik", "Lin", "Rel", "Jiang" , "Wang"  and
##' "JC" (Jaccard similarity coefficient) methods
##' @return result of graph.data.frame()
##' @noRd
build_emap_graph <- function(y, geneSets, color, cex_line, min_edge,
                             pair_sim  = NULL, method = NULL) {
  
  if (!is.numeric(min_edge) | min_edge < 0 | min_edge > 1) {
    stop('"min_edge" should be a number between 0 and 1.')
  }
  
  if (is.null(dim(y)) | nrow(y) == 1) {  # when just one node
    g <- graph.empty(0, directed=FALSE)
    g <- add_vertices(g, nv = 1)
    V(g)$name <- as.character(y$Description)
    V(g)$color <- "red"
    return(g)
  } else {
    w <- pair_sim[as.character(y$Description), as.character(y$Description)]
  }
  
  wd <- melt(w)
  wd <- wd[wd[,1] != wd[,2],]
  # remove NA
  wd <- wd[!is.na(wd[,3]),]
  if (method != "JC") {
    # map id to names
    wd[, 1] <- y[wd[, 1], "Description"]
    wd[, 2] <- y[wd[, 2], "Description"]
  }
  
  g <- graph.data.frame(wd[, -3], directed=FALSE)
  E(g)$width <- sqrt(wd[, 3] * 5) * cex_line
  
  # Use similarity as the weight(length) of an edge
  E(g)$weight <- wd[, 3]
  g <- delete.edges(g, E(g)[wd[, 3] < min_edge])
  idx <- unlist(sapply(V(g)$name, function(x) which(x == y$Description)))
  cnt <- sapply(geneSets[idx], length)
  V(g)$size <- cnt
  colVar <- y[idx, color]
  V(g)$color <- colVar
  return(g)
}



##' Get an iGraph object
##'
##' @param x Enrichment result.
##' @param y as.data.frame(x).
##' @param n Number of enriched terms to display.
##' @param color variable that used to color enriched terms, e.g. 'pvalue',
##' 'p.adjust' or 'qvalue'.
##' @param cex_line Scale of line width.
##' @param min_edge The minimum similarity threshold for whether 
##' two nodes are connected, should between 0 and 1, default value is 0.2.
##'
##' @return an iGraph object
##' @noRd
get_igraph <- function(x, y,  n, color, cex_line, min_edge){
  geneSets <- geneInCategory(x) ## use core gene for gsea result
  if (is.numeric(n)) {
    y <- y[1:n, ]
  } else {
    y <- y[match(n, y$Description),]
    n <- length(n)
  }
  
  if (n == 0) {
    stop("no enriched term found...")
  }
  
  g <- build_emap_graph(y = y, geneSets = geneSets, color = color,
                        cex_line = cex_line, min_edge = min_edge,
                        pair_sim = x@termsim, method = x@method)
}


##' Merge the compareClusterResult file
##'
##' @param yy A data.frame of enrichment result.
##'
##' @return a data.frame
##' @noRd
merge_compareClusterResult <- function(yy) {
  yy_union <- yy[!duplicated(yy$ID),]
  yy_ids <- lapply(split(yy, yy$ID), function(x) {
    ids <- unique(unlist(strsplit(x$geneID, "/")))
    cnt <- length(ids)
    list(ID=paste0(ids, collapse="/"), cnt=cnt)
  })
  
  ids <- vapply(yy_ids, function(x) x$ID, character(1))
  cnt <- vapply(yy_ids, function(x) x$cnt, numeric(1))
  
  yy_union$geneID <- ids[yy_union$ID]
  yy_union$Count <- cnt[yy_union$ID]
  yy_union$Cluster <- NULL
  yy_union
}

##' Get the an ggraph object
##'
##' @importFrom ggplot2 ylim
##' @param y A data.frame of enrichment result.
##' @param g An igraph object.
##' @param y_union A data.frame of enrichment result.
##' @param cex_category Numeric, scale of pie plot.
##' @param pie Proportion of clusters in the pie chart, one of 'equal' (default) or 'Count'.
##' @param layout Layout of the map.
##' @noRd
build_ggraph <- function(y, g, y_union, cex_category, pie, layout){
  ## when y just have one line
  if(is.null(dim(y)) | nrow(y) == 1) {
    title <- y$Cluster
    p <- ggraph(g, "tree") + geom_node_point(color="red", size=5 * cex_category) +
      geom_node_text(aes_(label=~name)) + theme_void() +
      labs(title=title)
    return(p)
  }
  
  if(is.null(dim(y_union)) | nrow(y_union) == 1) {
    ##return(ggraph(g) + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name)))
    p <- ggraph(g, "tree")
    ID_Cluster_mat <- prepare_pie_category(y, pie=pie)
    
    ID_Cluster_mat <- cbind(ID_Cluster_mat,1,1,0.1*cex_category)
    colnames(ID_Cluster_mat) <- c(colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],
                                  "x", "y", "radius")
    
    
    p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
                             cols=names(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],
                             color=NA)+
      xlim(-3,3) + ylim(-3,3) + coord_equal()+
      geom_node_text(aes_(label=~name), repel=TRUE) +
      theme_void()+labs(fill = "Cluster")
    return(p)
    
  }
  ggraph(g, layout=layout)
}

##' Keep selected category in enrichment result
##'
##' @param showCategory A number or a vector of terms. If it is a number, 
##' the first n terms will be displayed. If it is a vector of terms, 
##' the selected terms will be displayed.
##' @param x Enrichment result
##' @param split Separate result by 'category' variable.
##' @noRd
get_selected_category <- function(showCategory, x, split) {
  if (is.numeric(showCategory)) {
    y <- fortify(x, showCategory = showCategory,
                 includeAll = TRUE, split = split)
    
  } else {
    y <- as.data.frame(x)
    y <- y[y$Description %in% showCategory, ]
    y <- fortify(y, showCategory=NULL,
                 includeAll = TRUE, split = split)
  }
  y$Cluster <- sub("\n.*", "", y$Cluster)
  return(y)
}

##' Convert a list of gene IDs to igraph object.
##'
##'
##' @title Convert gene IDs to igraph object
##' @param inputList A list of gene IDs.
##' @return A igraph object.
##' @importFrom igraph graph.data.frame
##' @author Guangchuang Yu
##' @noRd
list2graph <- function(inputList) {
  x <- list2df(inputList)
  g <- graph.data.frame(x, directed=FALSE)
  return(g)
}

##' Convert a list of gene IDs to data.frame object.
##'
##'
##' @title Convert gene IDs to data.frame object
##' @param inputList A list of gene IDs
##' @return a data.frame object.
##' @noRd
list2df <- function(inputList) {
  # ldf <- lapply(1:length(inputList), function(i) {
  ldf <- lapply(seq_len(length(inputList)), function(i) {
    data.frame(categoryID=rep(names(inputList[i]),
                              length(inputList[[i]])),
               Gene=inputList[[i]])
  })
  
  do.call('rbind', ldf)
}

##' Get the location of group label
##'
##' @param pdata2 data of a ggraph object
##' @param label_format A numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' @return a data.frame object.
##' @noRd
get_label_location <- function(pdata2, label_format) {
  label_func <- default_labeller(label_format)
  if (is.function(label_format)) {
    label_func <- label_format
  }
  label_x <- stats::aggregate(x ~ color, pdata2, mean)
  label_y <- stats::aggregate(y ~ color, pdata2, mean)
  data.frame(x = label_x$x, y = label_y$y,
             label = label_func(label_x$color))
}

##' Add group label to a ggplot2 object
##'
##' @param repel a logical value, whether to correct the position of the label.
##' @param shadowtext a logical value, whether to use shadow font. 
##' @param p a ggplot2 object.
##' @param label_location a data.frame with the location of group label.
##' @param label_group a numeric value, default size of group label.
##' @param cex_label_group scale of group labels size.
##' @param ... additional parameters.
##' @return a ggplot2 object.
##' @noRd
add_group_label <- function(repel, shadowtext, p, label_location, 
                            label_group, cex_label_group, ...) {
  if (!repel) {
    if (shadowtext) {
      p <- p + geom_shadowtext(data = label_location,
                               aes_(x =~ x, y =~ y, label =~ label), colour = "black",
                               size = label_group * cex_label_group, bg.color = "white", bg.r = 0.1)
    } else {
      p <- p + geom_text(data = label_location,
                         aes_(x =~ x, y =~ y, label =~ label), colour = "black",
                         size = label_group * cex_label_group)
    }
    
    return(p)
  }
  
  if (shadowtext) {
    p <- p + ggrepel::geom_text_repel(data = label_location,
                                      aes_(x =~ x, y =~ y, label =~ label), colour = "black",
                                      size = label_group * cex_label_group, bg.color = "white", bg.r = 0.1,
                                      show.legend = FALSE, ...)
  } else {
    p <- p + ggrepel::geom_text_repel(data = label_location,
                                      aes_(x =~ x, y =~ y, label =~ label), colour = "black",
                                      size = label_group * cex_label_group, 
                                      show.legend = FALSE, ...)
  }
  return(p)   
}

##' Add node label to a ggplot2 object
##'
##' @param p a ggplot2 object.
##' @param data it is uesd as the `data` parameter of function `ggraph::geom_node_text`, a data.frame or NULL.
##' @param label_location a data.frame with the location of group label.
##' @param label_size_node a numeric value to indicate the font size of the node label.
##' @param cex_label_node a numeric value to indicate the scale of node label size.
##' @param shadowtext  a logical value, whether to use shadow font. 
##' @return a ggplot2 object.
##' @noRd
add_node_label <- function(p, data, label_size_node, cex_label_node, shadowtext) {
  if (shadowtext) {
    p <- p + geom_node_text(aes_(label=~name), data = data,
                            size = label_size_node * cex_label_node, bg.color = "white", repel=TRUE)
  } else {
    p <- p + geom_node_text(aes_(label=~name), data = data,
                            size = label_size_node * cex_label_node, repel=TRUE)
  }
  return(p)
}


# has_package <- function(pkg){
# if (!requireNamespace(pkg, quietly  = TRUE)) {
# stop(paste0(pkg, " needed for this function to work. Please install it."),
# call. = FALSE)
# }
# }



##' @method as.data.frame compareClusterResult
##' @export
as.data.frame.compareClusterResult <- function(x, ...) {
  as.data.frame(x@compareClusterResult, ...)
}


##' Prepare pie data for genes in cnetplot.
##' The function only works for compareClusterResult
##'
##' @importFrom DOSE geneID
##' @param y a data.frame converted from compareClusterResult
##' @return a data.frame
##' @noRd
prepare_pie_gene <- function(y) {
  gene_pie <- tibble::as_tibble(y[,c("Cluster", "Description", "geneID")])
  gene_pie$geneID <- strsplit(gene_pie$geneID, '/')
  gene_pie2 <- as.data.frame(tidyr::unnest(gene_pie, cols=geneID))
  gene_pie2 <- unique(gene_pie2)
  prepare_pie_data(gene_pie2, pie =  "equal", type = "gene")
}


##' Prepare pie data for categories in cnetplot/emapplot.
##' The function only works for compareClusterResult
##'
##' @param y a data.frame converted from compareClusterResult
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default)
##' or 'Count'
##' @return a data.frame
##' @noRd
prepare_pie_category <- function(y, pie = "equal") {
  pie <- match.arg(pie, c("equal", "count", "Count"))
  if (pie == "count") pie <- "Count"
  
  pie_data <- y[,c("Cluster", "Description", "Count")]
  pie_data[,"Description"] <- as.character(pie_data[,"Description"])
  prepare_pie_data(pie_data, pie = pie)
}




prepare_pie_data <- function(pie_data, pie = "equal",type = "category") {
  if(type == "category"){
    ID_unique <- unique(pie_data[,2])
  } else {
    ID_unique <- unique(pie_data[,3])
  }
  
  Cluster_unique <- unique(pie_data[,1])
  ID_Cluster_mat <- matrix(0, nrow = length(ID_unique), ncol = length(Cluster_unique))
  rownames(ID_Cluster_mat) <- ID_unique
  colnames(ID_Cluster_mat) <- Cluster_unique
  ID_Cluster_mat <- as.data.frame(ID_Cluster_mat, stringAsFactors = FALSE)
  if(pie == "Count") {
    for(i in seq_len(nrow(pie_data))) {
      ID_Cluster_mat[pie_data[i,2],pie_data[i,1]] <- pie_data[i,3]
    }
    for(kk in seq_len(ncol(ID_Cluster_mat))) {
      ID_Cluster_mat[,kk] <- as.numeric(ID_Cluster_mat[,kk])
    }
    return(ID_Cluster_mat)
  }
  for(i in seq_len(nrow(pie_data))) {
    if(type == "category"){
      ID_Cluster_mat[pie_data[i,2],pie_data[i,1]] <- 1
    } else {
      ID_Cluster_mat[pie_data[i,3],pie_data[i,1]] <- 1
    }
    
  }
  return(ID_Cluster_mat)
}


##' create color palette for continuous data
##'
##'
##' @title color_palette
##' @param colors colors of length >=2
##' @return color vector
##' @export
##' @examples
##' color_palette(c("red", "yellow", "green"))
##' @author guangchuang yu
color_palette <- function(colors) {
  # has_package("grDevices")
  grDevices::colorRampPalette(colors)(n = 299)
}


sig_palette <- color_palette(c("red", "yellow", "blue"))

heatmap_palette <- color_palette(c("red", "yellow", "green"))

overlap_ratio <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x,y)))
}

fc_readable <- function(x, foldChange = NULL) {
  if (is.null(foldChange))
    return(NULL)
  
  if(x@readable) {
    gid <- names(foldChange)
    if (is(x, 'gseaResult')) {
      ii <- gid %in% names(x@geneList)
    } else {
      ii <- gid %in% x@gene
    }
    gid[ii] <- x@gene2Symbol[gid[ii]]
    names(foldChange) <- gid
  }
  return(foldChange)
}

# fc_palette <- function(fc) {
# if (all(fc > 0, na.rm=TRUE)) {
# palette <- color_palette(c("blue", "red"))
# } else if (all(fc < 0, na.rm=TRUE)) {
# palette <- color_palette(c("green", "blue"))
# } else {
## palette <- color_palette(c("darkgreen", "#0AFF34", "#B3B3B3", "#FF6347", "red"))
# }
# return(palette)
# }

update_n <- function(x, showCategory) {
  if (!is.numeric(showCategory)) {
    return(showCategory)
  }
  
  ## geneSets <- geneInCategory(x) ## use core gene for gsea result
  n <- showCategory
  if (nrow(x) < n) {
    n <- nrow(x)
  }
  
  return(n)
}

extract_geneSets <- function(x, n) {
  n <- update_n(x, n)
  geneSets <- geneInCategory(x) ## use core gene for gsea result
  y <- as.data.frame(x)
  geneSets <- geneSets[y$ID]
  names(geneSets) <- y$Description
  if (is.numeric(n)) {
    return(geneSets[1:n])
  }
  return(geneSets[n]) ## if n is a vector of Description
}

##' Internal plot function for plotting compareClusterResult
##'
##'
##' @title plotting-clusterProfile
##' @param clProf.reshape.df data frame of compareCluster result
##' @param x x variable
##' @param type one of dot and bar
##' @param by one of percentage and count
##' @param title graph title
##' @param font.size graph font size
##' @param colorBy one of pvalue or p.adjust
##' @return ggplot object
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 %+%
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 scale_color_continuous
##' @importFrom ggplot2 guide_colorbar
##' @importFrom DOSE theme_dose
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
plotting.clusterProfile <- function(clProf.reshape.df,
                                    x = ~Cluster,
                                    type = "dot",
                                    colorBy = "p.adjust",
                                    by = "geneRatio",
                                    title="",
                                    font.size=12) {
  Description <- Percentage <- Count <- Cluster <- GeneRatio <- p.adjust <- pvalue <- NULL # to satisfy codetools
  if (type == "bar") {
    if (by == "percentage") {
      p <- ggplot(clProf.reshape.df,
                  aes(x=Description, y = Percentage, fill=Cluster))
    } else if (by == "count") {
      p <- ggplot(clProf.reshape.df,
                  aes(x=Description, y = Count, fill=Cluster))
    } else {
      
    }
    p <- p +
      geom_bar() +
      coord_flip()
  }
  # if (type == "dot") {
  #     if (by == "rowPercentage") {
  #         p <- ggplot(clProf.reshape.df,
  #                     aes_(x = x, y = ~Description, size = ~Percentage))
  #     } else if (by == "count") {
  #         p <- ggplot(clProf.reshape.df,
  #                     aes_(x = x, y = ~Description, size = ~Count))
  #     } else if (by == "geneRatio") {
  #         p <- ggplot(clProf.reshape.df,
  #                     aes_(x = x, y = ~Description, size = ~GeneRatio))
  #     } else {
  #         ## nothing here
  #     }
  #     p <- ggplot(clProf.reshape.df,
  #                 aes_(x = x, y = ~Description, size = by))
  #     if (any(colnames(clProf.reshape.df) == colorBy)) {
  #         p <- p +
  #             geom_point() +
  #             aes_string(color=colorBy) +
  #             scale_color_continuous(low="red", high="blue",
  #                                    guide=guide_colorbar(reverse=TRUE))
  #         ## scale_color_gradientn(guide=guide_colorbar(reverse=TRUE), colors = sig_palette)
  #     } else {
  #         p <- p + geom_point(colour="steelblue")
  #     }
  # }
  p <- p + xlab("") + ylab("") + ggtitle(title) +
    theme_dose(font.size)
  ## theme(axis.text.x = element_text(colour="black", size=font.size, vjust = 1)) +
  ##     theme(axis.text.y = element_text(colour="black",
  ##           size=font.size, hjust = 1)) +
  ##               ggtitle(title)+theme_bw()
  ## p <- p + theme(axis.text.x = element_text(angle=angle.axis.x,
  ##                    hjust=hjust.axis.x,
  ##                    vjust=vjust.axis.x))
  return(p)
}




##' Get the distance of the label
##'
##' @param dimension one of 1 and 2
##' @param label_location label_location
##' @noRd
get_label_diss <- function(dimension, label_location) {
  nn <- nrow(label_location)
  label_dis <- matrix(NA, nrow = nn, ncol = nn)
  colnames(label_dis) <- rownames(label_dis) <- label_location$label
  for (i in seq_len(nn - 1)) {
    for (j in (i + 1):nn) {
      label_dis[i ,j] <- label_location[i, dimension] -  label_location[j, dimension]
    }
  }
  label_diss <- reshape2::melt(label_dis)
  label_diss <- label_diss[label_diss[,1] != label_diss[,2], ]
  label_diss <- label_diss[!is.na(label_diss[,3]), ]
  label_diss[, 1] <- as.character(label_diss[, 1])
  label_diss[, 2] <- as.character(label_diss[, 2])
  return(label_diss)
}



# adjust_location <- function(label_location, x_adjust, y_adjust) {
# label_diss_x <- get_label_diss(1, label_location)
# label_diss_y <- get_label_diss(2, label_location)

# label_diss_large <- which(abs(label_diss_y[, 3]) < y_adjust) %>%
# intersect(which(label_diss_y[, 3] > 0)) %>%
# intersect(which(abs(label_diss_x[, 3]) < x_adjust))

# label_diss_small <- which(abs(label_diss_y[, 3]) < y_adjust) %>%
# intersect(which(label_diss_y[, 3] < 0)) %>%
# intersect(which(abs(label_diss_x[, 3]) < x_adjust))

# label_location[label_diss_y[label_diss_large, 1], 2] <- label_location[label_diss_y[label_diss_large, 2], 2] + y_adjust
# label_location[label_diss_y[label_diss_small, 1], 2] <- label_location[label_diss_y[label_diss_small, 2], 2] - y_adjust
# return(label_location)
# }


#' ep_str_wrap internal string wrapping function
#' @param string the string to be wrapped
#' @param width the maximum number of characters before wrapping to a new line
#' @noRd
ep_str_wrap <- function(string, width) {
  x <- gregexpr(' ', string)
  vapply(seq_along(x),
         FUN = function(i) {
           y <- x[[i]]
           n <- nchar(string[i])
           len <- (c(y,n) - c(0, y)) ## length + 1
           idx <- len > width
           j <- which(!idx)
           if (length(j) && max(j) == length(len)) {
             j <- j[-length(j)]
           }
           if (length(j)) {
             idx[j] <- len[j] + len[j+1] > width
           }
           idx <- idx[-length(idx)] ## length - 1
           start <- c(1, y[idx] + 1)
           end <- c(y[idx] - 1, n)
           words <- substring(string[i], start, end)
           paste0(words, collapse="\n")
         },
         FUN.VALUE = character(1)
  )
}

#' default_labeller
#'
#' default labeling function that uses the
#' internal string wrapping function `ep_str_wrap`
#' @noRd
default_labeller <- function(n) {
  function(str){
    str <- gsub("_", " ", str)
    ep_str_wrap(str, n)
  }
}

# x = ego1
# showCategory = 30;color = "p.adjust"; cex_line = 1;with_edge = TRUE;
# nWords = 4; nCluster = NULL; split = NULL; min_edge = 0.2;
# cex_label_group = 1; label_style = "shadowtext"; group_legend = FALSE;
# cex_category = 1; label_format = 30; repel = FALSE; shadowtext = TRUE;
# layout = "nicely"
emapplot_cluster.enrichResult <- function(x, showCategory = 30,
                                          color = "p.adjust", cex_line = 1,
                                          with_edge = TRUE,
                                          nWords = 4, nCluster = NULL,
                                          split = NULL, min_edge = 0.2,
                                          cex_label_group = 1, 
                                          label_style = "shadowtext", 
                                          group_legend = FALSE, cex_category = 1, 
                                          label_format = 30, repel = FALSE, 
                                          shadowtext = TRUE, layout = "nicely", ...){
  
  
  has_pairsim(x)
  label_group <- 3
  n <- update_n(x, showCategory)
  y <- as.data.frame(x)
  
  g <- get_igraph(x=x, y=y, n=n, color=color, cex_line=cex_line,
                  min_edge=min_edge)
  if(n == 1) {
    return(ggraph(g, "tree") + geom_node_point(color="red", size=5) +
             geom_node_text(aes_(label=~name)))
  }
  edgee <- igraph::get.edgelist(g)
  
  #edge_w <- E(g)$weight
  ## set.seed(123)
  #lw <- layout_with_fr(g, weights=edge_w)
  
  #p <- ggraph::ggraph(g, layout=lw)
  p <- ggraph(g, layout = layout)
  # cluster_label1 <- lapply(clusters, function(i){i[order(y[i, "pvalue"])[1]]})
  
  ## Using k-means clustering to group
  pdata2 <- p$data
  dat <- data.frame(x = pdata2$x, y = pdata2$y)
  colnames(pdata2)[5] <- "color2"
  
  if(is.null(nCluster)){
    pdata2$color <- kmeans(dat, floor(sqrt(nrow(dat))))$cluster
  } else {
    if(nCluster > nrow(dat)) nCluster <- nrow(dat)
    pdata2$color <- kmeans(dat, nCluster)$cluster
  }
  
  goid <- y$ID
  cluster_color <- unique(pdata2$color)
  clusters <- lapply(cluster_color, function(i){goid[which(pdata2$color == i)]})
  cluster_label <- sapply(cluster_color, get_wordcloud, pdata2 = pdata2,
                          nWords=nWords)
  names(cluster_label) <- cluster_color
  pdata2$color <- cluster_label[as.character(pdata2$color)]
  p$data <- pdata2
  ## Take the location of each group's center nodes as the location of the label
  # label_func <- default_labeller(label_format)
  # if (is.function(label_format)) {
  #     label_func <- label_format
  # }
  
  # label_x <- stats::aggregate(x ~ color, pdata2, mean)
  # label_y <- stats::aggregate(y ~ color, pdata2, mean)
  # label_location <- data.frame(x = label_x$x, y = label_y$y,
  #                              # label = label_x$color)
  #                              label = label_func(label_x$color))
  label_location <- get_label_location(pdata2, label_format)
  ## Adjust the label position up and down to avoid overlap
  # rownames(label_location) <- label_location$label
  # label_location <- adjust_location(label_location, x_adjust, y_adjust)
  ## use spread.labs
  # label_location$y <- TeachingDemos::spread.labs(x = label_location$y, mindiff = cex_label_group*y_adjust)
  show_legend <- c(group_legend, FALSE)
  names(show_legend) <- c("fill", "color")
  
  if (with_edge) {
    p <-  p +  ggraph::geom_edge_link(alpha = .8,
                                      aes_(width =~ I(width*cex_line)), colour='darkgrey')
  }
  
  if (label_style == "shadowtext") {
    p <- p + ggforce::geom_mark_ellipse(aes_(x =~ x, y =~ y, color =~ color,
                                             fill =~ color), show.legend = show_legend)
  } else {
    p <- p + ggforce::geom_mark_ellipse(aes_(x =~ x, y =~ y, color =~ color,
                                             fill =~ color, label =~ color), show.legend = show_legend)
  }
  
  if (group_legend) p <- p + scale_fill_discrete(name = "groups")
  p <- p + ggnewscale::new_scale_fill() +
    geom_point(shape = 21, aes_(x =~ x, y =~ y, fill =~ color2,
                                size =~ size)) +
    scale_size_continuous(name = "number of genes",
                          range = c(3, 8) * cex_category) +
    scale_fill_continuous(low = "red", high = "blue", name = color,
                          guide = guide_colorbar(reverse = TRUE)) + 
    theme(legend.title = element_text(size = 10),
          legend.text  = element_text(size = 10)) +
    theme(panel.background = element_blank())
  if (label_style == "ggforce") return(p)  
  add_group_label(repel = repel, shadowtext = shadowtext, p = p, 
                  label_location = label_location, label_group = label_group, 
                  cex_label_group = cex_label_group, ...)                 
}



