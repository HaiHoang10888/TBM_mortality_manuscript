
#=================================================================
#define function to filter genes( by expression and variance)
# get 75% of gene with higher experssion first and then filter by variance
#=================================================================
filter_genes <- function(expr, pct=0.75, apply_vst=FALSE){
  if(nrow(expr) == 0){
    stop("No expression data found.")
  }
  if(pct == 0){
    stop("Argument pct cannot be zero as that would remove all genes.")
  }
  expr <- expr_pct_filter(expr, pct)
  temp <- as.matrix(expr)
  rownames(temp) <- rownames(expr)
  colnames(temp) <- names(expr)
  expr <- temp
  expr_var <- matrixStats::rowVars(expr, na.rm=TRUE)
  expr <- expr[which(expr_var!=0),]
  if (apply_vst){
    expr <- vst(expr)
  }
  return(as.data.frame(expr))
}

# expr A data.frame containing expression values
# n_genes (Optional) Number of genes to be selected
# filter_pval P-value cutoff for gene selection

select_genes <- function(expr, n_genes, filter_pval=0.1){
  expr <- as.matrix(expr)
  expr_var <- matrixStats::rowVars(expr)
  names(expr_var) <- rownames(expr)
  
  mean_var <- mean(expr_var)
  var_var <- var(expr_var)
  
  ah <- mean_var^2/var_var + 2
  bh <- mean_var*(ah - 1)
  
  p <- sapply(expr_var, function(x) {
    ig <- pracma::gammainc(bh/x, ah)['uppinc']
    g <- gamma(ah)
    return(1 - ig/g)
  })
  
  names(p) <- gsub('.uppinc', '', names(p))
  
  if(!missing(n_genes)){
    filter_pval <- sort(p)[n_genes]
    names(filter_pval) <- NULL
  }
  
  selected <- which(p <= filter_pval)
  if(length(selected) == 0) warning("No gene left after filtering!")
  return(names(selected))
}

# Filter genes based on expression.
# pct percentage of most expressed genes to maintain

expr_pct_filter <- function(expr, pct=0.75){
  rows <- floor(pct*nrow(expr))
  val <- apply(expr, 1, mean)
  sel_rows <- order(val, decreasing=TRUE)[1:rows]
  expr <- expr[sel_rows, ]
  return(expr)
}