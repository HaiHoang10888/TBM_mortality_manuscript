tmodLimmaDecideTests
function (f, genes, lfc.thr = 0.5, pval.thr = 0.05, filter.unknown = FALSE, 
          adjust.method = "BH", mset = "all") 
{
  .check.limma(f)
  cn <- colnames(f$coefficients)
  if (is.null(cn)) 
    cn <- 1:ncol(f$coefficients)
  if (length(genes) == 1) {
    if (!genes %in% colnames(f$genes)) 
      stop(sprintf("column %s not in f$genes", genes))
    genes <- f$genes[, genes]
  }
  lfc <- f$coefficients
  pval <- f$p.value
  pval <- apply(pval, 2, function(x) p.adjust(x, method = adjust.method))
  tmodDecideTests(genes, lfc = lfc, pval = pval, lfc.thr = lfc.thr, 
                  pval.thr = pval.thr, labels = cn, filter.unknown = filter.unknown, 
                  mset = mset)
}

#=====================================
#Get Mset data
#=====================================
mset= "DC"
mset <- .getmodules2(NULL, mset, known.only = FALSE)

#=====================================
#Get Reactome pathway of MSET
#=====================================

library(ReactomePA)
set_DC <- "DC.M3.2"

genelist <- paste0("mset$MODULES2GENES$",set_DC)
genelist <- eval(parse( text=genelist ))
genelist <- unname(genelist)

genelist_entrez <- clusterProfiler::bitr(genelist,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F)
genelist_entrez <- subset(genelist_entrez,!is.na(ENTREZID))

x <- enrichPathway(gene=genelist_entrez$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
x_f <- filter(x, p.adjust < .05, qvalue < 0.2)
x_f@result$Description


paste(x_f@result$Description,collapse="; ")


tmodDecideTests
function (g, lfc = NULL, pval = NULL, lfc.thr = 0.5, pval.thr = 0.05, 
          labels = NULL, filter.unknown = FALSE, mset = "all") 
{
  mset <- .getmodules2(NULL, mset, known.only = filter.unknown)
  if (is.null(lfc) && is.null(pval)) {
    N <- sapply(mset$MODULES$ID, function(m) sum(g %in% mset$MODULES2GENES[[m]]))
    return(data.frame(ID = mset$MODULES$ID, N = N))
  }
  if ((!is.null(lfc) && is.null(dim(lfc))) || (!is.null(pval) && 
                                               is.null(dim(pval)))) {
    one.dim <- TRUE
  }
  else {
    one.dim <- FALSE
  }
  nr <- NULL
  nc <- NULL
  if (!is.null(pval)) {
    pval <- as.matrix(pval)
    if (!is.null(colnames(pval)) && is.null(labels)) 
      labels <- colnames(pval)
    nr <- nrow(pval)
    nc <- ncol(pval)
    if (is.null(lfc)) {
      lfc <- matrix(Inf, ncol = nc, nrow = nr)
    }
  }
  if (!is.null(lfc)) {
    lfc <- as.matrix(lfc)
    if (is.null(labels) && !is.null(colnames(lfc))) 
      labels <- colnames(lfc)
    if (is.null(nr)) 
      nr <- nrow(lfc)
    if (is.null(nc)) 
      nc <- ncol(lfc)
    if (is.null(pval)) {
      pval <- matrix(0, ncol = nc, nrow = nr)
    }
  }
  if (length(pval.thr) == 1) 
    pval.thr <- rep(pval.thr, nc)
  if (length(lfc.thr) == 1) 
    lfc.thr <- rep(lfc.thr, nc)
  if (length(pval.thr) != nc) 
    stop(sprintf("pval.thr must have length 1 or %d", 
                 nc))
  if (length(lfc.thr) != nc) 
    stop(sprintf("lfc.thr must have length 1 or %d", 
                 nc))
  if (is.null(labels)) 
    labels <- paste0("X.", 1:nc)
  if (!identical(dim(lfc), dim(pval)) || nrow(lfc) != length(g)) {
    stop("Dimensions of lfc and pval must be identical")
  }
  if (length(labels) != ncol(lfc)) {
    stop("Length of labels must be equal to ncol(lfc) and ncol(pval)")
  }
  lfc.a <- abs(lfc)
  signif.up <- sapply(1:nc, function(i) pval[, i] < pval.thr[i] & 
                        lfc[, i] > lfc.thr[i])
  signif.down <- sapply(1:nc, function(i) pval[, i] < pval.thr[i] & 
                          lfc[, i] < -lfc.thr[i])
  count.m <- function(m) {
    sel <- which(g %in% mset$MODULES2GENES[[m]])
    up <- colSums(signif.up[sel, , drop = F])
    down <- colSums(signif.down[sel, , drop = F])
    N <- length(sel)
    return(cbind(down = down, N = N - (up + down), up = up))
  }
  res <- lapply(mset$MODULES$ID, count.m)
  nmod <- length(res)
  res <- simplify2array(res)
  ret <- lapply(1:nc, function(i) {
    .x <- t(res[i, , ])
    rownames(.x) <- mset$MODULES$ID
    .x
  })
  names(ret) <- labels
  return(ret)
}



mset

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

