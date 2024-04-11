library(stringr)
gitdir <- paste0(strsplit(getwd(), split = "/github_repos", fixed=T)[[1]][1], "/github_repos/")

source(paste0(gitdir,"/R_functions/WGCNA_CEMitool_ORA_function_5May2022.R"))
source(paste0(gitdir,"/R_functions/WGCNA_CEMitool_Filter_genes_May2022.R"))
source(paste0(gitdir,"/R_functions/WGCNA_CEMitool_Module_functions_May2022.R"))


get_GS_WGCNA <- function(Expr,pheno,pheno_index,fit_type="lm"){
  y = as.data.frame(pheno[,pheno_index])
  names(y) = "y"
  n= ncol(Expr)
  data_output <- data.frame(Gene=rep(NA,n),GeneSignificance= rep(NA,n), GeneEffect=rep(NA,n))
  # i=1
  for (i in 1:ncol(Expr)){
    gene_name <- names(Expr)[i]
    data_fit <- data.frame(gene = Expr[,i], y=y)
    if(fit_type=="lm") {
      fit <- lm(y ~ gene, data=data_fit)
      data_output$Gene[i] <- gene_name
      data_output$GeneSignificance[i] <- summary(fit)$coefficients[2,4]
      data_output$GeneEffect[i] <- summary(fit)$coefficients[2,1]
    }
    if(fit_type=="glm") {
      fit <- glm(y ~ gene, data=data_fit, family = "binomial")
      data_output$Gene[i] <- gene_name
      data_output$GeneSignificance[i] <- summary(fit)$coefficients[2,4]
      data_output$GeneEffect[i] <- plogis(coef(fit))[2]
    }
  }
  return(data_output)
}

## get p value and HR
get_GSCox_WGCNA <- function(Expr,pheno,pheno_index,covar =NULL, evdays=97,fit_type = "cox"){
  names_y <- c("ttdeath","evdeath")
  if (!is.null(covar)) {
    covar_names <- c(paste0("cov",1:covar))
    names_y <- c(c("ttdeath","evdeath"),covar_names)
  }
  y = as.data.frame(pheno[,pheno_index])
  names(y) = names_y
  n= ncol(Expr)
  data_output <- data.frame(Gene=rep(NA,n),GeneSignificance= rep(NA,n), GeneEffect=rep(NA,n))
  # i=1
  for (i in 1:ncol(Expr)){
    gene_name <- names(Expr)[i]
    data_fit <- data.frame(gene = Expr[,i], y)
    if (is.null(covar)) {
      fit <- coxph(Surv(pmin(ttdeath,evdays),ifelse(ttdeath<=evdays,evdeath,0))~gene, data =data_fit)
    } else if (covar==1) {
      fit <- coxph(Surv(pmin(ttdeath,evdays),ifelse(ttdeath<=evdays,evdeath,0))~gene + cov1, data =data_fit)
    } else if (covar==2) {
      fit <- coxph(Surv(pmin(ttdeath,evdays),ifelse(ttdeath<=evdays,evdeath,0))~gene + cov1 + cov2, data =data_fit)
    } else if (covar==3) {
      fit <- coxph(Surv(pmin(ttdeath,evdays),ifelse(ttdeath<=evdays,evdeath,0))~gene + cov1 + cov2 + cov3, data =data_fit)
    }
    data_output$Gene[i] <- gene_name
    data_output$GeneSignificance[i] <- summary(fit)$coefficients[1,5]
    data_output$GeneEffect[i] <- round(summary(fit)$coefficients[1,2],2)
  }
  return(data_output)
}

## get p value and HR by gene high low
get_GSCox_level <- function(Expr,pheno,pheno_index,covar =NULL, evdays=97,fit_type = "cox"){
  names_y <- c("ttdeath","evdeath")
  if (!is.null(covar)) {
    covar_names <- c(paste0("cov",1:covar))
    names_y <- c(c("ttdeath","evdeath"),covar_names)
  }
  y = as.data.frame(pheno[,pheno_index])
  names(y) = names_y
  n= ncol(Expr)
  data_output <- data.frame(Gene=rep(NA,n),GeneSignificance= rep(NA,n), GeneEffect=rep(NA,n))
  # i=1
  for (i in 1:ncol(Expr)){
    gene_name <- names(Expr)[i]
    data_fit <- data.frame(gene = Expr[,i], y)
    data_fit <- data_fit %>% 
      mutate(gene1 = ifelse(gene> quantile(gene, 0.5, na.rm=T), 'high', 'low'),
             gene=factor(gene1,levels = c("low","high")))
    if (is.null(covar)) {
      fit <- coxph(Surv(pmin(ttdeath,evdays),ifelse(ttdeath<=evdays,evdeath,0))~gene, data =data_fit)
    } else if (covar==1) {
      fit <- coxph(Surv(pmin(ttdeath,evdays),ifelse(ttdeath<=evdays,evdeath,0))~gene + cov1, data =data_fit)
    } else if (covar==2) {
      fit <- coxph(Surv(pmin(ttdeath,evdays),ifelse(ttdeath<=evdays,evdeath,0))~gene + cov1 + cov2, data =data_fit)
    }
    data_output$Gene[i] <- gene_name
    data_output$GeneSignificance[i] <- summary(fit)$coefficients[1,5]
    data_output$GeneEffect[i] <- round(summary(fit)$coefficients[1,2],2)
  }
  return(data_output)
}



