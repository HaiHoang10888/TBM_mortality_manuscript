
## ======================= 1 ======================= function to do progress bar
f_bar <- function(iterator){
  pb <- txtProgressBar(min = 1, max = iterator - 1, style = 3)
  count <- 0
  function(...) {
    count <<- count + length(list(...)) - 1
    setTxtProgressBar(pb, count)
    flush.console()
    rbind(...) # this can feed into .combine option of foreach
  }
}


## ======================= 2 ======================= lasso cox selection
lasso_cox_boost <- function(x,y, dat_hubgenes, clini_vars =NULL ,evdays=97, n, B,
                            lambda.1se = TRUE,nfolds=10,alpha=1,standardize = T,
                            Multi=T){
  if (Multi==T) {
    Full_x = rbind(x[[1]]$data[dat_hubgenes$gene],
                   x[[2]]$data[dat_hubgenes$gene]) %>%
      as.matrix()
    if(!is.null(clini_vars)) {
      Full_y_var <- rbind(y[[1]]$data[clini_vars],
                          y[[2]]$data[clini_vars])
      Full_x <- cbind(Full_x,Full_y_var) %>%
        as.matrix()
    }
    
    Full_y = rbind(y[[1]]$data,
                   y[[2]]$data) %>%
      mutate(time=ifelse(ttdeath > evdays,evdays,ifelse(ttdeath==0,1,ttdeath)),
             status=ifelse(ttdeath>evdays,0,evdeath)) %>%
      dplyr::select(time,status) %>%
      as.matrix()
  } else {
    # which(dat_hubgenes$gene %nin% colnames(x))
    # class(x)
    Full_x = x[,dat_hubgenes$gene] %>%
      as.matrix()
    if(!is.null(clini_vars)) {
      Full_y_var <-  y[clini_vars]
      Full_x <- cbind(Full_x,Full_y_var) %>%
        as.matrix()
    }
    Full_y = y %>%
      mutate(time=ifelse(ttdeath > evdays,evdays,ifelse(ttdeath==0,1,ttdeath)),
             status=ifelse(ttdeath>evdays,0,evdeath)) %>%
      dplyr::select(time,status) %>%
      as.matrix()
  }

  if (!all(rownames(Full_x)==rownames(Full_y))) stop("rownames exp data differs rownames trait data")
  index <- sample(n, n*B, replace = T)
  dim(index) <- c(n,B) # matrix i=1
  temp <- foreach (i=1:B, .combine=f_bar(B), .packages=c("glmnet")) %dopar% {
    cat("\r", i, "of", B)
    boot.x <- Full_x[index[,i],]
    boot.y <- Full_y[index[,i],]
    cvfit <- cv.glmnet(boot.x, boot.y, family = "cox", type.measure = "C",
                       nfolds=nfolds, alpha = alpha, standardize = standardize)
    if (lambda.1se){
      best.lamb <- cvfit$lambda.1se
    }else{
      best.lamb <- cvfit$lambda.min
    }
    cof <- coef(cvfit, s = best.lamb)
    #return(cof)
  }
  #cof <- cof
  vars <- table(rownames(temp)[which(temp != 0)])
  name.vars <- names(vars)
  Tab1 <- matrix(100 * as.vector(vars) / B, dimnames = list(name.vars, "(%)")) %>%
    as.data.frame() %>%
    mutate(gene=rownames(.),freq=`(%)`) %>%
    dplyr::select(gene,freq) %>% 
    dplyr::arrange(freq,)
  return(Tab1)
}

## ======================= 3 ======================= group lasso logit regression
group_lasso_glm_boost <- function(x,y, dat_hubgenes, evdays=97, n, B,lambda.1se = TRUE,
                                  nfolds=5){
  Full_x = rbind(x[[1]]$data[dat_hubgenes$gene],
                 x[[2]]$data[dat_hubgenes$gene]) %>%
    as.matrix()
  
  Full_y = rbind(y[[1]]$data,
                 y[[2]]$data) %>%
    mutate(time=ifelse(ttdeath > evdays,evdays,ifelse(ttdeath==0,1,ttdeath)),
           status=ifelse(ttdeath>evdays,0,evdeath),
           event = ifelse(status==0,-1,1)) %>%
    dplyr::select(event) %>%
    as.matrix()
  if (!all(rownames(Full_x)==rownames(Full_y))) stop("rownames exp data differs rownames trait data")
  index <- sample(n, n*B, replace = T)
  dim(index) <- c(n,B) # matrix
  cv.group <- as.numeric(factor(dat_hubgenes$module, levels = c("blue","brown","red","cyan"),labels = c(1,2,3,4))) 
  # i=1
  temp <- foreach (i=1:B, .combine=f_bar(B), .packages=c("glmnet","gglasso")) %dopar% {
    cat("\r", i, "of", B)
    boot.x <- Full_x[index[,i],]
    boot.y <- Full_y[index[,i],]
    boot.x.design <- model.matrix(~.,as.data.frame(boot.x))[,-1] # exclude the intercept
    gr_cv <- cv.gglasso(x = boot.x.design,
                        y = boot.y,
                        group=cv.group,
                        loss="logit",
                        pred.loss="loss",
                        intercept = F,
                        nfolds=nfolds)
    if (lambda.1se){
      best.lamb <- gr_cv$lambda.1se
    }else{
      best.lamb <- gr_cv$lambda.min
    }
    cof <- coef(gr_cv, s = best.lamb)
    # return(cof)
  }
  #cof <- cof
  vars <- table(rownames(temp)[which(temp != 0)])
  name.vars <- names(vars)
  Tab1 <- matrix(100 * as.vector(vars) / B, dimnames = list(name.vars, "(%)")) %>%
    as.data.frame() %>%
    mutate(gene=rownames(.),freq=`(%)`) %>%
    dplyr::select(gene,freq) %>% 
    dplyr::arrange(freq)
  return(Tab1)
}

data_cox= NULL
## ======================= 4 ======================= function to draw KM for a gene
export_KM_genes <- function(data, gene, covar=NULL, evdays=97,output,name,x.text=15, y.text=0.2,
                            plot_export=T){
  pheno_index = c(which(names(data) %in% c("ttdeath","evdeath")),
                  which(names(data) %in% gene),
                  which(names(data) %in% covar))
  data <- data[,pheno_index]
  colnames(data)[1:3] <- c("ttdeath","evdeath","gene")
  mean_gene= round(mean(data$gene, na.rm = T),1)
  range_gene= range(data$gene, na.rm = T)
  range_gene = paste0(mean_gene," [",round(range_gene[1],1),";", round(range_gene[2],1),"]")
  data <- data %>% 
    mutate(gene1 = ifelse(gene> quantile(gene, 0.5, na.rm=T), 'high', 'low'),
           gene=factor(gene1,levels = c("low","high")))
  fit <- survfit(Surv(pmin(ttdeath,n),ifelse(ttdeath<=n,evdeath,0))~ gene,data = data)
  plot <- ggsurvplot(fit,data=data,conf.int=TRUE, pval=F, risk.table=TRUE)
  fit_HZ <- coxph(Surv(pmin(ttdeath,n),ifelse(ttdeath<=n,evdeath,0))~ gene,data = data)
  # str <- paste0("fit_HZ <- coxph(Surv(pmin(ttdeath,n),ifelse(ttdeath<=n,evdeath,0))~ ",
  #               paste(c("gene",covar), collapse = " + "), ",data =data)")
  # eval(parse( text=str ))
  Pval <- rstatix::p_round(summary(fit_HZ)$coefficients[1,5])
  HZval <- round(summary(fit_HZ)$coefficients[1,2],1)
  HZval_U <- round(summary(fit_HZ)$conf.int[1,4],1)
  HZval_L <- round(summary(fit_HZ)$conf.int[1,3],1)
  HZval <- paste0(HZval," [",HZval_L,"; ",HZval_U,"]")
  plot$plot <- plot$plot +
    ggplot2::annotate("text",
                      x = x.text, y = y.text, # x and y coordinates of the text
                      label = paste0("HZ = ",HZval," \n ", "p = ", Pval), size = 5)
  if (plot_export==T) {
    jpeg(file =  paste0(output,name,"_KM.jpeg"), wi = 7, 
         he = 7, units = "in",res=480)
    print(plot$plot)
    dev.off()
  }
  
  # 
  # ggsave(paste0(output,name,"_KM.jpeg"), plot,width = 7, height = 7,
  #        units = "in",dpi = 480)
  dat <- data.frame(gene=gene,range_gene=range_gene,HZ=HZval,Pval=Pval)
  return(dat)
}

## ======================= 5 ======================= function to calculate AUC for glm and 1000 boostrap for a model
glm_auc <- function(formula,data, B){
  fit <- glm(formula, data = data, family=binomial(link="logit")) # Fit model on orginal data
  #summary(fit)
  prob_orig <- predict(fit, data, type="response")# predict y from the above model
  pred_orig <- ROCR::prediction(prob_orig, data$evdeath_3M)
  auc_pef_orig <- ROCR::performance(pred_orig, measure="auc")
  auc_orig <- round(auc_pef_orig@y.values[[1]],3)# AUC orginal data
  n <- nrow(data)
  index <- sample(n, n*B, replace = T)
  dim(index) <- c(n,B)
  auc_boot <- NULL
  auc_boot_orig <- NULL
  auc_corrected <- NULL
  o <- NULL
  #data_auc <- NULL ### what is it used for?
  num_warn <- 0
  
  err <- NULL
  i=1
  for (i in 1:B){
    cat("\r", i, "of", B)#### what is this?
    flush.console() ### what is this
    boot.data <- data[index[,i],]
    tryCatch({
      g <- glm(formula, family=binomial(link="logit"), data=boot.data)
      prob <- predict(g, boot.data, type="response")
      pred <- ROCR::prediction(prob, boot.data$evdeath_3M)
      auc_pef <- ROCR::performance(pred, measure="auc")
      auc_boot[i] <- round(auc_pef@y.values[[1]],3)
      
      prob_boot <- predict(g, data, type="response")
      pred_boot <- ROCR::prediction(prob_boot, data$evdeath_3M)
      auc_pef_boot <- ROCR::performance(pred_boot, measure="auc")
      auc_boot_orig[i] <- round(auc_pef_boot@y.values[[1]],3) 
      o[i] <- auc_boot[i] - auc_boot_orig[i]
      auc_corrected[i] <- auc_orig - o[i] 
    },
    warning = function(warn)
    {
      auc_boot[i]<<- NA
      num_warn <<- num_warn + 1
      
    }
    )
  }
  median_cor_auc = median(auc_corrected, na.rm = TRUE)
  mean_cor_auc = mean(auc_corrected, na.rm = TRUE)
  ci.025 = quantile(auc_corrected, 0.025, na.rm = TRUE)
  ci.975 = quantile(auc_corrected, 0.975, na.rm = TRUE)
  num_actual = B - num_warn 
  tab_cor_auc <- data.frame(
    median_cor_auc = median_cor_auc,
    mean_cor_auc = mean_cor_auc,
    ci.025 = ci.025, ## should we convert it to value?
    ci.975 = ci.975,
    number_actual = B-num_warn ## 
  )
  auc_corrected <- auc_corrected ## can we remove it, depends?
  rownames(tab_cor_auc) <- "Value" ## ahihi
  return(list(tab_cor_auc, num_warn , auc_boot_orig, auc_boot, auc_corrected, o))
}


## ======================= 6 ======================= function to plot AUC of a model
Plot_AUC <- function(formula,data,corrected.AUC=NULL, model=NULL,
                     status = "status"){
  glm_m <- glm(formula, family=binomial(link="logit"), data=data)
  data$prob_glm_m <- predict(glm_m, data, type="response")
  colnames(data)[which(colnames(data)==status)] <- "status"
  results <- OptimalCutpoints::optimal.cutpoints(X = "prob_glm_m",
                               status = "status",
                               tag.healthy = "No",
                               direction = c("<"),
                               methods = "Youden",
                               data = data)
  Youden <- as.numeric(results$Youden$Global$optimal.cutoff[[1]][[1]])
  sen <- as.numeric(results$Youden$Global$optimal.cutoff$Se[[1]][[1]])
  spe <- as.numeric(results$Youden$Global$optimal.cutoff$Sp[[1]][[1]])
  roc3 <- data %>%  
    pROC::roc(status, prob_glm_m, plot = TRUE, legacy.axes = TRUE,  xlab = "1 - Specificity",
              ylab = "Sensitivity",
              percent = FALSE, col = "#377eb8", lwd = 5,print.auc = TRUE, print.auc.x = 0.45, ci = TRUE)# -> plot the ROC again 
  roc3 <- pROC::ci.sp(roc3, sensitivities=seq(0, 1, .01), boot.n=1000) #-> define CI based on ROC plot above
  plot(roc3, type="bars", col="#377eb822")
  title(main = paste("AUC of model",model), line = 2.3, cex.main = 1)
  legend(0.8, 0.2, legend = c(paste0("Optimum-corrected AUC: ",corrected.AUC),
                              paste0("Cut-off (Youden index): ", round(Youden,2)), 
                              paste0("Sensitivity: ",round(sen,2)), 
                              paste0("Specificity: ",round(spe,2))),
         cex=0.8,lwd=2)
}


## ======================= 6 ======================= function to do bootstrap calibration 

# formula = evdeath_3M ~ MCEMP1 + NELL2 + ZNF354C + PKD1 + AGE + TBMGRADE + Log_CSF_Lymp + Log_BL_Neu
# data=df_All
# B=10
# 
# test <- ca_bt_glm(formula=formula,data=df_All, B=10)
# summary(model1)
# offset(data$pred)
calibration_glm <- function(formula, data, B=1000, return_slope=T) {
  f <- glm (formula, data = data, family=binomial(link="logit"))
  data$pred <- predict(f, data)
  tryCatch({
    model1 <- glm(evdeath_3M ~ 1 + offset(pred), data = data, family = binomial)
    model2 <- glm(evdeath_3M ~ 1 + pred,  data = data, family = binomial)
    intlarge_origin <- model1$coefficients[[1]]
    int_origin <- model2$coefficients[[1]]
    slope_origin <- model2$coefficients[[2]]
  },
  warning = function(warn1) {
    print(i)
  }
  )
  n <- nrow(data)
  index <- sample(n, n*B, replace = T)
  dim(index) <- c(n,B)
  a1_boot <- NULL
  a2_boot <- NULL
  b1_boot <- NULL
  a1_boot_origin <- NULL
  a2_boot_origin <- NULL
  b1_boot_origin <- NULL
  opti_a1 <- NULL
  opti_a2 <- NULL
  opti_b1 <- NULL
  num_warn1 <- 0
  warn1 <- NULL
  num_warn2 <- 0
  warn2 <- NULL
  num_warn3 <- 0
  warn3 <- NULL
  for (i in 1:B){
    cat("\r", i, "of", B) 
    flush.console()
    boot.data <- data[index[,i],]
    tryCatch({
      g <- glm (formula, data = boot.data, family=binomial(link="logit"))
    },
    warning = function(warn2) {
      num_warn1 <<- num_warn1 + 1 
      warn1[i] <<- i
    }
    )
    tryCatch({
      data$pred1 <- predict(g, data)
      m1_origin <- glm(evdeath_3M ~ 1 + offset(pred1), data = data, family = binomial)
      a1_boot_origin[i] <- m1_origin$coefficients[[1]]
      m2_origin <- glm(evdeath_3M ~ 1 + pred1,  data = data, family = binomial)
      a2_boot_origin[i] <- m2_origin$coefficients[[1]]
      b1_boot_origin[i] <- m2_origin$coefficients[[2]]
    },
    warning = function(warn3) {
      print(i)
      num_warn2 <<- num_warn2 + 1 
      warn2[i] <<- i
    }
    )
    tryCatch({
      boot.data$pred <- predict(g, boot.data)
      m1 <- glm(evdeath_3M ~ 1 + offset(pred), data = boot.data, family = binomial)
      a1_boot[i] <- m1$coefficients[[1]]
      m2 <- glm(evdeath_3M ~ 1 + pred,  data = boot.data, family = binomial)
      a2_boot[i] <- m2$coefficients[[1]]
      b1_boot[i] <- m2$coefficients[[2]]
    },
    warning = function(warn4) {
      print(i)
      num_warn3 <<- num_warn3 + 1 
      warn3[i] <<- i
    }
    )
    opti_a1[i] =  a1_boot[i] - a1_boot_origin[i]
    opti_a2[i] =  a2_boot[i] - a2_boot_origin[i]
    opti_b1[i] =  b1_boot[i] - b1_boot_origin[i]
  }
  intlarge_correct = intlarge_origin - mean(opti_a1, na.rm = TRUE)
  int_correct = int_origin - mean(opti_a2, na.rm = TRUE)
  slope_correct = slope_origin - mean(opti_b1, na.rm = TRUE)
  tab_cal <- data.frame(
    cal_large_corrected = intlarge_correct,
    intercept_corrected = int_correct,
    slope_corrected = slope_correct,
    cal_large_orig = intlarge_origin,
    intercept_orig = int_origin,
    slope_orig = slope_origin
  )
  rownames(tab_cal) <- "Value"
  if (return_slope==T) {
    return(tab_cal)
  } else return(list(results = tab_cal, warn1 = num_warn1, warn2 = num_warn2, warn3 = num_warn3))
  #return(tab_cal)
  # return(list(tab_cal, num_warn1, num_warn2, num_warn3, a1_boot, a2_boot, b1_boot,
  #             a1_boot_origin, a2_boot_origin, b1_boot_origin, opti_a1, opti_a2, opti_b1))
  
}

boot632plus_glm <- function(formula, data, B){
  fit <- glm(formula, data = data, family = binomial(link = "logit"))
  pred_fit <- predict(fit, data, type="response")
  auc_orig <- round(performance(prediction(pred_fit, data$definite_tb), measure="auc")@y.values[[1]],4)
  outcome <- data$definite_tb
  brier_orig <- round(mean((pred_fit - outcome)^2),4) 
  
  
  n <- nrow(data)
  
  S <- matrix(integer(1), nrow=n, ncol=B)
  W <- matrix(TRUE, nrow=n, ncol=B)
  for(i in 1 : B) {
    S[, i] <- s <- sample(n, replace=TRUE)
    W[s, i] <- FALSE  
  }
  nomit <- drop(W %*% rep(1,ncol(W)))  #no. boot samples omitting each obs
  if(min(nomit) == 0)
    stop("not every observation omitted at least once ",
         "in bootstrap samples.\nRe--run with larger B")
  W <- apply(W / nomit, 2, sum) / n
  
  auc_train <- 0
  auc_test <- 0
  brier_train <- 0
  brier_test <- 0
  
  num_warn <- 0
  warn <- NULL
  
  for (i in 1:B){
    cat("\r", i, "of", B) 
    flush.console()
    train <- data[S[,i],]
    test <- data[-S[,i],]
    
    tryCatch({
      g <- glm(formula, family=binomial(link="logit"), data=train)
      pred_train <- predict(g, train, type="response")
      auc_trainj <- round(performance(prediction(pred_train, train$definite_tb), measure="auc")@y.values[[1]],4)
      outcome_train <- train$definite_tb
      brier_trainj<- mean((pred_train - outcome_train)^2) 
      
      pred_test <- predict(g, test, type="response")
      auc_testj <- round(performance(prediction(pred_test, test$definite_tb), measure="auc")@y.values[[1]],4)
      outcome_test <- test$definite_tb
      brier_testj<- mean((pred_test - outcome_test)^2)
      
      wt <- W[i]
      
      na <- is.na(auc_trainj + auc_testj + brier_trainj + brier_testj)
      auc_trainj[na] <- 0
      auc_train <- auc_train + auc_trainj
      auc_testj[na] <- 0
      auc_test <- auc_test + auc_testj*wt
      brier_trainj[na] <- 0
      brier_train <- brier_train + brier_trainj
      brier_testj[na] <- 0
      brier_test <- brier_test + brier_testj*wt
      
    },
    warning = function(warn) {
      num_warn <<- num_warn + 1 
      warn[i] <<- i
    }
    )
    
  }
  num_actual = B - num_warn
  warn <- warn
  
  auc_train <- round(auc_train/num_actual,4)
  auc_test <- round(auc_test,4)
  brier_train <- round(brier_train/num_actual,4)
  brier_test <- round(brier_test,4)
  
  tab_auc <- data.frame(
    name = c("AUC", "Brier score"),
    index_orig = c(auc_orig, brier_orig),
    training = c(auc_train, brier_train),
    test = c(auc_test, brier_test),
    optimism = c(round(0.632*(auc_orig - auc_test),4), round(0.632*(brier_orig - brier_test),4)),
    index_corrected = c(round(auc_orig - 0.632*(auc_orig - auc_test),4), round(brier_orig - 0.632*(brier_orig - brier_test),4)),
    number_actual = c(num_actual, num_actual)
  )
  return(list(tab_auc, warn))
}




## ======================= 7 ======================= function to draw bloxplot for a target genes
Boxplot_for_genes <- function(d,x, y,xlab=NULL,ylab=NULL,title=NULL,
                              scale_y_reverse=F,scale_y_reverse_lim=NULL,
                              label.p.x=NULL,label.p.y=NULL, jitter.w=0.2,jitter.h=0.2){
  colnames(d)[which(colnames(d)==x)] <- "x"
  colnames(d)[which(colnames(d)==y)] <- "y"
  p <- ggplot(d, aes(x = x, y = y)) +
    geom_boxplot(size=1,position=position_dodge(1),width=0.7,
                 fill="white",outlier.colour="red", 
                 outlier.shape=NA,outlier.size=3)+
    geom_jitter(aes(fill=x,colour=x),width=jitter.w, height=jitter.h,size=2,alpha = 0.5)+
    scale_color_manual(values=c("#de2d26","#3182bd","#cc4c02","#00441b","#252525"))+
    #theme_minimal()+
    theme_bw()+
    labs(x = xlab,y = ylab, title = title,subtitle = "")+
    theme(plot.title = element_text(size = 25),
          legend.position = "none",
          legend.text = element_text(size = 20) ,
          legend.title  = element_text(size = 20,color="black",face="bold"),
          axis.text = element_text(size = 25),
          axis.text.x=element_text(color="black",size=13,face="bold"),
          axis.title.x = element_text(size = 20,color="black",face="bold"),
          axis.text.y=element_text(color="black",size=20),
          axis.title.y = element_text(size = 20,color="black",face="bold"),
          axis.title= element_text(size = 30,color="black",face="bold"))+
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)), size = 7,label.x = label.p.x, label.y = label.p.y)
  #stat_compare_means(comparisons = my_comparisons,size = 5)+
  if (scale_y_reverse==T) {
    p <- p + scale_y_reverse(lim=scale_y_reverse_lim)
  }
  return(p)
}

#======================= 8 ======================= Function to draw module network
drawn_network_WGCNA <- function(TOM,gene_list, node_cutoff = NULL , gradient_color=NULL,
                                module_color = NULL,coexp_genes = NULL,PPI_genes = NULL,
                                nodesize=50, Hubsize=75, title = NULL, collection = NULL){
  adj <- TOM
  adj <- adj[gene_list,gene_list]
  diag(adj) <- NA
  #max(adj, na.rm = T)
  adj[adj < quantile(adj, probs = 0.3, na.rm = T)] = 0
  
  #adj[adj != 1] = 0
  color <- module_color 
  network <- graph.adjacency(adj, mode ="undirected", weighted=T)
  if(!is.null(coexp_genes)) {
    label_name <- V(network)$name
    label_name[which(label_name %nin% coexp_genes)] <- ""
    nodesize <- rep(nodesize,length(label_name))
    nodesize[which(label_name %in% coexp_genes)] <- Hubsize
    network <- set.vertex.attribute(network, "label_name", value=label_name)
    network <- set.vertex.attribute(network, "nodesize", value=nodesize)
  }
  # V(network)$name
  ## adding atribute of gradient color
  network <- set.vertex.attribute(network, "gradient_color", value=gradient_color)
  network <- set.vertex.attribute(network, "gradient_color", value=gradient_color)
  V(network)$gradient_color
  network <- simplify(network)  # removes self-loops
  V(network)$color <- color
  
  if (!is.null(node_cutoff)) {
    node_cutoff <- min(node_cutoff,length(degree(network)))
    degree_cutoff <- sort(as.vector(degree(network)),decreasing = T)[node_cutoff] 
    deg_filter <- degree(network) < degree_cutoff
    network <- igraph::delete.vertices(network, deg_filter)
  }
  
  E(network)$width <- E(network)$weight
  library(RCy3)
  cytoscapePing() ## open Cytoscape by cmd Cytoscape &
  createNetworkFromIgraph(network, title = title,
                          collection = paste0("Igraph Collection ", collection))
}


# coexp <- coexp_genes
# PPI <- PPI_genes
# 
# vertex.size[which(names(vertex.size) %in% coexp)] <- 20*2
# vertex.size[which(names(vertex.size) %in% PPI)] <- 20*2
# 
# vertex.label.color[which(names(vertex.label.color)%in% PPI)] <- "red"
# vertex.label.color[which(names(vertex.label.color)%in% coexp)] <- "red"
# 
# #jpeg(file = paste0("Plots_WGCNA_4041/",module,"_network_N.jpeg"), wi = 12, he = 9,units = 'in', res = 300)
# plot(network, edge.arrow.size = 0.2, vertex.size = vertex.size, label.dist = 10,
#      vertex.label.cex	=1.5,
#      vertex.label.color=vertex.label.color,layout=layout.fruchterman.reingold(network))

#======================= 9 ======================= Function to draw correlation heatmap
corrGeneClin <- function(df, h_vars,v_vars,outfile=NULL, cor_method='spearman', 
                         coef_lab=F, size_lab=2, symmetric=FALSE, hc.order = F){
  cor_method <- cor_method
  col=c('#3182bd', 'white', '#de2d26') ###  de2d26","#3182bd c('darkblue', 'white', 'firebrick')
  
  # pval matrix
  cor.pval.test <- function(x1, x2, ...) {
    r <- ncol(x1)
    n <- ncol(x2)
    p.mat<- matrix(NA, r, n)
    for (i in 1:r) {
      for (j in 1:n) {
        tmp <- cor.test(x1[, i], x2[, j], ...)
        p.mat[i, j] <- tmp$p.value
      }
    }
    rownames(p.mat) <-  colnames(x1)
    colnames(p.mat) <- colnames(x2)
    return(p.mat)
  }
  
  # Matrix 1
  which(v_vars %nin% colnames(df))
  mat1 <- cor(df[, h_vars], df[, v_vars],
              use='pairwise.complete.obs', method = cor_method)
  
  # matrix of the p-value of the correlation
  p.mat1 <- cor.pval.test(df[, h_vars], df[, v_vars],
                          method=cor_method, use='pairwise.complete.obs')
  if (symmetric==F) {
    if(coef_lab == T){
      plot1 <-
        ggcorrplot(mat1,
                   p.mat = p.mat1, insig = 'blank',
                   colors = col,
                   outline.color = 'black',
                   lab = T, lab_size = size_lab,
                   hc.order=hc.order)
    } else {
      plot1 <-
        ggcorrplot(mat1,
                   p.mat = p.mat1, insig = 'blank',
                   colors = col,
                   outline.color = 'black',
                   hc.order=hc.order)
    }
  }
  
  if (symmetric==T) {
    if(coef_lab == T){
      plot1 <-
        ggcorrplot(mat1,
                   p.mat = p.mat1, insig = 'blank',
                   colors = col,
                   type = 'lower',
                   outline.color = 'black',
                   show.diag = T, lab = T,
                   lab_size = size_lab, hc.order=hc.order) +
        scale_fill_gradient2(low=col[1], mid= col[2], high=col[3], na.value = 'gray50', limits=c(-1,1))+
        labs(fill="Corr")+
        theme(axis.text.x = element_blank())
    } else {
      plot1 <-
        ggcorrplot(mat1,
                   p.mat = p.mat1, insig = 'blank',
                   colors = col,
                   type = 'lower',
                   outline.color = 'black',
                   show.diag = T, hc.order=hc.order) +
        scale_fill_gradient2(low=col[1], mid= col[2], high=col[3], na.value = 'gray50', limits=c(-1,1))+
        labs(fill="Corr")+
        theme(axis.text.x = element_blank())
    }
  }
  
  
  return(plot1)
}
