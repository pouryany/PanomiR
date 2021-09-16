
# The function to generate a table of pathways and genes associations
# takes a list of pathways each element is a list of genes
# returns a table of pathway-gene members

# necessary to process pathway activity profiles.
# There would be a premade Pathway_Gene_Tab going out with the package.

Pathway_Gene_Tab <- function(path.address = F,
                             out.dir  = F){
    pathList  <- readRDS(path.address)
    pathList2 <- lapply(pathList,
                        function(X){clusterProfiler::bitr(X,"ENTREZID",
                                                          "ENSEMBL",
                                                          OrgDb =  org.Hs.eg.db)})
    
    temp       <- lapply(names(pathList2),
                         function(X){data.frame(Pathway = (X), 
                                                pathList2[[X]])}
    )
    temp       <- do.call(rbind,(temp))
    PathExpTab <- as_tibble(temp)
    saveRDS(PathExpTab,out.dir)
}



# The function to generate pathway summary statistics 
# Requires an expression matrix. 
# Pathway reference is pathGene Tab


Path_Summary <- function(exprs.mat, 
                         pathway.ref,
                         id, 
                         z.normalize=F,
                         method = F,
                         de.genes = NULL,
                         trim = 0, 
                         t.scores = NULL){
    
    # There is a confusion about the data format here. 
    # Make sure it is consistent.
    # Current version only works with ENSEMBL.
    exprs.mat  <- rownames_to_column(as.data.frame(exprs.mat), var = id)
    
    if (!is.null(de.genes)) {
      
        if(is.null(t.scores))
          stop("Provide tscores/pvalues")
      
        pathway.ref  <- dplyr::inner_join(pathway.ref,t.scores, 
                                          by = c("ENSEMBL"))
        
        pathway.ref  %<>%  dplyr::group_by(.,Pathway) %>%
            dplyr::filter(., abs(t) >= median(abs(t))) %>%
            dplyr::select(.,-t)  
    }
    
    
    if(method == "none"){
      
        exprs.mat  <- exprs.mat  %>% 
            mutate_if(.,is.numeric,function(X){X})
        
    } else if (method == "x"){
      
        exprs.mat      <- exprs.mat  %>% mutate_if(.,is.numeric,rank) %>%
            mutate_if(.,is.numeric,function(X){X}) 
        
    } else if (method == "x2"){
      
        exprs.mat      <- exprs.mat  %>% mutate_if(.,is.numeric,rank) %>%
            mutate_if(.,is.numeric,function(X){X * X})
        
    } else  {
        stop("invalid choice of summarization function")
    }
    
    
    PathExpTab    <- inner_join(pathway.ref,exprs.mat,by = id)
    
    PathExpTab    <- PathExpTab %>% group_by(., Pathway) %>% 
        dplyr::summarise_if(is.numeric,mean, na.rm=T,trim = trim)
    
    
    PathExp           <- as.data.frame(PathExpTab[,-1])
    rownames(PathExp) <- as.character(pull(PathExpTab[,1]))
    
    if(z.normalize){
        PathExp <- apply(PathExp, 2, function(X){(X - mean(X))/sd(X)})
    }
    
    return(PathExp)
}











# Takes two dataframes and return the correlation between their rows Cartesian

cor.pval.table <- function(df1, df2, df.name1 = "Pathway",
                           df.name2 = "miRNA", method = "pearson",qval = T,
                           pre.pairs = NULL){
    # Function for calculating correlation among rows of two dataframes.
    # Particular application here is to find correlation among a bunch of 
    # miRNAs and Pathway expression profiles. 
    
    if(!all(colnames(df1) == colnames(df2))){
        stop("Mismatch in colnames of the dataframes")
    }
    if(is.null(pre.pairs)){
        cor.data.frame <- expand.grid( rownames(df1),
                                       rownames(df2))
        colnames(cor.data.frame) <- c(df.name1,df.name2)
        cor.data.frame <- cor.data.frame %>% mutate(.,cor= NA, pval =NA) 
    }else if(!is.null(pre.pairs)){
        if(class(pre.pairs) != "data.frame")
            warning("The provided preset pairs is not a dataframe")
        
        if(ncol(pre.pairs) < 2)
            stop("Not enough columns in the preset pairs")
        if(!all(pre.pairs[,1] %in% rownames(df1)))
            stop("Nonexistent pairs. \n")
        if(!all(pre.pairs[,2] %in% rownames(df2)))
            stop("Nonexistent pairs. \n")
        
        cor.data.frame <- pre.pairs
        colnames(cor.data.frame) <- c(df.name1,df.name2)
    }
    # cor.data.frame <- expand.grid( rownames(df1),
    #                                rownames(df2))
    # colnames(cor.data.frame) <- c(df.name1,df.name2)
    cor.data.frame <- cor.data.frame %>% mutate(.,cor= NA, pval =NA)
    
    temp <- apply(cor.data.frame,1,function(X){
        path <- as.character(unlist(X[(df.name1)]))
        mir  <- as.character(unlist(X[(df.name2)]))
        pval <- cor.test(df1[path,],
                         df2[mir,],method = method)[["p.value"]]
        cor  <- cor.test(df1[path,],
                         df2[mir,],method = method)[["estimate"]]
        return(c(cor,pval))
    })
    cor.data.frame[,c(3,4)] <- t(temp)
    cor.data.frame <- cor.data.frame %>% 
        dplyr::mutate(., adj = p.adjust(pval,method = "fdr"))
    
    if(qval){
        adjs  <- qvalue::qvalue(as.vector(cor.data.frame$pval))
        adjs  <- adjs$qvalues
        cor.data.frame$adj <- adjs  
    }
    
    return(cor.data.frame)
}





cor.est.table  <- function(df1, df2, df.name1 = "Pathway",
                           df.name2 = "miRNA", method = "pearson"){
    # Function for calculating correlation among rows of two dataframes.
    # Particular application here is to find correlation among a bunch of 
    # miRNAs and Pathway expression profiles. 
    if(!all(colnames(df1) == colnames(df2))){
        stop("Mismatch in colnames of the dataframes")
    }
    
    cor.data.frame <- expand.grid( rownames(df1),
                                   rownames(df2))
    colnames(cor.data.frame) <- c(df.name1,df.name2)
    cor.data.frame <- cor.data.frame %>% mutate(.,cor= NA)
    
    temp           <- apply(cor.data.frame,1,function(X){
        path <- as.character(unlist(X[(df.name1)]))
        mir  <- as.character(unlist(X[(df.name2)]))
        cor <- cor(df1[path,], df2[mir,],method = method)
        return(cor)
    })
    cor.data.frame[,3] <- (temp)
    return(cor.data.frame)
}





# Wrapper function for correlation estimation

cor.table      <- function(df1, df2, df.name1 = "Pathway",df.name2 = "miRNA",
                           method = "pearson",pvals = T,qval = T, 
                           pre.pairs = NULL){
    if(pvals){
        return(cor.pval.table(df1 = df1,df2 = df2,method = method,
                              df.name1 = df.name1,df.name2 = df.name2,
                              qval = qval,pre.pairs = pre.pairs ))
    } else{
        return(cor.est.table(df1, df2, df.name1,df.name2, method))
    }
}






cor.dif  <- function(df1,df2,lab.vec, df.name1 = "Pathway",
                     df.name2 = "miRNA", method = "pearson",
                     verbose = F, as.df = F,pre.pairs = NULL){
    # df1 is the Pathway expression
    # df2 is the microRNA expression matrix
    
    unique <- unique(lab.vec)
    unique <- unique[order(unique)]
    if(length(unique) != 2){
        stop("The labels vector must contain two distinct types")
    }
    include <- ifelse(lab.vec == unique[1],T,F)
    
    cor.1   <- cor.table(df1 = df1[,include],df2 = df2[,include],
                         df.name1 = df.name1, df.name2 = df.name2,
                         method = method, pvals = F,pre.pairs = pre.pairs)
    cor.2   <- cor.table(df1 = df1[,!include], df2 = df2[,!include], 
                         df.name1 = df.name1, df.name2 = df.name2,
                         method = method, pvals = F,pre.pairs = pre.pairs)
    
    if(verbose){
        print(paste0("Calculating correlation difference between ", unique[1], 
                     " and ", unique[2]))
    }
    if(as.df){
        cor.1$cor.2  <-  cor.2$cor
        cor.1$difCor <-  cor.1$cor - cor.1$cor.2
        return(cor.1)
    } else {
        return(cor.1$cor - cor.2$cor) 
    }
    
}





cor.null <- function(df1,df2,df.name1 = "Pathway",df.name2 = "miRNA",lab.vec,
                     iter = 1000,cores =1, method = "pearson",pre.pairs = NULL,
                     verbose = F){
    
    perms <- mclapply(1:iter,function(X){
        set.seed(X)
        perm <- sample(1:length(lab.vec))
        perm <- lab.vec[perm]
        return(cor.dif(df1 = df1,df2 =  df2,lab.vec =  perm,df.name1 = df.name1,
                       df.name2 = df.name2,as.df = F,
                       verbose = verbose, method = method,
                       pre.pairs = pre.pairs))
    },mc.cores = cores)
    return(perms)
}

get.cor.dif.pval <- function(df1,df2,lab.vec, df.name1 = "Pathway",
                             df.name2 = "miRNA", method = "pearson",
                             pre.pairs = NULL,
                             iter = 1000, cores =1, verbose = F){
    
    nulls    <- cor.null(df1 = df1,df2 = df2, lab.vec = lab.vec,iter = iter,
                         cores = cores,verbose = verbose, method = method,
                         pre.pairs = pre.pairs)
    nulls    <- do.call(cbind,nulls)
    observed <- cor.dif(df1 = df1,df2 = df2,lab.vec = lab.vec,
                        df.name1 = df.name1, df.name2 = df.name2,
                        as.df = F,pre.pairs = pre.pairs,
                        method = method,verbose = verbose)
    hits.pos <- apply(nulls, 2, function(X){return(X > observed)})
    hits.neg <- apply(nulls, 2, function(X){return(X < observed)})
    
    hits     <- vector("numeric",length(observed))
    hits[observed > 0]     <- (rowSums(hits.pos)/iter)[observed > 0]  
    hits[observed <= 0]    <- (rowSums(hits.neg)/iter)[observed <= 0]  
    hits[hits == 0]        <- 1/iter
    
    if(is.null(pre.pairs)){
        cor.data.frame <- expand.grid(rownames(df1), rownames(df2))
    } else if(!is.null(pre.pairs)){
        if(!all(pre.pairs[,1] %in% rownames(df1)))
            stop("Nonexistent pairs. \n")
        if(!all(pre.pairs[,2] %in% rownames(df2)))
            stop("Nonexistent pairs. \n")
        cor.data.frame <- pre.pairs
    }
    
    
    colnames(cor.data.frame) <- c(df.name1,df.name2)
    cor.data.frame <- cor.data.frame %>% mutate(., pval =NA)
    cor.data.frame[,3] <- (hits)
    return(cor.data.frame)
    
}



###==========================================================
#-----
#-----
#----- The following are miRNA-Pathway pval aggregations
#-----
#-----
###==========================================================



pCut.fn <- function(enriches, pathways, is.selector, thresh=0.05){
  
  if (is.selector==T){
    enriches <- enriches %>% mutate(.,hit2=ifelse(pval < thresh,1,0))
  }
  
  temp.enrich <- enriches[enriches$y %in% pathways, ]
  selector <- temp.enrich %>% group_by(x)  %>% 
              dplyr::summarise(n = n(),k = sum(hit2)) %>% 
              arrange(.,x)
  
  if (is.selector==T){
    selector <- selector %>% filter(.,k > thresh * length(pathways))
    return(list('selector'=selector, 'enriches0'=enriches))
  } else {
    return(selector)
  }
}




AggInv.fn <- function(enriches, pathways, is.selector = TRUE, thresh=NULL){
  
  if (is.selector==T){
    
    enriches <- enriches %>% mutate(., ES2 =  qnorm(1 - pval))
    min.es   <- min(enriches$ES2[!is.infinite(enriches$ES2)])
    enriches <- enriches %>% mutate(., ES2 = ifelse(is.infinite(.$ES2),
                                                    min.es,
                                                    .$ES2)
                                    )
  }
  
  temp.enrich <- enriches[enriches$y %in% pathways, ]
  selector <- temp.enrich %>%
    group_by(x) %>%
    dplyr::summarise(n = n(), k = mean(ES2))%>%
    arrange(., x)
  
  if (is.selector==T){
    return(list('selector'=selector, 'enriches0'=enriches))
  } else {
    return(selector)
  }
}




AggLog.fn <- function(enriches, pathways, is.selector, thresh=0.1){
 
  enriches <- enriches %>% mutate(., ES =  -log(pval))
  if (is.selector==T){
    enriches <- enriches %>% mutate(., ES =  -log(pval))
  }
  
  temp.enrich <- enriches[enriches$y %in% pathways, ]
  selector <- temp.enrich %>%
    group_by(x)  %>%
    dplyr::summarise(n = n(), k = mean(ES))
  
  if (is.selector==T){
    selector <- selector %>% filter(.,k*n > thresh * length(pathways))
    return(list('selector'=selector, 'enriches0'=enriches))
  } else {
    return(selector)
  }
}



sumz.fn <- function(enriches, pathways, is.selector, thresh=NULL){
  enriches1 <- enriches %>% mutate(., pval =  ifelse(pval >= 0.999, 0.999, pval))
  enriches1 <- enriches1 %>% mutate(., pval =  ifelse(pval <= 1.0e-16, 1.0e-16, pval))
  
  temp.enrich <- enriches1[enriches1$y %in% pathways,]
  agg.p.tab <- vector()
  for (i in unique(temp.enrich$x)) {
    temp <- temp.enrich[temp.enrich$x == i, ]
    t.pval <- metap::sumz(temp$pval)
    #print(paste0(i, ": ", t.pval$p))
    agg.p.tab <- rbind(agg.p.tab, c(i, t.pval$p, nrow(temp)))    
  }
  
  selector <- tibble(
    "x" = agg.p.tab[, 1],
    "k" = signif(as.numeric(agg.p.tab[, 2]), 4),
    "n" =  agg.p.tab[, 3]
  )
  
  if (is.selector==T){
    return(list('selector'=selector, 'enriches0'=enriches))
  } else {
    return(selector)
  }
}





sumlog.fn <- function(enriches, pathways, is.selector, thresh=NULL){
  enriches1 <- enriches %>% mutate(., pval =  ifelse(pval >= 0.999, 0.999, pval))
  enriches1 <- enriches1 %>% mutate(., pval =  ifelse(pval <= 1.0e-16, 1.0e-16, pval))
  
  temp.enrich <- enriches1[enriches1$y %in% pathways,]
  agg.p.tab <- vector()
  for (i in unique(temp.enrich$x)) {
    temp <- temp.enrich[temp.enrich$x == i, ]
    t.pval <- metap::sumlog(temp$pval)
    #print(paste0(i, ": ", t.pval$p))
    agg.p.tab <- rbind(agg.p.tab, c(i, t.pval$p, nrow(temp)))    
  }
  
  selector <- tibble(
    "x" = agg.p.tab[, 1],
    "k" = signif(as.numeric(agg.p.tab[, 2]), 4),
    "n" =  agg.p.tab[, 3]
  )
  
  if (is.selector==T){
    return(list('selector'=selector, 'enriches0'=enriches))
  } else {
    return(selector)
  }
}



lancaster.fn <- function(enriches, pathways, is.selector, thresh=NULL){
  temp.enrich <- enriches[enriches$y %in% pathways, ]
  selector <- temp.enrich %>%
    group_by(x) %>%
    dplyr::summarise(n = n(), k = lancaster(pval,weight)/n())%>%
    arrange(., x)
  
  if (is.selector==T){
    return(list('selector'=selector, 'enriches0'=enriches))
  } else {
    return(selector)
  }
}



pCut.cover.fn <- function(selector, cover.name) {
  selector <- selector %>%
    dplyr::mutate(., !!cover.name := k/n)
  return(selector)
}

AggInv.cover.fn <- function(selector, cover.name) {
  selector <- selector %>%
    dplyr::mutate(., !!cover.name := k)
  return(selector)
}

AggLog.cover.fn <- AggInv.cover.fn

sumz.cover.fn <- AggInv.cover.fn

lancaster.cover.fn <- AggInv.cover.fn

#### Working

#' Outputs a table of sampling data(rows are miRNA and cols are samples)
#' 
#' @param enrich.null Enrichment dataset with x (miRNA), y (pathway) and pval (probability of observing x in pathway cluster).
#' @param selector Table with x(miRNA) in pathway cluster.
#' @param samp.rate Sampling rate.
#' @param fn Methodology function.
#' @param n_paths Number of pathways in pathway cluster.
#' @param sampling.data.file If file exists, load file. Else, perform random sampling
#' @param save.sampling If TRUE, data is saved.
#' @param jack.knife If TRUE, conduct sampling with one less pathway, used for jack knifing
#' @param num.cores number of cores used
#' @return Outputs of sampling data.



sampling.data.base2 <- function(enrich.null,
                               selector,
                               samp.rate,
                               fn,
                               n_paths,
                               sampling.data.file,
                               jack.knife=FALSE,
                               save.sampling,
                               num.cores=1){
  if(!all(hasName(selector,c("x"))))
    stop("The selector table needs a column x (miRNA name)")
  if(!all(hasName(enrich.null,c("x","y","pval"))))
    stop(paste0("The enrichment table needs a column x (miRNA name)",
                ",a column y (pathway name), and a pval column"))
  
  if (!file.exists(sampling.data.file)){
    all.paths   <- unique(enrich.null$y)
    #temp.n_paths <- n_paths
    
    if (jack.knife==TRUE) {temp.n_paths <- n_paths-1}
    samp.size.vec <- c(n_paths,n_paths-1,100,50)
    
    out.list <- list()
    for (temp.n_paths in samp.size.vec){
      
      temp   <- mclapply(1:(samp.rate), function(Y){
        set.seed(Y)        
        null.paths <- sample(all.paths,temp.n_paths,replace = F)
        sel.null   <- fn(enriches=enrich.null, pathways=null.paths, is.selector=F)
        return(sel.null$k)
        
      }, mc.cores = num.cores )
    # build null distribution of K 
    
    
    
    temp <- do.call(rbind, temp)
    temp <-  t(temp)
    
    rownames(temp) <- selector$x
    colnames(temp) <- sapply(1:(samp.rate), function(Y){paste0("sample_", Y)})
    samp.tag <- paste0("SampSize_",temp.n_paths)
    
    out.list[[samp.tag]] <- temp
    }
    if (save.sampling==TRUE) {
      saveRDS(out.list, file=sampling.data.file)
      print(paste0(sampling.data.file, " saved."))
    }
    
  }
  else{
    print(paste0("Skipping sampling, ", sampling.data.file, " exists."))
    out.list <- readRDS(sampling.data.file)
    print(paste0(sampling.data.file, " loaded."))
  }
  
  return(out.list)
  
}




#' Outputs a table with col x (miRNA), probability of observing k (depending on methodology) against a random distribution and cover of methodology
#' 
#' @param sampling.data Random distribution data
#' @param selector Table with x(miRNA) in pathway cluster and observed k (depending on methodology).
#' @param m method name.
#' @param cover.fn Cover of methodology function.
#' @return Outputs a new selector table with col x, pval and cover.

method.prob.base2 <- function(sampling.data, selector,m, n_paths = 100,cover.fn=NULL)
  {
  
  if(!all(hasName(selector,c("x","k"))))
    stop("The selector table needs a column x (miRNA name) and a column k (miRNA hits)")
  
  # obtain means and sds for distribution, assume CLT
  means <- rowMeans(sampling.data) 
  sds   <- apply(sampling.data, 1, sd)
  sds   <- sds *10/sqrt(n_paths)
  
  pval.name <- paste0(m, '_pval')
  cover.name <- paste0(m, '_cover')
  
  # obtain p-vals
  p_vals <- pnorm(selector$k, mean=means, sd=sds, lower.tail=FALSE)
  selector <- selector %>%
    dplyr::mutate(., !!pval.name := p_vals) %>%
    cover.fn(., cover.name) %>%
    dplyr::select(.,-c(k,n))
  return(selector)
}




#' Outputs a table with col x (miRNA), probability of observing k (depending on methodology) against a random distribution with jack-knifing of the pathway cluster (removing a pathway at a time)
#' 
#' @param selector Table with x(miRNA) in pathway cluster and observed k (depending on methodology).
#' @param pathways Pathways in pathway cluster.
#' @param enrich.null Enrichment dataset with x (miRNA), y (pathway) and pval (probability of observing x in pathway cluster).
#' @param fn Methodology function.
#' @param jack.knife.data Random distribution data with jack-knifing (i.e. one less pathway)
#' @param m method name
#' @param num.cores number of cores
#' @return Outputs a new selector table with col x, pval_jk
#' refer to 05-02-test.R


jack.knife.base2 <- function(selector, pathways, enrich.null, fn, jack.knife.data, m, num.cores=1){
  
  # obtain means and sds for distribution, assume CLT
  n_paths      <- length(pathways)
  sample.means <- rowMeans(jack.knife.data) 
  sample.sds   <- apply(jack.knife.data, 1, sd)
  sample.sds   <- sample.sds *10/sqrt(n_paths-1)
  
  # remove one pathway at a time and obtain K for each miRNA 
  temp1 <- mclapply(1:length(pathways), function(X){
    temp.pathways  <- pathways[-X]
    temp.selector   <- fn(enriches=enrich.null, pathways=temp.pathways, is.selector=F)
    # obtain p-values using the means and sds obtain above
    p_vals <- pnorm(temp.selector$k, mean=sample.means, sd=sample.sds, lower.tail = FALSE)
    return(p_vals)
  }, mc.cores = num.cores)
  
  
  # rows <- number of pathways; col <- number of miRNAs
  # contain p-values
  temp1 <- do.call(rbind, temp1)
  temp1 <- t(temp1)
  
  ## added from here
  
  jack.knife.name <- paste0(m, "_pval_jk")
  
  # obtain aggregate p-values
  means <- rowMeans(temp1)
  selector <- selector %>%
    dplyr::mutate(., !!jack.knife.name := means)
  
  return(selector)
}











# Function imported from https://github.com/pouryany/CovariateAnalysis/blob/master/R/getDesignMatrix.R

# Function to optain desing matrix (modified from covairates pipeline of Menachem Former)
getDesignMatrix <- function(covariatesDataFrame, Intercept = T, RELEVELS=list()) {
  
  ROWNAMES = rownames(covariatesDataFrame)
  COLNAMES = colnames(covariatesDataFrame)
  
  FACTOR_COVARIATE_NAMES <- names(covariatesDataFrame)[sapply(covariatesDataFrame,is.factor)]
  FACTOR_COVARIATE_NAMES = setdiff(FACTOR_COVARIATE_NAMES, FACTOR_COVARIATE_NAMES[!(FACTOR_COVARIATE_NAMES %in% colnames(covariatesDataFrame))])
  NUMERIC_COVARIATE_NAMES = setdiff(COLNAMES, FACTOR_COVARIATE_NAMES)
  
  # Ensure the factors are in fact of type factor, and the quantitative variables are numeric:
  covariatesDataFrame = as.data.frame( lapply(colnames(covariatesDataFrame), function(column) {if (column %in% FACTOR_COVARIATE_NAMES) {fac = as.factor(covariatesDataFrame[, column]); if (column %in% names(RELEVELS)) {fac = relevel(fac, ref=RELEVELS[[column]])}; return(fac)} else {return(as.numeric(covariatesDataFrame[, column]))}}) )
  rownames(covariatesDataFrame) = ROWNAMES
  colnames(covariatesDataFrame) = COLNAMES
  
  contra = NULL
  MAX_NUM_CATS = Inf
  catData = covariatesDataFrame[, FACTOR_COVARIATE_NAMES, drop=FALSE]
  if (ncol(catData) > 0) {
    numCats = sapply(colnames(catData), function(col) nlevels(factor(catData[, col])))
    EXCLUDE_CATEGORICAL_COLS = names(numCats)[numCats <= 1 | numCats > MAX_NUM_CATS]
    if (!is.null(EXCLUDE_CATEGORICAL_COLS) && length(EXCLUDE_CATEGORICAL_COLS) > 0) {
      warning(paste("Excluding categorical variables with less than 2", ifelse(is.infinite(MAX_NUM_CATS), "", paste(" or more than ", MAX_NUM_CATS, sep="")), " categories: ", paste(paste("'", EXCLUDE_CATEGORICAL_COLS, "'", sep=""), collapse=", "), sep=""))
      FACTOR_COVARIATE_NAMES = setdiff(FACTOR_COVARIATE_NAMES, EXCLUDE_CATEGORICAL_COLS)
      covariatesDataFrame = covariatesDataFrame[, !(colnames(covariatesDataFrame) %in% EXCLUDE_CATEGORICAL_COLS), drop=FALSE]
    }
    
    # Inspired by http://stackoverflow.com/questions/4560459/all-levels-of-a-factor-in-a-model-matrix-in-r
    #
    # And, already ensured above that covariatesDataFrame[, FACTOR_COVARIATE_NAMES] satisfies:
    # 1) fac is of type factor.
    # 2) fac is releveled as designated in RELEVELS.
    if (Intercept)
      contra = lapply(FACTOR_COVARIATE_NAMES, function(column) {fac = covariatesDataFrame[, column]; fac = contrasts(fac);})
    else
      contra = lapply(FACTOR_COVARIATE_NAMES, function(column) {fac = covariatesDataFrame[, column]; fac = contrasts(fac,contrasts=F);})
    names(contra) = FACTOR_COVARIATE_NAMES
  }
  
  # Inspired by http://stackoverflow.com/questions/5616210/model-matrix-with-na-action-null :
  current.na.action = getOption('na.action')
  # Model matrix will now include "NA":
  options(na.action='na.pass')
  
  if(Intercept)
    design = model.matrix(~ ., data=covariatesDataFrame, contrasts.arg=contra)
  else
    design = model.matrix(~ 0 + ., data=covariatesDataFrame, contrasts.arg=contra)
  
  rownames(design) = rownames(covariatesDataFrame)
  
  options(na.action=current.na.action)
  
  return(list(design=design, covariates=COLNAMES, factorsLevels=sapply(contra, colnames, simplify=FALSE), numericCovars=NUMERIC_COVARIATE_NAMES, covariatesDataFrame=covariatesDataFrame))
}




# Function imported from https://github.com/pouryany/CovariateAnalysis
# Modified from http://stackoverflow.com/questions/13088770/how-to-write-linearly-dependent-column-in-a-matrix-in-terms-of-linearly-independ
# Function to find linearly dependednt columns of a matrix
linColumnFinder <- function(mat){
  
  mat[is.na(mat)] = 0
  
  # If the matrix is full rank then we're done
  if(qr(mat)$rank == ncol(mat)){
    return(list(indepCols = seq(1,ncol(mat),1), relations = "Matrix is of full rank"))
  }
  
  m <- ncol(mat)
  # cols keeps track of which columns are linearly independent
  cols <- 1
  All.message <- c()
  for(i in seq(2, m)){
    ids <- c(cols, i)
    mymat <- mat[, ids]
    if(qr(mymat)$rank != length(ids)){
      # Regression the column of interest on the previous columns to figure out the relationship
      o <- lm(mat[,i] ~ as.matrix(mat[,cols]) + 0)
      # Construct the output message
      start <- paste0(colnames(mat)[i], " = ")
      # Which coefs are nonzero
      nz <- !(abs(coef(o)) <= .Machine$double.eps^0.5)
      tmp <- colnames(mat)[cols[nz]]
      vals <- paste(coef(o)[nz], tmp, sep = "*", collapse = " + ")
      message <- paste0(start, vals)      
      All.message <- c(All.message,message)
    } else {
      # If the matrix subset was of full rank
      # then the newest column in linearly independent
      # so add it to the cols list
      cols <- ids
    }
  }
  return(list(indepCols = cols, relations = All.message))
}
