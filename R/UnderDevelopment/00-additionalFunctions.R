





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



# Function imported from https://github.com/pouryany/CovariateAnalysis
# Modified from http://stackoverflow.com/questions/13088770/how-to-write-linearly-dependent-column-in-a-matrix-in-terms-of-linearly-independ
# Function to find linearly dependednt columns of a matrix



#' #'  The function calculate targeting score of miRNA w.r.t to a cluster 
#' #' of pathways via lancaster aggregation method.
#' #' @param enriches a table of miRNA pathway enrichments. Universe
#' #' @param pathways queried pathways. e.g. cluster pathways
#' #' @param is.selector internal argument
#' #' @param thresh internal argument
#' #' @return a  scoring of miRNAs in a cluster of pathways
#' #' @import metap
#' #' @import dplyr
#' lancaster.fn <- function(enriches, pathways, is.selector, thresh=NULL){
#'   temp.enrich <- enriches[enriches$y %in% pathways, ]
#'   selector <- temp.enrich %>%
#'     group_by(x) %>%
#'     dplyr::summarise(n = n(), k = lancaster(pval,weight)/n())%>%
#'     arrange(., x)
#'   
#'   if (is.selector==T){
#'     return(list('selector'=selector, 'enriches0'=enriches))
#'   } else {
#'     return(selector)
#'   }
#' }


