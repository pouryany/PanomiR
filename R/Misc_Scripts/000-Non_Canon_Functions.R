
#' Performs differential expression analysis for pathways using LIMMA package with gene counts and random pertubation of labels
#' 
#' @param genes.counts gene counts, rows refer to genes and columns to samples
#' @param pathways pathways table, containing pathway names and genes with id specified
#' @param covariates covariates/metadata file; rows matches the columns of gene.counts 
#' @param condition to be examined (tumor vs normal etc); must exist in covariates column
#' @param adjust.covars adjustment covariates like batch; if NULL, no adjustments performed
#' @param covariates.correction if T, performs covariates correction; requires **adjust.covars**; (limma)
#' @param quantile.norm if T, performs quantile normalization on pathway summary statistics; from *preprocess* package
#' @param out.dir output directory
#' @param save.RDS.name if not NULL, saves output as RDS using save name, if NULL, does not save output 
#' @param id id matching genes to pathways; rownames of gene.counts 
#' @param de.genes if not NULL, add t-scores to pathways summary statistics; filter by genes t-scores
#' @param min.path.size minimum pathway size 
#' @param method define method to use for pathway summary statistics; specifications in documentations 
#' @param trim filter pathways with mean less than trim threshold in pathway summary statistics  
#' @param genes.counts.log if T, log(genes.counts)   
#' @param samp.rate sampling rate  
#' @param pval.test probability of adjusted p-values of original DE pathways against the sampling data
#' @return list contain differentially expressed pathways as DEP, pathway summary statistics as pathwaySummaryStats and iterations of labels pertubation

performDPAwithLabelPertubation <- function(genes.counts, 
                                           pathways, 
                                           covariates, 
                                           condition,
                                           adjust.covars=NULL,
                                           covariates.correction=F,
                                           quantile.norm=F,
                                           out.dir='',
                                           save.RDS.name = NULL,
                                           id='ENSEMBL', 
                                           de.genes=NULL, 
                                           min.path.size=10, 
                                           method='x2', 
                                           trim=0.025,
                                           genes.counts.log = T,
                                           samp.rate=1000,
                                           pval.test=T)
{
    output.list <- list()
    
    # generate original output
    output <- DifferentialPathwayAnalysis(genes.counts=genes.counts, 
                                          pathways=pathways, 
                                          covariates=covariates, 
                                          condition=condition,
                                          adjust.covars=adjust.covars,
                                          covariates.correction=covariates.correction,
                                          quantile.norm=quantile.norm,
                                          out.dir=out.dir,
                                          save.RDS.name = NULL,
                                          id=id, 
                                          de.genes=de.genes, 
                                          min.path.size=min.path.size, 
                                          method=method, 
                                          trim=trim, 
                                          genes.counts.log = genes.counts.log)
    pathwaySummaryStats <- output$pathwaySummaryStats
    
    pvals <- c()
    
    for (i in 1:samp.rate){
        set.seed(i)
        listname <- paste0('ITERATION_',i)
        
        # perturb sample labels
        temp.covariates <- covariates
        set.seed(i)
        temp.covariates[,condition] <- sample(temp.covariates[,condition])
        
        conditions.type <- as.character(unique(covariates[,condition]))
        if (length(conditions.type) > 2)
            stop('Please compare only 2 conditions at once.')
        cond1 <- paste0(condition, conditions.type[1])
        cond2 <- paste0(condition, conditions.type[2])
        
        
        
        design.mat <- getDesignMatrix(temp.covariates[,c(condition,adjust.covars),drop = F], Intercept = F)
        design.mat$design <- design.mat$design[,linColumnFinder(design.mat$design)$indepCols]
        
        
        
        # limma DE analysis
        FIT <- lmFit(pathwaySummaryStats, design.mat$design)
        colnames(FIT$coefficients) <- gsub('-','_', colnames(FIT$coefficients))
        contrasts.name <- paste0(cond1, '-', cond2)
        contrast  <- makeContrasts(contrasts=c(contrasts.name),
                                   levels = colnames(FIT$coefficients))
        FIT.CONTR <- contrasts.fit(FIT, contrasts=contrast)
        FIT.CONTR <- eBayes(FIT.CONTR)
        tt <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)
        
        output.list[[listname]] <- tt
        
        # store the adjusted p-values based on the original DE Pathway analysis
        if (pval.test == T){
            tt <- tt[rownames(output$DEP),]
            pvals <- cbind(pvals, tt$adj.P.Val)
        }
    }
    
    if (pval.test == T){
        # derive the probaility of the adjusted p-values in the original DE Pathway analysis
        pvals.mean <- rowMeans(pvals)
        pvals.sd <- apply(pvals, 1, sd)
        pvals.prob <- pnorm(tr$adj.P.Val, mean=pvals.mean, sd=pvals.sd)
        output$DEP$perturb.labels.pvals.prob <- pvals.prob
    }
    
    output[['iterations']] <- output.list
    if (!is.null(save.RDS.name))
        saveRDS(output, paste0(out.dir, save.RDS.name))
    return(output)
}

#' Performs differential expression analysis for pathways using LIMMA package with gene counts and random pertubation of genes
#' 
#' @param genes.counts gene counts, rows refer to genes and columns to samples
#' @param pathways pathways table, containing pathway names and genes with id specified
#' @param covariates covariates/metadata file; rows matches the columns of gene.counts 
#' @param condition to be examined (tumor vs normal etc); must exist in covariates column
#' @param adjust.covars adjustment covariates like batch; if NULL, no adjustments performed
#' @param covariates.correction if T, performs covariates correction; requires **adjust.covars**; (limma)
#' @param quantile.norm if T, performs quantile normalization on pathway summary statistics; from *preprocess* package
#' @param out.dir output directory
#' @param save.RDS.name if not NULL, saves output as RDS using save name, if NULL, does not save output 
#' @param id id matching genes to pathways; rownames of gene.counts 
#' @param de.genes if not NULL, add t-scores to pathways summary statistics; filter by genes t-scores
#' @param min.path.size minimum pathway size 
#' @param method define method to use for pathway summary statistics; specifications in documentations 
#' @param trim filter pathways with mean less than trim threshold in pathway summary statistics  
#' @param genes.counts.log if T, log(genes.counts)   
#' @param samp.rate sampling rate  
#' @param pval.test probability of adjusted p-values of original DE pathways against the sampling data
#' @return list contain differentially expressed pathways as DEP, pathway summary statistics as pathwaySummaryStats and iterations of genes pertubation

performDPAwithGenePertubation <- function(genes.counts, 
                                          pathways, 
                                          covariates, 
                                          condition,
                                          adjust.covars=NULL,
                                          covariates.correction=F,
                                          quantile.norm=F,
                                          out.dir='',
                                          save.RDS.name = NULL,
                                          id='ENSEMBL', 
                                          de.genes=NULL, 
                                          min.path.size=10, 
                                          method='x2', 
                                          trim=0.025,
                                          genes.counts.log = T,
                                          samp.rate=1000,
                                          pval.test=T)
{
    output.list <- list()
    
    # generate original output
    output <- DifferentialPathwayAnalysis(genes.counts=genes.counts, 
                                          pathways=pathways, 
                                          covariates=covariates, 
                                          condition=condition,
                                          adjust.covars=adjust.covars,
                                          covariates.correction=covariates.correction,
                                          quantile.norm=quantile.norm,
                                          out.dir=out.dir,
                                          save.RDS.name = NULL,
                                          id=id, 
                                          de.genes=de.genes, 
                                          min.path.size=min.path.size, 
                                          method=method, 
                                          trim=trim, 
                                          genes.counts.log = genes.counts.log)
    pvals <- c()
    
    for (i in 1:samp.rate){
        set.seed(i)
        listname <- paste0('ITERATION_',i)
        
        # perturb genes
        temp.genes.counts <- gene.counts
        rownames(temp.genes.counts) <- sample(rownames(temp.genes.counts))
        
        # run DPA with perturbed genes
        tt <- DifferentialPathwayAnalysis(genes.counts=temp.genes.counts, 
                                          pathways=pathways, 
                                          covariates=covariates, 
                                          condition=condition,
                                          adjust.covars=adjust.covars,
                                          covariates.correction=covariates.correction,
                                          quantile.norm=quantile.norm,
                                          out.dir=out.dir,
                                          save.RDS.name = NULL,
                                          id=id, 
                                          de.genes=de.genes, 
                                          min.path.size=min.path.size, 
                                          method=method, 
                                          trim=trim, 
                                          genes.counts.log = genes.counts.log)
        
        tt <- tt$DEP
        output.list[[listname]] <- tt
        
        # store the adjusted p-values based on the original DE Pathway analysis
        if (pval.test == T){
            tt <- tt[rownames(output$DEP),]
            pvals <- cbind(pvals, tt$adj.P.Val)
        }
    }
    
    if (pval.test == T){
        # derive the probaility of the adjusted p-values in the original DE Pathway analysis
        pvals.mean <- rowMeans(pvals)
        pvals.sd <- apply(pvals, 1, sd)
        pvals.prob <- pnorm(tT$adj.P.Val, mean=pvals.mean, sd=pvals.sd)
        output$DEP$perturb.labels.pvals.prob <- pvals.prob
    }
    
    output[['iterations']] <- output.list
    if (!is.null(save.RDS.name))
        saveRDS(output, paste0(out.dir, save.RDS.name))
    return(output)
}






performDPAwithGenePertubation2 <- function(genes.counts, 
                                           pathways, 
                                           covariates, 
                                           condition,
                                           adjust.covars=NULL,
                                           covariates.correction=F,
                                           quantile.norm=F,
                                           out.dir='',
                                           save.RDS.name = NULL,
                                           id='ENSEMBL', 
                                           de.genes=NULL, 
                                           min.path.size=10, 
                                           method='x2', 
                                           trim=0.025,
                                           genes.counts.log = T,
                                           samp.rate=1000,
                                           pval.test=T,
                                           numcores = 1)
{
    output.list <- list()
    
    # generate original output
    output <- DifferentialPathwayAnalysis(genes.counts=genes.counts, 
                                          pathways=pathways, 
                                          covariates=covariates, 
                                          condition=condition,
                                          adjust.covars=adjust.covars,
                                          covariates.correction=covariates.correction,
                                          quantile.norm=quantile.norm,
                                          out.dir=out.dir,
                                          save.RDS.name = NULL,
                                          id=id, 
                                          de.genes=de.genes, 
                                          min.path.size=min.path.size, 
                                          method=method, 
                                          trim=trim, 
                                          genes.counts.log = genes.counts.log)
    pvals <- c()
    
    
    library(doParallel)
    library(foreach)
    library(parallel)
    
    registerDoParallel(numcores)
    
    
    
    par.paths <- foreach(i = 1:samp.rate) %dopar% {
        
        set.seed(i)
        listname <- paste0('ITERATION_',i)
        
        # perturb genes
        temp.genes.counts <- genes.counts
        rownames(temp.genes.counts) <- sample(rownames(temp.genes.counts))
        
        # run DPA with perturbed genes
        tt <- DifferentialPathwayAnalysis(genes.counts=temp.genes.counts, 
                                          pathways=pathways, 
                                          covariates=covariates, 
                                          condition=condition,
                                          adjust.covars=adjust.covars,
                                          covariates.correction=covariates.correction,
                                          quantile.norm=quantile.norm,
                                          out.dir=out.dir,
                                          save.RDS.name = NULL,
                                          id=id, 
                                          de.genes=de.genes, 
                                          min.path.size=min.path.size, 
                                          method=method, 
                                          trim=trim, 
                                          genes.counts.log = genes.counts.log)
        
        tt <- tt$DEP
        return(tt)
        
    }
    
    names(par.paths) <- 1:samp.rate
    
    #output[['iterations']] <- output.list
    
    return(par.paths)
}


# A variation of density plot where you could edit axis label sizes
# Not an important function I believe
pDens2 <- function (object, group = NULL, col = NULL, main = NULL, 
                    legend = "topleft", cex.lab = 1,
                    ...) 
{
    E <- as.matrix(object)
    narray <- ncol(E)
    if (is.null(group)) 
        group <- colnames(E)
    if (is.null(group)) 
        group <- 1:narray
    group <- as.factor(group)
    ngroup <- nlevels(group)
    if (is.null(col)) 
        col <- 1:ngroup
    col <- rep(col, length = ngroup)
    if (is.logical(legend)) {
        legend.position <- "topleft"
    }
    else {
        legend.position <- as.character(legend)
        legend <- TRUE
    }
    legend.position <- match.arg(legend.position, c("bottomright", 
                                                    "bottom", "bottomleft", "left", "topleft", "top", "topright", 
                                                    "right", "center"))
    arraycol <- group
    levels(arraycol) <- col
    arraycol <- as.vector(arraycol)
    npoint <- 512
    X <- Y <- matrix(0, npoint, narray)
    for (a in 1:ncol(E)) {
        d <- density(E[, a], n = npoint, na.rm = TRUE, ...)
        X[, a] <- d$x
        Y[, a] <- d$y
    }
    matplot(X, Y, xlab = "Intensity", ylab = "Density", main = main, 
            type = "l", col = arraycol, lwd = 2, lty = 1, cex.lab = cex.lab)
    if (legend && ngroup > 1) 
        legend(legend.position, lwd = 2, legend = levels(group), 
               col = col)
    invisible(list(X = X, Y = Y))
}

