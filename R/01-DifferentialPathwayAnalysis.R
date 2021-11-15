#' Differential Expression Analysis For Pathways
#' 
#' Performs differential expression analysis for pathways using LIMMA package 
#' with gene counts 
#' 
#' @param genes.counts Gene counts, rows refer to genes and columns to samples.
#' @param pathways Pathways table, containing pathway names and genes with id
#'   specified.
#' @param covariates Covariates/metadata file; rows matches the columns of
#'   gene.counts.
#' @param condition Condition to be examined (tumor vs normal etc); must exist 
#'   in covariates column.
#' @param adjust.covars Adjustment covariates like batch; if NULL,
#'   no adjustments performed.
#' @param covariates.correction If T, performs covariates detection and
#'   correction; requires **adjust.covars**; (limma).
#' @param quantile.norm If T, performs quantile normalization on pathway
#'   summary statistics; from *preprocess* package.
#' @param out.dir Output directory.
#' @param save.RDS.name If not NULL, saves output as RDS using save name,
#'   if NULL, does not save output.
#' @param id ID matching genes to pathways; rownames of gene.counts.
#' @param de.genes If not NULL, add t-scores to pathways summary statistics;
#'   filter by genes t-scores.
#' @param min.path.size Minimum pathway size.
#' @param method Define method to use for pathway summary statistics;
#'   specifications in documentations.
#' @param trim Filter pathways with mean less than trim threshold
#'   in pathway summary statistics.
#' @param genes.counts.log If T, log(genes.counts).
#' @param contrast.conds Provide a contrast expression to be used in Limma
#'   comparison. This is necessary if you have more than two levels in the
#'   condition covariate.
#' @return List containing differentially expressed pathways as DEP and pathway
#'   summary statistics as pathwaySummaryStats.
#' @export
DifferentialPathwayAnalysis <- function(genes.counts, 
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
                                        contrast.conds = NA)
{
  
  if (substring(out.dir, nchar(out.dir))!='/')
    out.dir <- paste0(out.dir, '/')
  if (!dir.exists(out.dir))
    stop('Output directory does not exist.')

  # select pathways with genes in the gene count data and a minimum pathway set size  
  pathways       <- as.data.frame(pathways)
  genes.pathways <- pathways[pathways[,id] %in% rownames(genes.counts),]
  
  genes.pathways %<>% group_by(.,Pathway) %>%
    dplyr::summarise(., n = n()) %>% 
    filter(., n >=min.path.size)
  
  pathways <- pathways %>% 
    filter(., Pathway %in% genes.pathways$Pathway)
  
  # use de genes as a filter
  if (!is.null(de.genes)){
    t.scores <- de.genes %>% dplyr::mutate(., !!id:=rownames(de.genes)) %>% dplyr::select(., c(ENSEMBL, t))
  } 
  
  # log gene counts if gene counts are not log-transformed yet; essential for path summary statistics
  if (genes.counts.log == T)
    genes.counts <- log(genes.counts)

  # generate pathway summary statistics
  pathwaySummaryStats <- Path_Summary(genes.counts,
                                   pathways,
                                   id = id,
                                   method = method,
                                   de.genes = de.genes,
                                   z.normalize = F,
                                   trim = trim,
                                   t.scores = t.scores)
  
  # filter pathways with na values in the pathway summary statistics and z-normalize the pathway summary statistics 
  pathwaySummaryStats <- pathwaySummaryStats[rowSums(is.na(pathwaySummaryStats))==0,]
  pathwaySummaryStats <- apply(pathwaySummaryStats, 2, function(X){(X - mean(X))/stats::sd(X)})
  
  # perform quantile normalization if needed
  # Add importing
  if (quantile.norm == T){
    pathways.names <- rownames(pathwaySummaryStats)
    pathwaySummaryStats <- preprocessCore::normalize.quantiles(pathwaySummaryStats)
    rownames(pathwaySummaryStats) <- pathways.names
    colnames(pathwaySummaryStats) <- rownames(covariates)
  }
  
  # set factors in covariates if needed
  ### Fix this part later
  # for (covs in colnames(covariates)){
  #   if(!is.numeric(covs))
  #   covariates[,covs] <- as.factor(covariates[,covs])
  # }

  # perform covariates correction and create design matrix for limma DE analysis;
  # if no adjust.covars are available, then design matrix only consider the condition in question.
  FIT.res <- NULL
  if (covariates.correction==T){
    stop('Under development. Please use Covariates.correction = F option')
    # if (is.null(adjust.covars))
    #   stop('Covariates correction requires adjusted covariates.')
    # output <- CovariatesCorrelationCorrection(pathwaySummaryStats,
    #                                           covariates,
    #                                           adjust.covars,
    #                                           min.iter=20,
    #                                           voom=F,
    #                                           out.dir='',
    #                                           save.RDS.name=NULL,
    #                                           plot=F)
    # 
    # adjust.covars <- unique(c(adjust.covars, output$primeVariables))
    # adjust.covars <- adjust.covars[adjust.covars!=condition]
    # design.mat    <- getDesignMatrix(covariates[,c(condition,adjust.covars),drop = F], Intercept = F)
    # design.mat$design <- design.mat$design[,linColumnFinder(design.mat$design)$indepCols]
  } else {
    if (is.null(adjust.covars)){
      conditions <- as.data.frame(covariates[,condition])
      colnames(conditions) <- condition
      rownames(conditions) <- rownames(covariates)
      design.mat <- getDesignMatrix(conditions, Intercept=F)
    }
    else {
      design.mat <- getDesignMatrix(covariates[,c(condition,adjust.covars),drop = F], Intercept = F)
      design.mat$design <- design.mat$design[,linColumnFinder(design.mat$design)$indepCols]
      
      # Getting pathway residuals
      res.des.mat        <- getDesignMatrix(covariates[,c(adjust.covars),
                                                       drop = F],
                                            Intercept = F)
      res.des.mat$design <- res.des.mat$design[,
                                               linColumnFinder(res.des.mat$design)$indepCols]
      FIT.res  <- lmFit(pathwaySummaryStats, res.des.mat$design)
      FIT.res  <- residuals.MArrayLM(FIT.res,pathwaySummaryStats)
    }
  }
  
  conditions.type <- as.character(unique(covariates[,condition]))
  if(is.na(contrast.conds)){
    if (length(conditions.type) > 2)
      stop('Please compare only 2 conditions at once.')
    cond1 <- paste0(condition, conditions.type[1])
    cond2 <- paste0(condition, conditions.type[2])
    contrasts.name <- paste0(cond1, '-', cond2) 
  }else{
    contrasts.name <- contrast.conds
  }
 
  # limma DE analysis
  FIT       <- lmFit(pathwaySummaryStats, design.mat$design)
  
  colnames(FIT$coefficients) <- gsub('-','_', colnames(FIT$coefficients))
  
  contrast  <- makeContrasts(contrasts=c(contrasts.name),
                             levels = colnames(FIT$coefficients))
  FIT.CONTR <- contrasts.fit(FIT, contrasts=contrast)
  FIT.CONTR <- eBayes(FIT.CONTR)
  tT        <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)
  #}
  tT$contrast  <- contrasts.name
  
  output       <- list('DEP'=tT, 
                       'pathwaySummaryStats'=pathwaySummaryStats,
                       "contrast" = contrasts.name,
                       "PathwayResiduals" = FIT.res)
  
  
  if (!is.null(save.RDS.name))
    saveRDS(output, paste0(out.dir, save.RDS.name))
  return(output)
}
