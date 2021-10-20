library(data.table)
library(limma)
library(dplyr)
library(stringr)
library(parallel)

#' Outputs a table of miRNA-pathway correlation
#'
#' @param pathway.summaries pathway summary statistics
#' @param pathway.clusters pathway clusters, obtained from MappingPathwaysClusters
#' @param mirna.counts miRNA count data
#' @param covariates covariates/metadata file; rows matches the columns of gene.counts 
#' @param condition to be examined (tumor vs normal etc); must exist in covariates column
#' @param prime.condition condition of interest; usually tumor/disease
#' @param compare.condition condition of interest; usually normal
#' @param adjust.covars adjustment covariates like batch; if NULL, no adjustments performed
#' @param covariates.correction if T, performs covariates correction; requires **adjust.covars**; (limma)
#' @param mirna.log if T, performs log-transformation of miRNA count data
#' @param cor.method correlation method
#' @param sig.cor.thresh significant correlation threshold for prime condition
#' @param num.cores number of cores available
#' @param top.clust top n clusters to perform function
#' @param save.RDS.name if T, saves RDS file with chosen name
#' @return table of miRNA and their associated adj pvalue for correlation

miRNAPathwayCorrelation <- function(pathway.summaries,
                                    pathway.clusters,
                                    mirna.counts,
                                    covariates,
                                    condition,
                                    prime.cond,
                                    compare.cond,
                                    adjust.covars = NULL,
                                    covariates.correction = F,
                                    mirna.log = F,
                                    cor.method = 'pearson',
                                    sig.cor.thresh = 0.05,
                                    num.cores = 1,
                                    top.clust = 2,
                                    save.RDS.name = NULL)
{
  output <- list()
  for (clustNo in 1:top.clust){
    clustName <- paste0('Cluster',clustNo)

    # select pathways within cluster
    sel.paths <- as.character(pathway.clusters[pathway.clusters$cluster == clustNo,]$Pathway)
    temp.pathway.summaries <- pathway.summaries[sel.paths, rownames(covariates)]
    mirna.counts <- mirna.counts[,rownames(covariates)]
    
    # log transform miRNA counts
    if (mirna.log == T)
      mirna.counts <- log(mirna.counts+3)
    
    # adjust for covariates in miRNA count data and pathway summary statistics
    if (!is.null(adjust.covars)){
      design.mat <- getDesignMatrix(covariates[,(adjust.covars),drop = F], Intercept = F)
      design.mat$design <- design.mat$design[,linColumnFinder(design.mat$design)$indepCols]
      FIT <- lmFit(mirna.counts, design.mat$design)
      res.mirna.counts <- residuals.MArrayLM(FIT, mirna.counts)
      FIT <- lmFit(temp.pathway.summaries, design.mat$design)
      res.pathway.summaries <- residuals.MArrayLM(FIT, temp.pathway.summaries)
    } else {
      res.mirna.counts <- mirna.counts
      res.pathway.summaries <- temp.pathway.summaries
    }
    
    # obtain samples for chosen prime and compare conditions
    prime.cond.samples <- rownames(covariates)[covariates[,condition] == prime.cond]
    compare.cond.samples <- rownames(covariates)[covariates[,condition] == compare.cond]
    
    # derive correlation for miRNA count data and pathway summary statistics
    prime.cond.cor <- cor.table(res.pathway.summaries[,prime.cond.samples], res.mirna.counts[,prime.cond.samples], method=cor.method)
    compare.cond.cor <- cor.table(res.pathway.summaries[,compare.cond.samples], res.mirna.counts[,compare.cond.samples], method=cor.method)
    
    # obtain miRNA-pathway with significant correlation in prime condition
    temp <- prime.cond.cor[prime.cond.cor$adj < sig.cor.thresh,]
    combined.cor <- inner_join(temp, compare.cond.cor, by=c('Pathway', 'miRNA'))
    combined.cor$adj.y <- p.adjust(combined.cor$pval.y, method='fdr')
    
    output[[clustName]] <- combined.cor 
  }
  if (!is.null(save.RDS.name)){
    saveRDS(output, save.RDS.name)
  }
  return(output)
}

