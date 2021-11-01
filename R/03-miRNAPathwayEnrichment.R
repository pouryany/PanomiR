<<<<<<< HEAD
# library(dplyr)
# library(parallel)
# library(clusterProfiler)
# library(org.Hs.eg.db)

#' Outputs enrichment probability of miRNAs based on pathway clusters
#'
#' @param mir.sets a table of miRNAs and a list of their interactions with genes in ENTREZ ID 
#' @param pathways.sets a table of pathways and a list of their interactions with genes in ENTREZ ID
#' @param gene.selection a table of genes with dtype; if not NULL, select only genes from given table
#' @param mir.selection a table of miRNA names; if not NULL, select only miRNAs from given table
#' @param from.id id of genes in genes.selection
#' @param to.id id of genes used in pcxn and pathways set
#' @param min.path.size filter out pathways with sets less than given value
#' @param num.cores number of cores available
#' @param out.dir output directory
#' @param save.RDS.name if not NULL, saves output as RDS using save name 
#' @import dplyr
#' @import parallel
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @return table of enrichment, each row contains mirna-pathway and its enrichment p-values
#' @export

miRNAPathwayEnrichment <- function(mir.sets, 
                                   pathways.sets, 
                                   genes.selection=NULL, 
                                   mir.selection=NULL, 
                                   from.id='ENSEMBL', 
                                   to.id='ENTREZID', 
                                   min.path.size=9, 
                                   num.cores=1, 
                                   out.dir='', 
                                   save.RDS.name=NULL)
  {
  
  
  if (substring(out.dir, nchar(out.dir))!='/')
    out.dir <- paste0(out.dir, '/')
  if (!dir.exists(out.dir))
    stop('Output directory does not exist.')
  
  # select pathways with minimum set size
  paths.sel  <- sapply(pathways.sets,length)
  pathways.sets   <- pathways.sets[paths.sel > min.path.size]
  paths.ref <- Reduce(union,pathways.sets)

  # select miRNAs with targets in pathways
  mir.sets <- lapply(mir.sets, function(X){X[X %in% paths.ref]})
  
  # select pathways with selected genes of interest
  if (!is.null(genes.selection)){
    gene.df     <- bitr(genes.selection, fromType = from.id,
                        toType = to.id,
                        OrgDb = org.Hs.eg.db)
    pathways.sets <- lapply(pathways.sets, function(X){
      X[X %in% gene.df[,c(to.id)]]})
    
    mir.sets <- lapply(mir.sets, function(X){
      X[X %in% gene.df[,c(to.id)]]})
    paths.ref <- Reduce(union,pathways.sets)
  }
  
  # select miRNAs of interest
  if (!is.null(mir.selection)){
    
    mir.sets <- mir.sets[names(mir.sets) %in% mir.selection]
    
  }
  
  sel.vec <- sapply(mir.sets,length)
  mir.sets <- mir.sets[sel.vec > min.path.size]
  iterator <- (merge(names(mir.sets),names(pathways.sets)))
  iterator <- iterator %>% mutate_all(.,as.character)
  all <- length(paths.ref)
  
  # find enrichment p-value of each miRNA target set and each pathway set
  enrichs   <- mclapply(1:nrow(iterator), 
                        function(Y){
                          X <- iterator[Y,]
                          q <- length(intersect(unlist(pathways.sets[X[[2]]]), 
                                                unlist(mir.sets[X[[1]]])))  
                          m <- length(unlist(mir.sets[X[[1]]]))
                          n <- all - m
                          k <- length(unlist(pathways.sets[X[[2]]]))
                          pval <- phyper(q-1,m,n,k,lower.tail = F,log.p = F)
                          return(c("pval" = pval,
                                   "Intersect"    = q, 
                                   "mirset_Size"  = m,
                                   "not_mirset"   = n, 
                                   "pathway_Size" = k,
                                   "ratio_in"     = q/(k-q),
                                   "ratio_out"    = (m-q)/n,
                                   "ratio_ratios" = (q/(k-q))/((m-q)/n) ))
                        }, mc.cores = num.cores)
  
  temp      <- do.call(rbind,enrichs)
  temp      <- as.data.frame(temp)
  iterator  <- cbind(iterator,temp)
  
  if (!is.null(save.RDS.name)){
    saveRDS(iterator, paste0(out.dir, save.RDS.name))
  }
  return (iterator)
}

=======
library(dplyr)
library(parallel)
library(clusterProfiler)
library(org.Hs.eg.db)

#' Outputs enrichment probability of miRNAs based on pathway clusters
#'
#' @param mir.sets a table of miRNAs and a list of their interactions with genes in ENTREZ ID 
#' @param pathways.sets a table of pathways and a list of their interactions with genes in ENTREZ ID
#' @param gene.selection a table of genes with dtype; if not NULL, select only genes from given table
#' @param mir.selection a table of miRNA names; if not NULL, select only miRNAs from given table
#' @param from.id id of genes in genes.selection
#' @param to.id id of genes used in pcxn and pathways set
#' @param min.path.size filter out pathways with sets less than given value
#' @param num.cores number of cores available
#' @param out.dir output directory
#' @param save.RDS.name if not NULL, saves output as RDS using save name 
#' @return table of enrichment, each row contains mirna-pathway and its enrichment p-values

miRNAPathwayEnrichment <- function(mir.sets, 
                                   pathways.sets, 
                                   genes.selection=NULL, 
                                   mir.selection=NULL, 
                                   from.id='ENSEMBL', 
                                   to.id='ENTREZID', 
                                   min.path.size=9, 
                                   num.cores=1, 
                                   out.dir='', 
                                   save.RDS.name=NULL)
  {
  
  
  if (substring(out.dir, nchar(out.dir))!='/')
    out.dir <- paste0(out.dir, '/')
  if (!dir.exists(out.dir))
    stop('Output directory does not exist.')
  
  # select pathways with minimum set size
  paths.sel  <- sapply(pathways.sets,length)
  pathways.sets   <- pathways.sets[paths.sel > min.path.size]
  paths.ref <- Reduce(union,pathways.sets)

  # select miRNAs with targets in pathways
  mir.sets <- lapply(mir.sets, function(X){X[X %in% paths.ref]})
  
  # select pathways with selected genes of interest
  if (!is.null(genes.selection)){
    gene.df     <- bitr(genes.selection, fromType = from.id,
                        toType = to.id,
                        OrgDb = org.Hs.eg.db)
    pathways.sets <- lapply(pathways.sets, function(X){
      X[X %in% gene.df[,c(to.id)]]})
    
    mir.sets <- lapply(mir.sets, function(X){
      X[X %in% gene.df[,c(to.id)]]})
    paths.ref <- Reduce(union,pathways.sets)
  }
  
  # select miRNAs of interest
  if (!is.null(mir.selection)){
    
    mir.sets <- mir.sets[names(mir.sets) %in% mir.selection]
    
  }
  
  sel.vec <- sapply(mir.sets,length)
  mir.sets <- mir.sets[sel.vec > min.path.size]
  iterator <- (merge(names(mir.sets),names(pathways.sets)))
  iterator <- iterator %>% mutate_all(.,as.character)
  all <- length(paths.ref)
  
  # find enrichment p-value of each miRNA target set and each pathway set
  enrichs   <- mclapply(1:nrow(iterator), 
                        function(Y){
                          X <- iterator[Y,]
                          q <- length(intersect(unlist(pathways.sets[X[[2]]]), 
                                                unlist(mir.sets[X[[1]]])))  
                          m <- length(unlist(mir.sets[X[[1]]]))
                          n <- all - m
                          k <- length(unlist(pathways.sets[X[[2]]]))
                          pval <- phyper(q-1,m,n,k,lower.tail = F,log.p = F)
                          return(c("pval" = pval,
                                   "Intersect"    = q, 
                                   "mirset_Size"  = m,
                                   "not_mirset"   = n, 
                                   "pathway_Size" = k,
                                   "ratio_in"     = q/(k-q),
                                   "ratio_out"    = (m-q)/n,
                                   "ratio_ratios" = (q/(k-q))/((m-q)/n) ))
                        }, mc.cores = num.cores)
  
  temp      <- do.call(rbind,enrichs)
  temp      <- as.data.frame(temp)
  iterator  <- cbind(iterator,temp)
  
  if (!is.null(save.RDS.name)){
    saveRDS(iterator, paste0(out.dir, save.RDS.name))
  }
  return (iterator)
}

>>>>>>> parent of 1392f08 (functions metadata fixed/@import calls)
