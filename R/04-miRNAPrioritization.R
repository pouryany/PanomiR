# library(parallel)
# library(metap)
# library(dplyr)

#' Outputs a table of miRNA (ordered) with respective p-values derived from method for prioritization
#'
#' @param enriches miRNA-pathway enrichment dataset, obtained from miRNAPathwayEnrichment
#' @param pathway.clusters pathway clusters, obtained from MappingPathwaysClusters
#' @param method a vector of methods (pCut, AggInv, AggLog, sumz, sumlog)
#' @param method.thresh a vector of methods threshold for each method in method, if NULL use default thresh values in method
#' @param mir.path.fdr.selection FDR threshold for mir-pathway adjusted p-values ; filter edges with adjusted p-values less than given threshold
#' @param top.clust top x clusters to perform miRNA prioritization on
#' @param samp.rate sampling rate for CLT
#' @param num.cores number of cores available
#' @param out.dir output directory 
#' @param data.dir data directory
#' @param save.csv if T, saves csv file for each cluster in top.clust in out.dir 
#' @param save.sampling if T, saves sampling data as RDS for each cluster in top.clust in data.dir 
#' @param save.jack.knife if T, saves jack-knifed sampling data as RDS for each cluster in top.clust in data.dir 
#' @param prefix prefix for all saved data 
#' @import dplyr
#' @import metap
#' @import dplyr
#' @import parallel
#' @return table of miRNA and p-values, each row contains a mirna and its associated p-values from the methods
#' @export

#### Alternative Prioritzation root
# Name to be updated
# Need to make sure it uses direct bootstrapping for small n instead of approximation

PrioritizeMicroRNA <- function(enriches0,
                                pathway.clusters,
                                method='AggInv',
                                method.thresh=NULL,
                                mir.path.fdr.thresh=0.25,
                                top.clust=2,
                                samp.rate=1000,
                                out.dir='',
                                data.dir='',
                                save.sampling=T,
                                run.jack.knife=T,
                                save.jack.knife=F,
                                num.cores=1,
                                save.csv=T,
                                prefix='')
{
  
  if (substring(out.dir, nchar(out.dir))!='/')
    out.dir <- paste0(out.dir, '/')
  if (!dir.exists(out.dir)){
    warning('Output directory does not exist.')
    dir.create(out.dir,recursive = TRUE)
  }
    
  
  
  if(out.dir=='/')
    out.dir=''
  
  if (substring(data.dir, nchar(data.dir))!='/')
    data.dir <- paste0(data.dir, '/')
  if (!dir.exists(data.dir))
    stop('Data directory does not exist.')
  
  if(data.dir=='/')
    data.dir=''
  
  output = list()
  
  # count miRNA-pathway enrichment with p-value less than threshold
  enriches   <- enriches0 %>% filter(., Intersect != 0)
  enriches   %<>% group_by(., y) %>% mutate(., path_fdr = p.adjust(pval, method = "fdr"))
  enriches   <- enriches %>% mutate(.,hit =ifelse(path_fdr < mir.path.fdr.thresh, 1, 0))
  
  # Need to fix this function to be able to work on specific clusters instead of all top clusters.
  for (clustNo in 1:top.clust){
    
    clustName <- paste0('Cluster', clustNo)
    
    print(paste0('Working on ', clustName, '.'))
    
    # select pathways in cluster
    pathways    <- as.character(pathway.clusters[pathway.clusters$cluster == clustNo, ]$Pathway)
    n_paths     <- length(pathways)
    
    # formulate number of miRNA-pathway enrichment with p-value less than threshold for each miRNA
    temp.enrich <- enriches[enriches$y %in% pathways, ]
    selector   <- temp.enrich %>% group_by(x)  %>%
      dplyr::summarise(., "cluster_hits" = sum(hit))
    
    # perform p-value aggregation based on methodlogy provided
    for (i in 1:length(method)){
      m <- method[i]
      
      print(paste0('Performing ', m, ' function.'))
      fn <- get(paste0(m, '.fn'))
      cover.fn <- get(paste0(m, '.cover.fn'))
      
      if (!is.null(method.thresh)){
        m.thresh <- method.thresh[i]
        temp <- fn(enriches=enriches0, pathways, is.selector=T, thresh=m.thresh)
      } else{
        temp <- fn(enriches=enriches0, pathways, is.selector = T)
      }
      
      m.selector <- temp$selector
      m.enriches0 <- temp$enriches0
      
      if (nrow(m.selector) < 3){
        print(paste0('Skipping ', m, ' function due to low number of miRNA after filter'))
        next
      }
      
      enrich.null <- m.enriches0 %>% filter(.,x %in% m.selector$x)
      
      sampling.data.dir <- paste0(data.dir, prefix, 'Sampling_Data/')
      
      if (save.sampling==TRUE){
        if(!dir.exists(sampling.data.dir))dir.create(sampling.data.dir,recursive = T)
      }
      
      

      sampling.data.filename  <- paste0(prefix, m, "_", samp.rate, "_samples.RDS")
      sampling.data.file <- paste0(sampling.data.dir, '/', sampling.data.filename)
     
      
      
       if(m %in% c("AggInv","AggLog")){
         # perform sampling
         sampling.data <- samplingDataBase(enrich.null,
                                              m.selector,
                                              samp.rate,
                                              fn,
                                              n_paths,
                                              sampling.data.file,
                                              jack.knife=FALSE,
                                              save.sampling = save.sampling,
                                              num.cores = 8)
         
         m.selector    <- methodProbBase(sampling.data = sampling.data[[paste0("SampSize_",100)]],
                                            selector      = m.selector,
                                            m             = m,
                                            n_paths       = n_paths,
                                            cover.fn      = cover.fn)
         
         
         
       } else{
         names(m.selector)[2:3] <- paste0(m, "_",names(m.selector)[2:3])
      }
     
      
      print(paste0(m, " Method Done"))
      
      # perform jack-knife
      if (run.jack.knife ==T){
        # jack.knife.data.filename  <- paste0(prefix, m, "_", samp.rate, "_samples_jack_knifed.RDS")
        # jack.knife.data.file <- paste0(sampling.data.dir, '/', jack.knife.data.filename)
        # 
        # jack.knife.data <- sampling.data.base(enrich.null, m.selector, samp.rate, fn, n_paths, 
        #                                       sampling.data.file=jack.knife.data.file, jack.knife=TRUE, save.sampling = save.jack.knife)
        # 
        
        sampling.data <- samplingDataBase(enrich.null,
                                             m.selector,
                                             samp.rate,
                                             fn,
                                             n_paths,
                                             sampling.data.file,
                                             jack.knife=FALSE,
                                             save.sampling = save.sampling,
                                             num.cores = 8)
        m.selector <- jackKnifeBase(selector    = m.selector,
                                       pathways    = pathways,
                                       enrich.null = enrich.null,
                                       fn = fn,
                                       jack.knife.data = 
                                         sampling.data[[paste0("SampSize_",100)]],
                                       m  = m, 
                                       num.cores = 8)
        
        print(paste0(m, " JackKnifing Method Done!"))
        
        m.selector <- m.selector[,c(1,3,2,4)]
      } else {
        m.selector <- m.selector[,c(1,3,2)]
      }
      
      selector <- merge(selector, m.selector, all=T)
      
      met      <- paste0(m,"_pval")
      selector <- selector %>% dplyr::arrange(.,!!sym(met))
      met2     <- paste0(m,"_fdr")
      selector <- selector %>% dplyr::mutate(.,
                                      !!met2 := p.adjust(p = !!sym(met)
                                                         ,method = "fdr")
                                      )
    }
    
    if (save.csv==T){
      save.name <- paste0(prefix, samp.rate, '_samples_clustNo_', clustNo, '.csv')
      print(paste0(save.name, ' saved!'))
      write.csv(selector, paste0(out.dir, save.name))
    }
    
    output[[clustName]] <- selector
  }
  
  return(output)
}




