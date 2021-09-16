library(limma)
library(edgeR)
library(org.Hs.eg.db)
library(gplots)
library(RColorBrewer)
library(pROC)
library(aggregation)
source('DifferentialPathwayAnalysis.R')
source('MappingPathwaysClusters.R')
source('miRNAPrioritization.R')

generate_simulated_data <- function(genes.counts,
                                    covariates,
                                    condition,
                                    out.dir = '',
                                    log = T,
                                    plot = T,
                                    seed = 42,
                                    plot.genes = 500,
                                    prefix = '',
                                    display.genes = NULL)
{
  if (substring(out.dir, nchar(out.dir))!='/')
    out.dir <- paste0(out.dir, '/')
  if (!dir.exists(out.dir))
    stop('Output directory does not exist.')
  
  fig.dir <- paste0(out.dir, 'Figures/')
  
  if (!dir.exists(fig.dir))
    dir.create(fig.dir, recursive = T)
  
  fig.dir <- paste0(fig.dir, prefix)
  
  set.seed(seed)
  if (log == T)
    genes.counts <- log(genes.counts+3)
  geneMeans <- rowMeans(genes.counts)
  geneSds <- apply(genes.counts, 1, sd)
  simulated.genes.counts <- (genes.counts - geneMeans)/geneSds
  
  simulated.genes.counts <- simulated.genes.counts[,sample(colnames(simulated.genes.counts))]
  condition_data <- covariates[,condition]
  covariates <- covariates[colnames(simulated.genes.counts),]
  covariates[,condition] <- condition_data
  
  if (plot == T){
    if (is.null(display.genes))
      display.genes <- sample(1:nrow(genes.counts), plot.genes)
    fig.name <- paste0(fig.dir,'original_genes_counts.jpg')
    jpeg(fig.name)
    display.ori <- genes.counts[display.genes,]
    display.ori <- ((display.ori-min(display.ori))*2/(max(display.ori)-min(display.ori))-1)*6
    heatmap.2(display.ori, col=bluered(100), dendrogram='none', 
              Colv=NA, Rowv=NA, trace="none", colRow=NA, labCo=NA, labRow = NA)
    dev.off()
    
    fig.name <- paste0(fig.dir,'simulated_genes_counts.jpg')
    jpeg(fig.name)
    display.simulated <- simulated.genes.counts[display.genes,]
    heatmap.2(display.simulated, col=bluered(100), dendrogram='none', 
              Colv=NA, Rowv=NA, trace="none", colRow=NA,labCol=NA, labRow = NA)
    dev.off()
  }
  
  return(list('simulated_genes_counts'=simulated.genes.counts,'simulated_covariates'=covariates, 'display_genes'=display.genes))
}

induce.miRNA.DE <- function(simulated.genes.counts,
                            mir.sets,
                            cond1.samples,
                            cond2.samples,
                            out.dir = '',
                            mirs.interest = NULL,
                            mir.selection = NULL,
                            seed = 42,
                            mir.set.size = 50,
                            plot = T,
                            mapping = F,
                            from.id = 'ENTREZID',
                            to.id = 'ENSEMBL',
                            min.target.size = 0,
                            max.target.size = NULL,
                            fixed.targets = 1,
                            random.targets = F,
                            effect = 'BOTH',
                            set.both.ratio = 0.5,
                            uniform.effect = T,
                            effect.mean = 0.5,
                            effect.sd = 0,
                            prefix = '',
                            display.genes = NULL,
                            plot.genes = 500)
{
  if (substring(out.dir, nchar(out.dir))!='/')
    out.dir <- paste0(out.dir, '/')
  if (!dir.exists(out.dir))
    stop('Output directory does not exist.')
  
  fig.dir <- paste0(out.dir, 'Figures/')
  
  if (!dir.exists(fig.dir))
    dir.create(fig.dir, recursive = T)
  
  fig.dir <- paste0(fig.dir, prefix)
  
  set.seed(seed)
  
  if (is.null(mirs.interest)){
    print('USING RANDOM miRNAs')
    if (!is.null(mir.selection))
      mir.sets <- mir.sets[names(mir.sets) %in% mir.selection]
    if (!is.null(max.target.size))
      mir.sets <- mir.sets[sapply(mir.sets, length)<max.target.size]
    mir.sets <- mir.sets[sapply(mir.sets, length)>min.target.size]
    mirs.interest <- sample(names(mir.sets), mir.set.size)
  } else
    print('USING CUSTOM miRNAs')
  
  for (i in 1:mir.set.size){
    set.seed(seed+i)
    mir.interest <- mirs.interest[i]
    targets <- unlist(mir.sets[mir.interest])
    if (mapping == T){
      targets <- mapIds(org.Hs.eg.db, keys = targets, keytype = from.id, column = to.id)
      targets <- targets[!is.na(targets)]
    }
    targets <- targets[targets %in% rownames(simulated.genes.counts)]
    
    if (random.targets == T){
      fixed.targets <- sample(1:10, 1)/10
    }
    targets <- sample(targets, floor(length(targets)*fixed.targets))
    
    effects <- rnorm(length(targets), mean = effect.mean, sd = effect.sd)
    
    if (is.null(set.both.ratio)){
      both.ratio <- sample(1:10, 1)/10
    } else{
      both.ratio <- set.both.ratio
    }
    
    # some of the mir targets are up-regulated and some are down-regulated
    if (effect == 'BOTH' & uniform.effect == F){
      up.targets <- sample(targets, floor(both.ratio * length(targets)))
      down.targets <- targets[!targets %in% up.targets]
      
      effects <- rnorm(length(up.targets), mean = effect.mean, sd = effect.sd)
      
      simulated.genes.counts[up.targets,cond1.samples] <- simulated.genes.counts[up.targets,cond1.samples] + effects
      simulated.genes.counts[up.targets,cond2.samples] <- simulated.genes.counts[up.targets,cond2.samples] - effects
      
      effects <- rnorm(length(down.targets), mean = effect.mean, sd = effect.sd)
      
      simulated.genes.counts[down.targets,cond1.samples] <- simulated.genes.counts[down.targets,cond1.samples] - effects
      simulated.genes.counts[down.targets,cond2.samples] <- simulated.genes.counts[down.targets,cond2.samples] + effects
    } 
    # the first 50% of miRNAs will have up-regulatory effect and the rest have down-regulatory effects
    else if (effect == 'BOTH' & uniform.effect == T){
      if (i <= floor(mir.set.size)*both.ratio) {
        simulated.genes.counts[targets,cond1.samples] <- simulated.genes.counts[targets,cond1.samples] + effects
        simulated.genes.counts[targets,cond2.samples] <- simulated.genes.counts[targets,cond2.samples] - effects
      } else{
        simulated.genes.counts[targets,cond1.samples] <- simulated.genes.counts[targets,cond1.samples] - effects
        simulated.genes.counts[targets,cond2.samples] <- simulated.genes.counts[targets,cond2.samples] + effects
      }
    } 
    # all the mir targets have up-regulatory effect
    else if (effect == 'UP'){
      simulated.genes.counts[targets,cond1.samples] <- simulated.genes.counts[targets,cond1.samples] + effects
      simulated.genes.counts[targets,cond2.samples] <- simulated.genes.counts[targets,cond2.samples] - effects
    }
    # all the mir targets have down-regulatory effect
    else if (effect == 'DOWN'){
      simulated.genes.counts[targets,cond1.samples] <- simulated.genes.counts[targets,cond1.samples] - effects
      simulated.genes.counts[targets,cond2.samples] <- simulated.genes.counts[targets,cond2.samples] + effects
    }
    else {
      stop('Choose only BOTH, UP or DOWN for effects.')
    }
  }
  
  if (plot == T){
    if (is.null(display.genes))
      display.genes <- sample(1:nrow(genes.counts), plot.genes)
    fig.name <- paste0(fig.dir,'simulated_genes_counts_DE.jpg')
    jpeg(fig.name)
    display.simulated <- simulated.genes.counts[display.genes,]
    heatmap.2(display.simulated, col=bluered(100), dendrogram='none', 
              Colv=NA, Rowv=NA, trace="none", colRow=NA,labCol=NA, labRow = NA)
    dev.off()
  }
  return(list('simulated_genes_counts'=simulated.genes.counts,'mirs_interest'=mirs.interest))
}

ori.genes.counts <- readRDS('../Data/TCGA-LIHC-GENE-COUNTS-NEW.RDS')
ori.covariates <- read.csv('../Data/TCGA-LIHC-COV.csv', row.names = 1)

mir.sets <- readRDS('../Data/NORMALIZED_MIRSETS_ENSEMBL.RDS')
LIHC.mirs <- rownames(readRDS('../Data/TCGA-LIHC-miRNA-RPM.RDS'))

output <- generate_simulated_data(genes.counts = ori.genes.counts,
                                  covariates = ori.covariates,
                                  out.dir = '../Simulated/NON_DE/',
                                  condition = 'shortLetterCode'
                                  )
# 
# saveRDS(output, '../Simulated/NON_DE/SIMULATED.RDS')

# output <- readRDS('../Simulated/NON_DE/SIMULATED.RDS')

simulated.genes.counts <- output$simulated_genes_counts
simulated.cov <- output$simulated_covariates
display.genes <- output$display_genes

TP.samples <- rownames(simulated.cov)[simulated.cov$shortLetterCode=='TP']
NT.samples <- rownames(simulated.cov)[simulated.cov$shortLetterCode=='NT']

type <- simulated.cov$shortLetterCode
batch <- simulated.cov$plate
design <- model.matrix(~0+type+batch)
gene.min <- apply(simulated.genes.counts, 1, min)
logcpm <- edgeR::cpm(simulated.genes.counts-gene.min)
logcpm <- limma::voom(logcpm, design)
contrast.mat <- limma::makeContrasts(contrasts = 'typeTP-typeNT',
                                     levels = design)
fit <- lmFit(logcpm, design)
fit <- contrasts.fit(fit, contrast.mat)
fit <- eBayes(fit, trend=T)
tt <- topTable(fit, coef=1, number=Inf)
sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.5)

rm(type, batch, design, gene.min, logcpm, contrast.mat, fit, tt)

pathways <- readRDS('../../Data/preprocessed/MSigDBPathGeneTab.RDS')
condition <- 'shortLetterCode'
pcxn <- readRDS('../../Data/GeneSets/improved_PCxN_MSigDB.RDS')
enriches0 <- readRDS('../Data/LIHCGenesLIHCMirsENRICHMENT.RDS')
out.dir <- '../Simulated/BOTH/'
method <- 'AggInv'

####

sampling.output <- list()

set.seed(1)
set.seed(42)
for (i in 1:100){
  print(paste0('Running iteration ', i))
  temp.output <- list()
  prefix <- paste0('new_min100max500_mean05_sd01_iteration', i, '_')
  plot <- F
  if (i%%10 == 0)
    plot <- T
  mir.set.size <- sample(40:70, 1)
  temp.output2 <- induce.miRNA.DE(simulated.genes.counts = output$simulated_genes_counts,
                                  mir.sets,
                                  cond1.samples = TP.samples,
                                  cond2.samples = NT.samples,
                                  mir.selection = LIHC.mirs,
                                  out.dir = out.dir,
                                  mirs.interest = NULL,
                                  min.target.size = 100,
                                  max.target.size = 500,
                                  mir.set.size = mir.set.size,
                                  seed = i,
                                  plot = plot,
                                  random.targets = T,
                                  effect = 'BOTH',
                                  set.both.ratio = NULL,
                                  uniform.effect = F,
                                  effect.mean = 0.5,
                                  effect.sd = 0.1,
                                  prefix = prefix,
                                  display.genes = display.genes)
  
  simulated.genes.counts <- temp.output2$simulated_genes_counts
  temp.output$mirs_interest <- temp.output2$mirs_interest
  
  type <- simulated.cov$shortLetterCode
  batch <- simulated.cov$plate
  design <- model.matrix(~0+type+batch)
  gene.min <- apply(simulated.genes.counts, 1, min)
  logcpm <- edgeR::cpm(simulated.genes.counts-gene.min)
  logcpm <- limma::voom(logcpm, design)
  contrast.mat <- limma::makeContrasts(contrasts = 'typeTP-typeNT',
                                       levels = design)
  fit <- lmFit(logcpm, design)
  fit <- contrasts.fit(fit, contrast.mat)
  fit <- eBayes(fit, trend=T)
  tt <- topTable(fit, coef=1, number=Inf)
  
  temp.output$No_of_DE_genes <- sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.5)
  temp.output$No_of_different_genes_counts <- sum(rowMeans(temp.output2$simulated_genes_counts != output$simulated_genes_counts))
  
  rm(type, batch, design, gene.min, logcpm, contrast.mat, fit, tt)
  
  temp.de.paths <- DifferentialPathwayAnalysis(simulated.genes.counts,
                                               pathways,
                                               covariates = simulated.cov,
                                               condition,
                                               adjust.covars='plate',
                                               log = F)
  
  temp.de.paths <- temp.de.paths$DEP
  
  temp.pathway.clusters <- MappingPathwaysClusters(pcxn, 
                                                   temp.de.paths,
                                                   out.dir = out.dir,
                                                   plot=plot,
                                                   subplot = F,
                                                   top.paths = 100,
                                                   prefix =prefix)
  
  temp.output3 <- miRNAPrioritization(enriches0,
                                      temp.pathway.clusters,
                                      method,
                                      out.dir=out.dir,
                                      data.dir=out.dir,
                                      samp.rate=1000,
                                      run.jack.knife = F,
                                      save.sampling = F,
                                      save.jack.knife=F,
                                      prefix = prefix,
                                      save.csv=F)
  
  temp.output$prioritized.mirna <- temp.output3
  sampling.output[[i]] <- temp.output
}

sampling.output <- sampling.output1
temp.df <- data.frame()

for (j in 1:length(sampling.output)){
  mirs.interest <- sampling.output[[j]]$mirs_interest
  clust1.mirna <- sampling.output[[j]]$prioritized.mirna$Cluster1
  clust2.mirna <- sampling.output[[j]]$prioritized.mirna$Cluster2

  clust1.mirna$adj.p <- p.adjust(clust1.mirna$AggInv_pval, method='fdr')
  clust1.mirna <- clust1.mirna[order(clust1.mirna$adj.p),]
  clust2.mirna$adj.p <- p.adjust(clust2.mirna$AggInv_pval, method='fdr')
  clust2.mirna <- clust2.mirna[order(clust2.mirna$adj.p),]

  cutoffs <- c(0.05)

  for (cutoff in cutoffs){
    clust1.top <- sum(clust1.mirna$adj.p < cutoff)
    clust2.top <- sum(clust2.mirna$adj.p < cutoff)
    # clust1.top <- 100
    # clust2.top <- 100
    
    clust1.prioritized <- clust1.mirna$x[1:clust1.top]
    clust2.prioritized <- clust2.mirna$x[1:clust2.top]
    
    clust1.recovered <- clust1.prioritized[clust1.prioritized %in% mirs.interest]
    clust2.recovered <- clust2.prioritized[clust2.prioritized %in% mirs.interest]
    
    clust1.fisher <- phyper(q = length(clust1.recovered),
                            m = length(mirs.interest),
                            n = nrow(clust1.mirna) - length(mirs.interest),
                            k = length(clust1.prioritized),
                            lower.tail = F,log.p = F)
    
    clust2.fisher <- phyper(q = length(clust2.recovered),
                            m = length(mirs.interest),
                            n = nrow(clust2.mirna) - length(mirs.interest),
                            k = length(clust2.prioritized),
                            lower.tail = F,log.p = F)
    
    combined.mirna <- unique(c(clust1.mirna$x, clust2.mirna$x))
    combined.prioritized <- unique(c(clust1.prioritized, clust2.prioritized))
    combined.recovered <- unique(c(clust1.recovered, clust2.recovered))
    combined.fisher <- phyper(q = length(combined.recovered),
                              m = length(mirs.interest),
                              n = length(combined.mirna) - length(mirs.interest),
                              k = length(combined.prioritized),
                              lower.tail = F,log.p = F)
    
    df.input <- c(clust1.top/nrow(clust1.mirna)*100, length(clust1.recovered)/length(mirs.interest)*100, clust1.fisher,
                  clust2.top/nrow(clust2.mirna)*100, length(clust2.recovered)/length(mirs.interest)*100, clust2.fisher,
                  length(combined.prioritized), length(combined.recovered)/length(mirs.interest)*100, combined.fisher)
    # df.input <- c(clust1.top/nrow(clust1.mirna), length(clust1.recovered)/length(mirs.interest), clust1.fisher,
    #               clust2.top/nrow(clust2.mirna), length(clust2.recovered)/length(mirs.interest), clust2.fisher,
    #               length(combined.prioritized), length(combined.recovered)/length(mirs.interest), combined.fisher)
    temp.df <- rbind(temp.df, df.input)
    colnames(temp.df) <- c('clust1.top', 'clust1.recovered', 'clust1.fisher',
                           'clust2.top', 'clust2.recovered', 'clust2.fisher',
                           'combined.top', 'combined.recovered', 'combined.fisher')
  }
}

df.means <- apply(temp.df, 2, mean)
df.sds <- apply(temp.df, 2, sd)
df.means
df.sds

