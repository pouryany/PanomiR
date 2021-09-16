library(dplyr)
library(org.Hs.eg.db)
library(pROC)

iter <- c('x_50_LIHCGenesLIHCMirs_1000_samples_clustNo_',
          'x2_all_LIHCGenesLIHCMirs_1000_samples_clustNo_',
          '20_UP_20_DOWN_x2_1000_samples_clustNo_',
          '30_UP_30_DOWN_x2_1000_samples_clustNo_',
          '40_UP_40_DOWN_x2_1000_samples_clustNo_',
          'x2_all_top100_1000_samples_clustNo_',
          'x2_all_top100_connected_1000_samples_clustNo_',
          'x2_all_cor05_1000_samples_clustNo_',
          'x2_all_cor05_connected_1000_samples_clustNo_')

names <- c('x_50',
           'x2_all',
           '20_up_down',
           '30_up_down',
           '40_up_down',
           'top100',
           'top100_connected',
           'cor05',
           'cor05_connected'
           )

# iter <- c('x2_all_BRCAGenesBRCAMirs_1000_samples_clustNo_')
# names <- c('x2_all')

# cutoffs <- c(0.05, 0.01, 0.005, 0.001, 0.00005, 0.00001)
cutoffs <- c(100)

output <- list()

DE.mirna.data <- readRDS(paste0('../Data/TCGA-LIHC-DEmiRNA.RDS'))
# sig.DE.mirna <- rownames(DE.mirna.data[DE.mirna.data$adj.P.Val<0.05 & abs(DE.mirna.data$logFC) > 0.5,])
sig.DE.mirna <- rownames(DE.mirna.data)[1:30]

for (cutoff in cutoffs){
  # cutoff <- cutoffs[1]
  
  sig.P.values <- c()
  
  listname <- paste0('cutoff', cutoff)
  listname <- gsub('\\.*','', listname)
  
  clust <- 2
  
  for (i in 1:clust){
    cat(paste0('Cluster ', i, '\n'))
    for (j in 1:length(iter)){
      name <- names[j]
      filename <- paste0('../Output/Prioritization/',iter[j], i,'.csv')
      if (!file.exists(filename))
        next
      clust.mirna <- read.csv(filename, row.names = 1)
      clust.mirna$adj.p <- p.adjust(clust.mirna$AggInv_pval, method='fdr')
      clust.mirna <- clust.mirna[order(clust.mirna$adj.p),]
      # top <- nrow(clust.mirna[clust.mirna$adj.p < cutoff,])
      top <- 30
      cat(paste0(filename, ' loaded.\n'))
      cat(paste0('ROWS: ', nrow(clust.mirna), ' TOP: ', top, '\n'))
      # roc_obj <- roc(response=ifelse(clust.mirna[1:top,]$x %in% sig.DE.mirna, 'DE', 'non-DE'), predictor=clust.mirna[1:top,]$adj.p, direction = '>', levels=c('non-DE', 'DE'))
      # cat(paste0(name, ', AUC: ', auc(roc_obj), '\n'))
      fisher.p <- fisher(clust.mirna, DE.mirna=sig.DE.mirna, top, name)
      if (fisher.p < 0.05)
        sig.P.values <- c(sig.P.values, paste0(name, ' Cluster ', i, ' FISHER P-val: ', fisher.p, ' TOP: ', top))
      # mwu.p <- mann.whitney(clust.mirna, DE.mirna.data, top, name)
      # if (mwu.p < 0.05)
      #   sig.P.values <- c(sig.P.values, paste0(name, ' Cluster ', i, ' MWU P-val: ', mwu.p, ' TOP: ', top))
    }
  }
  if (is.null(sig.P.values))
    sig.P.values <- c('Empty')
  output[[listname]] <- sig.P.values
}

## FISHER EXACT TEST

fisher <- function(clust.mirna, DE.mirna, top=50, name=''){
  clust.top <- as.character(clust.mirna[1:top,]$x)
  q <- sum(clust.top %in% DE.mirna)
  m <- length(DE.mirna)
  n <- nrow(clust.mirna) - m
  k <- length(clust.top)
    
  pval <- phyper(q,m,n,k,lower.tail = F,log.p = F)
  cat(paste0(name, ' FISHER: ',pval, ', intersect: ', q, '\n'))
  if (pval == 0 & q < 10)
    pval <- 1
  return (pval)
}

## MANN WHITNEY TEST

mann.whitney <- function(clust.mirna, DE.mirna.data, top=50, name=''){
  clust.mirna <- clust.mirna[1:top,]
  clust.mirna.in.DE <- DE.mirna.data[rownames(DE.mirna.data) %in% clust.mirna$x,]
  other.mirna <- DE.mirna.data[!rownames(DE.mirna.data) %in% clust.mirna$x,]
  res <- wilcox.test(abs(clust.mirna.in.DE$t), abs(other.mirna$t), alternative='greater')
  cat(paste0(name, ' MWU: ', res$p.value, ', # in DE: ', nrow(clust.mirna.in.DE), ', # not in DE: ', nrow(other.mirna), '\n'))
  return (res$p.value)
}


