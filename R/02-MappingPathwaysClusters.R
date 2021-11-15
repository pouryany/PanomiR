#' Outputs a table with pathways and their respective clusters
#'
#' @param pcxn pathways network (edge list of pathways)
#' @param de.paths differential expressed pathways, obtained from
#'   *DifferentialPathwayAnalysis*
#' @param clust.fn clustering algorithm
#' @param edge.fdr.thresh FDR threshold for pathway-pathway adjusted p-values;
#'   filter edges with adjusted p-values less than given threshold
#' @param cor.thresh correlation threshold for pathway-pathway correlation;
#'   filter pathways with correlation less than given threshold
#' @param path.fdr.thresh FDR threshold for DE pathways adjusted p-values;
#'   filter pathways with adjusted p-values less than given threshold
#' @param top.paths  use only top x paths; if NULL, use all paths
#' @param seed set seed
#' @param plot if T, store graph plot in Figures directory of plots
#' @param subplot if T, store inidividual clusters plots and connected plots
#'   in Figures directory of plots
#' @param top.clusts plot figures for top x clusters
#' @param out.dir output directory
#' @param save.csv.name if not NULL, saves output as csv using save name 
#' @param prefix add prefix to plots
#' @param weighted True if you wish to include correlation weights in clustering
#' @return a table with each row containing a pathway and its respective cluster
#' @export

MappingPathwaysClusters <- function(pcxn,
                            de.paths,
                            clust.fn=NULL,
                            edge.fdr.thresh=0.05,
                            cor.thresh=0.316,
                            path.fdr.thresh=0.05,
                            top.paths=200,
                            seed=2,
                            plot=T,
                            subplot=T,
                            top.clusts=2,
                            prefix='',
                            out.dir='',
                            save.csv.name=NULL,
                            weighted = F)
  {
  
  if (substring(out.dir, nchar(out.dir))!='/')
    out.dir <- paste0(out.dir, '/')
  if (!dir.exists(out.dir))
    stop('Output directory does not exist.')
  
  fig.dir <- paste0(out.dir, 'Figures/')
  
  if (!dir.exists(fig.dir))
    dir.create(fig.dir, recursive = T)
  
  fig.dir <- paste0(fig.dir, prefix)
  
  # filter pcxn edges with less than edge fdr threshold and correlation threshold
  pcxn$p.Adjust <- stats::p.adjust(pcxn$p.value, method='fdr')
  res <- pcxn[pcxn$p.Adjust < edge.fdr.thresh,]
  res <- res[abs(res$PathCor) > cor.thresh,]
  all.paths  <- union(unique(res$Pathway.B),unique(res$Pathway.A))
  net <- res[,c(1,2,4)]
  if(weighted){
    net$weight <- abs(net$PathCor)
  }
  net <- as.matrix(net)
  net <- graph_from_data_frame(net, directed = F)
  
  set_vertex_attr(net,"col",value = "red")
  
  # filter de pathways with less than path fdr threshold
  de.paths <- de.paths[de.paths$adj.P.Val < path.fdr.thresh,]
  if (nrow(de.paths) < 10){
    stop('Please use a more lenient threshold for differentially expressed pathways adjusted p-values.')
  }
  
  # choose top n pathways
  if (!is.null(top.paths)){
    de.paths <- de.paths[1:top.paths,]
  }
  de.paths1 <- de.paths
  de.paths <- de.paths[which(rownames(de.paths) %in% all.paths),]
  
  V(net)$shape <- ifelse(V(net)$name %in% de.paths, "square","circle")
  sub.mods <- induced_subgraph(net,rownames(de.paths))

  # choose clustering function for nodes
  set.seed(seed)
  if (is.null(clust.fn)){
    clusts <- cluster_edge_betweenness(sub.mods)
  } else{
    clusts <- clust.fn(sub.mods)
  }
  
  # setting up plot
  zz <- as.factor(clusts$membership)
  zz <- forcats::fct_infreq(zz,ordered = NA)
  levels(zz) <- 1:length(unique(zz))
  clusts$membership <- as.numeric(zz)
  
  cols <- brewer.pal(8,"Set2")
  E(sub.mods)$color <- ifelse(as.numeric(E(sub.mods)$PathCor)>0,"#E41A1C","#377EB8")
  shape.inds <- rownames(de.paths1) %in% V(sub.mods)$name
  shape.inds <- de.paths1[shape.inds,c(1,2)]
  shape.inds<- shape.inds[V(sub.mods)$name,]
  
  shape.dirs <- as.numeric(shape.inds$logFC > 0)
  names(shape.dirs) <- rownames(shape.inds)
  shape.inds <- shape.inds$logFC > 0
  
  V(sub.mods)$shape <- ifelse(shape.inds,"square","circle")
  
  if (plot==T){
    legend_cats <- data.frame(attr = c("Up-regulated","Down-regulated"),
                              shape = unique(V(sub.mods)$shape))
    node.cols   <- cols[clusts$membership]
    small.clust <- which(table(clusts$membership) <= table(clusts$membership)[5],useNames = T)
    node.cols[clusts$membership %in% small.clust] <- NA
    
    
    pdf(paste0(fig.dir,"PCxNCorGraph.pdf"),
        width = 18,height = 11)
    set.seed(seed)
    plot(sub.mods, vertex.size = 5,vertex.label =NA,
         vertex.color = node.cols)
    legend(x = "bottomleft",      ## position, also takes x,y coordinates
           legend = legend_cats$attr,
           pch =  c(0,1),
           bty = "n",cex=1.6)
    legend(x = "topleft", legend=c("Positive Cor", "Negative Cor"),
           col=c("#E41A1C","#377EB8"), lty=1, lwd =2, cex = 1.6 , bty = "n")
    
    dev.off()
  }
  
  val.tab     <- as.data.frame(ends((sub.mods),E(sub.mods)))
  val.tab$cor <- ifelse(as.numeric(E(sub.mods)$PathCor) > 0,1,0)
  val.tab$V2  <- as.character(val.tab$V2)
  val.tab$V1  <-  as.character(val.tab$V1)

  V(sub.mods)$shape <- ifelse(shape.inds,"square","circle")
  
  
  if(subplot==T){
    for(k in 1:top.clusts){
      
      keep   <- which((clusts$membership) ==k)
      sub.mods2 <- induced_subgraph(sub.mods,keep)
      if(length(V(sub.mods2)) < 2) next
      paths.out  <- V(sub.mods)$name
      paths.out  <- as.data.frame(cbind("Pathway" =paths.out,"cluster"=clusts$membership))
      
      pdf(paste0(fig.dir,"PCxNCorGraph_","Cluster_",k,".pdf"))
      
      set.seed(seed+1)
      plot(sub.mods2, edge.width= 1.3, vertex.size = 5,vertex.label =NA, 
           vertex.color = cols[clusts$membership[keep]],
           legend = T,
           layout = layout.fruchterman.reingold)
      legend(x = "bottomleft",      ## position, also takes x,y coordinates
             legend = legend_cats$attr,
             pch =  c(0,1),
             bty = "n",cex=1.4)
      legend(x = "topleft", legend=c("Positive Cor", "Negative Cor"),
             col=c("#E41A1C","#377EB8"), lty=1, lwd =2, cex = 1.4 , bty = "n")
      
      dev.off()
    }
  }
  

  '%notin%' <- Negate('%in%')
  
  remove <- which(table(clusts$membership) < 4)
  remove <- which((clusts$membership %notin% remove))
  
  sub.mods2 <- induced_subgraph(sub.mods,remove)
  
  paths.out  <- V(sub.mods)$name
  paths.out  <- as.data.frame(cbind("Pathway" =paths.out,"cluster"=clusts$membership))
  
  if(!is.null(save.csv.name))
    utils::write.csv(paths.out,paste0(out.dir,save.csv.name))
  
  if(subplot==T){
    
    pdf(paste0(fig.dir,"ConnectedPathways_PCxNCorGraph.pdf"),
        width = 18,height = 11)
    set.seed(seed-1)
    plot(sub.mods2, edge.width= 1.3, vertex.size = 5,vertex.label =NA, 
         vertex.color = cols[clusts$membership[remove]],
         legend = T,
         layout = layout_components)
    legend(x = "bottomright",   y=200,   ## position, also takes x,y coordinates
           legend = legend_cats$attr,
           pch =  c(0,1),
           bty = "n",cex=1.4)
    legend(x = "topleft", legend=c("Positive Cor", "Negative Cor"),
           col=c("#E41A1C","#377EB8"), lty=1, lwd =2, cex = 1.4 , bty = "n")
    dev.off()
  }
  
  return(paths.out)
}
