#' Outputs a table with pathways and their respective clusters
#'
#' @param pcxn pathways network (edge list of pathways)
#' @param dePathways differential expressed pathways, obtained from
#'   *DifferentialPathwayAnalysis*
#' @param clusteringFunction clustering algorithm
#' @param edgeFDR FDR threshold for pathway-pathway adjusted p-values;
#'   filter edges with adjusted p-values less than given threshold
#' @param correlationCutOff cut-off threshold for pathway-pathway correlation;
#'   filter pathways with correlation less than given threshold
#' @param pathwayFDR FDR threshold for DE pathways adjusted p-values;
#'   filter pathways with adjusted p-values less than given threshold
#' @param topPathways  use only top x paths; if NULL, use all paths
#' @param seed set seed
#' @param plot if TRUE, store graph plot in Figures directory of plots
#' @param subplot if TRUE, store inidividual clusters plots and connected plots
#'   in Figures directory of plots
#' @param topClusters plot figures for top x clusters
#' @param outDir output directory
#' @param saveNameCSV if not NULL, saves output as csv using save name
#' @param prefix add prefix to plots
#' @param weighted True if you wish to include correlation weights in clustering
#' @return a table with each row containing a pathway and its respective cluster
#' @export
mappingPathwaysClusters <- function(pcxn,
                                    dePathways,
                                    clusteringFunction = NULL,
                                    edgeFDR = 0.05,
                                    correlationCutOff = 0.316,
                                    pathwayFDR = 0.05,
                                    topPathways = 200,
                                    seed = 2,
                                    plot = TRUE,
                                    subplot = TRUE,
                                    topClusters = 2,
                                    prefix = "",
                                    outDir = "",
                                    saveNameCSV = NULL,
                                    weighted = FALSE) {
    if (substring(outDir, nchar(outDir)) != "/") {
        outDir <- paste0(outDir, "/")
    }
    if (!dir.exists(outDir)) {
        stop("Output directory does not exist.")
    }
    
    figDir <- paste0(outDir, "Figures/")
    
    if (!dir.exists(figDir)) {
        dir.create(figDir, recursive = TRUE)
    }
    
    figDir <- paste0(figDir, prefix)
    
    # filter pcxn edges with edge fdr threshold and correlation threshold
    pcxn$p.Adjust <- stats::p.adjust(pcxn$p.value, method = "fdr")
    res <- pcxn[pcxn$p.Adjust < edgeFDR, ]
    res <- res[abs(res$PathCor) > correlationCutOff, ]
    allPathways <- union(unique(res$Pathway.B), unique(res$Pathway.A))
    net <- res[, c(1, 2, 4)]
    if (weighted) {
        net$weight <- abs(net$PathCor)
    }
    net <- as.matrix(net)
    net <- igraph::graph_from_data_frame(net, directed = FALSE)
    
    igraph::set_vertex_attr(net, "col", value = "red")
    
    # filter de pathways with less than path fdr threshold
    dePathways <- dePathways[dePathways$adj.P.Val < pathwayFDR, ]
    if (nrow(dePathways) < 10) {
      stop(
        "Please use a more lenient threshold for pathway adjusted p-values."
      )
    }

    # choose top n pathways
    if (!is.null(topPathways)) {
        dePathways <- dePathways[seq_len(topPathways), ]
    }
    dePathways1 <- dePathways
    dePathways  <- dePathways[which(rownames(dePathways) %in% allPathways), ]
    
    igraph::V(net)$shape <- ifelse(igraph::V(net)$name %in% dePathways,
                                   "square",
                                   "circle")
    subNet <- igraph::induced_subgraph(net, rownames(dePathways))
    
    # choose clustering function for nodes
    set.seed(seed)
    if (is.null(clusteringFunction)) {
        clusts <- igraph::cluster_edge_betweenness(subNet)
    } else {
        clusts <- clusteringFunction(subNet)
    }
    
    # setting up plot
    zz <- as.factor(clusts$membership)
    zz <- forcats::fct_infreq(zz, ordered = NA)
    levels(zz) <- seq_along(unique(zz))
    clusts$membership <- as.numeric(zz)
    
    cols <- RColorBrewer::brewer.pal(8, "Set2")
    
    igraph::E(subNet)$color <-
        ifelse(as.numeric(igraph::E(subNet)$PathCor) > 0, "#E41A1C", "#377EB8")
    
    shapeIndex <- rownames(dePathways1) %in% igraph::V(subNet)$name
    shapeIndex <- dePathways1[shapeIndex, c(1, 2)]
    shapeIndex <- shapeIndex[igraph::V(subNet)$name, ]
    
    shapeDirs <- as.numeric(shapeIndex$logFC > 0)
    names(shapeDirs) <- rownames(shapeIndex)
    shapeIndex <- shapeIndex$logFC > 0
    
    igraph::V(subNet)$shape <- ifelse(shapeIndex, "square", "circle")
    
    if (plot == TRUE) {
        legend_cats <- data.frame(
            attr = c("Up-regulated", "Down-regulated"),
            shape = unique(igraph::V(subNet)$shape)
        )
        nodeColors <- cols[clusts$membership]
        
        small.clust <-
            which(table(clusts$membership) <= table(clusts$membership)[5],
                  useNames = TRUE)
        
        nodeColors[clusts$membership %in% small.clust] <- NA

        grDevices::pdf(
            paste0(figDir, "PCxNCorGraph.pdf"),
            width = 18, height = 11
        )
        set.seed(seed)
        plot(
            subNet,
            vertex.size = 5, vertex.label = NA,
            vertex.color = nodeColors
        )
        graphics::legend(
            x = "bottomleft", # position, also takes x,y coordinates
            legend = legend_cats$attr,
            pch = c(0, 1),
            bty = "n",
            cex = 1.6
        )
        graphics::legend(
            x = "topleft",
            legend = c("Positive Cor", "Negative Cor"),
            col = c("#E41A1C", "#377EB8"),
            lty = 1,
            lwd = 2,
            cex = 1.6,
            bty = "n"
        )
        grDevices::dev.off()
    }
    
    valTab <- as.data.frame(igraph::ends((subNet), igraph::E(subNet)))
    valTab$cor <- ifelse(as.numeric(igraph::E(subNet)$PathCor) > 0, 1, 0)
    valTab$V2 <- as.character(valTab$V2)
    valTab$V1 <- as.character(valTab$V1)
    
    igraph::V(subNet)$shape <- ifelse(shapeIndex, "square", "circle")

    if (subplot == TRUE) {
        for (k in seq_len(topClusters)) {
            keep     <- which((clusts$membership) == k)
            subNet2 <- igraph::induced_subgraph(subNet, keep)
            
            if (length(igraph::V(subNet2)) < 2) next
            
            pathsOut <- igraph::V(subNet)$name
            
            pathsOut <-
                as.data.frame(cbind("Pathway" = pathsOut,
                                    "cluster" = clusts$membership))
            
            grDevices::pdf(
                paste0(figDir, "PCxNCorGraph_", "Cluster_", k, ".pdf")
            )
            
            set.seed(seed + 1)
            plot(subNet2,
                 edge.width = 1.3, vertex.size = 5, vertex.label = NA,
                 vertex.color = cols[clusts$membership[keep]],
                 legend = TRUE,
                 layout = igraph::layout.fruchterman.reingold
            )
            graphics::legend(
                x = "bottomleft", # position, also takes x,y coordinates
                legend = legend_cats$attr,
                pch = c(0, 1),
                bty = "n",
                cex = 1.4
            )
            graphics::legend(
                x = "topleft",
                legend = c("Positive Cor", "Negative Cor"),
                col = c("#E41A1C", "#377EB8"),
                lty = 1,
                lwd = 2,
                cex = 1.4,
                bty = "n"
            )
            grDevices::dev.off()
        }
    }
    
    "%notin%" <- Negate("%in%")
    
    remove <- which(table(clusts$membership) < 4)
    remove <- which((clusts$membership %notin% remove))
    
    subNet2 <- igraph::induced_subgraph(subNet, remove)
    
    pathsOut <- igraph::V(subNet)$name
    
    pathsOut <-
        as.data.frame(cbind("Pathway" = pathsOut,
                            "cluster" = clusts$membership))
    
    if (!is.null(saveNameCSV)) {
        utils::write.csv(pathsOut, paste0(outDir, saveNameCSV))
    }
    
    if (subplot == TRUE) {
        grDevices::pdf(paste0(figDir, "ConnectedPathways_PCxNCorGraph.pdf"),
                       width = 18, height = 11
        )
        set.seed(seed - 1)
        plot(subNet2,
             edge.width = 1.3, vertex.size = 5, vertex.label = NA,
             vertex.color = cols[clusts$membership[remove]],
             legend = TRUE,
             layout = igraph::layout_components
        )
        graphics::legend(
            x = "bottomright",
            y = 200, ## position, also takes x,y coordinates
            legend = legend_cats$attr,
            pch = c(0, 1),
            bty = "n",
            cex = 1.4
        )
        graphics::legend(
            x = "topleft",
            legend = c("Positive Cor", "Negative Cor"),
            col = c("#E41A1C", "#377EB8"),
            lty = 1,
            lwd = 2,
            cex = 1.4,
            bty = "n"
        )
        grDevices::dev.off()
    }
    return(pathsOut)
}
