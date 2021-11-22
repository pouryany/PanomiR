#' Plots clusters of pathways with associated directionality.
#'
#' @param subNet pathways network (edge list of pathways)
#' @param subplot if TRUE, store individual clusters plots and connected plots
#'   in Figures directory of plots
#' @param topClusters plot figures for top x clusters
#' @param outDir output directory
#' @param prefix add prefix to plots
#' @return a set of plots for DE-PCXN and subclusters
#' @export
clusterPlot <- function(subNet,
                        plot = TRUE,
                        subplot = FALSE,
                        topClusters = 2,
                        prefix = "",
                        outDir = "") {
    if (substring(outDir, nchar(outDir)) != "/") {
        outDir <- paste0(outDir, "/")
    }
    if (!dir.exists(outDir)) {
        stop("Output directory does not exist.")
    }
    
    figDir <- outDir

    legend_cats <- data.frame(
        attr = c("Up-regulated", "Down-regulated"),
        shape = unique(igraph::V(subNet)$shape)
    )
    cols <- RColorBrewer::brewer.pal(8, "Set2")
    clustMems  <- igraph::V(subNet)$cluster
    nodeColors <- cols[clustMems]
    
    small.clust <-
        which(table(clustMems) <= table(clustMems)[5],
              useNames = TRUE)
    
    nodeColors[clustMems %in% small.clust] <- NA
    
    grDevices::pdf(
        paste0(figDir, "PCxNCorGraph.pdf"),
        width = 18, height = 11
    )
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

    if (subplot == TRUE) {
        for (k in seq_len(topClusters)) {
            keep     <- which((clustMems) == k)
            subNet2 <- igraph::induced_subgraph(subNet, keep)
            
            if (length(igraph::V(subNet2)) < 2) next
      
            grDevices::pdf(
                paste0(figDir, "PCxNCorGraph_", "Cluster_", k, ".pdf")
            )
            
            plot(subNet2,
                 edge.width = 1.3, vertex.size = 5, vertex.label = NA,
                 vertex.color = cols[clustMems[keep]],
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
     
     remove <- which(table(clustMems) < 4)
     remove <- which((clustMems %notin% remove))
    # 
     subNet2 <- igraph::induced_subgraph(subNet, remove)
  
    if (subplot == TRUE) {
        grDevices::pdf(paste0(figDir, "ConnectedPathways_PCxNCorGraph.pdf"),
                       width = 18, height = 11
        )
        plot(subNet2,
             edge.width = 1.3, vertex.size = 5, vertex.label = NA,
             vertex.color = cols[clustMems[remove]],
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
    return(NULL)
}
