#' Plots clusters of pathways with associated directionality.
#'
#' @param subNet pathways network (edge list of pathways)
#' @param subplot if TRUE, store individual clusters plots and connected plots
#'   in Figures directory of plots
#' @param topClusters plot figures for top x clusters
#' @param outDir output directory
#' @param prefix add prefix to plots
#' @param plotSave saves the plot if set true. Otherwise display
#' @return a set of plots for DE-PCXN and subclusters
#' @examples
#' data(miniTestsPanomiR)
#' clusterPlot(miniTestsPanomiR$miniPathClusts$DE_PCXN, plotSave = FALSE)
#' @export
clusterPlot <- function(subNet, subplot = FALSE, topClusters = 2, prefix = "",
                        outDir = ".", plotSave = TRUE) {
    if (!dir.exists(outDir)) {
        stop("Output directory does not exist.")
        }
    figDir <- paste0(outDir, "/", prefix)
    legendCats <- data.frame(attr = c("Up-regulated", "Down-regulated"),
                shape = unique(igraph::V(subNet)$shape))

    cols <- RColorBrewer::brewer.pal(8, "Set2")
    clstMems  <- igraph::V(subNet)$cluster
    nodeColors <- cols[clstMems]
    smallClust <- which(table(clstMems) <= table(clstMems)[5], useNames = TRUE)

    nodeColors[clstMems %in% smallClust] <- NA

    .clusterPlotHelper(plotSave, figDir, subNet, nodeColors, legendCats)

    .clusterSubPlotHelper(subplot, topClusters, clstMems, subNet,
                        figDir, cols, legendCats)

    "%notin%" <- Negate("%in%")
    remove <- which(table(clstMems) < 4)
    remove <- which((clstMems %notin% remove))
    subNet2 <- igraph::induced_subgraph(subNet, remove)

    if (subplot == TRUE) {
        grDevices::pdf(paste0(figDir, "ConnectedPathways_PCxNCorGraph.pdf"),
            width = 18, height = 11)
        plot(subNet2, edge.width = 1.3, vertex.size = 5, vertex.label = NA,
            vertex.color = cols[clstMems[remove]], legend = TRUE,
            layout = igraph::layout_components)
        graphics::legend(x = "bottomright", y = 200, legend = legendCats$attr,
            pch = c(0, 1), bty = "n", cex = 1.4)
        graphics::legend(x = "topleft",
            legend = c("Positive Cor", "Negative Cor"),
            col = c("#E41A1C", "#377EB8"), lty = 1, lwd = 2, cex = 1.4,
            bty = "n")
        grDevices::dev.off()
    }
}
