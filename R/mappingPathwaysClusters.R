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
#' @param plotOut if TRUE, store graph plot in Figures directory of plots
#' @param subplot if TRUE, store individual clusters plots and connected plots
#'   in Figures directory of plots
#' @param topClusters plot figures for top x clusters
#' @param outDir output directory
#' @param saveNameCSV if not NULL, saves output as csv using save name
#' @param prefix add prefix to plots
#' @param weighted True if you wish to include correlation weights in clustering
#' @return a list where the first item is a table with each row containing
#'   a pathway and its respective cluster. The second item is an igraph object.
#' @export
mappingPathwaysClusters <- function(pcxn, dePathways, clusteringFunction = NULL,
                            edgeFDR = 0.05, correlationCutOff = 0.316,
                            pathwayFDR = 0.05, topPathways = 200,
                            plotOut = TRUE, subplot = TRUE, topClusters = 2,
                            prefix = "", outDir = ".", saveNameCSV = NULL,
                            weighted = FALSE) {
    outDir <- paste0(outDir, "/")
    if (!dir.exists(outDir)) stop("Output directory does not exist.")
    figDir <- paste0(outDir, "Figures/")
    if (!dir.exists(figDir) && plotOut) dir.create(figDir, recursive = TRUE)

    tempNet     <- pcxnToNet(pcxn, edgeFDR, correlationCutOff, weighted)
    net         <- tempNet$net
    allPathways <- tempNet$allPathways
    dePathways <- dePathways[dePathways$adj.P.Val < pathwayFDR, ]
    if (nrow(dePathways) < 10) stop("Relax pathway adjusted p-values.")

    if (!is.null(topPathways)) {
        dePathways <- dePathways[seq_len(topPathways), ]
    }
    dePathways1 <- dePathways
    dePathways  <- dePathways[which(rownames(dePathways) %in% allPathways), ]

    igraph::V(net)$shape <- ifelse(igraph::V(net)$name %in% dePathways,
                                "square", "circle")
    subNet <- igraph::induced_subgraph(net, rownames(dePathways))
    # choose clustering function for nodes
    if (is.null(clusteringFunction)) {
        clusts <- igraph::cluster_edge_betweenness(subNet)
    } else {
        clusts <- utils::getFromNamespace(clusteringFunction, "igraph")(subNet)
    }
    clusts <- .levelFixerMapper(clusts)
    subNet <- .shapeColNet(subNet, dePathways1)
    igraph::V(subNet)$cluster <- clusts$membership
    pathsOut <- igraph::V(subNet)$name
    pathsOut <- as.data.frame(cbind("Pathway" = pathsOut,
                                    "cluster" = clusts$membership))
    if (plotOut == TRUE) {
        clusterPlot(subNet = subNet, subplot = subplot,
                topClusters = topClusters, outDir = figDir, prefix = prefix)
    }
    if (!is.null(saveNameCSV)) {
        utils::write.csv(pathsOut, paste0(outDir, saveNameCSV))
    }
    return(list("Clustering" = pathsOut, "DE_PCXN" = subNet,
                "Cluter_method" = eval(clusteringFunction)))
}
