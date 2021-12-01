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
#' @param plot if TRUE, store graph plot in Figures directory of plots
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
mappingPathwaysClusters <- function(pcxn,
                                    dePathways,
                                    clusteringFunction = NULL,
                                    edgeFDR = 0.05,
                                    correlationCutOff = 0.316,
                                    pathwayFDR = 0.05,
                                    topPathways = 200,
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

    # filter pcxn edges with edge fdr threshold and correlation threshold
    pcxn$pAdjust <- stats::p.adjust(pcxn$p.value, method = "fdr")
    res <- pcxn[pcxn$pAdjust < edgeFDR, ]
    res <- res[abs(res$PathCor) > correlationCutOff, ]
    allPathways <- union(unique(res$Pathway.B), unique(res$Pathway.A))
    net <- res[, c("Pathway.A", "Pathway.B", "PathCor")]
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
                                "square", "circle")

    subNet <- igraph::induced_subgraph(net, rownames(dePathways))

    # choose clustering function for nodes
    if (is.null(clusteringFunction)) {
        clusts <- igraph::cluster_edge_betweenness(subNet)
    } else {
        clusts <- utils::getFromNamespace(clusteringFunction, "igraph")(subNet)
    }

    # setting up plot
    zz <- as.factor(clusts$membership)
    zz <- forcats::fct_infreq(zz, ordered = NA)
    levels(zz) <- seq_along(unique(zz))
    clusts$membership <- as.numeric(zz)

    igraph::E(subNet)$color <-
        ifelse(as.numeric(igraph::E(subNet)$PathCor) > 0, "#E41A1C", "#377EB8")

    shapeIndex <- rownames(dePathways1) %in% igraph::V(subNet)$name
    shapeIndex <- dePathways1[shapeIndex, c(1, 2)]
    shapeIndex <- shapeIndex[igraph::V(subNet)$name, ]

    shapeDirs <- as.numeric(shapeIndex$logFC > 0)
    names(shapeDirs) <- rownames(shapeIndex)
    shapeIndex <- shapeIndex$logFC > 0

    igraph::V(subNet)$shape   <- ifelse(shapeIndex, "square", "circle")
    igraph::V(subNet)$cluster <- clusts$membership
    pathsOut <- igraph::V(subNet)$name
    pathsOut <- as.data.frame(cbind(
        "Pathway" = pathsOut,
        "cluster" = clusts$membership
    ))

    if (plot == TRUE) {
        clusterPlot(subNet = subNet, subplot = subplot,
                topClusters = topClusters, outDir = figDir, prefix = prefix)
    }

    if (!is.null(saveNameCSV)) {
        utils::write.csv(pathsOut, paste0(outDir, saveNameCSV))
    }

    return(list("Clustering" = pathsOut, "DE-PCXN" = subNet,
                "Cluter_method" = eval(clusteringFunction)))
}
