#' Prioritize miRNA
#'
#' Outputs a table of miRNA (ordered) with respective p-values derived from
#'   method for prioritization
#'
#' @param enriches0 miRNA-pathway enrichment dataset obtained from
#'   miRNAPathwayEnrichment.
#' @param pathway.clusters Pathway clusters, obtained from
#'   MappingPathwaysClusters.
#' @param method Vector of methods (pCut, AggInv, AggLog, sumz, sumlog).
#' @param method.thresh Vector of methods threshold for each method in method,
#'   if NULL use default thresh values in method.
#' @param mir.path.fdr.thresh for calculating miRNA-pathway hits in the input
#'   cluster based on standard enrichment analysis.
#' @param top.clust Top x clusters to perform miRNA prioritization on.
#' @param samp.rate Sampling rate for CLT.
#' @param num.cores Number of CPU cores to use, must be at least one.
#' @param out.dir Output directory.
#' @param data.dir Data directory.
#' @param save.csv If T, saves CSV file for each cluster in top.clust in
#'   out.dir.
#' @param save.sampling If T, saves sampling data as RDS for each cluster in
#'   top.clust in data.dir.
#' @param run.jack.knife If True, jacknifing will be performed.
#' @param save.jack.knife If T, saves jack-knifed sampling data as RDS for each
#'   cluster in top.clust in data.dir.
#' @param prefix Prefix for all saved data.
#' @return Table of miRNA and p-values, each row contains a miRNA and its
#'   associated p-values from the methods.
#' @export
PrioritizeMicroRNA <- function(enriches0,
                               pathway.clusters,
                               method = "AggInv",
                               method.thresh = NULL,
                               mir.path.fdr.thresh = 0.25,
                               top.clust = 2,
                               samp.rate = 1000,
                               out.dir = "",
                               data.dir = "",
                               save.sampling = T,
                               run.jack.knife = T,
                               save.jack.knife = F,
                               num.cores = 1,
                               save.csv = T,
                               prefix = "") {
  if (substring(out.dir, nchar(out.dir)) != "/") {
    out.dir <- paste0(out.dir, "/")
  }
  if (!dir.exists(out.dir)) {
    warning("Output directory does not exist.")
    dir.create(out.dir, recursive = TRUE)
  }

  if (out.dir == "/") {
    out.dir <- ""
  }

  if (substring(data.dir, nchar(data.dir)) != "/") {
    data.dir <- paste0(data.dir, "/")
  }
  if (!dir.exists(data.dir)) {
    stop("Data directory does not exist.")
  }

  if (data.dir == "/") {
    data.dir <- ""
  }

  output <- list()

  # count miRNA-pathway enrichment with p-value less than threshold
  enriches <- enriches0 %>% dplyr::filter(., Intersect != 0)
  enriches %<>% dplyr::group_by(., y) %>% dplyr::mutate(., path_fdr = stats::p.adjust(pval, method = "fdr"))
  enriches <- enriches %>% dplyr::mutate(., hit = ifelse(path_fdr < mir.path.fdr.thresh, 1, 0))

  # Need to fix this function to be able to work on specific clusters instead of all top clusters.
  for (clustNo in 1:top.clust) {
    clustName <- paste0("Cluster", clustNo)

    print(paste0("Working on ", clustName, "."))

    # select pathways in cluster
    pathways <- as.character(pathway.clusters[pathway.clusters$cluster == clustNo, ]$Pathway)
    n_paths <- length(pathways)

    # formulate number of miRNA-pathway enrichment with p-value less than threshold for each miRNA
    temp.enrich <- enriches[enriches$y %in% pathways, ]
    selector <- temp.enrich %>%
      dplyr::group_by(x) %>%
      dplyr::summarise(., "cluster_hits" = sum(hit))

    # perform p-value aggregation based on methodlogy provided
    for (i in 1:length(method)) {
      m <- method[i]

      print(paste0("Performing ", m, " function."))
      fn <- get(paste0(m, ".fn"))
      cover.fn <- get(paste0(m, ".cover.fn"))

      if (!is.null(method.thresh)) {
        m.thresh <- method.thresh[i]
        temp <- fn(enriches = enriches0, pathways, is.selector = T, thresh = m.thresh)
      } else {
        temp <- fn(enriches = enriches0, pathways, is.selector = T)
      }

      m.selector <- temp$selector
      m.enriches0 <- temp$enriches0

      if (nrow(m.selector) < 3) {
        print(paste0("Skipping ", m, " function due to low number of miRNA after filter"))
        next
      }

      enrich.null <- m.enriches0 %>% dplyr::filter(., x %in% m.selector$x)

      sampling.data.dir <- paste0(data.dir, prefix, "Sampling_Data/")

      if (save.sampling == TRUE) {
        if (!dir.exists(sampling.data.dir)) dir.create(sampling.data.dir, recursive = T)
      }

      sampling.data.filename <- paste0(prefix, m, "_", samp.rate, "_samples.RDS")
      sampling.data.file <- paste0(sampling.data.dir, "/", sampling.data.filename)

      if (m %in% c("AggInv", "AggLog")) {
        # perform sampling
        sampling.data <- samplingDataBase(enrich.null,
          m.selector,
          samp.rate,
          fn,
          n_paths,
          sampling.data.file,
          jack.knife = FALSE,
          save.sampling = save.sampling,
          num.cores = 8
        )

        m.selector <- methodProbBase(
          sampling.data = sampling.data[[paste0("SampSize_", 100)]],
          selector = m.selector,
          m = m,
          n_paths = n_paths,
          cover.fn = cover.fn
        )
      } else {
        names(m.selector)[2:3] <- paste0(m, "_", names(m.selector)[2:3])
      }


      print(paste0(m, " Method Done"))

      # perform jack-knife
      if (run.jack.knife == T) {
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
          jack.knife = FALSE,
          save.sampling = save.sampling,
          num.cores = 8
        )
        m.selector <- jackKnifeBase(
          selector = m.selector,
          pathways = pathways,
          enrich.null = enrich.null,
          fn = fn,
          jack.knife.data =
            sampling.data[[paste0("SampSize_", 100)]],
          m = m,
          num.cores = 8
        )

        print(paste0(m, " JackKnifing Method Done!"))

        m.selector <- m.selector[, c(1, 3, 2, 4)]
      } else {
        m.selector <- m.selector[, c(1, 3, 2)]
      }

      selector <- merge(selector, m.selector, all = T)

      met <- paste0(m, "_pval")
      selector <- selector %>% dplyr::arrange(., !!rlang::sym(met))
      met2 <- paste0(m, "_fdr")
      selector <- selector %>% dplyr::mutate(
        .,
        !!met2 := stats::p.adjust(
          p = !!rlang::sym(met),
          method = "fdr"
        )
      )
    }

    if (save.csv == T) {
      save.name <- paste0(prefix, samp.rate, "_samples_clustNo_", clustNo, ".csv")
      print(paste0(save.name, " saved!"))
      utils::write.csv(selector, paste0(out.dir, save.name))
    }

    output[[clustName]] <- selector
  }

  return(output)
}
