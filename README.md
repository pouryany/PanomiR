<!-- badges: start -->

[![R-CMD-check](https://github.com/pouryany/PanomiR/workflows/R-CMD-check/badge.svg)](https://github.com/pouryany/PanomiR/actions)
[![lint](https://github.com/pouryany/PanomiR/workflows/lint/badge.svg)](https://github.com/pouryany/PanomiR/actions)
<!-- badges: end -->

## Introduction

PanomiR is a package for pathway and microRNA Analysis of gene
expression data. This document provides details about how to install and
utilize various functionality in PanomiR.

For questions, comments, and other queries, contact <pouryany@gmail.com>

## Installation

PanomiR can be accessed via Bioconductor. To install, start R (version
\> “4.1) and run the following code.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("PanomiR")
```

You can also install the latest development version of PanomiR using
GitHub.

``` r
devtools::install_github("pouryany/PanomiR")
```

## 1. Pathway summarization

PanomiR can generate pathway activity profiles given a gene expression
dataset and a list of pathways. This section uses an example dataset
from the Cancer Genome Atlas (TCGA) Liver Hepatocellular Carcinoma
(LIHC) dataset to generate Pathway summary statistics.

``` r
library(PanomiR)

# Downloading example LIHC expression dataset
LIHC_url  <- url(paste0("https://github.com/pouryany/PanomiR_paper",
                        "/raw/main/data/TCGA_LIHC.RDS"))
TCGA_LIHC <- readRDS(LIHC_url)

# Downloading example LIHC covariates dataset
LIHC_cov_url  <- url(paste0("https://github.com/pouryany/PanomiR_paper",
                        "/raw/main/data/covariates_TCGA_LIHC.RDS"))
cov_TCGA_LIHC <- readRDS(LIHC_cov_url)

# Pathway reference from the PanomiR package
data("path_gene_table")

# Generating pathway summary statistics 

summaries <- pathwaySummary(TCGA_LIHC,path_gene_table,method = "x2")

head(summaries)[,1:2]
```

    ##                                                       TCGA-2V-A95S-01A-11R-A37K-07
    ## Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS                                  107100510
    ## Pathway.KEGG_CITRATE_CYCLE_TCA_CYCLE                                     137540732
    ## Pathway.KEGG_PENTOSE_PHOSPHATE_PATHWAY                                   105321837
    ## Pathway.KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS                    115609938
    ## Pathway.KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM                              78860500
    ## Pathway.KEGG_GALACTOSE_METABOLISM                                         88480592
    ##                                                       TCGA-2Y-A9GS-01A-12R-A38B-07
    ## Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS                                  131945603
    ## Pathway.KEGG_CITRATE_CYCLE_TCA_CYCLE                                     149435696
    ## Pathway.KEGG_PENTOSE_PHOSPHATE_PATHWAY                                   135163359
    ## Pathway.KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS                    162386997
    ## Pathway.KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM                              92557981
    ## Pathway.KEGG_GALACTOSE_METABOLISM                                        103639806

## 2. Differential Pathway activation

Once you generate the pathway activity profiles, as discussed in the
last section, there are several analysis that you can perform. We have
bundled some of the most important ones into standalone functions. Here,
we describe differential pathway activation profiling, which is
examining differenes in pathway activity profiles in user-determined
conditions.

``` r
output0 <- differentialPathwayAnalysis(geneCounts = TCGA_LIHC,
                                       pathways =  path_gene_table,
                                       covariates = cov_TCGA_LIHC, 
                                       condition = 'shortLetterCode',
                                       adjustCovars ='plate')

de.paths <- output0$DEP

head(de.paths)
```

    ##                                                           logFC   AveExpr
    ## Pathway.REACTOME_NUCLEAR_SIGNALING_BY_ERBB4          -0.4006321 -0.708800
    ## Pathway.KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION -0.4031298 -2.381353
    ## Pathway.KEGG_JAK_STAT_SIGNALING_PATHWAY              -0.3663317 -1.076826
    ## Pathway.REACTOME_CLASS_A1_RHODOPSIN_LIKE_RECEPTORS   -0.4057902 -2.098614
    ## Pathway.REACTOME_GPCR_LIGAND_BINDING                 -0.3692202 -1.998562
    ## Pathway.REACTOME_HDL_MEDIATED_LIPID_TRANSPORT        -0.9826705  1.460851
    ##                                                              t      P.Value
    ## Pathway.REACTOME_NUCLEAR_SIGNALING_BY_ERBB4          -13.26122 1.633566e-33
    ## Pathway.KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION -13.24757 1.853574e-33
    ## Pathway.KEGG_JAK_STAT_SIGNALING_PATHWAY              -12.58752 7.847204e-31
    ## Pathway.REACTOME_CLASS_A1_RHODOPSIN_LIKE_RECEPTORS   -12.30811 9.768330e-30
    ## Pathway.REACTOME_GPCR_LIGAND_BINDING                 -12.24499 1.720745e-29
    ## Pathway.REACTOME_HDL_MEDIATED_LIPID_TRANSPORT        -11.89976 3.721198e-28
    ##                                                         adj.P.Val        B
    ## Pathway.REACTOME_NUCLEAR_SIGNALING_BY_ERBB4          1.130680e-30 65.36049
    ## Pathway.KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION 1.130680e-30 65.23564
    ## Pathway.KEGG_JAK_STAT_SIGNALING_PATHWAY              3.191196e-28 59.26021
    ## Pathway.REACTOME_CLASS_A1_RHODOPSIN_LIKE_RECEPTORS   2.979341e-27 56.76963
    ## Pathway.REACTOME_GPCR_LIGAND_BINDING                 4.198618e-27 56.21046
    ## Pathway.REACTOME_HDL_MEDIATED_LIPID_TRANSPORT        7.566436e-26 53.17513
    ##                                                                                 contrast
    ## Pathway.REACTOME_NUCLEAR_SIGNALING_BY_ERBB4          shortLetterCodeTP-shortLetterCodeNT
    ## Pathway.KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION shortLetterCodeTP-shortLetterCodeNT
    ## Pathway.KEGG_JAK_STAT_SIGNALING_PATHWAY              shortLetterCodeTP-shortLetterCodeNT
    ## Pathway.REACTOME_CLASS_A1_RHODOPSIN_LIKE_RECEPTORS   shortLetterCodeTP-shortLetterCodeNT
    ## Pathway.REACTOME_GPCR_LIGAND_BINDING                 shortLetterCodeTP-shortLetterCodeNT
    ## Pathway.REACTOME_HDL_MEDIATED_LIPID_TRANSPORT        shortLetterCodeTP-shortLetterCodeNT

## 3. Finding clusters of pathways

PanomiR provides a function to find groups coordinated differentially
activated pathways based on the PCxN methodology.

``` r
# using an updated version of pcxn 
pcxn_url  <- url(paste0("https://github.com/pouryany/PanomiR_paper",
                        "/raw/main/data/pcxn_panomir.RDS"))
pcxn_net  <- readRDS(pcxn_url)

set.seed(2)
pathwayClustsLIHC <- mappingPathwaysClusters(pcxn = pcxn_net, 
                            dePathways = de.paths[1:300,],
                            topPathways = 200,
                            outDir=".",
                            plotOut = FALSE,
                            subplot = FALSE,
                            prefix='',
                            clusteringFunction = "cluster_louvain",
                            correlationCutOff = 0.1)


head(pathwayClustsLIHC$Clustering)
```

    ##                                        Pathway cluster
    ## 1           Pathway.KEGG_TRYPTOPHAN_METABOLISM       6
    ## 2         Pathway.KEGG_BETA_ALANINE_METABOLISM       3
    ## 3        Pathway.KEGG_LINOLEIC_ACID_METABOLISM       3
    ## 4              Pathway.KEGG_RETINOL_METABOLISM       3
    ## 5 Pathway.KEGG_DRUG_METABOLISM_CYTOCHROME_P450       3
    ## 6                 Pathway.KEGG_DNA_REPLICATION       2

## 4. Prioritizing miRNAs per cluster of pathways.

PanomiR identifies miRNAs that target clusters of pathways, as defined
in the last section. In order to this, you would need a reference table
of miRNA-Pathway association score (enrichment). We recommend using a
customized miRNA-Pathway association table, tailored to your
experimental data. Here, we provide a pre-processed table for LIHC table
and in the next section, we will explain how to generate the customized
tables.

Note that in the example below, we use a sampling rate of 50, the
recommended rate is between 500-1000. Also, we set the saveSampling
argument to FALSE. This argument should be set to TRUE if you wish to
save your sampling and check for different outputs from the clustering
algorithms or pathway thresholds.

``` r
# using an updated version of pcxn 
enrch_url    <- url(paste0("https://github.com/pouryany/PanomiR_paper",
                        "/raw/main/data/LIHC_ENRICHMENT_TargetScan03.RDS"))
tableEnrich  <- readRDS(enrch_url)

set.seed(1)
output2 <- prioritizeMicroRNA(enriches0 = tableEnrich,
                              pathClust = pathwayClustsLIHC$Clustering,
                              topClust = 1,
                              sampRate = 50, 
                              method = c("aggInv"),
                              outDir = "Output/",
                              dataDir = "outData/",
                              saveSampling = FALSE,
                              runJackKnife = FALSE,
                              numCores = 1,
                              prefix = "outmiR",
                              saveCSV = FALSE)
```

    ## [1] "Working on Cluster1."
    ## [1] "Performing aggInv function."
    ## [1] "aggInv Method Done"

``` r
head(output2$Cluster1)
```

    ##                  x cluster_hits aggInv_cover  aggInv_pval   aggInv_fdr
    ## 1  hsa-miR-371a-5p            0  -0.98581204 2.166514e-46 9.402670e-44
    ## 2 hsa-miR-505-3p.2            0  -1.16947734 6.777526e-31 1.470723e-28
    ## 3   hsa-miR-556-5p            0  -0.89339753 6.940239e-28 1.004021e-25
    ## 4  hsa-miR-1298-5p            0  -1.94996850 8.025049e-27 8.707179e-25
    ## 5     hsa-miR-1278            1  -0.40782882 2.378328e-26 2.064388e-24
    ## 6   hsa-miR-325-3p            3  -0.08758886 2.755195e-25 1.992925e-23

## 5. miRNA-Pathway enrichment tables

PanomiR best performs on tissue/experiment-customized datasets. In order
to do this, you need to create a customized enrichment table. You can
simply do so by using the pathway and miRNA list that we have provided
as a part of the package. simply, plug in the name of the genes present
(expressed) in your experiment in the following code

``` r
# using an updated version of pcxn 
data("msigdb_c2")
data("targetScan_03")


customeTableEnrich <- miRNAPathwayEnrichment(mirSets = targetScan_03,
                                              pathwaySets = msigdb_c2,
                                              geneSelection = yourGenes,
                                              mirSelection = yourMicroRNAs,
                                              fromID = "ENSEMBL",
                                              toID = "ENTREZID",
                                              minPathSize = 9,
                                              numCores = 1,
                                              outDir = ".",
                                              saveOutName = NULL)
```

In the above section, the field `fromID` denotes the gene representation
format of your input list. Here is a quick example that runs fast.

``` r
# using an updated version of pcxn 
data("msigdb_c2")
data("targetScan_03")

tempEnrich <-miRNAPathwayEnrichment(targetScan_03[1:20],msigdb_c2[1:20])

head(tempEnrich)
```

    ##                                 x                                       y
    ## 1                hsa-miR-101-3p.1 Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS
    ## 2                    hsa-miR-1193 Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS
    ## 3                    hsa-miR-1197 Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS
    ## 4                hsa-miR-124-3p.1 Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS
    ## 5 hsa-miR-124-3p.2/hsa-miR-506-3p Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS
    ## 6                 hsa-miR-1252-5p Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS
    ##        pval Intersect mirset_Size not_mirset pathway_Size   ratio_in  ratio_out
    ## 1 0.6041595         1          10        695           62 0.01639344 0.01294964
    ## 2 0.3528596         2          14        691           62 0.03333333 0.01736614
    ## 3 1.0000000         0          15        690           62 0.00000000 0.02173913
    ## 4 0.4540258         3          28        677           62 0.05084746 0.03692762
    ## 5 0.8722208         1          22        683           62 0.01639344 0.03074671
    ## 6 1.0000000         0          16        689           62 0.00000000 0.02322206
    ##   ratio_ratios
    ## 1    1.2659381
    ## 2    1.9194444
    ## 3    0.0000000
    ## 4    1.3769492
    ## 5    0.5331772
    ## 6    0.0000000

## Session info

``` r
sessionInfo()
```

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Mojave 10.14.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] PanomiR_0.99.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] igraph_1.2.8     knitr_1.36       magrittr_2.0.1   tidyselect_1.1.1
    ##  [5] R6_2.5.1         rlang_0.4.12     fastmap_1.1.0    fansi_0.5.0     
    ##  [9] stringr_1.4.0    dplyr_1.0.7      tools_4.1.2      parallel_4.1.2  
    ## [13] xfun_0.28        utf8_1.2.2       withr_2.4.3      htmltools_0.5.2 
    ## [17] ellipsis_0.3.2   yaml_2.2.1       digest_0.6.28    tibble_3.1.6    
    ## [21] lifecycle_1.0.1  crayon_1.4.2     purrr_0.3.4      vctrs_0.3.8     
    ## [25] glue_1.5.0       evaluate_0.14    rmarkdown_2.11   limma_3.50.0    
    ## [29] stringi_1.7.5    compiler_4.1.2   pillar_1.6.4     forcats_0.5.1   
    ## [33] generics_0.1.1   pkgconfig_2.0.3
