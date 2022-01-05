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

PanomiR can be accessed via Bioconductor. To install, start R and run
the following code.

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

## Overview

PanomiR is a pipeline to prioritize disease-associated miRNAs based on
activity of disease-associated pathways. The input datasets for PanomiR
are (a) a gene expression disease dataset along with covariates, (b) a
background collection of pathways/genesets, and (c) a collection of
miRNAs containing gene targets.

The general workflow of PanomiR is (a) generation of pathway summary
statistics from gene expression data, (b) detection of differentially
activated pathways, (c) finding coherent groups, or clusters, of
differentially activated pathways, and (d) detecting miRNAs targeting
each group of pathways.

Individual steps of the workflow can be used in isolation to carry out
different analyses. The following sections outline each step and
material needed to execute PanomiR.

## 1. Pathway summarization

PanomiR can generate pathway activity profiles given a gene expression
dataset and a list of pathways.Pathway summaries are numbers that
represent the overall activity of genes that belong to each pathway.
These numbers are calculated based on a methodology previously described
in Altschuler et al. 2013 and Joachim and Altschuler et al. 2018 (2013;
2018). Briefly, genes in each sample are ranked by their expression
values and then pathway summaries are calculated as the average
rank-squared of genes within a pathway. The summaries are then center
and scaled (zNormalized) across samples.

The default list of background pathways in PanomiR is formatted into a
table (`data("path_gene_table")`). The table is based on canonical
pathways collection of Molecular Signatures Database (MSigDB) V6.2 and
it contains annotated pathways from a variety of sources excluding KEGG
dataset (2011).

\*\* Users interested in using other pathway/geneset backgrounds, such
as newer versions of MSigDB or KEGG, should refer to the
[appendix](#geneset) of this manual.

This section uses a reduced example dataset from The Cancer Genome Atlas
(TCGA) Liver Hepatocellular Carcinoma (LIHC) dataset to generate Pathway
summary statistics (2017).

``` r
library(PanomiR)

# Pathway reference from the PanomiR package
data("path_gene_table")
data("miniTestsPanomiR")
# Generating pathway summary statistics 

summaries <- pathwaySummary(miniTestsPanomiR$mini_LIHC_Exp,
                            path_gene_table, method = "x2",
                            zNormalize = TRUE)

head(summaries)[,1:2]
```

    ##                         TCGA-BC-A10S-01A-22R-A131-07
    ## BIOCARTA_41BB_PATHWAY                     -0.1506216
    ## BIOCARTA_ACE2_PATHWAY                     -0.5676447
    ## BIOCARTA_ACH_PATHWAY                      -0.3211747
    ## BIOCARTA_ACTINY_PATHWAY                    1.4363526
    ## BIOCARTA_AGPCR_PATHWAY                    -0.1948523
    ## BIOCARTA_AGR_PATHWAY                       0.6802993
    ##                         TCGA-BC-4073-01B-02R-A131-07
    ## BIOCARTA_41BB_PATHWAY                     -0.1269436
    ## BIOCARTA_ACE2_PATHWAY                     -0.8327436
    ## BIOCARTA_ACH_PATHWAY                      -0.4390042
    ## BIOCARTA_ACTINY_PATHWAY                    1.4975456
    ## BIOCARTA_AGPCR_PATHWAY                    -0.2499193
    ## BIOCARTA_AGR_PATHWAY                       0.5420588

## 2. Differential Pathway activation

Once you generate the pathway activity profiles, as discussed in the
last section, there are several analysis that you can perform. We have
bundled some of the most important ones into standalone functions. Here,
we describe differential pathway activation profiling, which is
examining differences in pathway activity profiles in user-determined
conditions.

At this stage you need to provide a pathway-gene association table, an
expression dataset, and a covariates table. You need to specity what
covariates you would like to contrast. You also need to provide a
contrast, as formatted in limma. If the contrast is not provided, the
function assumes the first two levels of the provided contrast
covariate. **Note:** make sure the contrast covariate is formatted as
factor.

``` r
output0 <- differentialPathwayAnalysis(
                        geneCounts = miniTestsPanomiR$mini_LIHC_Exp,
                        pathways =  path_gene_table,
                        covariates = miniTestsPanomiR$mini_LIHC_Cov,
                        condition = 'shortLetterCode')

de.paths <- output0$DEP

head(de.paths)
```

    ##                                                                                                                           logFC
    ## REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                                                                           -0.9159376
    ## BIOCARTA_AKT_PATHWAY                                                                                                 -0.5744103
    ## PID_IL5_PATHWAY                                                                                                      -0.6219876
    ## REACTOME_ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING                                                         0.4899259
    ## REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS -1.0070977
    ## BIOCARTA_PTC1_PATHWAY                                                                                                 0.6625938
    ##                                                                                                                         AveExpr
    ## REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                                                                            0.3044281
    ## BIOCARTA_AKT_PATHWAY                                                                                                  0.3123897
    ## PID_IL5_PATHWAY                                                                                                       0.4240432
    ## REACTOME_ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING                                                         1.3942628
    ## REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS  1.2536069
    ## BIOCARTA_PTC1_PATHWAY                                                                                                -1.5149008
    ##                                                                                                                               t
    ## REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                                                                           -10.404966
    ## BIOCARTA_AKT_PATHWAY                                                                                                  -6.770069
    ## PID_IL5_PATHWAY                                                                                                       -6.255756
    ## REACTOME_ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING                                                          6.096421
    ## REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS  -5.905357
    ## BIOCARTA_PTC1_PATHWAY                                                                                                  5.717801
    ##                                                                                                                           P.Value
    ## REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                                                                           1.942463e-06
    ## BIOCARTA_AKT_PATHWAY                                                                                                 6.903010e-05
    ## PID_IL5_PATHWAY                                                                                                      1.276971e-04
    ## REACTOME_ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING                                                        1.555564e-04
    ## REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS 1.979562e-04
    ## BIOCARTA_PTC1_PATHWAY                                                                                                2.519967e-04
    ##                                                                                                                        adj.P.Val
    ## REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                                                                           0.002012391
    ## BIOCARTA_AKT_PATHWAY                                                                                                 0.035757593
    ## PID_IL5_PATHWAY                                                                                                      0.040289104
    ## REACTOME_ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING                                                        0.040289104
    ## REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS 0.041016532
    ## BIOCARTA_PTC1_PATHWAY                                                                                                0.041612753
    ##                                                                                                                             B
    ## REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                                                                           5.240095
    ## BIOCARTA_AKT_PATHWAY                                                                                                 2.126311
    ## PID_IL5_PATHWAY                                                                                                      1.550780
    ## REACTOME_ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING                                                        1.364364
    ## REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS 1.135618
    ## BIOCARTA_PTC1_PATHWAY                                                                                                0.905467
    ##                                                                                                                                                 contrast
    ## REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                                                                           shortLetterCodeTP-shortLetterCodeNT
    ## BIOCARTA_AKT_PATHWAY                                                                                                 shortLetterCodeTP-shortLetterCodeNT
    ## PID_IL5_PATHWAY                                                                                                      shortLetterCodeTP-shortLetterCodeNT
    ## REACTOME_ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING                                                        shortLetterCodeTP-shortLetterCodeNT
    ## REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS shortLetterCodeTP-shortLetterCodeNT
    ## BIOCARTA_PTC1_PATHWAY                                                                                                shortLetterCodeTP-shortLetterCodeNT

## 3. Finding clusters of pathways

PanomiR provides a function to find groups coordinated differentially
activated pathways based on a pathway co-expression network (PCxN)
previously described in Pita-Juárez and Altschuler et al (2018).
Briefly, PCxN is a network where nodes are pathways and links are
co-expression between the nodes. It is formatted into a table were rows
represent edges. The edges of PCxN are marked by two numbers, 1- a
correlation co-efficient and 2- a significance adjusted p-value.
Cut-offs for both of these numbers can be manually set using PanomiR
functions. See function manuals for more info.

PCxN and its associated genesets are already released and can be
accessed through [pcxn Bioconductor
package](http://bioconductor.org/packages/release/bioc/html/pcxn.html)
and [pcxnData
package](http://bioconductor.org/packages/release/data/experiment/html/pcxnData.html).

Here we have provided a small version of PCxN for tutorial purposes. A
more recent version of PCxN based on MSigDB V6.2 is available through
the data repository accompanying PanomiR manuscript, which can be found
[here](https://github.com/pouryany/PanomiR_paper).

``` r
# using an updated version of pcxn 

set.seed(2)
pathwayClustsLIHC <- mappingPathwaysClusters(
                            pcxn = miniTestsPanomiR$miniPCXN, 
                            dePathways = de.paths[1:300,],
                            topPathways = 200,
                            outDir=".",
                            plot = FALSE,
                            subplot = FALSE,
                            prefix='',
                            clusteringFunction = "cluster_louvain",
                            correlationCutOff = 0.1)


head(pathwayClustsLIHC$Clustering)
```

    ##                      Pathway cluster
    ## 1       BIOCARTA_NO1_PATHWAY       1
    ## 2       BIOCARTA_AKT_PATHWAY       1
    ## 3       BIOCARTA_ALK_PATHWAY       1
    ## 4     BIOCARTA_RANKL_PATHWAY       1
    ## 5       BIOCARTA_MCM_PATHWAY       3
    ## 6 BIOCARTA_CELLCYCLE_PATHWAY       3

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
set.seed(1)
output2 <- prioritizeMicroRNA(enriches0 = miniTestsPanomiR$miniEnrich,
                    pathClust = miniTestsPanomiR$miniPathClusts$Clustering,
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

    ##                                 x cluster_hits aggInv_cover  aggInv_pval
    ## 1                hsa-miR-101-3p.2            6   -1.9566603 0.0001216703
    ## 2                hsa-miR-101-3p.1            4   -0.3395771 0.0006214715
    ## 3 hsa-miR-124-3p.2/hsa-miR-506-3p            7   -0.2357761 0.0008599272
    ## 4                 hsa-miR-1247-5p            4   -1.6599230 0.0021625662
    ## 5                 hsa-miR-1249-3p            1   -2.4578993 0.0042061415
    ## 6                 hsa-miR-1252-5p            4   -0.7572036 0.0050836835
    ##    aggInv_fdr
    ## 1 0.002433406
    ## 2 0.005732848
    ## 3 0.005732848
    ## 4 0.010812831
    ## 5 0.016824566
    ## 6 0.016945612

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

tempEnrich <-miRNAPathwayEnrichment(targetScan_03[1:30],msigdb_c2[1:30])

head(tempEnrich)
```

    ##                                 x                     y      pval Intersect
    ## 1     hsa-miR-103a-3p/hsa-miR-107 BIOCARTA_RELA_PATHWAY 0.4476283         1
    ## 2                hsa-miR-124-3p.1 BIOCARTA_RELA_PATHWAY 0.1342510         2
    ## 3 hsa-miR-124-3p.2/hsa-miR-506-3p BIOCARTA_RELA_PATHWAY 0.4476283         1
    ## 4                 hsa-miR-1252-5p BIOCARTA_RELA_PATHWAY 1.0000000         0
    ## 5   hsa-miR-1271-5p/hsa-miR-96-5p BIOCARTA_RELA_PATHWAY 1.0000000         0
    ## 6     hsa-miR-103a-3p/hsa-miR-107  BIOCARTA_NO1_PATHWAY 0.7143586         1
    ##   mirset_Size not_mirset pathway_Size   ratio_in  ratio_out ratio_ratios
    ## 1          14        378           16 0.06666667 0.03439153    1.9384615
    ## 2          16        376           16 0.14285714 0.03723404    3.8367347
    ## 3          14        378           16 0.06666667 0.03439153    1.9384615
    ## 4          10        382           16 0.00000000 0.02617801    0.0000000
    ## 5          17        375           16 0.00000000 0.04533333    0.0000000
    ## 6          14        378           33 0.03125000 0.03439153    0.9086538

## Session info

``` r
sessionInfo()
```

    ## R Under development (unstable) (2021-12-03 r81290)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Mojave 10.14.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] PanomiR_0.99.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] igraph_1.2.10    knitr_1.37       magrittr_2.0.1   tidyselect_1.1.1
    ##  [5] R6_2.5.1         rlang_0.4.12     fastmap_1.1.0    fansi_0.5.0     
    ##  [9] stringr_1.4.0    dplyr_1.0.7      tools_4.2.0      parallel_4.2.0  
    ## [13] xfun_0.29        utf8_1.2.2       DBI_1.1.2        withr_2.4.3     
    ## [17] htmltools_0.5.2  ellipsis_0.3.2   assertthat_0.2.1 yaml_2.2.1      
    ## [21] digest_0.6.29    tibble_3.1.6     lifecycle_1.0.1  crayon_1.4.2    
    ## [25] purrr_0.3.4      vctrs_0.3.8      glue_1.6.0       evaluate_0.14   
    ## [29] rmarkdown_2.11   limma_3.51.2     stringi_1.7.6    compiler_4.2.0  
    ## [33] pillar_1.6.4     forcats_0.5.1    generics_0.1.1   pkgconfig_2.0.3

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-ally2017comprehensive" class="csl-entry">

Ally, Adrian, Miruna Balasundaram, Rebecca Carlsen, Eric Chuah, Amanda
Clarke, Noreen Dhalla, Robert A Holt, et al. 2017. “Comprehensive and
Integrative Genomic Characterization of Hepatocellular Carcinoma.”
*Cell* 169 (7): 1327–41.

</div>

<div id="ref-altschuler2013pathprinting" class="csl-entry">

Altschuler, Gabriel M, Oliver Hofmann, Irina Kalatskaya, Rebecca Payne,
Shannan J Ho Sui, Uma Saxena, Andrei V Krivtsov, et al. 2013.
“Pathprinting: An Integrative Approach to Understand the Functional
Basis of Disease.” *Genome Medicine* 5 (7): 1–13.

</div>

<div id="ref-joachim2018relative" class="csl-entry">

Joachim, Rose B, Gabriel M Altschuler, John N Hutchinson, Hector R Wong,
Winston A Hide, and Lester Kobzik. 2018. “The Relative Resistance of
Children to Sepsis Mortality: From Pathways to Drug Candidates.”
*Molecular Systems Biology* 14 (5): e7998.

</div>

<div id="ref-liberzon2011molecular" class="csl-entry">

Liberzon, Arthur, Aravind Subramanian, Reid Pinchback, Helga
Thorvaldsdóttir, Pablo Tamayo, and Jill P Mesirov. 2011. “Molecular
Signatures Database (MSigDB) 3.0.” *Bioinformatics* 27 (12): 1739–40.

</div>

<div id="ref-pita2018pathway" class="csl-entry">

Pita-Juárez, Yered, Gabriel Altschuler, Sokratis Kariotis, Wenbin Wei,
Katjuša Koler, Claire Green, Rudolph E Tanzi, and Winston Hide. 2018.
“The Pathway Coexpression Network: Revealing Pathway Relationships.”
*PLoS Computational Biology* 14 (3): e1006042.

</div>

</div>
