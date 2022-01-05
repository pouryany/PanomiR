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

    ##                                                       TCGA-BC-A10S-01A-22R-A131-07
    ## Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS                                 1.31081432
    ## Pathway.KEGG_CITRATE_CYCLE_TCA_CYCLE                                    2.12764679
    ## Pathway.KEGG_PENTOSE_PHOSPHATE_PATHWAY                                  0.93814543
    ## Pathway.KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS                   1.34561508
    ## Pathway.KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM                           -0.02061692
    ## Pathway.KEGG_GALACTOSE_METABOLISM                                      -0.02148356
    ##                                                       TCGA-BC-4073-01B-02R-A131-07
    ## Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS                                  1.1074424
    ## Pathway.KEGG_CITRATE_CYCLE_TCA_CYCLE                                     1.7690134
    ## Pathway.KEGG_PENTOSE_PHOSPHATE_PATHWAY                                   1.2218070
    ## Pathway.KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS                    1.5658403
    ## Pathway.KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM                             0.2796341
    ## Pathway.KEGG_GALACTOSE_METABOLISM                                        0.4472539

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

    ##                                                                                                                                   logFC
    ## Pathway.REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                                                                           -0.8868271
    ## Pathway.KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION                                                                         -0.6062216
    ## Pathway.BIOCARTA_AKT_PATHWAY                                                                                                 -0.5485422
    ## Pathway.REACTOME_ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING                                                         0.5203503
    ## Pathway.PID_IL5_PATHWAY                                                                                                      -0.5943319
    ## Pathway.REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS -0.9673573
    ##                                                                                                                                 AveExpr
    ## Pathway.REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                                                                            0.2979081
    ## Pathway.KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION                                                                         -2.2091556
    ## Pathway.BIOCARTA_AKT_PATHWAY                                                                                                  0.3064294
    ## Pathway.REACTOME_ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING                                                         1.3834223
    ## Pathway.PID_IL5_PATHWAY                                                                                                       0.4171509
    ## Pathway.REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS  1.2390324
    ##                                                                                                                                       t
    ## Pathway.REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                                                                           -10.200349
    ## Pathway.KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION                                                                          -7.948648
    ## Pathway.BIOCARTA_AKT_PATHWAY                                                                                                  -6.797107
    ## Pathway.REACTOME_ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING                                                          6.356020
    ## Pathway.PID_IL5_PATHWAY                                                                                                       -6.174177
    ## Pathway.REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS  -5.705278
    ##                                                                                                                                   P.Value
    ## Pathway.REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                                                                           2.937157e-06
    ## Pathway.KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION                                                                         2.274107e-05
    ## Pathway.BIOCARTA_AKT_PATHWAY                                                                                                 7.774231e-05
    ## Pathway.REACTOME_ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING                                                        1.295512e-04
    ## Pathway.PID_IL5_PATHWAY                                                                                                      1.610341e-04
    ## Pathway.REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS 2.878368e-04
    ##                                                                                                                                adj.P.Val
    ## Pathway.REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                                                                           0.003583332
    ## Pathway.KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION                                                                         0.013872055
    ## Pathway.BIOCARTA_AKT_PATHWAY                                                                                                 0.031615208
    ## Pathway.REACTOME_ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING                                                        0.039292310
    ## Pathway.PID_IL5_PATHWAY                                                                                                      0.039292310
    ## Pathway.REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS 0.043435623
    ##                                                                                                                                      B
    ## Pathway.REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                                                                           4.8824453
    ## Pathway.KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION                                                                         3.1444207
    ## Pathway.BIOCARTA_AKT_PATHWAY                                                                                                 2.0285870
    ## Pathway.REACTOME_ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING                                                        1.5525264
    ## Pathway.PID_IL5_PATHWAY                                                                                                      1.3478228
    ## Pathway.REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS 0.7963496
    ##                                                                                                                                                         contrast
    ## Pathway.REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING                                                                           shortLetterCodeTP-shortLetterCodeNT
    ## Pathway.KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION                                                                         shortLetterCodeTP-shortLetterCodeNT
    ## Pathway.BIOCARTA_AKT_PATHWAY                                                                                                 shortLetterCodeTP-shortLetterCodeNT
    ## Pathway.REACTOME_ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING                                                        shortLetterCodeTP-shortLetterCodeNT
    ## Pathway.PID_IL5_PATHWAY                                                                                                      shortLetterCodeTP-shortLetterCodeNT
    ## Pathway.REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS shortLetterCodeTP-shortLetterCodeNT

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

    ##                                                      Pathway cluster
    ## 1                         Pathway.KEGG_TRYPTOPHAN_METABOLISM       5
    ## 2                       Pathway.KEGG_BETA_ALANINE_METABOLISM       6
    ## 3 Pathway.KEGG_GLYCOSPHINGOLIPID_BIOSYNTHESIS_GANGLIO_SERIES       7
    ## 4                 Pathway.KEGG_DRUG_METABOLISM_OTHER_ENZYMES       8
    ## 5                               Pathway.KEGG_DNA_REPLICATION       1
    ## 6                                    Pathway.KEGG_PROTEASOME       1

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

    ##                             x cluster_hits aggInv_cover  aggInv_pval
    ## 1                hsa-miR-1179            2    -1.041824 3.866337e-08
    ## 2                hsa-miR-1197            1    -0.319530 7.653476e-06
    ## 3             hsa-miR-1224-5p            1    -2.454211 2.310664e-02
    ## 4    hsa-miR-1-3p/hsa-miR-206            0    -1.463553 2.819640e-02
    ## 5             hsa-miR-1251-5p            1    -2.286925 1.274129e-01
    ## 6 hsa-miR-103a-3p/hsa-miR-107            1    -1.187462 1.316088e-01
    ##     aggInv_fdr
    ## 1 7.732674e-07
    ## 2 7.653476e-05
    ## 3 1.409820e-01
    ## 4 1.409820e-01
    ## 5 4.386960e-01
    ## 6 4.386960e-01

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

    ##                             x                                       y      pval
    ## 1    hsa-miR-1-3p/hsa-miR-206 Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS 0.1901934
    ## 2            hsa-miR-101-3p.1 Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS 0.5697729
    ## 3 hsa-miR-103a-3p/hsa-miR-107 Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS 1.0000000
    ## 4             hsa-miR-1185-5p Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS 0.1629978
    ## 5                hsa-miR-1193 Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS 0.3583626
    ## 6                hsa-miR-1197 Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS 1.0000000
    ##   Intersect mirset_Size not_mirset pathway_Size   ratio_in   ratio_out
    ## 1         2          11        834           62 0.03333333 0.010791367
    ## 2         1          11        834           62 0.01639344 0.011990408
    ## 3         0          10        835           62 0.00000000 0.011976048
    ## 4         2          10        835           62 0.03333333 0.009580838
    ## 5         2          17        828           62 0.03333333 0.018115942
    ## 6         0          18        827           62 0.00000000 0.021765417
    ##   ratio_ratios
    ## 1     3.088889
    ## 2     1.367213
    ## 3     0.000000
    ## 4     3.479167
    ## 5     1.840000
    ## 6     0.000000

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
