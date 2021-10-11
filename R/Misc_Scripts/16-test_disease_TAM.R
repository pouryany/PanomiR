library(dplyr)
TAM.results1 <- read.delim("~/Downloads/TAM_cluster1.txt")

TAM.results1 <- TAM.results1 %>%
                    filter(.,Category == "Disease") %>%
                    arrange(.,FDR,desc(Count))


TAM.results2 <- read.delim("~/Downloads/TAM_cluster2.txt")

TAM.results2 <- TAM.results2 %>%
    filter(.,Category == "Disease") %>%
    arrange(.,FDR,desc(Count))




TAM.results3 <- read.delim("~/Downloads/TAM_cluster3.txt")

TAM.results3 <- TAM.results3 %>%
    filter(.,Category == "Disease") %>%
    arrange(.,FDR,desc(Count))




temp <- read.csv("../test_cases/LIHC/Output/cluster_louvain_Prioritization_Tarbase/x2_LIHCGene_Tarbase1000_samples_clustNo_1.csv")

temp <- temp$x
set.seed(1)
sample(temp,size = 25,replace = F)
