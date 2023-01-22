# Unsupervised Hierarchical Clustering
# David Chen
# 01/20/2022

rm(list=ls())
library(matrixStats)

UTIL_DIR <- "~/repos/davechen/hf_ml/" #TODO: MASK
source(paste0(UTIL_DIR,"utils_gen.R"))
source(paste0(UTIL_DIR,"utils_hm.R"))

load(paste0(DIR,"230120_objects.RData"))

## Reusable sample annotations for heatmaps:
SAMP_ANNOT <- meta[ , c("HeartFailure","Subtype","Sex","AgeGroup","PrincipalGroup")]
rownames(SAMP_ANNOT) <- meta$GEO_accession

## Inter-dataset variance distribution:
sampSds <- rowSds(expr)
plot(sort(sampSds), main="Inter-sample SD")
summary(sampSds)
sum(sampSds >= 0.98998)

binaryHcTree <- wrapper_hm(expr[sampSds >= 0.98998,], 2, SAMP_ANNOT)

## Extract hierarchical clusters:
resClust <- getClusters(binaryHcTree, 2, "GEO_accession")
meta <- merge(meta, resClust, by="GEO_accession")
meta$Cluster <- paste0("Cluster", meta$Cluster)

## Test overlap of Expression Cluster 1 w/ PCGroup B:
cTabConcord <- table(
  isClust1 = meta$Cluster=="Cluster1",
  isGroupB = meta$PrincipalGroup=="B"
) 
cTabConcord
round( 100*prop.table(cTabConcord, 2), 1 )
mcnemar.test(cTabConcord)
