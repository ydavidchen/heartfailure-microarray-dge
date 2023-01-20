# Unsupervised Hierarchical Clustering
# 01/19/2022

rm(list=ls())
library(matrixStats)
library(pheatmap)
source("<PREFIXED_MASKED>/r_utils.R")

CL_PARAMS <- c("euclidean", "ward.D")

ANN_COLOR <- list(
  AgeGroup = c(Over50="black", Under50="lightgray"),
  Sex = c(Male="black", Female="lightgray"),
  HeartFailure = c(Yes="black", No="lightgray"),
  Subtype = c(Ischemic="black", Dilated="lightgray"),
  PrincipalGroup = c(A="black", B="lightgray"),
  HasRepeat  = c(Yes="black", No="lightgray")
)

load(paste0(DIR,"230115_objects.RData"))

## Reusable sample annotations for heatmaps:
samp_annot <- meta[ , c("HeartFailure","Subtype","Sex","AgeGroup","PrincipalGroup")]
samp_annot$Subtype[samp_annot$Subtype=="(control)"] <- NA #blank positions in tracking bar
rownames(samp_annot) <- meta$GEO_accession

wrapper_hm <- function(mat, plot=TRUE) {
  #'@description Wrapper to reduce code duplication & increase consistency
  #'@details Uses objects from .GlobalEnv directly
  #'@return hclust tree object
  phObj <- pheatmap(
    mat,
    cutree_cols = 2,
    clustering_distance_cols = CL_PARAMS[1],
    clustering_distance_rows = CL_PARAMS[1],
    clustering_method = CL_PARAMS[2], 
    annotation_col = samp_annot,
    annotation_colors = ANN_COLOR,
    color = HEAT_COLS,
    fontsize = 12,
    show_rownames = FALSE,
    show_colnames = FALSE,
    border_color = NA
  )
  return(phObj$tree_col)
}

# ------------------- Unsupervised Hierarchical Clustering -------------------
## Inter-dataset variance distribution:
sampSds <- rowSds(expr)
plot(sort(sampSds), main="Inter-sample SD")
summary(sampSds)
sum(sampSds >= 0.9905226)

binaryHcTree <- wrapper_hm(expr[sampSds >= 0.9905226, ], FALSE)

## Test hierarchical cluster membership:
resClust <- getClusters(binaryHcTree, 2, "GEO_accession")
meta <- merge(meta, resClust, by="GEO_accession")
meta$Cluster <- paste0("Cluster", meta$Cluster)

cTabConcord <- table(PCA=meta$PrincipalGroup, HC=meta$Cluster) 
cTabConcord <- cTabConcord[ , c(2,1)]
cTabConcord
prop.table(cTabConcord, 2)
mcnemar.test(cTabConcord)

caret::confusionMatrix(
  factor(meta$PrincipalGroup == "A"), 
  factor(meta$Cluster == "Cluster2")
) #confirms above

# ------------------- Visualize Loci Identified by LIMMA Supervised Analysis -------------------
LIMMA_HF <- read.csv(paste0(DIR,"results/limma_binaryhf.csv"), row.names=1, stringsAsFactors=FALSE)

quantile(LIMMA_HF$absFc, 0.99)
diffClustTreeHF <- wrapper_hm(expr[rownames(LIMMA_HF)[LIMMA_HF$absFc>=2], ], FALSE)
resClustHF <- getClusters(diffClustTreeHF, 2, "GEO_accession")
colnames(resClustHF)[2] <- "SupervisedHF"
meta <- merge(meta, resClustHF, by="GEO_accession")

cTabSup <- table(Cluster=meta$SupervisedHF, Fail=meta$HeartFailure)
cTabSup <- cTabSup[c(2,1), ]
cTabSup
prop.table(cTabSup, 2)
mcnemar.test(cTabSup)
caret::confusionMatrix(
  factor(meta$SupervisedHF==2),
  factor(meta$isHF)
) #verify pval

resClustSub <- wrapper_hm(expr[c("IPO5","RABEP1","PCYOX1","MTSS1"), ], FALSE)
