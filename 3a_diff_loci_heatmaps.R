# Hierarchical Clustering of LIMMA Differential Loci
# David Chen
# 01/20/2022

rm(list=ls())

UTIL_DIR <- "*************PATH MASKED*************"
source(paste0(UTIL_DIR,"utils_gen.R"))
source(paste0(UTIL_DIR,"utils_hm.R"))

FC_THRESH <- 1.5
NUM_HC <- 2+1 #considers 1+ sub-clusters
LIMMA_HF <- read.csv(paste0(DIR,"results/limma_binaryhf.csv"), row.names=1, stringsAsFactors=FALSE)

load(paste0(DIR,"230120_objects.RData"))

## Refine heatmap annotation objects:
DIFFGENE_ANNOT <- LIMMA_HF[ , "Direction", drop=FALSE]
DIFFGENE_ANNOT$Direction <- gsub("ed", "", DIFFGENE_ANNOT$Direction, ignore.case=FALSE)

ANN_COLOR[["Direction"]] <- c(`Higher in Fail`="black", `Lower in Fail`="lightgray")
ANN_COLOR[["Subtype"]]  <- ANN_COLOR[["Subtype"]][1:2]

SAMP_ANNOT <- meta[ , c("Subtype","HeartFailure","Sex","AgeGroup")]
SAMP_ANNOT$Subtype[SAMP_ANNOT$Subtype=="(control)"] <- NA
rownames(SAMP_ANNOT) <- meta$GEO_accession

quantile(LIMMA_HF$absFc, 0.99)

diffClustTreeHF <- wrapper_hm(
  expr[rownames(LIMMA_HF)[LIMMA_HF$absFc>=FC_THRESH & LIMMA_HF$adj.P.Val<0.05], ], 
  NUM_HC, 
  annot_samps = SAMP_ANNOT,
  annot_gene = DIFFGENE_ANNOT
)

## Extract cluster membership:
resClustHF <- getClusters(diffClustTreeHF, NUM_HC, "GEO_accession")
colnames(resClustHF)[2] <- "LIMMAclusters"

meta <- merge(meta, resClustHF, by="GEO_accession")
table(meta$HeartFailure, meta$LIMMAclusters)

## Test 1: Concordance of Supervised Cluster w/ Status as paired-design hypothesis
cTabLimma <- table(isCluster23=meta$LIMMAclusters>1, Fail=meta$isHF)
cTabLimma
prop.table(cTabLimma, 2)
mcnemar.test(cTabLimma)

## Test 2: Subgroup (i.e. within Clusters 2 & 3) associations as independent hypothesis
metaSub <- subset(meta, LIMMAclusters>1 & !is.na(Subtype) & isHF)
cTabSubtype <- table(Cluster=metaSub$LIMMAclusters, Subtype=metaSub$Subtype)
cTabSubtype
prop.table(cTabSubtype, 2)
fisher.test(cTabSubtype)
