# Binary LIMMA Models
# David Chen
# 01/20/2023

rm(list=ls())
library(limma)
library(EnhancedVolcano)

UTIL_DIR <- "~/repos/davechen/hf_ml/" #TODO: MASK
source(paste0(UTIL_DIR,"utils_gen.R"))

load(paste0(DIR,"230120_objects.RData"))
meta$Status <- ifelse(meta$isHF, "Failed", "Control")
meta$Subtype[meta$Subtype=="(control)"] <- NA
stopifnot(identical(colnames(expr), meta$GEO_accession))

wrapper_limma <- function(design, dat, outcomeName, nonRefGroup,
                          sampBlock=NULL, dupCor=0,
                          fcThresh=1.5, vXlim=c(-2,2), vYlim=c(0,50), 
                          labelSize=5, pointSize=2) {
  #'@description Reusable wrapper to run robust Bonferroni-adjusted Robust LIMMA for maximal consistency
  #'@param fcThresh Absolute fold change threshold. Fold increase = x*fc, fold decrease = x/fc
  #'@param outcomeCol Group field + non-referent group name
  #'@return DataFrame of LIMMA results with EnhancedVolcano plotted
  rawP <- 0.05 / nrow(dat) #Bonferroni
  lfc <- log2(fcThresh)
  
  fit <- lmFit(dat, design, method="robust", block=sampBlock, correlation=dupCor)
  fit <- eBayes(fit, robust=TRUE)
  
  res_lm <- topTable(fit, coef=paste0(outcomeName,nonRefGroup), number=Inf, sort.by="p", adjust.method="bonferroni")
  print(paste("Number of significant:", sum(res_lm$adj.P.Val<0.05)))
  
  res_lm$foldSigned <- 2^res_lm$logFC
  res_lm$absFc[res_lm$logFC>0] <- res_lm$fold[res_lm$logFC>0]
  res_lm$absFc[res_lm$logFC<0] <- 1/res_lm$fold[res_lm$logFC<0]
  
  res_lm$Direction <- NA
  res_lm$Direction[res_lm$adj.P.Val<0.05 & res_lm$logFC > lfc] <- paste("Higher in", nonRefGroup)
  res_lm$Direction[res_lm$adj.P.Val<0.05 & res_lm$logFC < -lfc] <- paste("Lower in", nonRefGroup)
  print( table(res_lm$Direction) )
  
  keyvals <- rep("dimgray", nrow(res_lm))
  keyvals[grepl("Higher",res_lm$Direction)] <- "orange"
  keyvals[grepl("Lower",res_lm$Direction)] <- "skyblue"
  keyvals[is.na(res_lm$Direction)] <- "(n.s.)"
  names(keyvals) <- res_lm$Direction
  
  vPlt <- EnhancedVolcano(
    res_lm, 
    rownames(res_lm),
    "logFC", 
    "P.Value",
    title = NULL,
    subtitle = NULL,
    FCcutoff = lfc, 
    pCutoff = rawP,
    pointSize = pointSize,
    labSize = labelSize,
    legendLabSize = 17,
    xlim = vXlim,
    ylim = vYlim,
    colAlpha = 0.5,
    colCustom = keyvals
  )
  print(vPlt)
  return(res_lm)
}

# -------------------------- Primary Endpoint: HF vs. control --------------------------
design1 <- model.matrix( ~ Status + Age + Sex, data=meta)

dupcor1 <- duplicateCorrelation(expr, design1, block=meta$RepeatSubject)$consensus.correlation
dupcor1 #0.2613018

# png("~/Downloads/Fig1C.png",10,6.5,"in",res=300)
resHf <- wrapper_limma(
  design1, expr, 
  "Status", "Failed", 
  meta$RepeatSubject, dupcor1,
  1.5, c(-1.5,1.5), c(0,55), 
  labelSize =12
)
# dev.off()

# write.table(rownames(expr), file=paste0(DIR, "results/BackgroundSet.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
# write.csv(resHf, file=paste0(DIR,"results/limma_binaryhf.csv"), row.names=TRUE, quote=FALSE)

# -------------------------- Subgroup Analysis: Dilated vs. Ischemic --------------------------
## Restrict samples:
metaSub <- subset(meta, isHF & !is.na(Subtype))

exprSub <- expr[ , colnames(expr) %in% metaSub$GEO_accession]
stopifnot(identical(colnames(exprSub), metaSub$GEO_accession)) #if not, match

## LIMMA model setup:
design2 <- model.matrix( ~ Subtype + Age + Sex, data=metaSub)

dupcor2 <- duplicateCorrelation(exprSub, design2, block=metaSub$RepeatSubject)$consensus.correlation
dupcor2 #0.2633701

# png("~/Downloads/Fig2B.png",8,5,"in",res=300)
resSub <- wrapper_limma(
  design2, exprSub,
  "Subtype", "Ischemic", 
  metaSub$RepeatSubject, dupcor2,
  1.25, c(-0.7,0.7), c(0,7.75), 
  labelSize = 5,
  pointSize = c(rep(3,6), rep(1, 5553-6))
)
# dev.off()

# write.csv(resSub, file=paste0(DIR,"results/limma_subtypes.csv"), row.names=TRUE, quote=FALSE)
# save(list=c("design1","design2","dupcor1","dupcor2","resHf","resSub"), file=paste0(DIR,"results/limmaDgeIO.RData"), compress=TRUE)
