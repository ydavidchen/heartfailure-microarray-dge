# Binary LIMMA Models
# David Chen

rm(list=ls())
library(limma)
library(EnhancedVolcano)
source("<PREFIXED_MASKED>/r_utils.R")

wrapper_limma <- function(design, dat, outcomeName, nonRefGroup,
                          sampBlock=NULL, dupCor=0,
                          fcThresh=1.5, vXlim=c(-2,2), vYlim=c(0,50)) {
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
  keyvals[grepl("Lower",res_lm$Direction)] <- "navy"
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
    # hline = rawP,
    # cutoffLineType = "blank",
    xlim = vXlim,
    ylim = vYlim,
    colAlpha = 0.5,
    colCustom = keyvals
  )
  print(vPlt)
  return(res_lm)
}

# -------------------------- Primary Endpoint: HF vs. control --------------------------
load(paste0(DIR,"230115_objects.RData"))
# write.table(rownames(expr), file=paste0(DIR, "results/BackgroundSet.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)

meta$Status <- ifelse(meta$isHF, "Failed", "Control")
meta$Subtype[meta$Subtype=="(control)"] <- "Control"
stopifnot(identical(colnames(expr), meta$GEO_accession))

design1 <- model.matrix( ~ Status + Age + Sex, data=meta)

dupcor1 <- duplicateCorrelation(expr, design1, block=meta$RepeatSubject)$consensus.correlation
dupcor1 #0.262561

resHf <- wrapper_limma(
  design1, expr, 
  "Status", "Failed", 
  meta$RepeatSubject, dupcor1,
  1.5, c(-1.5,1.5), c(0,55)
)

# write.csv(resHf, file=paste0(DIR,"results/limma_binaryhf.csv"), row.names=TRUE, quote=FALSE)

# -------------------------- Subgroup Analysis: Dilated vs. Ischemic --------------------------
metaSub <- subset(meta, isHF)
exprSub <- expr[ , colnames(expr) %in% metaSub$GEO_accession]
stopifnot(identical(colnames(exprSub), metaSub$GEO_accession)) #if not, match

design2 <- model.matrix( ~ Subtype + Age + Sex, data=metaSub)

dupcor2 <- duplicateCorrelation(exprSub, design2, block=metaSub$RepeatSubject)$consensus.correlation
dupcor2 #0.2136145

resSub <- wrapper_limma(
  design2, exprSub,
  "Subtype", "Ischemic", 
  metaSub$RepeatSubject, dupcor2,
  1.25, c(-0.75,0.75), c(0,7)
)

# write.csv(resSub, file=paste0(DIR,"results/limma_subtypes.csv"), row.names=TRUE, quote=FALSE)
# save(list=c("design1","design2","dupcor1","dupcor2","resHf","resSub"), file=paste0(DIR,"results/limmaDgeIO.RData"), compress=TRUE)
