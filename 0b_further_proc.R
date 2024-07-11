# Further Processing of Joined Datasets for Downstream Analyses
# Date created: 10/27/2022; Updated 01/05/2022
# Copyright (C) 2023-4 Y. David Chen

rm(list=ls())
library(matrixStats)
library(sva)

DIR <- "*************PATH MASKED*************"
STRNA <- c("", "nan", "NaN", "NAN", "N/A", "n/a", "---")

## Load Python3-processed data:
meta <- read.csv(paste0(DIR,"microarray_samples.csv"), stringsAsFactors=FALSE, na.strings=STRNA)

meta$isHF <- as.logical(toupper(meta$isHF))
meta$HeartFailure <- ifelse(meta$isHF, "Yes", "No")
meta$HeartFailure <- factor(meta$HeartFailure, levels=c("Yes","No"))

meta$Subtype <- gsub("ICM", "Ischemic", meta$Subtype, ignore.case=FALSE)
meta$Subtype <- gsub("DCM", "Dilated", meta$Subtype, ignore.case=FALSE)
table(meta$Subtype)

meta$AgeGroup <- ifelse(meta$Age > 50, "Over50", "Under50")

## Label repeated measurements, if any:
meta$HasRepeat <- ! is.na(meta$RepeatSubject)
table(meta$Subtype, meta$HasRepeat)
meta$RepeatSubject[! meta$HasRepeat] <- meta$GEO_accession[! meta$HasRepeat]
sum(duplicated(meta$RepeatSubject) | duplicated(meta$RepeatSubject, fromLast=TRUE))

## Expression data:
expr <- read.csv(paste0(DIR,"merged_microarrays_for_dge.csv"), row.names="ID_REF", na.strings=STRNA)
stopifnot(! anyDuplicated(rownames(expr)))
stopifnot(identical(meta$GEO_accession, colnames(expr)))

## ComBat:
MODCOMBAT <- model.matrix( ~ 1, data=meta)
combatExpr <- ComBat(dat=expr, batch=meta$Dataset, mod=MODCOMBAT, par.prior=TRUE)

## Additional processing:
expr <- combatExpr
expr <- scale(expr)
quantile(expr, c(0.000001, 0.999999))
range(expr)

## Execute PCA & normalize PC1-2:
pca <- prcomp(t(expr), center=FALSE, scale.=FALSE) #centering & var-scaling already performed manually
summary(pca)$importance[c(2,3), 1:7]

resPca <- as.data.frame(scale(pca$x[ , c(1,2)]))
resPca$GEO_accession <- rownames(resPca)
resPca$PrincipalGroup <- ifelse(resPca$PC1 > resPca$PC2, "A", "B") #surrogate variable of global expression
colnames(resPca) <- gsub("PC", "Dimension", colnames(resPca), ignore.case=FALSE)
meta <- merge(meta, resPca, by="GEO_accession")

## Data visualization: 
ggplot(meta, aes(Dimension2, Dimension1, color=Dataset)) +
  geom_point(size=5, alpha=0.5) +
  GTHEME

qplot(expr, xlab="Value", ylab="Freq", xlim=c(-7,7)) + GTHEME

## Export:
save(list=c("expr","meta", "pca"), file=paste0(DIR,"230120_objects.RData"), compress=TRUE) 
write.csv(expr, paste0(DIR,"230120_expression.csv"), row.names=TRUE, quote=FALSE)
write.csv(meta, paste0(DIR,"230120_clin.csv"), row.names=FALSE, quote=FALSE)
