# Global Expression Analysis
# Created: 11/02/2022; Updated: 01/20/2022

rm(list=ls())
library(ggplot2)
library(lme4)

GCOLSC <- scale_color_brewer(palette="Set1", direction=-1)
GABLINE <- geom_abline(slope=1, intercept=0, linetype="dashed")
GLIMS <- c(-2.82, 2.82)
GBREAKS_AGE <- seq(10, 80, 10)

UTIL_DIR <- "*************PATH MASKED*************"
source(paste0(UTIL_DIR,"utils_gen.R"))

load(paste0(DIR,"230120_objects.RData"))
summary(pca)$importance[c(2,3), 1:2]

# ------------------------------ Statistical Tests ------------------------------
## Univariate un-adjusted:
cTabRegion <- table(meta$HeartFailure, meta$PrincipalGroup)
cTabRegion
prop.table(cTabRegion, 2)
uRegionFit <- fisher.test(cTabRegion)
uRegionFit
uRegionFit$p.value #exact pval

## Multivariate Linear Mixed-effect:
mRegionFit <- glmer(
  factor(PrincipalGroup=="A") ~ isHF + Age + Sex + (1|RepeatSubject), 
  data = meta, 
  family = binomial
)

resRegionFit <- exp(cbind(
  aOR = fixef(mRegionFit),
  confint(mRegionFit, parm="beta_", method="Wald")
))
cbind(resRegionFit, pval=summary(mRegionFit)$coefficient[ , 4])

## Data Visualization
pltPc <- subset(meta, ! is.na(Subtype))
pltPc$`Heart Failure` <- pltPc$Subtype 
pltPc$`Heart Failure`[pltPc$Subtype=="(control)"] <- "No (Healthy)"
pltPc$`Heart Failure`[pltPc$Subtype!="(control)"] <- paste0("Yes (", pltPc$`Heart Failure`[pltPc$Subtype!="(control)"], ")")
table(pltPc$`Heart Failure`) #verify

ggplot(pltPc, aes(Dimension2, Dimension1, color=`Heart Failure`, shape=Sex, size=Age, alpha=Age)) +
  geom_point() +
  scale_size_continuous(breaks=GBREAKS_AGE, range=c(1, 8)) +
  scale_alpha_continuous(breaks=GBREAKS_AGE, range=c(0.2,0.9)) +
  xlim(GLIMS) + ylim(GLIMS) +
  annotate("text", -2, 2, label="A", size=20, color="gray", alpha=0.5) +
  annotate("text", 1.7, -1.7, label="B", size=20, color="gray", alpha=0.5) +
  GABLINE + GCOLSC + GTHEME
