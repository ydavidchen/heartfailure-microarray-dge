# Global Expression Analysis
# Created: 11/02/2022; Updated: 01/16/2022

rm(list=ls())
library(ggplot2)
library(lme4)
source("<PREFIXED_MASKED>/r_utils.R")

GCOLSC <- scale_color_brewer(palette="Set1")
GABLINE <- geom_abline(slope=1, intercept=0, linetype="dashed")
GLIMS <- c(-2.82, 2.82)
GBREAKS_AGE <- seq(10, 80, 10)

load(paste0(DIR,"230115_objects.RData"))
summary(pca)$importance[c(2,3), 1:2]

# ------------------------------ Data Visualization ------------------------------
qplot(expr, xlab="Value", ylab="Freq", xlim=c(-7,7)) + GTHEME

## 0. SUPPLEMENTAL FIG: Quality of batch-effect & processing
ggplot(meta, aes(Dimension2, Dimension1, color=Dataset)) +
  geom_point(size=5, alpha=0.5) +
  GTHEME

## 1. Univariate test of global expression vs. disease status:
ggplot(meta, aes(Dimension2, Dimension1, color=HeartFailure, shape=Subtype)) +
  geom_point(size=5, alpha=0.5) +
  xlim(GLIMS) + ylim(GLIMS) +
  annotate("text", -2, 2, label="A", size=16, color="gray", alpha=0.75) +
  annotate("text", 1.7, -1.7, label="B", size=16, color="gray", alpha=0.75) +
  GABLINE + GCOLSC + GTHEME

## 2. Principal Group vs. Age Group:
ggplot(meta, aes(Dimension2, Dimension1, color=AgeGroup)) +
  geom_point(size=2.5, alpha=0.75) +
  scale_color_brewer(palette="Accent") +
  xlim(GLIMS) + ylim(GLIMS) +
  annotate("text", -2, 2, label="A", size=16, color="gray", alpha=0.75) +
  annotate("text", 1.7, -1.7, label="B", size=16, color="gray", alpha=0.75) +
  GABLINE + GTHEME

## 3. Principal Group vs. Sex:
ggplot(meta, aes(Dimension2, Dimension1, color=Sex)) +
  geom_point(size=2.5, alpha=0.5) +
  scale_color_brewer(palette="Dark2") +
  xlim(GLIMS) + ylim(GLIMS) +
  annotate("text", -2, 2, label="A", size=16, color="gray", alpha=0.75) +
  annotate("text", 1.7, -1.7, label="B", size=16, color="gray", alpha=0.75) +
  GABLINE + GTHEME

## MAIN FIGURE: Mapping Age & Stratifying on Sex
ggplot(meta, aes(Dimension2, Dimension1, shape=Subtype, color=HeartFailure, size=Age, alpha=Age)) +
  geom_point() +
  scale_size_continuous(breaks=GBREAKS_AGE, range=c(0.6, 4.8)) +
  scale_alpha_continuous(breaks=GBREAKS_AGE, range=c(0.2,0.9)) +
  facet_grid(~ Sex) +
  xlim(GLIMS) + ylim(GLIMS) +
  annotate("text", -2, 2, label="A", size=17, color="dimgray", alpha=0.75) +
  annotate("text", 1.7, -1.7, label="B", size=17, color="dimgray", alpha=0.75) +
  GABLINE + GCOLSC + GTHEME

# ------------------------------ Statistical Tests ------------------------------
## Univariate un-adjusted:
cTabRegion <- table(meta$PrincipalGroup, meta$HeartFailure)
cTabRegion
prop.table(cTabRegion)
uRegionFit <- fisher.test(cTabRegion)
uRegionFit

c(uRegionFit$estimate, cilb=uRegionFit$conf.int, pval=uRegionFit$p.value) #exact pval

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
