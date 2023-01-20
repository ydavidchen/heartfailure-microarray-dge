# Differential Candidate Genes
# David Chen

rm(list=ls())
library(ggplot2)
library(reshape2)
library(lme4)
library(lmtest)
source("<PREFIXED_MASKED>/r_utils.R")

PKEY <- "GEO_accession"

GTHEME <- theme_bw() +
  theme(axis.text=element_text(size=17,color="black"), 
        axis.title.x = element_blank(),
        axis.title.y=element_text(size=22,color="black"),
        strip.text=element_text(size=22,color="black"),
        legend.text=element_text(size=17,color="black"), legend.title=element_text(size=22,color="black"),
        legend.position="top",
        panel.spacing=unit(0.5, "lines"))

load(paste0(DIR,"230115_objects.RData"))

## Main figure: differential subtype gene boxplot
pltSubtypes <- as.data.frame(t(expr[c("IPO5","RABEP1","PCYOX1","MTSS1"), ]))
pltSubtypes <- resetIdx(pltSubtypes, PKEY)
pltSubtypes <- merge(pltSubtypes, meta[,c(PKEY,"Subtype","AgeGroup","Sex","RepeatSubject")], by=PKEY)
pltSubtypes <- subset(pltSubtypes, Subtype != "(control)")

## Data Visualization
pltSubtypes <- melt(pltSubtypes, variable.name="Gene", value.name="Expression")

ggplot(pltSubtypes, aes(Subtype, Expression)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.25, width=0.25) +
  scale_color_brewer(palette="Dark2") +
  facet_wrap( ~ Gene, nrow=1) +
  GTHEME
