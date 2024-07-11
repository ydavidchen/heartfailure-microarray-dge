# Differential Candidate Genes
# Last Update: 01/20/2022
# Copyright (C) 2023-4 Y. David Chen

rm(list=ls())
library(ggplot2)
library(reshape2)
library(lme4)
library(lmtest)

UTIL_DIR <- "*************PATH MASKED*************"
source(paste0(UTIL_DIR,"utils_gen.R"))
PKEY <- "GEO_accession"

DIFF_GENES <- c("EGR1","MTSS1","PCYOX1",
                "CLTB","ITGA7","IPO5")

GTHEME <- theme_bw() +
  theme(axis.text.x=element_text(size=20,color="black"), axis.text.y=element_text(size=15,color="black"),
        axis.title.x=element_blank(), axis.title.y=element_text(size=20,color="black"),
        strip.text=element_text(size=20,color="black",face="italic"),
        panel.spacing=unit(0.5, "lines"))

load(paste0(DIR,"230120_objects.RData"))

## Main figure: differential subtype gene boxplot
pltSubtypes <- as.data.frame(t(expr[DIFF_GENES, ]))
pltSubtypes <- resetIdx(pltSubtypes, PKEY)
pltSubtypes <- merge(pltSubtypes, meta[,c(PKEY,"Subtype","AgeGroup","Sex","RepeatSubject")], by=PKEY)
pltSubtypes <- subset(pltSubtypes, Subtype != "(control)" & !is.na(Subtype))

## Data Visualization
pltSubtypes <- melt(pltSubtypes, variable.name="Gene", value.name="Expression")

# png("~/Downloads/Fig2C.png", 10, 5, "in", res=300)
ggplot(pltSubtypes, aes(Subtype, Expression)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.25, width=0.25) +
  scale_color_brewer(palette="Dark2") +
  facet_wrap( ~ Gene, nrow=2, scales="free_y") +
  GTHEME
# dev.off()
