# Summary of Study Population
# Last Update: 01/20/2023
# Copyright (C) 2023-4 Y. David Chen

rm(list=ls())
library(tableone)
VARS <- c("Age","Sex","HasRepeat","Subtype")

UTIL_DIR <- "*************PATH MASKED*************"
source(paste0(UTIL_DIR,"utils_gen.R"))

load(paste0(DIR,"230120_objects.RData"))
summary(pca)$importance[c(2,3), 1:2]

t1 <- CreateTableOne(strata="HeartFailure", vars=VARS, data=meta, test=FALSE)
print(t1, showAllLevels=TRUE)

supp_t2 <- CreateTableOne(strata="Dataset", vars=VARS, data=meta, test=FALSE)
print(supp_t2, showAllLevels=TRUE)
