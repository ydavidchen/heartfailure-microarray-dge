# Table One Summary
# David Chen

rm(list=ls())
library(tableone)
source("<PREFIXED_MASKED>/r_utils.R")

load(paste0(DIR,"230115_objects.RData"))
summary(pca)$importance[c(2,3), 1:2]

t1 <- CreateTableOne(strata="Subtype", vars=c("Age","Sex","HasRepeat"), data=meta, includeNA=TRUE, test=FALSE)
print(t1, showAllLevels=TRUE)

t1_alt <- CreateTableOne(strata="HeartFailure", vars=c("Age","Sex","HasRepeat"), data=meta, test=FALSE)
print(t1_alt, showAllLevels=TRUE)

ts2 <- CreateTableOne(strata="Dataset", vars=c("Age","Sex","HasRepeat"), data=meta, test=FALSE)
print(ts2, showAllLevels=TRUE)
