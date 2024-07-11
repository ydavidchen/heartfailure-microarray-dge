# General Utility Module
# Copyright (C) 2023-4 Y. David Chen

library(ggplot2)

## Constants: 
DIR <- "*************PATH MASKED*************"

GTHEME <- theme_bw() +
  theme(axis.text=element_text(size=15,color="black"), axis.title=element_text(size=20,color="black"),
        strip.text=element_text(size=20,color="black"),
        legend.text=element_text(size=15,color="black"), legend.title=element_text(size=20,color="black"),
        legend.position="right", legend.box="vertical", legend.title.align=0,
        panel.spacing=unit(3, "lines"))

## Helpers:
resetIdx <- function(df, idxname="Gene") {
  #'@description Motivated by python Pandas' <DataFrame>.reset_index method
  df[idxname] <- rownames(df)
  rownames(df) <- NULL
  return(df[ , c(ncol(df),2:ncol(df)-1)])
}

getClusters <- function(hcObj, num_cl=3, rowname=NULL) {
  #'@param hcObj Allows user to directly use existing hclust object e.g. <pheatmap>$tree_col
  res_df <- data.frame(cutree(hcObj, k=num_cl))
  colnames(res_df) <- "Cluster"
  if(! is.null(rowname)) {
    res_df[rowname] <- rownames(res_df)
    rownames(res_df) <- NULL
    res_df <- res_df[ , c(2,1)]
  }
  return(res_df)
}
