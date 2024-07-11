# Reusable Objects for Heat Maps
# Copyright (C) 2023-4 Y. David Chen

library(pheatmap)

CL_PARAMS <- c("euclidean", "ward.D")

HEAT_COLS <- colorRampPalette(c("blue","white","red"))(1024)

ANN_COLOR <- list(
  AgeGroup = c(Over50="black", Under50="lightgray"),
  Sex = c(Male="black", Female="lightgray"),
  HeartFailure = c(Yes="black", No="lightgray"),
  Subtype = c(Ischemic="black", Dilated="lightgray", `(control)`="white"),
  PrincipalGroup = c(A="black", B="lightgray"),
  HasRepeat  = c(Yes="black", No="lightgray")
)

wrapper_hm <- function(mat, num_clusters=1, annot_samps=NULL, num_genes_cut=1, annot_gene=NULL) {
  #'@description Wrapper to reduce code duplication & increase consistency
  #'@return hclust tree object
  phObj <- pheatmap(
    mat,
    cutree_cols = num_clusters,
    cutree_rows = num_genes_cut,
    clustering_distance_cols = CL_PARAMS[1],
    clustering_distance_rows = CL_PARAMS[1],
    clustering_method = CL_PARAMS[2], 
    annotation_col = annot_samps,
    annotation_row = annot_gene,
    annotation_colors = ANN_COLOR,
    color = HEAT_COLS,
    fontsize = 15,
    show_rownames = FALSE,
    show_colnames = FALSE,
    border_color = NA
  )
  return(phObj$tree_col)
}

