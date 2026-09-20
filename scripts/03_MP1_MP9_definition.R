mp_scores <- sce.all_int@meta.data[, grep("^MP[1-9]$", colnames(sce.all_int@meta.data))]
norm_scores <- as.data.frame(scale(mp_scores, center = TRUE, scale = FALSE))
assign_mp <- function(cell_scores) {
  sorted <- sort(cell_scores, decreasing = TRUE)
  top_score <- sorted[1]
  second_score <- sorted[2]
  
  if (second_score <= 0.85 * top_score) {
    return(names(sorted)[1])  
  } else {
    return("Unresolved")
  }
}
sce.all_int@meta.data$cluster <- apply(norm_scores, 1, assign_mp)
sce.all_int@meta.data$cluster <- factor(
  sce.all_int@meta.data$cluster,
  levels = c(paste0("MP", 1:9), "Unresolved")
)
DimPlot(sce.all_int,label = T)
Idents(sce.all_int)<-sce.all_int$cluster
table(sce.all_int$cluster)
DimPlot(sce.all_int, reduction = "umap_MP", group.by = "cluster") + theme(aspect.ratio = 1)