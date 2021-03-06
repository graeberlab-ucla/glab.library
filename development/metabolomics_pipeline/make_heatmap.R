#heat.color = normal
#cluster_samples = TRUE
#matrix = MID
#samples = info
#width = 12
#height = 14

make_heatmap2 <- function (matrix, samples = samples, heat.color = normal, cluster_samples = TRUE, 
          width = NA, height = NA) 
{
  if (gsub("(.)*_|[0-9]", "", colnames(matrix))[1] == "Exp") {
    ext = "Relative Amounts"
  } else if (gsub("(.)*_|[0-9]", "", colnames(matrix))[1] == "MID") {
    ext = "MIDs"
  } else if (gsub("(.)*_|[0-9]", "", colnames(matrix))[1] == "FC") {
    ext = "Fractional Contribution"
  } else if (gsub("(.)*_|[0-9]", "", colnames(matrix))[1] == "Labeled") {
    ext = "Percent Labeled"
  }
  if (exists("samples") == F) 
    samples <- info
  ann <- select(samples, Condition, Cell.Number) %>% as.data.frame()
  rownames(ann) <- colnames(matrix)
  ann$Cell.Number <- rep(1, 6)
  if (exists("Norv") == T) {
    if (all(is.na(Norv))) 
      Norv <- rep(1, length(Norv))
    ann$Norvaline <- Norv
  } else ann$Norvaline <- 1
  ann_colors = list(Condition = colors[1:length(unique(gsub("_(.)*", 
                                                            "", colnames(matrix))))], Norvaline = c("white", "purple"), 
                    Cell.Number = c("white", "green"))
  names(ann_colors[["Condition"]]) <- unique(gsub("_(.)*", 
                                                  "", colnames(matrix)))
  matrix[is.na(matrix)] <- 0
  heatmap.title = paste(Title, "-Heatmap-", ext, ".pdf", sep = "")
  if (cluster_samples == F) 
    heatmap.title = paste(Title, "-Heatmap-", ext, "-Unclustered", 
                          ".pdf", sep = "")
  if (is.na(width) & is.na(height)) {
    pheatmap::pheatmap(matrix, cluster_row = T, cluster_col = cluster_samples, 
                       clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", 
                       color = colorRampPalette(heat.color)(100), border_color = "black", 
                       scale = "row", cellwidth = 20, cellheight = 10, 
                       annotation = ann, annotation_colors = ann_colors, 
                       show_colnames = F, main = paste(Title, ext, sep = "-"), 
                       filename = heatmap.title)
  }else {
    ### Added code to remove rows with zero variance
    zero_variance <- apply(matrix, 1, var)
    zero_variance_vec <- zero_variance[zero_variance == 0]
    matrix_filt <- matrix[!rownames(matrix) %in% names(zero_variance_vec),]
    matrix_remove <- matrix[rownames(matrix) %in% names(zero_variance_vec),]
    ### Changed matrix to matrix_filt
    pheatmap::pheatmap(matrix_filt, cluster_row = T, cluster_col = cluster_samples
                       ,clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
                       color = colorRampPalette(heat.color)(100), border_color = "black", 
                       scale = "row", cellwidth = 20, cellheight = 10, 
                       annotation = ann, annotation_colors = ann_colors, 
                       show_colnames = F, main = paste(Title, ext, sep = "-"), 
                       filename = heatmap.title, width = width, height = height)
    
    #pheatmap::pheatmap(matrix, cluster_row = T, cluster_col = cluster_samples, 
                       #clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", 
                       #color = colorRampPalette(heat.color)(100), border_color = "black", 
                       #scale = "row", cellwidth = 20, cellheight = 10, 
                       #annotation = ann, annotation_colors = ann_colors, 
                       #show_colnames = F, main = paste(Title, ext, sep = "-"), 
                       #filename = heatmap.title, width = width, height = height)
  }
  dev.off()
}
