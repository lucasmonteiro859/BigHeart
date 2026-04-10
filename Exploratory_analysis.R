library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

make_df <- function(matrix) {
  matrix_df <- data.frame(Gene = rownames(matrix), matrix, check.names = FALSE)
  df_plot <- melt(matrix_df, id.vars = "Gene")
  df_plot <- df_plot[df_plot$value > 0, ]
  colnames(df_plot) <- c("Gene", "Sample", "Counts")
  return(df_plot)
}

plot_density <- function(df_plot) { 
ggplot(df_plot, aes(x = Counts, color = Sample)) +
  geom_density(size = 1, alpha = 0.6) +
  labs(title = paste("Density of log2(counts + 1) for samples"),
       x = "log2(count + 1)",
       y = "Density",
       color = "Sample") +
  theme_minimal() +
  theme(legend.position = "right")
}

plot_box  <- function(df_plot) { 
ggplot(df_plot, aes(x = Sample, y = Counts)) +
  geom_boxplot(fill = "lightblue", outlier.size = 0.5) +
  labs(title = paste("Boxplot of log2(counts + 1) for samples"),
       x = "Sample",
       y = "log2(count + 1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # rotate labels for multiple samples
}

plot_PCA  <- function(matrix, condition) { 
  pca_input <- t(matrix)
  pca_results <- prcomp(pca_input, center = TRUE, scale. = TRUE)
  var_explained <- summary(pca_results)$importance[2, ]
  pca_df <- data.frame(
    Sample = rownames(pca_input),
    PC1 = pca_results$x[, "PC1"],
    PC2 = pca_results$x[, "PC2"],
    Group = condition  # make sure names match samples
  )
  
  plt <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5, size = 3) +
    labs(
      title = "PCA of log2(normalized counts + 1)",
      x = paste0("PC1 (", round(var_explained["PC1"]*100, 1), "%)"),
      y = paste0("PC2 (", round(var_explained["PC2"]*100, 1), "%)")
    ) +
    theme_minimal() +
    scale_color_brewer(palette = "Set1")
  
  print(plt)
  
  lib_size <- colSums(matrix)            
  genes_detected <- colSums(matrix > 0)  
  
  # Combine with PCA results
  pca_qc_df <- data.frame(
    Sample = pca_df$Sample,
    PC1 = pca_df$PC1,
    PC2 = pca_df$PC2,
    Group = pca_df$Group,
    lib_size = lib_size[pca_df$Sample],
    genes_detected = genes_detected[pca_df$Sample]
  )
  cor_table <- cor(pca_qc_df[, c("PC1", "PC2", "lib_size", "genes_detected")])
  return(cor_table)
}

plot_sample_heatmap <- function(log_counts, condition_vector) {
  sampleDists <- dist(t(log_counts))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- colnames(log_counts)
  colnames(sampleDistMatrix) <- colnames(log_counts)
  
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  annotation_col <- data.frame(Condition = factor(condition_vector))
  rownames(annotation_col) <- colnames(log_counts)
  
  pheatmap(
    sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors,
    annotation_col = annotation_col,
    main = "Normalized Disease Counts",
    display_numbers = FALSE,
    fontsize_row = 10,
    fontsize_col = 10
  )
}


df_plot_hs <- make_df(counts_matrix_hs)
plot_density(df_plot_hs)
plot_box(df_plot_hs)
cor_table_hs <- plot_PCA(counts_matrix_hs, condition_hs)
plot_sample_heatmap(counts_matrix_hs, condition_hs)


df_plot_ss <- make_df(counts_matrix_ss)
plot_density(df_plot_ss)
plot_box(df_plot_ss)
cor_table_ss <- plot_PCA(counts_matrix_ss, condition_ss)
plot_sample_heatmap(counts_matrix_ss, condition_ss)

