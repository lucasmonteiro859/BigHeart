library(DESeq2)
library(apeglm)
library(ashr)
library(readxl)
library(limma)
library(edgeR)
library(ggplot2)


RBP_df <- read_excel("data/Supplementary_Table_6_-_Gene_Expression_RBPs_R.xlsx")
RBP_ids <- RBP_df$gene_id
RBP_ids <- RBP_ids[RBP_ids != "ENSG00000167281"] # Giving problems

dds_hs <- DESeq(dds_hs, test = "Wald")
# resultsNames(dds_hs)

plotDispEsts(dds_hs)
vsd <- vst(dds_hs, blind = FALSE)
plotPCA(vsd)
library(pheatmap)
sampleDists <- dist(t(assay(vsd)))
pheatmap(as.matrix(sampleDists))

## HCM_vs_EA  
res_HCM_vs_EA   <- results(dds_hs, contrast=c("condition","HCM","EA"))
#res_HCM_vs_EA_shrunk <- lfcShrink(dds_hs, contrast=c("condition","HCM","EA"), type="ashr")
res_HCM_vs_EA_shrunk <- lfcShrink(dds_hs, coef="condition_HCM_vs_EA", type="apeglm")

summary(res_HCM_vs_EA_shrunk)
subset(res_HCM_vs_EA_shrunk, padj < 0.05 & abs(log2FoldChange) > 1)
subset(res_HCM_vs_EA_shrunk[RBP_ids, ], padj < 0.05 & abs(log2FoldChange) > 1)


## HCM_vs_Ctrl
res_HCM_vs_Ctrl <- results(dds_hs, contrast=c("condition","HCM","Ctrl"))
res_HCM_vs_Ctrl_shrunk <- lfcShrink(dds_hs, coef="condition_HCM_vs_Ctrl", type="apeglm")

summary(res_HCM_vs_Ctrl_shrunk)
subset(res_HCM_vs_Ctrl_shrunk, padj < 0.05 & abs(log2FoldChange) > 1)
sub1 <- subset(res_HCM_vs_Ctrl_shrunk[RBP_ids, ], padj < 0.05 & abs(log2FoldChange) > 1)
sub2 <- subset(RBP_df, padj_DCM_vs_Ctrl < 0.05 & abs(log2FC_DCM_vs_Ctrl) > 1 )
sub1
sub2
sub2$gene_id %in% rownames(sub1)

## 
res_EA_vs_Ctrl  <- results(dds_hs, contrast=c("condition","EA","Ctrl"))
res_EA_vs_Ctrl_shrunk <- lfcShrink(dds_hs, coef="condition_EA_vs_Ctrl", type="apeglm")

summary(res_EA_vs_Ctrl_shrunk)
subset(res_EA_vs_Ctrl_shrunk, padj < 0.05 & abs(log2FoldChange) > 1)
sub1 <- subset(res_EA_vs_Ctrl_shrunk[RBP_ids, ], padj < 0.05 & abs(log2FoldChange) > 1)
sub2 <- subset(RBP_df, padj_ICM_vs_Ctrl < 0.05 & abs(log2FC_ICM_vs_Ctrl) > 1 )
sub1
sub2
sub2$gene_id %in% rownames(sub1)



dds_ss <- DESeq(dds_ss, test = "Wald")

res_Sham_vs_Band   <- results(dds_ss, contrast=c("condition","Sham","Band"))
res_Sham_vs_Band_shrunk <- lfcShrink(dds_ss, coef="condition_Sham_vs_Band", type="apeglm")

summary(res_Sham_vs_Band_shrunk)
subset(res_Sham_vs_Band_shrunk, padj < 0.05 & abs(log2FoldChange) > 1)



### Limma
group <- factor(condition_hs)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

dge <- DGEList(counts = counts_matrix_hs, group = group)

# Filter lowly expressed genes (recommended for limma)
keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize
dge <- calcNormFactors(dge)

v <- voom(dge, design, plot = TRUE)
voom_mat <- v$E   # log2 CPM (normalized + precision weights)
weights <- v$weights

df_voom <- data.frame(
  Gene = rownames(voom_mat),
  voom_mat,
  check.names = FALSE
)

df_voom <- reshape2::melt(df_voom, id.vars = "Gene")
colnames(df_voom) <- c("Gene", "Sample", "Expression")

ggplot(df_voom, aes(Expression, color = Sample)) +
  geom_density(size = 1, alpha = 0.6) +
  theme_minimal() +
  labs(title = "voom logCPM density")

ggplot(df_voom, aes(Sample, Expression)) +
  geom_boxplot(outlier.size = 0.3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "voom-normalized expression")

pca <- prcomp(t(voom_mat))

pca_df <- data.frame(
  Sample = colnames(voom_mat),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Group = group
)

var_exp <- pca$sdev^2 / sum(pca$sdev^2)

ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCA (voom logCPM)",
    x = paste0("PC1 (", round(var_exp[1]*100,1), "%)"),
    y = paste0("PC2 (", round(var_exp[2]*100,1), "%)")
  )

library(pheatmap)

sample_dist <- dist(t(voom_mat))
sample_mat <- as.matrix(sample_dist)

annotation <- data.frame(Group = group)
rownames(annotation) <- colnames(voom_mat)

pheatmap(sample_mat,
         annotation_col = annotation,
         main = "Sample distances (voom)")


fit <- lmFit(v, design)
fit <- eBayes(fit)

plotSA(fit, main = "Mean–variance trend (voom)")

contr.matrix <- makeContrasts(
  HCM_vs_EA   = HCM - EA,
  #HCM_vs_Ctrl = HCM - Ctrl,
  #EA_vs_Ctrl  = EA - Ctrl,
  levels = design
)

fit2 <- contrasts.fit(fit, contr.matrix)
fit2 <- eBayes(fit2)

plotSA(fit2, main = "Mean–variance trend (voom)")

res_HCM_vs_EA_limma <- topTable(fit2, coef = "HCM_vs_EA", number = Inf)
deg_HCM_vs_EA_limma <- subset(res_HCM_vs_EA_limma, adj.P.Val < 0.05 & abs(logFC) > 1)
RBP_ids %in% rownames(deg_HCM_vs_EA_limma)

res_HCM_vs_Ctrl_limma <- topTable(fit2, coef = "HCM_vs_Ctrl", number = Inf)
deg_HCM_vs_Ctrl_limma <- subset(res_HCM_vs_Ctrl_limma, adj.P.Val < 0.05 & abs(logFC) > 1)
RBP_ids %in% rownames(deg_HCM_vs_Ctrl_limma)

res_EA_vs_Ctrl_limma <- topTable(fit2, coef = "EA_vs_Ctrl", number = Inf)
deg_EA_vs_Ctrl_limma <- subset(res_EA_vs_Ctrl_limma, adj.P.Val < 0.05 & abs(logFC) > 1)
RBP_ids %in% rownames(deg_EA_vs_Ctrl_limma)




group_hs <- factor(condition_hs)
design_hs <- model.matrix(~0 + group_hs)
colnames(design_hs) <- levels(group_hs)

dge_hs <- DGEList(counts = counts_matrix_hs, group = group_hs)

# Filter lowly expressed genes (recommended for limma)
keep <- filterByExpr(dge_hs, design_hs)
dge_hs <- dge_hs[keep, , keep.lib.sizes = FALSE]

# Normalize
dge_hs <- calcNormFactors(dge_hs)

v <- voom(dge_hs, design_hs, plot = TRUE)


group_ss <- factor(condition_ss)
design_ss <- model.matrix(~ 0 + group_ss)
colnames(design_ss) <- levels(group_ss)

dge_ss <- DGEList(counts = counts_matrix_ss, group = group_ss)

# Filter lowly expressed genes (recommended for limma)
keep <- filterByExpr(dge_ss, design_ss)
dge_ss <- dge_ss[keep, , keep.lib.sizes = FALSE]

dge_ss <- calcNormFactors(dge_ss)

v_ss <- voom(dge_ss, design_ss)

fit_ss <- lmFit(v_ss, design_ss)
fit_ss <- eBayes(fit_ss)

plotSA(fit_ss, main = "Mean–variance trend (voom)")
       
contrast_ss <- makeContrasts(
  Sham_vs_Band = Sham - Band,
  levels = design_ss
)

fit2_ss <- contrasts.fit(fit_ss, contrast_ss)
fit2_ss <- eBayes(fit2_ss)

res_Sham_vs_Band_limma <- topTable(fit2_ss, coef = "Sham_vs_Band", number = Inf)
deg_HCM_vs_EA_limma <- subset(res_Sham_vs_Band_limma, adj.P.Val < 0.05 & abs(logFC) > 1)
RBP_ids %in% rownames(deg_HCM_vs_EA_limma)

