library(DESeq2)
library(apeglm)
library(ashr)

dds_hs <- DESeq(dds_hs, test = "Wald")
# resultsNames(dds_hs)

plotDispEsts(dds_hs)
vsd <- vst(dds_hs, blind = FALSE)
plotPCA(vsd)
library(pheatmap)
sampleDists <- dist(t(assay(vsd)))
pheatmap(as.matrix(sampleDists))

res_HCM_vs_EA   <- results(dds_hs, contrast=c("condition","HCM","EA"))
summary(res_HCM_vs_EA)
sum(res_HCM_vs_EA$padj < 0.05, na.rm = TRUE)
# res_HCM_vs_EA_shrunk <- lfcShrink(dds_hs, contrast=c("condition","HCM","EA"), type="ashr")
res_HCM_vs_EA_shrunk <- lfcShrink(dds_hs, coef="condition_HCM_vs_EA", type="apeglm")
deg_HCM_vs_EA   <- subset(res_HCM_vs_EA_shrunk, padj < 0.05 & abs(log2FoldChange) > 1)

res_HCM_vs_Ctrl <- results(dds_hs, contrast=c("condition","HCM","Ctrl"))
summary(res_HCM_vs_Ctrl)
sum(res_HCM_vs_Ctrl$padj < 0.05, na.rm = TRUE)
res_HCM_vs_Ctrl_shrunk <- lfcShrink(dds_hs, coef="condition_HCM_vs_Ctrl", type="apeglm")
deg_HCM_vs_Ctrl <- subset(res_HCM_vs_Ctrl_shrunk, padj < 0.05 & abs(log2FoldChange) > 1)

res_EA_vs_Ctrl  <- results(dds_hs, contrast=c("condition","EA","Ctrl"))
summary(res_EA_vs_Ctrl)
sum(res_EA_vs_Ctrl$padj < 0.05, na.rm = TRUE)
res_EA_vs_Ctrl_shrunk <- lfcShrink(dds_hs, coef="condition_EA_vs_Ctrl", type="apeglm")
deg_EA_vs_Ctrl  <- subset(res_EA_vs_Ctrl_shrunk, padj < 0.05 & abs(log2FoldChange) > 1)



dds_ss <- DESeq(dds_ss, test = "Wald")

res_Sham_vs_Band   <- results(dds_ss, contrast=c("condition","Sham","Band"))
summary(res_Sham_vs_Band)
sum(res_Sham_vs_Band$padj < 0.05, na.rm = TRUE)
res_HCM_vs_EA_shrunk <- lfcShrink(dds_ss, coef="condition_Sham_vs_Band", type="apeglm")
deg_HCM_vs_EA   <- subset(res_HCM_vs_EA_shrunk, padj < 0.05 & abs(log2FoldChange) > 1)



