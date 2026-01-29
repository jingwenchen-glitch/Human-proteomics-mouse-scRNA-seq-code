# 清除系统环境变量，加载R包：
rm(list=ls())
options(stringsAsFactors = F) 
setwd("/Users/jingwenchen/Desktop/Ph.D/Jin Lab/scRNA-seq/Results/")

# 加载必要包 
library(Seurat)
library(data.table)
library(dplyr)
library(patchwork)
library(harmony)
library(ggplot2)
library(clustree)
library(cowplot)
BiocManager::install("lisi")

# 读取文件
combined <- readRDS("QC_after_filtered_merged_seurat.rds")
names(combined@assays$RNA@layers)
combined[["RNA"]]$counts 
LayerData(combined, assay = "RNA", layer = "counts")

# 使用JoinLayers函数对layers进行合并
combined <- JoinLayers(combined)
combined
# 查看combined内部的一些信息，以此来检查数据是否完整
dim(combined[["RNA"]]$counts )
as.data.frame(combined@assays$RNA$counts[1:10, 1:2])
head(combined@meta.data, 10)
table(combined$orig.ident) 
length(combined$orig.ident)

# 标准化 + 高变基因识别
combined <- NormalizeData(combined, 
                          normalization.method = "LogNormalize",
                          scale.factor = 1e4) 
combined <- FindVariableFeatures(combined)
p1 <- VariableFeaturePlot(combined) 
p1
# 数据归一化
combined <- ScaleData(combined)
# PCA线性降维
combined <- RunPCA(combined, features = VariableFeatures(object = combined))
# 可视化PCA结果
VizDimLoadings(combined, dims = 1:2, reduction = "pca")
DimPlot(combined, reduction = "pca") 
DimHeatmap(combined, dims = 1:12, cells = 500, balanced = TRUE)
# 运行 UMAP & TSNE 进行非线性降维可视化
combined <- RunUMAP(combined, dims = 1:20)
DimPlot(combined, reduction = "umap", label=F, group.by = "orig.ident") 
combined <- RunTSNE(combined, dims = 1:20)
DimPlot(combined, reduction = "tsne",label=F ) 

# 下游分析
# 构建最近邻图
combined <- FindNeighbors(combined, dims = 1:20)
# 进行聚类，resolution 可调节分辨率
#设置不同的分辨率，观察分群效果(选择哪一个？)
combined.all <- combined
for (res in c(0.1, 0.15, 0.3, 0.5, 0.7, 1)) {
  combined.all <- FindClusters(combined.all, #graph.name = "CCA_snn", 
                               resolution = res, algorithm = 1)
}
colnames(combined.all@meta.data)
apply(combined.all@meta.data[,grep("RNA_snn",colnames(combined.all@meta.data))],2,table)

p1_dim<-plot_grid(ncol = 3, DimPlot(combined.all, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
                    ggtitle("louvain_0.1"), DimPlot(combined.all, reduction = "umap", group.by = "RNA_snn_res.0.15") + 
                    ggtitle("louvain_0.15"), DimPlot(combined.all, reduction = "umap", group.by = "RNA_snn_res.0.3") + 
                    ggtitle("louvain_0.3"))
p1_dim
ggsave(plot=p1_dim, filename="Cluster/Dimplot_diff_resolution_low.png",width = 14, dpi = 300)

p2_dim <- plot_grid(ncol = 3, DimPlot(combined.all, reduction = "umap", group.by = "RNA_snn_res.0.5") + 
                      ggtitle("louvain_0.5"), DimPlot(combined.all, reduction = "umap", group.by = "RNA_snn_res.0.7") + 
                      ggtitle("louvain_0.7"), DimPlot(combined.all, reduction = "umap", group.by = "RNA_snn_res.1") + 
                      ggtitle("louvain_1"))
p2_dim
ggsave(plot=p2_dim, filename="Cluster/Dimplot_diff_resolution_high.png",width = 14, dpi = 300)

p3_tree <- clustree(combined.all@meta.data, prefix = "RNA_snn_res.")
p3_tree
ggsave(plot=p3_tree, filename="Cluster/Tree_diff_resolution.png", dpi = 300)
table(combined.all@active.ident) 

# 将 RNA_snn_res.0.15 设为 active.ident（默认聚类标签）
combined.all <- SetIdent(combined.all, value = "RNA_snn_res.0.15")
# 同时修改 meta.data 中的 seurat_clusters，使它等于 RNA_snn_res.0.15
combined.all$seurat_clusters <- Idents(combined.all)
# 保存降维聚类后的 Seurat 对象
saveRDS(list(combined = combined, combined_all = combined.all), file = "reduction_combined.rds")

devtools::install_github("immunogenomics/lisi")
library(lisi)
emb <- Embeddings(combined, "pca")                     # 或 "umap"
meta <- combined@meta.data[, c("sample","seurat_clusters")]
lisi_res <- compute_lisi(emb, meta, c("sample"))  # iLISI by sample
summary(lisi_res$sample)
