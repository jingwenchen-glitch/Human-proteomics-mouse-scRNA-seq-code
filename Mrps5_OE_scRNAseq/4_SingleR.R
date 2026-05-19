# 清除系统环境变量，加载R包：
rm(list=ls())
options(stringsAsFactors = F) 
setwd("/Users/jingwenchen/Desktop/Ph.D/Jin Lab/scRNA-seq/Results/")

# SingleR预测细胞类型
# 安装并加载必要包
if (!requireNamespace("SingleR", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("SingleR")
}
if (!requireNamespace("celldex", quietly = TRUE)) {
  BiocManager::install("celldex")
}

library(SingleR)
library(celldex)
library(Seurat)
library(SeuratObject)
library(SummarizedExperiment)
library(ggplot2)
library(patchwork)
library(cowplot)
library(clustree)
library(dplyr)

# 读取文件
combined.int <- readRDS("combined_with_proportions.rds")

# 载入参考数据集
ref <- celldex::ImmGenData()  # 适用于小鼠数据

# 准备 Seurat 对象中的表达矩阵
# 提取归一化后的表达矩阵（建议使用 LogNormalize 后的结果）
seurat_matrix <- GetAssayData(combined.int, assay = "RNA", layer = "data")
# 提取元数据用于后续赋值
seurat_meta <- combined.int@meta.data
# 运行 SingleR 进行注释
pred <- SingleR(test = seurat_matrix, 
                ref = ref, 
                labels = ref$label.main)  # 使用主标签（大类）
# 将预测结果添加回 Seurat 对象
combined.int$SingleR_label <- pred$labels
# 可视化看看预测结果
DimPlot(combined.int, 
        group.by = "SingleR_label", 
        repel = TRUE,
        label.size = 5) +
  theme(plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text = element_text(size = 12, face = "bold", color = "black")
  )

# 保存为 PDF, PNG
ggsave(filename = "SingleR/SingleR_Annotation_UMAP.pdf", width = 6, height = 4)
ggsave(filename = "SingleR/SingleR_Annotation_UMAP.png", width = 6, height = 4, dpi = 300)
# 保存SingleR后的 Seurat 对象
saveRDS(combined.int, file = "SingleR_reduction.rds")

combined.int <- readRDS("SingleR_reduction.rds")

# 基于SingleR的结果，对样本按照细胞大类进行再分群分析
table(combined.int$SingleR_label)  # 确认注释中有 T cells
# 提取标注为 T cells 的细胞
tcell_subset <- subset(combined.int, subset = SingleR_label == "T cells")

# 查看 cluster 5 的线粒体比例
VlnPlot(tcell_subset, features = "percent.mt", group.by = "seurat_clusters")+
  ggtitle("percent.mt") +
  ylab("Percentage (%)") +  # 设置Y轴标签  
  theme(
    axis.text.x = element_text(size = 12, face = "bold", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold")) + NoLegend()
ggsave("SingleR/percent_mt_by_cluster.pdf", width = 5, height = 4)
ggsave("SingleR/percent_mt_by_cluster.png", width = 5, height = 4, dpi = 300)

# 检查 nFeature/nCount 情况
VlnPlot(tcell_subset, features = c("nFeature_RNA", "nCount_RNA"), group.by = "seurat_clusters") & 
  theme(
    axis.text.x = element_text(size = 12, face = "bold", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_blank(),
    plot.title = element_text(size = 14, face = "bold"))
ggsave("SingleR/nFeature_nCount_by_cluster.pdf", width = 5, height = 4)
ggsave("SingleR/nFeature_nCount_by_cluster.png", width = 5, height = 4, dpi = 300)

# 过滤cluster 5
# 保留 cluster 0~4 的细胞
valid_cells <- WhichCells(tcell_subset, idents = c("0", "1", "2", "3", "4"))
tcell_clean <- subset(tcell_subset, cells = valid_cells)
dim(tcell_clean)

# 重新构建seurat对象并降维聚类
# 提取表达矩阵，构建新对象（可选：指定 min.cells / min.features）
tcell_counts <- GetAssayData(tcell_clean, layer = "counts")
tcell_meta <- tcell_clean@meta.data

# 创建新 Seurat 对象
tcell_clean <- CreateSeuratObject(counts = tcell_counts, meta.data = tcell_meta)

# 清理旧的 metadata
tcell_clean@meta.data <- tcell_clean@meta.data %>%
  dplyr::select(-any_of(c(
    # 聚类相关
    "seurat_clusters", 
    "cluster_id", 
    "cluster_with_prop",
    "RNA_snn_res.0.1", "RNA_snn_res.0.15", "RNA_snn_res.0.3", 
    "RNA_snn_res.0.5", "RNA_snn_res.0.7", "RNA_snn_res.1",
    "scDblFinder_score", "scDblFinder_class",
    "SingleR_label"
  )))

# 标准化 + 高变基因识别
tcell_clean <- NormalizeData(tcell_clean, 
                          normalization.method = "LogNormalize",
                          scale.factor = 1e4) 
tcell_clean <- FindVariableFeatures(tcell_clean)
p1 <- VariableFeaturePlot(tcell_clean) 
p1
# 数据归一化
tcell_clean <- ScaleData(tcell_clean)
# PCA线性降维
tcell_clean <- RunPCA(tcell_clean, features = VariableFeatures(object = tcell_clean))
# 可视化PCA结果
VizDimLoadings(tcell_clean, dims = 1:2, reduction = "pca")
DimPlot(tcell_clean, reduction = "pca") 
DimHeatmap(tcell_clean, dims = 1:12, cells = 500, balanced = TRUE)
# 运行 UMAP & TSNE 进行非线性降维可视化
tcell_clean <- RunUMAP(tcell_clean, dims = 1:20)
DimPlot(tcell_clean, reduction = "umap", label=F ) 
tcell_clean <- RunTSNE(tcell_clean, dims = 1:20)
DimPlot(tcell_clean, reduction = "tsne", label=F ) 

# 下游分析
# 构建最近邻图
tcell_clean <- FindNeighbors(tcell_clean, dims = 1:20)
# 进行聚类，resolution 可调节分辨率
for (res in c(0.1, 0.175, 0.3, 0.5, 0.7, 1)) {
  tcell_clean <- FindClusters(tcell_clean, #graph.name = "CCA_snn", 
                               resolution = res, algorithm = 1)
}
colnames(tcell_clean@meta.data)
apply(tcell_clean@meta.data[, grep("RNA_snn",colnames(tcell_clean@meta.data))], 2, table)

p1_dim<-plot_grid(ncol = 3, DimPlot(tcell_clean, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
                    ggtitle("louvain_0.1"), DimPlot(tcell_clean, reduction = "umap", group.by = "RNA_snn_res.0.175") + 
                    ggtitle("louvain_0.175"), DimPlot(tcell_clean, reduction = "umap", group.by = "RNA_snn_res.0.3") + 
                    ggtitle("louvain_0.3"))
p1_dim
ggsave(plot=p1_dim, filename="tcell_clean/Dimplot_diff_resolution_low.png", width = 14, dpi = 300)

p2_dim <- plot_grid(ncol = 3, DimPlot(tcell_clean, reduction = "umap", group.by = "RNA_snn_res.0.5") + 
                      ggtitle("louvain_0.5"), DimPlot(tcell_clean, reduction = "umap", group.by = "RNA_snn_res.0.7") + 
                      ggtitle("louvain_0.7"), DimPlot(tcell_clean, reduction = "umap", group.by = "RNA_snn_res.1") + 
                      ggtitle("louvain_1"))
p2_dim
ggsave(plot=p2_dim, filename="tcell_clean/Dimplot_diff_resolution_high.png", width = 14, dpi = 300)

p3_tree <- clustree(tcell_clean@meta.data, prefix = "RNA_snn_res.")
p3_tree
ggsave(plot=p3_tree, filename="tcell_clean/Tree_diff_resolution.png", dpi = 300)
table(tcell_clean@active.ident) 

# 将 RNA_snn_res.0.175 设为 active.ident（默认聚类标签）
tcell_clean <- SetIdent(tcell_clean, value = "RNA_snn_res.0.175")
# 同时修改 meta.data 中的 seurat_clusters，使它等于 RNA_snn_res.0.175
tcell_clean$seurat_clusters <- Idents(tcell_clean)
# 保存降维聚类后的 Seurat 对象
saveRDS(tcell_clean, "tcell_clean.rds")