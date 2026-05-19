# 清除系统环境变量，加载R包：
rm(list=ls())
options(stringsAsFactors = F) 
setwd("/Users/jingwenchen/Desktop/Ph.D/Jin Lab/scRNA-seq/Results/")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scDblFinder")
BiocManager::install("Seurat")
BiocManager::install("SeuratDisk")

# 加载必要包
library(Seurat)
library(data.table)
library(dplyr)
library(scDblFinder)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
# 读取 CellRanger 输出
NC.data <- Read10X("/Users/jingwenchen/Desktop/Ph.D/Jin Lab/scRNA-seq/Summary/1_Cellranger_result/NC/filtered_feature_bc_matrix/")
OE.data <- Read10X("/Users/jingwenchen/Desktop/Ph.D/Jin Lab/scRNA-seq/Summary/1_Cellranger_result/OE/filtered_feature_bc_matrix/")

# 创建 Seurat 对象 
NC <- CreateSeuratObject(NC.data, project = "NC", min.cells = 3, min.features = 200)
OE <- CreateSeuratObject(OE.data, project = "OE", min.cells = 3, min.features = 200)

# 封装函数
process_sample <- function(seurat_obj, sample_name, qc_only = FALSE, prefix = "QC/QC_filter") {
  n_cells_before <- ncol(seurat_obj)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  # 若只看QC，提前返回
  if (qc_only) {
    message("Only QC plots generated, sample not filtered.")
    return(seurat_obj)
  }
  # 基础质控
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  seurat_obj <- subset(seurat_obj, 
                       subset = nFeature_RNA > 500 & 
                         nFeature_RNA < 7000 &
                         percent.mt < 25)
  message("Filtered ", n_cells_before - ncol(seurat_obj), " cells from ", sample_name)
  # 前处理
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 30)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  # 多细胞识别
  sce <- as.SingleCellExperiment(seurat_obj)
  sce <- scDblFinder(sce, dbr=0.05, clusters = NULL)
  seurat_obj$scDblFinder_score <- sce$scDblFinder.score
  seurat_obj$scDblFinder_class <- sce$scDblFinder.class
  # 可视化
  p1 <- VlnPlot(seurat_obj, features = "scDblFinder_score", group.by = "scDblFinder_class") +
    ggtitle(paste0(sample_name, ": scDblFinder Score")) & theme(
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 14, angle = 0, hjust = 0.5, face = "bold"),
      legend.text = element_text(size = 12, face = "bold")) + NoLegend()
  p2 <- DimPlot(seurat_obj, group.by = "scDblFinder_class", label = TRUE) +
    ggtitle(paste0(sample_name, ": Doublet Classification")) & theme(
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 14, angle = 0, hjust = 0.5, face = "bold"),
      legend.text = element_text(size = 12, face = "bold"))
  p_combined <- p1 + p2 + plot_layout(ncol = 2)
  # 保存图
  ggsave(paste0(prefix, sample_name, "_doublet_detection.pdf"), plot = p_combined, width = 8, height = 6)
  ggsave(paste0(prefix, sample_name, "_doublet_detection.png"), plot = p_combined, width = 8, height = 6, dpi = 300)

  # 过滤掉多细胞
  seurat_obj <- subset(seurat_obj, subset = scDblFinder_class == "singlet")
  cat("Removed", n_cells_before - ncol(seurat_obj), "doublets from", sample_name, "\n")
  seurat_obj$sample <- sample_name
  return(seurat_obj)
}
# 仅查看 NC 样本质控前的数据分布
process_sample(NC, sample_name = "NC", qc_only = TRUE)
# 仅查看 OE 样本质控前的数据分布
process_sample(OE, sample_name = "OE", qc_only = TRUE)
# 可视化质控前数据
plot_qc_merge <- function(seurat_list, sample_names, prefix = "QC/QC_before_filter") {
  for (i in seq_along(seurat_list)) {
    seurat_list[[i]]$sample <- sample_names[i]
    # 确保 percent.mt 存在
    if (!"percent.mt" %in% colnames(seurat_list[[i]]@meta.data)) {
      seurat_list[[i]][["percent.mt"]] <- PercentageFeatureSet(seurat_list[[i]], pattern = "^mt-")
    }
  }
  merged_obj <- merge(seurat_list[[1]], y = seurat_list[-1])
  p3 <- VlnPlot(merged_obj,
                group.by = "orig.ident",
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                pt.size = 0.1, ncol = 3) + NoLegend() & theme(
                  axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5, face = "bold"),
                  axis.text.y = element_text(size = 12, face = "bold"),
                  axis.title = element_blank())
  ggsave(paste0(prefix, "_violin.pdf"), plot = p3, width = 5, height = 4)
  ggsave(paste0(prefix, "_violin.png"), plot = p3, width = 5, height = 4, dpi = 300)
  p4 <- FeatureScatter(merged_obj, "nCount_RNA", "nFeature_RNA",
                       group.by = "orig.ident", pt.size = 0.5) & theme(
                         axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5, face = "bold"),
                         axis.text.y = element_text(size = 12, face = "bold"),
                         axis.title = element_text(size = 14, face = "bold"),
                         legend.title = element_text(size = 14, face = "bold"),
                         legend.text = element_text(size = 12, face = "bold")
                       )
  ggsave(paste0(prefix, "_scatter.pdf"), plot = p4, width = 5, height = 4)
  ggsave(paste0(prefix, "_scatter.png"), plot = p4, width = 5, height = 4, dpi = 300)
}
plot_qc_merge(list(NC, OE), c("NC", "OE"))

# 分别处理两个样本
NC_clean <- process_sample(NC, "NC")
OE_clean <- process_sample(OE, "OE")

# 可视化质控后数据
plot_qc_after_filter <- function(seurat_list, sample_names, prefix = "QC/QC_after_filter") {
  # 添加分组信息
  for (i in seq_along(seurat_list)) {
    seurat_list[[i]]$sample <- sample_names[i]
  }
  # 合并数据
  merged_obj <- merge(seurat_list[[1]], y = seurat_list[-1])
  # 小提琴图（nFeature_RNA, nCount_RNA, percent.mt）
  p5 <- VlnPlot(merged_obj,
                group.by = "orig.ident",
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                pt.size = 0.1, ncol = 3) + NoLegend() & theme(
                  axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5, face = "bold"),
                  axis.text.y = element_text(size = 12, face = "bold"),
                  axis.title = element_blank())
  ggsave(paste0(prefix, "_violin.pdf"), plot = p5, width = 5, height = 4)
  ggsave(paste0(prefix, "_violin.png"), plot = p5, width = 5, height = 4, dpi = 300)
  p6 <- FeatureScatter(merged_obj, "nCount_RNA", "nFeature_RNA",
                       group.by = "orig.ident", pt.size = 0.5) & theme(
                         axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5, face = "bold"),
                         axis.text.y = element_text(size = 12, face = "bold"),
                         axis.title = element_text(size = 14, face = "bold"),
                         legend.title = element_text(size = 14, face = "bold"),
                         legend.text = element_text(size = 12, face = "bold")
                       )
  ggsave(paste0(prefix, "_scatter.pdf"), plot = p6, width = 5, height = 4)
  ggsave(paste0(prefix, "_scatter.png"), plot = p6, width = 5, height = 4, dpi = 300)
}
plot_qc_after_filter(seurat_list = list(NC_clean, OE_clean), sample_names = c("NC", "OE"))
# 合并成一个整体对象
combined_clean <- merge(NC_clean, y = OE_clean, add.cell.ids = c("NC", "OE"), project = "Combined_clean")
# 保存合并后的 Seurat 对象
saveRDS(combined_clean, file = "QC_after_filtered_merged_seurat.rds")
