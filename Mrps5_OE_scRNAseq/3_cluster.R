# 清除系统环境变量，加载R包：
rm(list=ls())
options(stringsAsFactors = F) 
setwd("/Users/jingwenchen/Desktop/Ph.D/Jin Lab/scRNA-seq/Results/")

library(Seurat)
library(Matrix)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(patchwork)

# 读取文件
objs <- readRDS("reduction_combined.rds")
combined.int <- objs$combined_all
table(combined.int@active.ident) 
colnames(combined.int@meta.data) 

# 基于每个 cluster 的平均表达谱，计算 cluster 间的相关性矩阵，并绘制热图
# 使用 AggregateExpression 获取 pseudo-bulk 表达矩阵（基于 cluster）
pseudo_bulk <- AggregateExpression(combined.int, group.by = "seurat_clusters", assay = "RNA", slot = "data")
expr_mat <- pseudo_bulk$RNA  # 提取 RNA assay 的表达矩阵
# 确保是 matrix 格式
expr_mat <- as.matrix(expr_mat)
# 计算相关性矩阵
cor_mat <- cor(expr_mat, method = "pearson")
# 去掉行名和列名中的 "g"
rownames(cor_mat) <- sub("^g", "", rownames(cor_mat))
colnames(cor_mat) <- sub("^g", "", colnames(cor_mat))

# 绘制热图
# 设置色条范围与标签
breaks <- seq(-1, 1, length.out = 100)
colors <- colorRampPalette(c("darkred", "white", "darkblue"))(length(breaks) - 0.5)
# 保存为 PDF
# PDF
pdf("Cluster/Cluster_Correlation_Heatmap.pdf", width = 7, height = 6)
pheatmap(cor_mat,
         color = colors,
         breaks = breaks,
         display_numbers = TRUE,
         number_color = "white",
         fontsize_number = 8,
         main = "Cluster Correlation (Pearson)",
         border_color = NA,
         angle_col = 0)   # 列名竖直显示
dev.off()

# PNG
png("Cluster/Cluster_Correlation_Heatmap.png", width = 7, height = 6, units = "in", res = 300)
pheatmap(cor_mat,
         color = colors,
         breaks = breaks,
         display_numbers = TRUE,
         number_color = "white",
         fontsize_number = 8,
         main = "Cluster Correlation (Pearson)",
         border_color = NA,
         angle_col = 0)   # 列名竖直显示
dev.off()
dev.new()

# 统计每个 cluster 在不同 sample 中的细胞数量和百分比
# 提取 cluster 和 sample 信息
df <- data.frame(
  cluster = combined.int$seurat_clusters,
  sample = combined.int$sample
)
# 按 cluster 和 sample 计数
cluster_sample_stats <- df %>%
  group_by(cluster, sample) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(percent = round(100 * cell_count / sum(cell_count), 2))
write.csv(cluster_sample_stats, "Cluster/cluster_sample_stats_counts_and_percentage.csv")

# 可视化为堆叠柱状图
ggplot(cluster_sample_stats, aes(x = factor(cluster), y = percent, fill = sample)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Cluster", y = "Proportion") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        legend.title = element_text(size = 14, face = "bold", color = "black"),
        legend.text = element_text(size = 12, face = "bold", color = "black")
        )
# 保存图像
ggsave("Cluster/Cluster_Sample_Composition_Barplot.pdf", width = 6, height = 4)
ggsave("Cluster/Cluster_Sample_Composition_Barplot.png", width = 6, height = 4, dpi = 300)

# 每个 sample 内部各 cluster 占比
df %>%
  group_by(sample, cluster) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(percent = round(100 * cell_count / sum(cell_count), 2)) -> sample_cluster_stats
write.csv(sample_cluster_stats, "Cluster/sample_cluster_stats_percentage.csv")

# 绘制堆积柱状图（cluster 占比，按 sample 分组）
ggplot(sample_cluster_stats, aes(x = sample, y = percent, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_minimal(base_size = 14) +
  labs(x = "Sample", y = "Percentage (%)") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        legend.title = element_text(size = 14, face = "bold", color = "black"),
        legend.text = element_text(size = 12, face = "bold", color = "black")
  )

ggsave("Cluster/Cluster_Percentage_By_Sample.pdf", width = 6, height = 4)
ggsave("Cluster/Cluster_Percentage_By_Sample.png", width = 6, height = 4, dpi = 300)

# 计算聚类比例并修改图例标签
# 计算每个聚类的比例
clusters <- combined.int@meta.data$RNA_snn_res.0.15
cluster_props <- prop.table(table(clusters)) * 100
cluster_props <- round(cluster_props, 2)  # 保留2位小数
# 创建带比例的新标签
new_labels <- paste0(names(cluster_props), " (", cluster_props, "%)")
# 创建用于UMAP点标签的列（只包含cluster编号）
combined.int$cluster_id <- as.character(clusters)
# 创建用于图例的列（包含比例）
combined.int$cluster_with_prop <- factor(clusters, 
                                         levels = names(cluster_props),
                                         labels = new_labels)
# 绘制UMAP图
p <- DimPlot(combined.int, 
             reduction = "umap",
             group.by = "cluster_id",  # 使用只含编号的分组
             label = TRUE,             # 显示点标签（cluster编号）
             label.size = 6,           # 标签大小
             repel = TRUE) +           # 避免标签重叠
  theme(plot.title = element_blank())
p
# 修改图例标签为带比例的形式
# 获取当前的颜色映射
color_mapping <- levels(combined.int$cluster_with_prop)
names(color_mapping) <- levels(combined.int$cluster_id)

# 添加带比例的新图例
p <- p + scale_color_discrete(name = "Clusters",
                              labels = color_mapping,
                              guide = guide_legend(override.aes = list(size = 3)))+       
  theme(plot.title = element_blank(),
        legend.title = element_text(size = 14, face = "bold", color = "black"),
        legend.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text = element_text(size = 12, face = "bold", color = "black")
  )

# 移除点上的颜色图例指示（只保留图例框）
p$layers[[1]]$aes_params$show.legend <- FALSE
print(p)
ggsave("Cluster/UMAP_with_proportions.pdf", p, width = 6, height = 4)
ggsave("Cluster/UMAP_with_proportions.png", p, width = 6, height = 4, dpi = 300)

saveRDS(combined.int, "combined_with_proportions.rds")

# UMAP 按 sample 着色
p_sample <- DimPlot(combined.int, group.by = "sample", label = FALSE) +
  theme(plot.title = element_blank(),
        legend.title = element_text(size = 14, face = "bold", color = "black"),
        legend.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text = element_text(size = 12, face = "bold", color = "black")
  )
# 按 cluster 分组绘制 UMAP
p_cluster <- DimPlot(combined.int, group.by = "seurat_clusters", label = TRUE, label.size = 6, repel = TRUE) +
  theme(plot.title = element_blank(),
        legend.title = element_text(size = 14, face = "bold", color = "black"),
        legend.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text = element_text(size = 12, face = "bold", color = "black")
  )
p_merge <- p_sample + p_cluster + plot_layout(ncol = 2)
p_merge
ggsave("Cluster/UMAP_sample_vs_cluster.pdf", p_merge, width = 8, height = 4)
ggsave("Cluster/UMAP_sample_vs_cluster.png", p_merge, width = 8, height = 4, dpi = 300)

library(scales)
library(grid)
# NC vs OE 按 cluster 分组绘制 UMAP
# 确保 cluster 设置为主身份
Idents(combined.int) <- "seurat_clusters"
# 准备每个 sample 的细胞数用于 facet label
cell_counts <- table(combined.int$orig.ident)
label_df <- data.frame(
  sample = names(cell_counts),
  count = as.numeric(cell_counts),
  label = paste0(names(cell_counts), "\n(", comma(as.numeric(cell_counts)), " cells)")
)
# 构建主图
p1 <- DimPlot(combined.int, group.by = "seurat_clusters", split.by = "sample", label = TRUE, label.size = 6) +
  theme(plot.title = element_blank())
p1 <- p1 + scale_color_discrete(name = "Clusters",
                                    labels = color_mapping,
                                    guide = guide_legend(override.aes = list(size = 3)))+
  facet_wrap(~sample, labeller = labeller(sample = setNames(label_df$label, label_df$sample)))+
  theme(plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text = element_text(size = 12, face = "bold", color = "black")
  )

# 移除点上的颜色图例指示（只保留图例框）
p1$layers[[1]]$aes_params$show.legend <- FALSE
print(p1)
# 保存图像
ggsave("Cluster/UMAP_clusters_split_by_sample.pdf", plot = p1, width = 9.5, height = 4)
ggsave("Cluster/UMAP_clusters_split_by_sample.png", plot = p1, width = 9.5, height = 4, dpi = 300)

# 绘制NC样本UMAP
# 提取 NC 子集
nc_obj <- subset(combined.int, subset = sample == "NC")
# UMAP 1：按 sample 分组，不按 cluster 分类
p_nc_sample <- DimPlot(nc_obj, group.by = "sample") +
  ggtitle("UMAP - NC (by Sample)") +
  theme(plot.title = element_text(hjust = 0.5))
# UMAP 2：按 seurat_clusters 分组（按 cluster 分类）
p_nc_cluster <- DimPlot(nc_obj, group.by = "seurat_clusters", label = TRUE, label.size = 3) +
  ggtitle("UMAP - NC (by Cluster)") +
  theme(plot.title = element_text(hjust = 0.5))
p_nc_merge <- p_nc_sample + p_nc_cluster + plot_layout(ncol = 2)
p_nc_merge
ggsave("Cluster/UMAP_NC.pdf", plot = p_nc_merge, width = 12, height = 5)
ggsave("Cluster/UMAP_NC.png", plot = p_nc_merge, width = 12, height = 5, dpi = 300)

# 绘制OE样本UMAP
# 提取 OE 子集
oe_obj <- subset(combined.int, subset = sample == "OE")
# UMAP 1：按 sample 分组，不按 cluster 分类
p_oe_sample <- DimPlot(oe_obj, group.by = "sample") +
  ggtitle("UMAP - OE (by Sample)") +
  theme(plot.title = element_text(hjust = 0.5))
# UMAP 2：按 seurat_clusters 分组（按 cluster 分类）
p_oe_cluster <- DimPlot(oe_obj, group.by = "seurat_clusters", label = TRUE, label.size = 3) +
  ggtitle("UMAP - OE (by Cluster)") +
  theme(plot.title = element_text(hjust = 0.5))
p_oe_merge <- p_oe_sample + p_oe_cluster + plot_layout(ncol = 2)
p_oe_merge
ggsave("Cluster/UMAP_OE.pdf", plot = p_oe_merge, width = 12, height = 5)
ggsave("Cluster/UMAP_OE.png", plot = p_nc_merge, width = 12, height = 5, dpi = 300)

# 设置 cluster 身份
Idents(combined.int) <- "seurat_clusters"
# 获取所有 cluster 编号
clusters <- levels(combined.int)
# 分别获取 NC 和 OE 的细胞名
nc_cells <- colnames(combined.int)[combined.int$sample == "NC"]
oe_cells <- colnames(combined.int)[combined.int$sample == "OE"]
# 提取子对象
nc_obj1 <- subset(combined.int, cells = nc_cells)
oe_obj1 <- subset(combined.int, cells = oe_cells)
# 开始循环每个 cluster
for (cl in clusters) {
  # 当前 cluster 的细胞
  cl_cells <- WhichCells(combined.int, idents = cl)
  # 只在当前样本中高亮的细胞
  nc_highlight <- intersect(cl_cells, nc_cells)
  oe_highlight <- intersect(cl_cells, oe_cells)
  # 分别绘图
  p_nc <- DimPlot(nc_obj1,
                  cells.highlight = nc_highlight,
                  cols = "lightgrey",
                  cols.highlight = "#E64B35FF",
                  pt.size = 1) +
    ggtitle(paste0("NC - Cluster ", cl)) +
    NoLegend()
  p_oe <- DimPlot(oe_obj1,
                  cells.highlight = oe_highlight,
                  cols = "lightgrey",
                  cols.highlight = "#4DBBD5FF",
                  pt.size = 1) +
    ggtitle(paste0("OE - Cluster ", cl)) +
    NoLegend()
  # 拼图 & 保存
  p_combined <- p_nc + p_oe + plot_layout(ncol = 2)
  ggsave(paste0("Cluster/Cluster", cl, "_Highlight_NC_OE.pdf"), plot = p_combined, width = 8, height = 4)
  ggsave(paste0("Cluster/Cluster", cl, "_Highlight_NC_OE.png"), plot = p_combined, width = 8, height = 4, dpi = 300)
}