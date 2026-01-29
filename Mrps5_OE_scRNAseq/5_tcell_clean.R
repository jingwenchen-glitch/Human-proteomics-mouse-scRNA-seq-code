# 清除系统环境变量，加载R包：
rm(list=ls())
options(stringsAsFactors = F) 
setwd("/Users/jingwenchen/Desktop/Ph.D/Jin Lab/scRNA-seq/Results/")

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
tcell_clean <- readRDS("tcell_clean.rds")

# 基于每个 cluster 的平均表达谱，计算 cluster 间的相关性矩阵，并绘制热图
library(Matrix)
library(ggplot2)
library(pheatmap)
# 1. 使用 AggregateExpression 获取 pseudo-bulk 表达矩阵（基于 cluster）
pseudo_bulk <- AggregateExpression(tcell_clean, group.by = "seurat_clusters", assay = "RNA", slot = "data")
expr_mat <- pseudo_bulk$RNA  # 提取 RNA assay 的表达矩阵
# 2. 确保是 matrix 格式
expr_mat <- as.matrix(expr_mat)
# 3. 计算相关性矩阵
cor_mat <- cor(expr_mat, method = "pearson")
# 去掉行名和列名中的 "g"
rownames(cor_mat) <- sub("^g", "", rownames(cor_mat))
colnames(cor_mat) <- sub("^g", "", colnames(cor_mat))
# 4. 绘制热图
# 设置色条范围与标签
breaks <- seq(-1, 1, length.out = 100)
colors <- colorRampPalette(c("darkred", "white", "darkblue"))(length(breaks) - 0.5)
# 保存为 PDF
# PDF
pdf("tcell_clean/Cluster_Correlation_Heatmap.pdf", width = 7, height = 6)
pheatmap(cor_mat,
         color = colors,
         breaks = breaks,
         display_numbers = TRUE,
         number_color = "white",
         fontsize_number = 8,
         main = "Cluster Correlation (Pearson)",
         border_color = NA,
         angle_col = 0)
dev.off()

# PNG
png("tcell_clean/Cluster_Correlation_Heatmap.png", width = 7, height = 6, units = "in", res = 300)
pheatmap(cor_mat,
         color = colors,
         breaks = breaks,
         display_numbers = TRUE,
         number_color = "white",
         fontsize_number = 8,
         main = "Cluster Correlation (Pearson)",
         border_color = NA,
         angle_col = 0)
dev.off()
dev.new()

# 统计每个 cluster 在不同 sample 中的细胞数量和百分比
# 提取 cluster 和 sample 信息
df <- data.frame(
  cluster = tcell_clean$seurat_clusters,
  sample = tcell_clean$sample
)
# 按 cluster 和 sample 计数
cluster_sample_stats <- df %>%
  group_by(cluster, sample) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(percent = round(100 * cell_count / sum(cell_count), 2))
write.csv(cluster_sample_stats, "tcell_clean/cluster_sample_stats_counts_and_percentage.csv")
# 可视化为堆叠柱状图
ggplot(cluster_sample_stats, aes(x = factor(cluster), y = percent, fill = sample)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_hline(yintercept = 0, color = "black", size = 0.8) +  # 在 y=0 加黑色轴线
  labs(x = "Cluster", y = "Percentage") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) + 
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold", color = "black")
  )
# 保存图像
ggsave("tcell_clean/Cluster_Sample_Composition_Barplot.pdf", width = 6, height = 4)
ggsave("tcell_clean/Cluster_Sample_Composition_Barplot.png", width = 6, height = 4, dpi = 300)

# 每个 sample 内部各 cluster 占比
df %>%
  group_by(sample, cluster) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(percent = round(100 * cell_count / sum(cell_count), 2)) -> sample_cluster_stats
write.csv(sample_cluster_stats, "tcell_clean/sample_cluster_stats_percentage.csv")
# 绘制堆积柱状图（cluster 占比，按 sample 分组）
p2 <- ggplot(sample_cluster_stats, aes(x = sample, y = percent, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, color = "white", size = 0.8) +  # 手动画轴线
  theme_minimal() +
  labs(y = "Percentage (%)") +
  scale_x_discrete(expand = expansion(add = 0.3))+
  scale_y_continuous(position = "right", expand = expansion(mult = c(0, 0.02))) +
  theme(
    panel.grid.major = element_blank(),    # 移除主网格线
    panel.grid.minor = element_blank(),    # 移除次网格线
    axis.line.y.right = element_line(size = 0.5, color = "white"),  # 显示右侧 y 轴线
    plot.title = element_text(hjust = 0.5),
    axis.ticks.y = element_line(size = 0.5, color = "black", lineend = "butt"),  # y 轴刻度线改为短横线
    axis.ticks.length.y = unit(3, "pt"),   # 控制刻度线长度
    axis.ticks.x = element_blank(),        # 去除 x 轴刻度线
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, angle = 0, face = "bold", color = "black", hjust = 0.5),
    axis.title.y = element_text(size = 14, angle = 0, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
 ) + NoLegend()
p2
ggsave("tcell_clean/Cluster_Percentage_By_Sample.pdf", p2, width = 4, height = 4)
ggsave("tcell_clean/Cluster_Percentage_By_Sample.png", p2, width = 4, height = 4, dpi = 300)

# 1. 计算聚类比例并修改图例标签
# 计算每个聚类的比例
clusters <- tcell_clean@meta.data$RNA_snn_res.0.175
cluster_props <- prop.table(table(clusters)) * 100
cluster_props <- round(cluster_props, 2)  # 保留2位小数
# 创建带比例的新标签
new_labels <- paste0(names(cluster_props), " (", cluster_props, "%)")
# 创建用于UMAP点标签的列（只包含cluster编号）
tcell_clean$cluster_id <- as.character(clusters)
# 创建用于图例的列（包含比例）
tcell_clean$cluster_with_prop <- factor(clusters, 
                                         levels = names(cluster_props),
                                         labels = new_labels)
# 绘制UMAP图
p <- DimPlot(tcell_clean, 
             reduction = "umap",
             group.by = "cluster_id",  # 使用只含编号的分组
             label = TRUE,             # 显示点标签（cluster编号）
             label.size = 6,           # 标签大小
             repel = TRUE) +           # 避免标签重叠
  theme(plot.title = element_blank())
p
# 修改图例标签为带比例的形式
# 获取当前的颜色映射
color_mapping <- levels(tcell_clean$cluster_with_prop)
names(color_mapping) <- levels(tcell_clean$cluster_id)

# 添加带比例的新图例
p <- p + scale_color_discrete(name = "Cluster",
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
ggsave("tcell_clean/UMAP_with_proportions.pdf", p, width = 6, height = 4)
ggsave("tcell_clean/UMAP_with_proportions.png", p, width = 6, height = 4, dpi = 300)

# UMAP 按 sample 着色
p_sample <- DimPlot(tcell_clean, group.by = "sample", label = FALSE) +
  theme(plot.title = element_blank(),
        legend.title = element_text(size = 14, face = "bold", color = "black"),
        legend.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text = element_text(size = 12, face = "bold", color = "black"))
# 按 cluster 分组绘制 UMAP
p_cluster <- DimPlot(tcell_clean, group.by = "seurat_clusters", label = TRUE, label.size = 6, repel = TRUE) +
  theme(plot.title = element_blank(),
        legend.title = element_text(size = 14, face = "bold", color = "black"),
        legend.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text = element_text(size = 12, face = "bold", color = "black")
  )
p_merge <- p_sample + p_cluster + plot_layout(ncol = 2)
p_merge
ggsave("tcell_clean/UMAP_sample_vs_cluster.pdf", p_merge, width = 8, height = 4)
ggsave("tcell_clean/UMAP_sample_vs_cluster.png", p_merge, width = 8, height = 4, dpi = 300)

# NC vs OE 按 cluster 分组绘制 UMAP
# 确保 cluster 设置为主身份
Idents(tcell_clean) <- "seurat_clusters"
cell_counts <- table(tcell_clean$orig.ident)
label_df <- data.frame(
  sample = names(cell_counts),
  count = as.numeric(cell_counts),
  label = paste0(names(cell_counts), "\n(", comma(as.numeric(cell_counts)), " cells)")
)
# 构建主图
p1 <- DimPlot(tcell_clean, group.by = "seurat_clusters", split.by = "sample", label = TRUE, label.size = 6) +
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
ggsave("tcell_clean/UMAP_clusters_split_by_sample.pdf", plot = p1, width = 9.5, height = 4)
ggsave("tcell_clean/UMAP_clusters_split_by_sample.png", plot = p1, width = 9.5, height = 4, dpi = 300)

saveRDS(tcell_clean, file = "tcell_clean_with_propotion.rds")

tcell_clean <- readRDS("tcell_clean_with_propotion.rds")

# 挖掘marker基因
up_markers <- FindAllMarkers(tcell_clean, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
# 保存为 CSV
write.csv(up_markers, "tcell_clean/up_markers.csv", row.names = FALSE) 
# 保存查找 marker 基因后的 Seurat 对象
saveRDS(up_markers, file = "up_markers.rds")

FeaturePlot(tcell_clean, features = c("Cd8a","Sell", "Tcf7" , "Ccr7", "Il7r", "Gzmb", "Gzmk",
                                      "Tox", "Tnf","Ifng","Il2", "Mki67"))

# dotplot
# 设置默认分组（如按 cluster 显示）
# 指定感兴趣基因
Idents(tcell_clean) <- "celltype"
genes_to_plot <- c("Cd3e", "Cd8a", 
                   "Lef1", "Sell", "Il7r", "Tcf7", "Foxo1", "Cxcr3",
                   "Mki67", "Top2a", "Pcna", "Cdc25b",
                   "Foxp3", "Il2ra", "Lrrc32", 
                   "Cd69", "Itgae", "Selplg")

DotPlot(tcell_clean, features = genes_to_plot) +
  scale_color_gradient(low = "white", high = "blue") +  # 表达量颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Gene", y = "Cluster")

# cluster命名
celltype <- c(
  "0" = "Tcm-like",
  "1" = "Proliferating_1", 
  "2" = "Treg", 
  "3" = "Trm-like", 
  "4" = "Proliferating_2"
)

# 应用到 Seurat 对象中：
# 设置当前 cluster 为 identity
tcell_clean <- SetIdent(tcell_clean, value = "seurat_clusters")

# 添加注释
tcell_clean$celltype <- plyr::mapvalues(
  tcell_clean$seurat_clusters,
  from = names(celltype),
  to = celltype
)
saveRDS(tcell_clean, file = "tcell_clean_celltype.rds")
tcell_clean <- readRDS("tcell_clean_celltype.rds")
# 1. 计算聚类比例并修改图例标签
# 计算每个聚类的比例
clusters <- tcell_clean@meta.data$RNA_snn_res.0.175
cluster_props <- prop.table(table(clusters)) * 100
cluster_props <- round(cluster_props, 2)  # 保留2位小数
# 创建带比例的新标签
new_labels <- paste0(names(cluster_props), " (", cluster_props, "%)")
# 创建用于UMAP点标签的列（只包含cluster编号）
tcell_clean$cluster_id <- as.character(clusters)
# 创建用于图例的列（包含比例）
tcell_clean$cluster_with_prop <- factor(clusters, 
                                        levels = names(cluster_props),
                                        labels = new_labels)
# 绘制UMAP图
p <- DimPlot(tcell_clean, 
             reduction = "umap",
             group.by = "cluster_id",  # 使用只含编号的分组
             label = TRUE,             # 显示点标签（cluster编号）
             label.size = 6,           # 标签大小
             repel = TRUE) +           # 避免标签重叠
  theme(plot.title = element_blank())
p
# 修改图例标签为带比例的形式
# 获取当前的颜色映射
color_mapping <- levels(tcell_clean$cluster_with_prop)
names(color_mapping) <- levels(tcell_clean$cluster_id)

# 添加带比例的新图例
p <- p + scale_color_discrete(name = "Cluster",
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
ggsave("tcell_clean/UMAP_with_proportions.pdf", p, width = 6, height = 4)
ggsave("tcell_clean/UMAP_with_proportions.png", p, width = 6, height = 4, dpi = 300)

# 生成 DotPlot 并提取数据
library(scales)
library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)

# 创建主点图
p <- DotPlot(tcell_clean, features = genes_to_plot)

# 提取 DotPlot 数据
dot_data <- p$data %>%
  mutate(
    gene = factor(features.plot, levels = genes_to_plot),                  # x 轴基因顺序
    cluster = factor(id, levels = rev(levels(tcell_clean)))               # y 轴 cluster 顺序
  )

# DotPlot 主图
# 确定 cluster 顺序（从上到下）
cluster_levels <- rev(levels(tcell_clean))  # 根据你的 DotPlot 顺序调整

# DotPlot 主图
p_main <- ggplot(dot_data, aes(x = gene, y = factor(id, levels = cluster_levels))) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient(low = "white", high = "blue") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, face = "bold.italic",color = "black", angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    axis.title = element_blank(),
    # 图例优化设置
    legend.position = "right",  # 图例放在右侧
    legend.box = "vertical",    # 垂直排列图例
    legend.box.just = "center",   # 图例左对齐
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # 减少边距
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(size = 12, face = "bold", color = "black"),
    legend.spacing = unit(1.0, "cm"),  # 图例项间距
    legend.key.height = unit(0.8, "cm"),  # 颜色图例高度
    legend.key.width = unit(0.5, "cm")    # 颜色图例宽度
  ) +
  labs(color = "log2FC",
       size = "Percent\nExpressed") +   # 添加换行符
  # 单独设置颜色图例
  guides(
    color = guide_colorbar(
      title.position = "top",   # 标题在顶部
      title.hjust = 0.5,        # 标题居中
      barwidth = unit(0.5, "cm"),  # 颜色条宽度
      barheight = unit(3, "cm")    # 颜色条高度（垂直）
    ),
    size = guide_legend(
      title.position = "top",   # 标题在顶部
      title.hjust = 0.5,        # 标题居中
      override.aes = list(color = "black")  # 确保点图例可见
    )
  )

# 左侧 cluster 注释条
cluster_colors <- hue_pal()(length(levels(tcell_clean)))  # 自动生成 cluster 颜色
names(cluster_colors) <- levels(tcell_clean)

p_anno <- ggplot(dot_data, aes(x = 1, y = factor(id, levels = cluster_levels), fill = id)) +
  geom_tile(color = "white", linewidth = 0.2, width = 0.2, height = 1) +
  scale_fill_manual(values = cluster_colors) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none",
        plot.margin = margin(r = 1, l = 0))
p_anno <- p_anno + coord_fixed(ratio = 1)  
# 水平组合（注释条在左，DotPlot 在右）
p_combined <- p_anno | p_main +
  plot_layout(widths = c(0.01, 0))  # 左侧注释条宽度比主图小

# 输出
print(p_combined)
ggsave("tcell_clean/DotPlot_with_annotation.pdf", width = 8, height = 4)
ggsave("tcell_clean/DotPlot_with_annotation.png", width = 8, height = 4, dpi = 300)


library(ggplot2)
library(scales)
library(grid)

# 1. 创建带比例的新标签
celltype_counts <- table(tcell_clean$celltype)
celltype_props <- round(prop.table(celltype_counts) * 100, 2)
new_labels <- paste0(names(celltype_props), " (", celltype_props, "%)")
tcell_clean$celltype_with_prop <- factor(tcell_clean$celltype,
                                         levels = names(celltype_props),
                                         labels = new_labels)
Idents(tcell_clean) <- "celltype_with_prop"

# 2. 准备每个 sample 的细胞数用于 facet_label
cell_counts <- table(tcell_clean$orig.ident)
label_df <- data.frame(
  sample = names(cell_counts),
  count = as.numeric(cell_counts),
  label = paste0(names(cell_counts), "\n(", comma(as.numeric(cell_counts)), " cells)")
)

# 3. 构建主图
p1 <- DimPlot(
  tcell_clean,
  split.by = "sample",          # 按 sample 分图
  label = FALSE,
  repel = TRUE
) +
  facet_wrap(~sample, labeller = labeller(sample = setNames(label_df$label, label_df$sample))) +
  scale_color_discrete(
    labels = setNames(levels(tcell_clean$celltype_with_prop), levels(tcell_clean$celltype)),
    guide = guide_legend(
      override.aes = list(size = 3),
      title.position = "top",
      title.theme = element_text(angle = 0, face = "bold", size = 14)
    )
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),       # 去除所有网格线
    axis.line = element_line(),         # 保留坐标轴线
    axis.ticks = element_line(),        # 保留刻度线
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(size = 14, face = "bold", color = "black"),
    legend.key.size = unit(1.0, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
    legend.box.spacing = unit(10, "pt"),
    axis.text.x = element_text(size = 12, face = "bold", color = "black", hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold", color = "black", hjust = 1),
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black")
  )

# 4. 显示
print(p1)

# 保存图像
ggsave("tcell_clean/Annotated_T_Cell_clean_by_sample.pdf", plot = p1, width = 9.5, height = 4)
ggsave("tcell_clean/Annotated_T_Cell_clean_by_sample.png", plot = p1, width = 9.5, height = 4, dpi = 300)

# 移除 p1 的 legend
p1 <- p1 + NoLegend()
p1
p3 <- p1 + p2 + 
  plot_layout(ncol = 2, widths = c(4, 1), guides = "collect") 
p3
ggsave("tcell_clean/Annotated_T_Cell_cluster_percentage.pdf", plot = p3, width = 7, height = 3.5)
ggsave("tcell_clean/Annotated_T_Cell_cluster_percentage.png", plot = p3, width = 7, height = 3.5, dpi = 300)

# 绘制UMAP图
p <- DimPlot(tcell_clean, 
             reduction = "umap",
             group.by = "celltype",  # 使用只含编号的分组
             label = FALSE,             # 显示点标签（cluster编号）
             label.size = 6,           # 标签大小
             repel = TRUE) +           # 避免标签重叠
  theme(plot.title = element_blank())
p <- DimPlot(tcell_clean,
             reduction = "umap",
             group.by = "celltype",
             label = FALSE)  # 不使用默认标签
p <- p + NoLegend()
p

# 提取坐标
coords <- Embeddings(tcell_clean, "umap") |> as.data.frame()
coords$celltype <- tcell_clean$celltype

# 取每个 celltype 的中心点（用于标签定位）
centers <- aggregate(coords[,1:2], by=list(coords$celltype), FUN=mean)
colnames(centers) <- c("celltype", "UMAP_1", "UMAP_2")

# 在原 DimPlot 上添加加粗标签
p + 
  ggrepel::geom_text_repel(
    data = centers,
    aes(x = UMAP_1, y = UMAP_2, label = celltype),
    fontface = "bold",
    size = 5
  )
# 修改图例标签为带比例的形式
# 获取当前的颜色映射
color_mapping <- levels(tcell_clean$celltype_with_prop)
names(color_mapping) <- levels(tcell_clean$celltype)

# 添加带比例的新图例
p <- p + scale_color_discrete(labels = color_mapping,
                              guide = guide_legend(override.aes = list(size = 3)))+       
  theme(plot.title = element_blank(),
        legend.key.size = unit(1.0, "cm"),
        legend.spacing.y = unit(0.3, "cm"),
        legend.box.spacing = unit(10, "pt"),
        legend.text = element_text(size = 14, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text = element_text(size = 12, face = "bold", color = "black")
  ) 

# 移除点上的颜色图例指示（只保留图例框）
p$layers[[1]]$aes_params$show.legend <- FALSE
print(p)
ggsave("tcell_clean/nolegend_UMAP.pdf", plot = p, width = 5, height = 4.5)
ggsave("tcell_clean/nolegend_UMAP.png", plot = p, width = 5, height = 4.5, dpi = 300)
ggsave("tcell_clean/UMAP_with_celltype_proportions.pdf", p, width = 7, height = 4.2)
ggsave("tcell_clean/UMAP_with_celltype_proportions.png", p, width = 7, height = 4.2, dpi = 300)

# 统计每个 cluster 在不同 sample 中的细胞数量和百分比
# 提取 cluster 和 sample 信息
df <- data.frame(
  cluster = tcell_clean$celltype,
  sample = tcell_clean$sample
)
# 按 cluster 和 sample 计数
cluster_sample_stats <- df %>%
  group_by(cluster, sample) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(percent = round(100 * cell_count / sum(cell_count), 2))
# 可视化为堆叠柱状图
ggplot(cluster_sample_stats, aes(x = factor(cluster), y = percent, fill = sample)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_hline(yintercept = 0, color = "black", size = 0.8) +
  labs(y = "Percentage") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = c("NC" = "#EA3323", "OE" = "#0000f5")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x = element_text(size = 14, face = "bold", color = "black", angle = 45, hjust = 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold", color = "black")
  )

# 保存图像
ggsave("tcell_clean/Celltype_Sample_Composition_Barplot.pdf", width = 6, height = 5)
ggsave("tcell_clean/Celltype_Sample_Composition_Barplot.png", width = 6, height = 5, dpi = 300)

## 感兴趣基因在不同样本的表达量
# 箱线图
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(scales)
# 感兴趣基因
genes_of_interest <- c("Mrps5",  
                       "Il7r", "Sell", "Tcf7", "Mki67", 
                       "Eomes", "Gzmb", "Prf1", 
                       "Pdcd1")

# 提取表达矩阵
expr_data <- FetchData(tcell_clean, vars = c("orig.ident", genes_of_interest)) %>%
  pivot_longer(cols = all_of(genes_of_interest),
               names_to = "Gene", values_to = "Expression") %>%
  dplyr::rename(Sample = orig.ident) #  %>%
# filter(Expression > 0)

# 将 Sample 设置为 factor
expr_data$Sample <- factor(expr_data$Sample, levels = c("NC", "OE"))
expr_data$Gene <- factor(expr_data$Gene, levels = genes_of_interest)

# 统计表达均值和表达比例
df_bar <- expr_data %>%
  group_by(Sample, Gene) %>%
  summarise(
    avg.exp = mean(Expression),
    #  avg.exp = ifelse(sum(Expression > 0) == 0, 0, mean(Expression[Expression > 0])),
    pct.exp = round(mean(Expression > 0) * 100, 2),
    .groups = "drop"
  )

# 计算每个基因的最大表达量，用于设置y轴上限
max_exp <- df_bar %>%
  group_by(Gene) %>%
  summarise(max_exp = max(avg.exp) * 1.1)  

# 合并到主数据
df_bar <- df_bar %>%
  left_join(max_exp, by = "Gene")

# 统计每个基因的 p 值
df_pval <- expr_data %>%
  dplyr::filter(Expression > 0) %>%
  group_by(Gene) %>%
  summarise(
    p_val = tryCatch(
      wilcox.test(Expression ~ Sample, exact = FALSE)$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_label = case_when(
      is.na(p_val) ~ "NA",
      p_val < 0.0001 ~ "****",
      p_val < 0.001 ~ "***",
      p_val < 0.01 ~ "**",
      p_val < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# 合并进 df_bar
df_bar <- df_bar %>%
  left_join(df_pval, by = "Gene")

# 重新计算标签位置数据框，确保包含 max_exp 列
df_label_pos <- df_bar %>%
  dplyr::filter(Sample == "OE") %>%
  dplyr::select(Gene, max_exp, p_label)  # 确保包含 max_exp

# 创建绘图
ggplot(df_bar, aes(x = Sample, y = avg.exp, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.6) +
  #  geom_text(aes(label = paste0(pct.exp, "%")), vjust = -0.5,size = 3) +
  facet_wrap(~ Gene, scales = "free_y", nrow = 2) +
  # 修复：确保使用数值计算位置
  geom_text(data = df_label_pos, 
            aes(x = 1.5, y = max_exp * 0.91, label = p_label),  # x=1.5表示两样本中间位置
            inherit.aes = FALSE, size = 7, fontface = "bold", color = "black") +
  scale_fill_manual(values = c("NC" = "#EA3323", "OE" = "#0000F5")) +
  theme_classic() +
  labs(y = "Average Expression") +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05)),
    labels = label_number(accuracy = 0.01)   # 保留两位小数
  ) + 
  theme(
    strip.text = element_text(size = 14, face = "bold.italic", color = "black"),
    axis.text.x = element_text(size = 14, face = "bold", angle = 0, hjust = 0.5, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black", vjust = 2.0),
    legend.position = "none",
    panel.spacing = unit(0.2, "lines")  # 最小分面间距
  )
ggsave("tcell_clean/interested_genes_bar_plot.pdf", width = 8, height = 5)
ggsave("tcell_clean/interested_genes_bar_plot.png", width = 8, height = 5, dpi = 300)

# 设置 cluster 为分组变量
Idents(tcell_clean) <- "seurat_clusters"

# 提取表达矩阵（含 cluster 和基因表达）
expr_data_1 <- FetchData(tcell_clean, vars = c("seurat_clusters", genes_of_interest)) %>%
  pivot_longer(cols = all_of(genes_of_interest),
               names_to = "Gene", values_to = "Expression") %>%
  dplyr::rename(Cluster = seurat_clusters)

# 将 Cluster 设置为 factor，按顺序排序
expr_data_1$Cluster <- factor(expr_data_1$Cluster, levels = sort(unique(expr_data_1$Cluster)))

# 绘图：小提琴图 + jitter
library(ggplot2)
ggplot(expr_data_1, aes(x = Cluster, y = Expression, fill = Cluster)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_jitter(width = 0.2, size = 0.1, alpha = 0.2) +
  facet_wrap(~ Gene, scales = "free_y", nrow = 2) +
  theme_classic() +
  labs(y = "Expression level", fill = "Cluster") +
  theme(
    strip.text = element_text(size = 14),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.title.x = element_blank()
  ) + NoLegend()
ggsave("Plot/interested_genes_violin_plot_cluster.pdf", width = 10, height = 6)
ggsave("Plot/interested_genes_violin_plot_cluster.png", width = 10, height = 6, dpi = 300)
