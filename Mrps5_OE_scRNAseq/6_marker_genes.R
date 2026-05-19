# 清除系统环境变量，加载R包：
rm(list=ls())
options(stringsAsFactors = F) 
setwd("/Users/jingwenchen/Desktop/Ph.D/Jin Lab/scRNA-seq/Results/")

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)  # 小鼠用 org.Mm.eg.db，人类用 org.Hs.eg.db

# 读取文件
tcell_clean <- readRDS("tcell_clean_celltype.rds")

# 挖掘marker基因
up_markers <- FindAllMarkers(tcell_clean, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
# 保存为 CSV
write.csv(up_markers, "Cluster_markers/up_markers.csv", row.names = FALSE) 
# 保存查找 marker 基因后的 Seurat 对象
saveRDS(up_markers, file = "up_markers.rds")

# 读取文件
cluster_markers <- readRDS("up_markers.rds")

# 查看每个cluster前几个marker
cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# 查看 top marker：
top10 <- cluster_markers %>% filter(p_val <= 0.01, p_val_adj < 0.05, avg_log2FC >= 0.26) %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
View(top10)
write.csv(top10, "Cluster_markers/top10_Upreulated_markers.csv") 
# 进一步过滤
cluster_markers.filtered <- cluster_markers %>%
  filter(p_val <= 0.01, avg_log2FC >= 0.26)

# 统计每个 cluster 上调 marker 基因数量
cluster_markers_counts <- cluster_markers.filtered %>%
  group_by(cluster) %>%
  summarise(cluster_markers_counts = n()) %>%
  arrange(as.numeric(as.character(cluster)))
# 保存为CSV
write.csv(cluster_markers.filtered, file = "Cluster_markers/Filtered_MarkerGene_List.csv", row.names = FALSE)
write.csv(cluster_markers_counts, file = "Cluster_markers/Upregulated_Markers_Barplot_colored.csv", row.names = FALSE)

# 绘图
# 各亚群上调基因数量统计柱状图
library(ggplot2)
library(dplyr)
# 获取cluster levels（顺序）
cluster_levels <- levels(Idents(tcell_clean))
# 生成DimPlot，提取颜色和对应cluster（通过ggplot_build）
p <- DimPlot(tcell_clean, group.by = "seurat_clusters", label = FALSE)
p
plot_data <- ggplot_build(p)$data[[1]]
# ggplot_build提取的颜色顺序，是按照点顺序的，不一定和cluster_levels顺序一一对应
# 用 cluster_levels 生成对应颜色，用 scales::hue_pal() 生成颜色
library(scales)
cluster_colors <- hue_pal()(length(cluster_levels))
names(cluster_colors) <- cluster_levels
# 确保 cluster_markers_counts$cluster 是因子且levels和cluster_levels一致
cluster_markers_counts$cluster <- factor(cluster_markers_counts$cluster, levels = cluster_levels)
# 绘图，使用手动颜色映射
ggplot(cluster_markers_counts, aes(x = cluster, y = cluster_markers_counts, fill = cluster)) +
  geom_col() +
  scale_fill_manual(values = cluster_colors) +
  labs(x = "Cluster", y = "Number of Marker genes") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
# 保存为pdf, png
ggsave("Cluster_markers/Upregulated_Markers_Barplot_colored.pdf", width = 6, height = 4)
ggsave("Cluster_markers/Upregulated_Markers_Barplot_colored.png", width = 6, height = 4, dpi = 300)

# 绘制火山图
# 添加分组标签
library(ggplot2)
library(dplyr)
library(ggrepel)
plot_cluster_volcano <- function(cluster_markers.filtered, cluster_colors,
                                 logfc_cutoff = 0.26,
                                 adj_pval_cutoff = 0.05,
                                 output_prefix = "Volcano_Cluster_ColorBand_Filled",
                                 width = 6, height = 4, dpi = 300) {
  # 处理 cluster 顺序
  cluster_levels_num <- sort(as.numeric(as.character(unique(cluster_markers.filtered$cluster))))
  cluster_levels <- as.character(cluster_levels_num)
  cluster_markers.filtered$cluster <- factor(cluster_markers.filtered$cluster, levels = cluster_levels)
  # 标记 p 值组
  volcano_df <- cluster_markers.filtered %>%
    mutate(pval_group = ifelse(p_val_adj < adj_pval_cutoff, "adjust Pvalue < 0.01", "adjust Pvalue >= 0.01"))%>%
    mutate(pval_group = factor(pval_group, levels = c("adjust Pvalue < 0.01", "adjust Pvalue >= 0.01")))
  # Top5 上调基因（每个 cluster）
  top_genes <- volcano_df %>%
    filter(pval_group == "adjust Pvalue < 0.01") %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 5) %>%
    ungroup()
  # 底部颜色条
  colorbar_df <- data.frame(
    cluster = factor(names(cluster_colors), levels = cluster_levels),
    fill_color = cluster_colors[cluster_levels]
  )
  ymin_fill <- -logfc_cutoff
  ymax_fill <- logfc_cutoff
  # 绘图
  p <- ggplot() +
    geom_rect(data = colorbar_df,
              aes(xmin = as.numeric(cluster) - 0.45,
                  xmax = as.numeric(cluster) + 0.45),
              ymin = ymin_fill, ymax = ymax_fill,
              fill = colorbar_df$fill_color,
              inherit.aes = FALSE) +
    # 添加 jitter 点（按显著性分组着色）
    geom_jitter(data = volcano_df,
                aes(x = cluster, y = avg_log2FC, color = pval_group),
                width = 0.35, size = 0.6, alpha = 0.7) +
    geom_text_repel(data = top_genes,
                    aes(x = cluster, y = avg_log2FC, label = gene),
                    size = 3, max.overlaps = Inf, segment.color = "grey50") +
    scale_color_manual(values = c("adjust Pvalue < 0.01" = "red", "adjust Pvalue >= 0.01" = "blue")) +
    theme_minimal(base_size = 14) +
    labs(title = "log2FC of Marker Genes by Cluster",
         x = "Cluster", y = "log2 Fold Change") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          legend.title = element_blank())
  # 保存图像（PDF & PNG）
  ggsave(paste0(output_prefix, ".pdf"), plot = p, width = width, height = height)
  ggsave(paste0(output_prefix, ".png"), plot = p, width = width, height = height, dpi = dpi)
  message("Volcano plot saved to: ", output_prefix, ".pdf / .png")
  return(p)
}
plot_cluster_volcano(
  cluster_markers.filtered = cluster_markers.filtered,
  cluster_colors = cluster_colors,
  output_prefix = "Cluster_markers/Volcano_MarkerCluster"
)

# 绘制差异倍数 vs 表达比例差值的火山图
library(ggplot2)
library(dplyr)
library(ggrepel)
plot_pctdiff_logfc <- function(cluster_markers.filtered,
                               cluster_colors,
                               output_prefix = "Volcano_pctFC_Top5_Colored",
                               width = 10, height = 4, dpi = 300) {
  # 添加表达比例差值列
  cluster_markers.filtered <- cluster_markers.filtered %>%
    mutate(pct.diff = pct.1 - pct.2) %>%
    filter(avg_log2FC >= 0.26, p_val_adj < 0.05)  # 只展示上调基因
  # 找出每个 cluster log2FC 最大的前 5 个基因
  top_genes <- cluster_markers.filtered %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 5) %>%
    mutate(highlight = TRUE) %>%
    ungroup()
  # 添加 segment 颜色列
  top_genes$segment_color <- cluster_colors[as.character(top_genes$cluster)]
  # 合并到 volcano_df
  volcano_df <- cluster_markers.filtered %>%
    left_join(top_genes %>% dplyr::select(gene, cluster, highlight, segment_color), 
              by = c("gene", "cluster")) %>%
    mutate(highlight = ifelse(is.na(highlight), FALSE, highlight),
           cluster = factor(cluster, levels = names(cluster_colors)))
  # 计算统一坐标轴范围
  x_range <- range(volcano_df$pct.diff, na.rm = TRUE)
  y_range <- range(volcano_df$avg_log2FC, na.rm = TRUE)
  # 绘图
  p <- ggplot(volcano_df, aes(x = pct.diff, y = avg_log2FC)) +
    # 灰色背景点
    geom_point(data = subset(volcano_df, !highlight),
               color = "lightgrey", size = 0.5, alpha = 0.5) +
    # 高亮点
    geom_point(data = subset(volcano_df, highlight),
               aes(color = cluster), size = 1.2, alpha = 0.9) +
    # gene label
    geom_text_repel(data = subset(volcano_df, highlight),
                    aes(label = gene, color = cluster),
                    segment.color = subset(volcano_df, highlight)$segment_color,
                    size = 3,
                    max.overlaps = Inf,
                    box.padding = 0.5,
                    point.padding = 0.3,
                    min.segment.length = 0) +
    scale_color_manual(values = cluster_colors) +
    coord_cartesian(xlim = x_range, ylim = y_range) +
    facet_wrap(~ cluster, scales = "fixed", ncol = 7) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Expression Proportion Difference",
      y = "Average log2 Fold Change",
      color = "Cluster") +
    NoLegend()
  # 保存
  ggsave(paste0(output_prefix, ".pdf"), plot = p, width = width, height = height)
  ggsave(paste0(output_prefix, ".png"), plot = p, width = width, height = height, dpi = dpi)
  message("Volcano pct plot saved to: ", output_prefix, ".pdf / .png")
  return(p)
}
plot_pctdiff_logfc(cluster_markers.filtered = cluster_markers.filtered,
                   cluster_colors = cluster_colors,
                   output_prefix = "Cluster_markers/Volcano_pctFC_Top5_Colored")

# 展示每个细胞群体的marker基因，在所有细胞群体的细胞中对应的平均表达量表格
# 提取所有 marker 基因名
library(Seurat)
library(dplyr)
# 所有 marker 基因（只保留正向差异的）
all_marker_genes <- cluster_markers.filtered %>%
  filter(p_val <= 0.01, avg_log2FC >= 0.26) %>%
  pull(gene) %>%
  unique()
# 计算每个cluster中这些基因的平均表达量
avg_expr_all <- AverageExpression(
  tcell_clean,
  features = all_marker_genes,
  group.by = "seurat_clusters",
  slot = "data"  # log-normalized expression
)$RNA
# avg_expr_all 是矩阵，行为基因，列为cluster
head(avg_expr_all)
# 保存为 CSV 文件
write.csv(avg_expr_all, file = "Cluster_markers/up_MarkerGenes_AverageExpression.csv", row.names = TRUE)

# 绘制热图
library(Seurat)
library(dplyr)
library(ggplot2)
# 筛选每个 cluster 的 top10 marker 基因（默认按 log2FC）
top10_markers <- cluster_markers.filtered %>%
  filter(avg_log2FC >= 0.26, p_val <= 0.01, p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  pull(gene) %>%
  unique()

# 绘制热图
DoHeatmap(tcell_clean,
          features = top10_markers,
          group.by = "seurat_clusters") +
  theme(axis.text.y = element_text(size = 6))  # 可根据需要调字体大小
ggsave("Cluster_markers/Top10_MarkerGenes_Heatmap.pdf", width = 8, height = 10)
ggsave("Cluster_markers/Top10_MarkerGenes_Heatmap.png", width = 8, height = 10, dpi = 300)

# 使用 pheatmap() 绘制 Marker 基因热图
library(dplyr)
library(scales)

top10_markers <- cluster_markers.filtered %>%
  filter(avg_log2FC >= 0.26, p_val <= 0.01, p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  pull(gene) %>%
  unique()

# 提取表达矩阵（log-normalized）
expr_mat <- GetAssayData(tcell_clean, slot = "data")[top10_markers, ]
# 获取每个细胞所属 cluster，并排序
cell_clusters <- Idents(tcell_clean)
cell_clusters <- factor(cell_clusters, levels = sort(as.numeric(levels(cell_clusters))))
# 按 cluster 对细胞列排序
sorted_cells <- names(sort(cell_clusters))
expr_mat <- expr_mat[, sorted_cells]
# 构建列注释
annotation_col <- data.frame(Cluster = cell_clusters[sorted_cells])
rownames(annotation_col) <- sorted_cells
# 设置自定义颜色：淡蓝 → 白 → 淡红
# custom_colors <- colorRampPalette(c("#ADD8E6", "white", "#F4A7B9"))(100)
custom_colors <- colorRampPalette(c("#EB52F7", "black", "yellow"))(100)
# 设置颜色分布范围（-2 到 2 之间线性映射）
breaks <- seq(-2, 2, length.out = 101)  # 必须比 colors 多 1 个值
# Cluster 注释颜色（可自定义或用 Seurat 提取的 cluster_colors）
cluster_palette <- hue_pal()(length(levels(cell_clusters)))
names(cluster_palette) <- levels(cell_clusters)
ann_colors <- list(Cluster = cluster_palette)
# 绘制热图
library(pheatmap)
# PDF 保存
pdf("Cluster_markers/Top10_MarkerGenes_Pheatmap.pdf", width = 8, height = 6)
pheatmap(expr_mat,
         scale = "row",  
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         show_rownames = TRUE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color = custom_colors,
         breaks = breaks,
         fontsize_row = 6,
         border_color = NA)
dev.off()

# PNG 保存
png("Cluster_markers/Top10_MarkerGenes_Pheatmap.png", width = 8, height = 6, units = "in", res = 300)
pheatmap(expr_mat,
         scale = "row",  
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         show_rownames = TRUE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color = custom_colors,
         breaks = breaks,
         fontsize_row = 6,
         border_color = NA)
dev.off()
dev.new()

# 细胞群体Top10 Marker基因的小提琴图
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
# 设置输出目录
output_dir <- "Cluster_markers/Violin_Plots"
dir.create(output_dir, showWarnings = FALSE)
clusters <- unique(cluster_markers.filtered$cluster)
for (cl in clusters) {
  # 取Top10基因
  genes <- cluster_markers.filtered %>%
    filter(cluster == cl, p_val_adj < 0.05) %>%
    slice_max(order_by = avg_log2FC, n = 10) %>%
    pull(gene) %>%
    head(10)
  # 单个基因小提琴图集合
  single_gene_plots <- lapply(genes, function(gene) {
    VlnPlot(
      tcell_clean,
      features = gene,
      group.by = "seurat_clusters",
      pt.size = 0.5
    ) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = gene) +
      NoLegend()
  })
  # 2行5列排列
  p_all <- wrap_plots(single_gene_plots, ncol = 5)
  # 保存到指定文件夹中
  ggsave(filename = file.path(output_dir, paste0("Cluster", cl, "_Top10_Marker_Violin.pdf")), plot = p_all, width = 14,height = 8)
  ggsave(filename = file.path(output_dir, paste0("Cluster", cl, "_Top10_Marker_Violin.png")), plot = p_all, width = 14,height = 8, dpi = 300)
}

# UMAP展示每个cluster的Top10 Marker基因
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
plot_cluster_marker_umap <- function(seurat_obj,
                                     cluster_markers.filtered,
                                     top_n = 10,
                                     output_dir = "Cluster_markers/UMAPs",
                                     ncol = 5) {
  dir.create(output_dir, showWarnings = FALSE)
  # 确保log2FC >= 0.26
  cluster_markers.filtered <- cluster_markers.filtered %>% filter(avg_log2FC >= 0.26)
  # 提取 cluster 列表
  clusters <- unique(cluster_markers.filtered$cluster)
  for (clust in clusters) {
    # 当前 cluster 的 top N 基因
    top_genes <- cluster_markers.filtered %>%
      filter(cluster == clust, p_val_adj < 0.05) %>%
      slice_max(order_by = avg_log2FC, n = top_n) %>%
      pull(gene)%>%
      head(10)
    # 生成 FeaturePlot 列表
    gene_plots <- lapply(top_genes, function(gene) {
      FeaturePlot(seurat_obj, features = gene, reduction = "umap") +
        ggtitle(gene) +
        theme_minimal(base_size = 10) +
        theme(plot.title = element_text(hjust = 0.5))
    })
    # 拼图并保存为 PDF
    p_combined <- wrap_plots(gene_plots, ncol = ncol)
    ggsave(filename = file.path(output_dir, paste0("Cluster_", clust, "_Markers.pdf")), plot = p_combined, width = 3 * ncol, height = 2 * ceiling(length(top_genes) / ncol))
    ggsave(filename = file.path(output_dir, paste0("Cluster_", clust, "_Markers.png")), plot = p_combined, width = 3 * ncol, height = 2 * ceiling(length(top_genes) / ncol), dpi = 300)
  }
}
plot_cluster_marker_umap(
  seurat_obj = tcell_clean,
  cluster_markers.filtered = cluster_markers.filtered,
  top_n = 10,
  output_dir = "Cluster_markers/UMAPs",
  ncol = 5  # 每行5图
)

# 每个 cluster 的 Top10 marker 基因气泡图
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
# 创建输出文件夹
output_dir <- "Cluster_markers/Bubble_Plots"
dir.create(output_dir, showWarnings = FALSE)
clusters <- unique(cluster_markers.filtered$cluster)
for (cl in clusters) {
  # 获取 Top10 marker 基因
  top_genes <- cluster_markers.filtered %>%
    filter(cluster == cl, p_val_adj < 0.05) %>%
    slice_max(order_by = avg_log2FC, n = 10) %>%
    pull(gene) %>%
    head(10)
  # 获取表达矩阵（log-normalized）
  expr <- FetchData(tcell_clean, vars = top_genes)
  # 获取 cluster 信息
  meta <- data.frame(
    cluster = Idents(tcell_clean),
    cell = colnames(tcell_clean)
  )
  # 合并表达和 meta 数据，计算所有 cluster 的表达水平
  plot_data <- meta %>%
    left_join(expr %>% mutate(cell = rownames(expr)), by = "cell") %>%
    pivot_longer(cols = all_of(top_genes), names_to = "gene", values_to = "expression") %>%
    group_by(cluster, gene) %>%
    summarise(
      avg_exp = mean(expression),
      pct_exp = mean(expression > 0),
      .groups = "drop"
    )
  # 强制 cluster 为因子，统一 y 轴顺序
  plot_data$cluster <- factor(plot_data$cluster, levels = sort(unique(plot_data$cluster)))
  # 气泡图
  p <- ggplot(plot_data, aes(x = gene, y = cluster)) +
    geom_point(aes(size = pct_exp, color = avg_exp)) +
    scale_size(range = c(2, 8)) +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = paste0("Cluster ", cl, " Top10 Marker Genes (All Clusters)"),
      x = "Gene",
      y = "Cluster",
      size = "Percent Expressed",
      color = "Average Expression"
    )
  # 保存图像
  ggsave(filename = file.path(output_dir, paste0("Cluster", cl, "_Top10_Marker_Bubble_AllClusters.pdf")), plot = p, width = 10, height = 6)
  ggsave(filename = file.path(output_dir, paste0("Cluster", cl, "_Top10_Marker_Bubble_AllClusters.png")), plot = p, width = 10, height = 6, dpi = 300)
}

# 每个 cluster 的 Top10 marker 基因ridge图
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggridges)
library(tidyr)
plot_cluster_ridge <- function(seurat_obj, cluster_markers.filtered, cluster_colors, output_dir = "Cluster_markers/RidgePlots", width = 8, height = 10, dpi = 300) {
  dir.create(output_dir, showWarnings = FALSE)
  clusters <- unique(cluster_markers.filtered$cluster)
  for (cl in clusters) {
    # 取当前cluster的Top10基因
    genes <- cluster_markers.filtered %>%
      filter(cluster == cl, p_val_adj < 0.05) %>%
      slice_max(order_by = avg_log2FC, n = 10) %>%
      pull(gene) %>%
      head(10)
    # 获取所有细胞这些基因的表达矩阵
    expr_df <- FetchData(seurat_obj, vars = genes)
    expr_df$cell <- rownames(expr_df)
    # 获取细胞所属cluster信息
    meta <- data.frame(
      cell = colnames(seurat_obj),
      cluster = as.character(Idents(seurat_obj))
    )
    # 合并表达数据和meta信息
    plot_data <- meta %>%
      left_join(expr_df, by = "cell") %>%
      pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression")
    # cluster转因子，确保顺序和颜色对应
    plot_data$cluster <- factor(plot_data$cluster, levels = names(cluster_colors))
    # 画山脊图，x轴表达量，y轴cluster，fill按cluster填充
    p <- ggplot(plot_data, aes(x = expression, y = cluster, fill = cluster)) +
      geom_density_ridges(alpha = 0.8, scale = 3) +
      scale_fill_manual(values = cluster_colors) +
      theme_minimal(base_size = 14) +
      labs(title = paste0("Cluster ", cl, " Top10 Marker Genes Ridge Plot"),
           x = "Expression",
           y = "Cell Cluster") +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      facet_wrap(~ gene, nrow = 2, scales = "free_y")
    
    # 保存PDF
    ggsave(filename = file.path(output_dir, paste0("Cluster", cl, "_Top10_Marker_RidgePlot.pdf")),plot = p, width = width, height = height)
    ggsave(filename = file.path(output_dir, paste0("Cluster", cl, "_Top10_Marker_RidgePlot.png")),plot = p, width = width, height = height, dpi = dpi)
  }
}
# 运行
plot_cluster_ridge(seurat_obj = tcell_clean, cluster_markers.filtered = cluster_markers.filtered, cluster_colors = cluster_colors)
