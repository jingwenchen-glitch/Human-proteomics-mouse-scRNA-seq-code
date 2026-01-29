# 清除系统环境变量，加载R包：
rm(list=ls())
options(stringsAsFactors = F) 
setwd("/Users/jingwenchen/Desktop/Ph.D/Jin Lab/Proteomics/Y_O")

# 安装并加载所需的R包
library(impute)
library(tidyverse)
library(limma)
library(readr)
library(dplyr)
library(tibble)

# 读取数据 --------------------------------------------------------
df <- read_tsv("MS_identified_information.txt")

# 过滤策略 --------------------------------------------------------
# 判断是否存在你指定的列（避免列名错误）
names(df)

# 核心过滤逻辑（修复括号错误）
df_filtered <- df %>%
  filter(
    `Unique peptides` >= 2,  # 至少2个唯一肽段
    # 每组至少2个样本有值（Intensity A/B 列）
    rowSums(!is.na(dplyr::select(., starts_with("A_")))) >= 2,
    rowSums(!is.na(dplyr::select(., starts_with("B_")))) >= 2
  )

# 3. 提取定量表达矩阵 ------------------------------------------------
quant_data <- df_filtered %>%
  dplyr::select(`Gene name`, A_1, A_2, A_3, B_1, B_2, B_3)

# 转换为矩阵格式（用于limma）
expr_mat <- quant_data %>%
  column_to_rownames("Gene name") %>%
  as.matrix()
# 剩余缺失值用kNN填补：
expr_mat_imputed <- impute.knn(expr_mat)$data

# log2转换 ---------------------------------------------------------
log2_expr <- log2(expr_mat_imputed)

# 绘制箱线图 -------------------------------------------------------
boxplot(log2_expr,
        col = rep(c("blue", "red"), each = 3),
        main = "Boxplot of log2 Intensities",
        xlab = "Samples", ylab = "log2(Intensity)")

# 差异分析（使用limma - 适合小样本）
# 创建样本分组信息
group <- factor(c("Y", "Y", "Y", "O", "O", "O"))

# 构建设计矩阵
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# 创建 contrasts（比较O vs Y）
# 指定比较组
contrast_matrix <- makeContrasts(O_vs_Y = O - Y, levels = design)

# 运行 limma 差异分析
# 建立线性模型并拟合
fit <- lmFit(log2_expr, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 提取差异结果
diff_result <- topTable(fit2, coef = "O_vs_Y", number = Inf, sort.by = "P")%>%
  rownames_to_column("Gene") %>%
  dplyr::select(Gene, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(log2FC = logFC,
         p.value = P.Value,
         adj.p.val = adj.P.Val)

# 标记显著差异蛋白
diff_result <- diff_result %>%
  mutate(regulation = case_when(
    adj.p.val < 0.05 & log2FC > 0.58 ~ "Up",
    adj.p.val < 0.05 & log2FC < -0.58 ~ "Down",
    TRUE ~ "NotSig"
  ))

# 设置阈值（logFC > 0.58 & adj.p.val < 0.05）
deg_genes <- diff_result %>%
  filter(abs(log2FC) > 0.58, adj.p.val < 0.05)

# 保存差异结果
write.csv(diff_result, "Protein_DE_results.csv", row.names = FALSE)
write.csv(deg_genes, "Protein_DE_filtered.csv", row.names = FALSE)

# 计算 replicate-null 分布（同组内部差异）
library(tidyverse)
# helper: 计算组内所有两两差值
pairwise_diffs <- function(mat, samples){
  combn(samples, 2, function(x){
    mat[, x[1]] - mat[, x[2]]
  }) %>% as.data.frame()
}

# 组内差异（null distribution）
Y_samples <- colnames(log2_expr)[1:3]
O_samples <- colnames(log2_expr)[4:6]

Y_null <- pairwise_diffs(log2_expr, Y_samples) %>% as.matrix() %>% as.vector()
O_null <- pairwise_diffs(log2_expr, O_samples) %>% as.matrix() %>% as.vector()
null_vec <- c(Y_null, O_null)

# 计算组间差异（between-group）
cross_diffs <- combn(c(Y_samples, O_samples), 2, simplify = FALSE) %>%
  keep(~ (.x[1] %in% Y_samples & .x[2] %in% O_samples) |
         (.x[2] %in% Y_samples & .x[1] %in% O_samples)) %>%
  map(~ as.numeric(log2_expr[, .x[1]] - log2_expr[, .x[2]])) %>%
  unlist()

# 合并并可视化分布
df_plot <- tibble(
  value = c(null_vec, cross_diffs),
  type = c(rep("null comparison", length(null_vec)),
           rep("true comparison", length(cross_diffs)))
)

ggplot(df_plot, aes(x = value, color = type)) +
  geom_density(size = 1.5) +
  xlab("log2 difference") + ylab("Density") +
  scale_color_manual(values = c("null comparison" = "blue", 
                                "true comparison" = "red")) +
  theme_minimal() +
  theme(
    legend.position = "top",   # 图例放在顶部
    legend.direction = "horizontal",  # 横向排列
    legend.key.height = unit(0.2, "cm"),
    legend.key.width  = unit(0.5, "cm"),  # 横向多给一点宽度
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    axis.text = element_text(size = 12, face = "bold", color = "black")
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))   # 横向排成一行

ggsave("density_comparison.pdf", width = 5, height = 4)
ggsave("density_comparison.png", width = 5, height = 4, dpi = 300)

ks.test(null_vec, cross_diffs)  # Kolmogorov-Smirnov test

library(ggpointdensity)
library(ggplot2)
library(hexbin)
library(patchwork)
# 把所有比较组合放入列表
pairs <- list(
  c("A_1","A_2"),
  c("B_1","B_2"),
  c("A_1","B_1")
)
# 找每个 hexbin 的最大计数
max_count <- max(sapply(pairs, function(p){
  h <- hexbin(log2_expr[,p[1]], log2_expr[,p[2]], xbins=100)
  max(h@count)
}))
# 统一 scale_fill_viridis_c 的 limits
cor_val <- cor(log2_expr[, "A_1"], log2_expr[, "A_2"], use = "pairwise.complete.obs")
p1 <- ggplot(as.data.frame(log2_expr), aes(x = .data[["A_1"]], y = .data[["A_2"]])) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(option = "plasma", limits=c(1, max_count)) +
  annotate("text", x = min(log2_expr[, "A_1"]), 
           y = max(log2_expr[, "A_2"]), 
           label = paste0("r = ", round(cor_val, 2)),
           hjust = 0, vjust = 0, size = 5, color = "black", fontface = "bold") +
  theme_classic() +
  labs(x = "Y_1 log2FC", y = "Y_2 log2FC") + 
  theme(plot.title = element_blank(),
        legend.title = element_text(size = 14, face = "bold", color = "black"),
        legend.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text = element_text(size = 12, face = "bold", color = "black"))

cor_val <- cor(log2_expr[, "B_1"], log2_expr[, "B_2"], use = "pairwise.complete.obs")
p2 <- ggplot(as.data.frame(log2_expr), aes(x = .data[["B_1"]], y = .data[["B_2"]])) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(option = "plasma", limits=c(1, max_count)) +
  annotate("text", x = min(log2_expr[, "B_1"]), 
           y = max(log2_expr[, "B_2"]), 
           label = paste0("r = ", round(cor_val, 2)),
           hjust = 0, vjust = 0, size = 5, color = "black", fontface = "bold") +
  theme_classic() +
  labs(x = "O_1 log2FC", y = "O_2 log2FC") +
  theme(plot.title = element_blank(),
        legend.title = element_text(size = 14, face = "bold", color = "black"),
        legend.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text = element_text(size = 12, face = "bold", color = "black"))

cor_val <- cor(log2_expr[, "A_1"], log2_expr[, "B_1"], use = "pairwise.complete.obs")
p3 <- ggplot(as.data.frame(log2_expr), aes(x = .data[["A_1"]], y = .data[["B_1"]])) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(option = "plasma", limits=c(1, max_count)) +
  annotate("text", x = min(log2_expr[, "A_1"]), 
           y = max(log2_expr[, "B_1"]), 
           label = paste0("r = ", round(cor_val, 2)),
           hjust = 0, vjust = 0, size = 5, color = "black", fontface = "bold") +
  theme_classic() +
  labs(x = "Y_1 log2FC", y = "O_1 log2FC") + 
  theme(plot.title = element_blank(),
        legend.title = element_text(size = 14, face = "bold", color = "black"),
        legend.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text = element_text(size = 12, face = "bold", color = "black"))

(p1 | p2 | p3) + plot_layout(guides = "collect") & theme(legend.position="right")
ggsave("replicate_comparison_hexbin.pdf", width = 12, height = 4)    
ggsave("replicate_comparison_hexbin.png", width = 12, height = 4, dpi = 300)                                    

# 样本相关矩阵（Correlation heatmap）
library(pheatmap)
# 计算相关矩阵
cor_mat <- cor(log2_expr, use = "pairwise.complete.obs", method = "pearson")

# 分组注释（用 group 信息）
annotation <- data.frame(
  Group = factor(c(rep("Y", 3), rep("O", 3)))
)
rownames(annotation) <- colnames(log2_expr)

# 绘制热图
pheatmap(cor_mat,
         annotation_col = annotation,
         annotation_row = annotation,
         display_numbers = TRUE,  # 显示相关系数
         number_format = "%.2f",
         main = "Sample-Sample Pearson Correlation")

library(ggplot2)
library(dplyr)
library(impute)

# PCA分析
# 1. 准备表达矩阵
# log2_expr 已经是 log2 转换并填补缺失值
expr_all <- log2_expr

# 2. 转置矩阵（行=样本，列=蛋白）
expr_all_t <- t(expr_all)

# 3. 运行 PCA
pca_all <- prcomp(expr_all_t, scale. = TRUE)

# 4. 提取 PCA 数据用于绘图
pca_all_df <- as.data.frame(pca_all$x) %>%
  mutate(Group = group)

# 5. 提取方差贡献率
percentVar_all <- round(summary(pca_all)$importance[2, 1:2] * 100, 1)

# 6. 绘制 PCA 图函数
hulls <- pca_all_df %>% 
  group_by(Group) %>%
  slice(chull(PC1, PC2)) %>%
  ungroup()
ggplot(pca_all_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  geom_point(size = 4) +
  geom_polygon(data = hulls, aes(x = PC1, y = PC2, fill = Group), alpha = 0.5, color = NA) +
  labs(
    x = paste0("PC1 (", percentVar_all[1], "%)"),
    y = paste0("PC2 (", percentVar_all[2], "%)"),
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("Y" = "blue", "O" = "red")) +
  scale_fill_manual(values = c("Y" = "blue", "O" = "red")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        panel.grid = element_blank(),                   # 去掉网格
        axis.line = element_line(color = "black"),      # 坐标轴线
        axis.ticks.length = unit(0.25, "cm"),          # 刻度线长度
        axis.ticks = element_line(color = "black"),    # 刻度线颜色
        legend.title = element_text(size = 14, face = "bold", color = "black"),
        legend.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text = element_text(size = 12, face = "bold", color = "black"))

# 7. 保存图片
ggsave("PCA.pdf", width = 5, height = 4)
ggsave("PCA.png", width = 5, height = 4, dpi = 300)

library(ggplot2)
library(ggrepel)
# 火山图（可视化差异）
# 按 log2FC 排序提取 top10
top_up <- diff_result %>% 
  filter(regulation == "Up") %>%
  arrange(desc(log2FC)) %>%
  slice_head(n = 20)
top_down <- diff_result %>%
  filter(regulation == "Down") %>%
  arrange(log2FC) %>%
  slice_head(n = 20)
top_labels <- bind_rows(top_up, top_down)

# 3. 火山图绘制 + 标注
ggplot(diff_result, aes(x = log2FC, y = -log10(adj.p.val), color = regulation)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NotSig" = "grey")) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = top_labels,
                  aes(label = Gene),
                  size = 4, box.padding = 0.5, max.overlaps = 100) +
  theme_minimal(base_size = 12) +
  labs(x = "log2 Fold Change", y = "-log10 adjusted P-value")
ggsave("deg_volcano.pdf", width = 8, height = 6)
ggsave("deg_volcano.png", width = 8, height = 6, dpi = 300)

library(pheatmap)
# 差异蛋白热图
# 设置阈值（logFC > 0.58 & adj.p.val < 0.05）
deg_genes <- diff_result %>%
  filter(abs(log2FC) > 0.58, adj.p.val < 0.05) %>%
  pull(Gene)
sig_matrix <- log2_expr[rownames(log2_expr) %in% deg_genes, ]

# 创建带有样本分组的列注释
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(sig_matrix)  

ann_colors <- list(Group = c(Y = "skyblue", O = "salmon"))

pdf("deg_pheatmap.pdf", width = 5, height = 8)
pheatmap(sig_matrix,
         scale = "row",
         name = "log2FC",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         show_colnames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize = 14,
         treeheight_row = 20,   # ← 控制行树高度，默认50
         treeheight_col = 20,   # ← 控制列树高度，默认50
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

png("deg_pheatmap.png", width = 5, height = 8, units = "in", res = 300)
pheatmap(sig_matrix,
         scale = "row",
         name = "log2FC",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         show_colnames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize = 14,
         treeheight_row = 20,   # ← 控制行树高度，默认50
         treeheight_col = 20,   # ← 控制列树高度，默认50
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()
dev.new()

# 功能富集分析
library(clusterProfiler)
library(org.Hs.eg.db) 

# 差异基因列表
# 上调和下调 DEG 筛选
up_genes <- diff_result %>%
  filter(log2FC > 0.58, adj.p.val < 0.05) %>%
  pull(Gene)
down_genes <- diff_result %>%
  filter(log2FC < -0.58, adj.p.val < 0.05) %>%
  pull(Gene)
deg_genes <- diff_result %>%
  filter(regulation != "NotSig") %>%
  pull(Gene)
# 获取背景基因（所有分析基因）
background_genes <- df$`Gene name`
# SYMBOL 转换为 ENTREZID
up_gene_df <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
down_gene_df <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
deg_gene_df <- bitr(deg_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
background_genes_df <- bitr(background_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# 富集分析
# GO
go_up <- enrichGO(gene = up_gene_df$ENTREZID, universe = background_genes_df$ENTREZID, OrgDb = org.Hs.eg.db, ont = "ALL",
                  pAdjustMethod = "BH", pvalueCutoff = 0.3, qvalueCutoff = 0.3, readable = TRUE)
go_down <- enrichGO(gene = down_gene_df$ENTREZID, universe = background_genes_df$ENTREZID, OrgDb = org.Hs.eg.db, ont = "ALL",
                    pAdjustMethod = "BH", pvalueCutoff = 0.3, qvalueCutoff = 0.3, readable = TRUE)
go_deg <- enrichGO(gene = deg_gene_df$ENTREZID, universe = background_genes_df$ENTREZID, OrgDb = org.Hs.eg.db, ont = "ALL",
                   pAdjustMethod = "BH", pvalueCutoff = 0.3, qvalueCutoff = 0.3, readable = TRUE)
write.csv(as.data.frame(go_up), "GO_result/GO_Up.csv", row.names = FALSE)
write.csv(as.data.frame(go_down), "GO_result/GO_Down.csv", row.names = FALSE)
write.csv(as.data.frame(go_deg), "GO_result/GO_Deg.csv", row.names = FALSE)

# KEGG
#kegg富集分析
kegg_up <- enrichKEGG(gene = up_gene_df$ENTREZID, universe = background_genes_df$ENTREZID, keyType = "kegg", organism = "hsa", 
                   pvalueCutoff =0.99, pAdjustMethod = "BH", qvalueCutoff =0.99)
kegg_up <- setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_down <- enrichKEGG(gene = down_gene_df$ENTREZID, universe = background_genes_df$ENTREZID, keyType = "kegg", organism = "hsa", 
                      pvalueCutoff =0.99, pAdjustMethod = "BH", qvalueCutoff =0.99)
kegg_down <- setReadable(kegg_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_deg <- enrichKEGG(gene = deg_gene_df$ENTREZID, universe = background_genes_df$ENTREZID, keyType = "kegg", organism = "hsa", 
                        pvalueCutoff =0.99, pAdjustMethod = "BH", qvalueCutoff =0.99)
kegg_deg <- setReadable(kegg_deg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write.csv(as.data.frame(kegg_up), "KEGG_result/KEGG_Up.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_down), "KEGG_result/KEGG_Down.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_deg), "KEGG_result/KEGG_Deg.csv", row.names = FALSE)

# 可视化
barplot(go_down, 
        drop = TRUE, 
        showCategory =10,
        label_format = 60,
        font.size = 12,
        split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GO_result/go_down_barplot.pdf", width =8, height = 6)
ggsave("GO_result/go_down_barplot.png", width =8, height = 6, dpi = 300)
dotplot(go_down,
        showCategory = 10, 
        label_format = 60,
        font.size = 12,
        split="ONTOLOGY",
        orderBy = "GeneRatio") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GO_result/go_down_dotplot.pdf", width =8, height = 6)
ggsave("GO_result/go_down_dotplot.png", width =8, height = 6, dpi = 300)

# 从 go_down 提取 geneID 并展开
# 转换为 data.frame
go_df <- as.data.frame(go_down)

# 筛选前 10 通路
go_data_top20 <- go_df %>%
  group_by(ONTOLOGY) %>%
  arrange(desc(Count)) %>%  
  slice_head(n = 10) %>% 
  mutate(Description = str_wrap(Description, width = 40)) %>%
  ungroup()

go_data_top20 <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "CC")) %>%  # 只保留 BP 和 CC
  group_by(ONTOLOGY) %>%
  arrange(desc(Count)) %>%  
  slice_head(n = 10) %>%
  mutate(Description = str_wrap(Description, width = 30)) %>%
  ungroup()

library(ggplot2)
library(ggforce)
library(stringr)
# 可视化---彗星图
ggplot(go_data_top20) +
  geom_link(aes(
    x = 0,
    xend = Count - 0.1,
    y = Description,
    yend = Description,
    color = -log10(p.adjust),
    alpha = after_stat(index),
    linewidth = after_stat(index)
  ),
  n = 200,
  show.legend = TRUE) + 
  guides(
    color = guide_colorbar(
      barwidth = 0.5,
      label.position = "left",
      title.position = "top",
      title.hjust = 0.5
    ),   # 保留color图例
    alpha = "none",                # 去掉alpha图例
    linewidth = "none"             # 去掉linewidth图例
  ) +
  geom_point(aes(x = Count,
                 y = Description),
             size = 6,
             shape = 21,
             fill = "white",   # 固定白色填充，空心圆效果
             color = "gray",
             stroke = 1) +
  geom_text(aes(x = Count,
                y = Description,
                label = Count),
            size = 3) +
  
  scale_color_gradient(low = "gray", high = "blue",
                       name = "-log10(adj.p)") +
  
  scale_y_discrete(position = "right") +
  xlim(0, max(go_data_top20$Count, na.rm = TRUE) + 2) +
  facet_wrap(~ONTOLOGY,
             scales="free_y",
             ncol=1,
             strip.position="left")+
  theme_bw() +
  theme(
    legend.position = "left",
    legend.justification = c(1, 0.5),
    legend.box.spacing = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = 12, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black", face = "bold"),
    axis.title.x = element_text(size = 12, color = "black", face = "bold") 
  ) +
  guides(color = guide_colorbar(barwidth = 0.6,
                                label.position = "left",
                                title.position = "top",
                                title.hjust = 0.5),  # Center the title
         fill="none")+
  xlab("Count") +
  ylab("GO Description")
ggsave("GO_result/GO_Comet.pdf", width = 8, height = 8)
ggsave("GO_result/GO_Comet.png", width = 8, height = 8, dpi = 300)

# 展开 geneID 列
go_long_top20 <- go_data_top20 %>%
  mutate(geneID = strsplit(as.character(geneID), "/")) %>%
  unnest(geneID)%>%
  left_join(diff_result, by = c("geneID" = "Gene"))  # 将 foldChange 列合并进来

# 用ggplot画heatplot
ggplot(go_long_top20, aes(x = geneID, y = Description, fill = log2FC)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("blue", "grey90"), limits = c(-2, 0), name = "log2FC",
                       oob = scales::squish) +
  facet_wrap(~ONTOLOGY,
             scales="free_y",
             ncol=1,
             strip.position="left")+
  scale_y_discrete(position = "right") +
  coord_cartesian(clip = "off", expand = FALSE) +   # 固定坐标系，不随标签扩张
  theme_minimal() +
  theme(
    plot.margin = margin(3, 20, 3, 3),  # 给右侧标签预留空间（第二个参数）
    axis.title.x = element_blank(),
    axis.text.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12, hjust = 0, color = "black", face = "bold"),  # ← 左对齐
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 12, color = "black", face = "bold", angle = 0),
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 12, color = "black", face = "bold"),
    legend.position = "right",
    legend.justification = c(1, 0.5),
    legend.box.spacing = unit(0.1, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.margin = margin(0, 0, 0, 0)
  )

ggsave("GO_result/GO_heatplot.pdf", width = 7, height = 5)
ggsave("GO_result/GO_heatplot.png", width = 7, height = 5, dpi = 300)

library(enrichplot)
#计算相似性
go_pt <- pairwise_termsim(go_down)
treeplot(go_pt, showCategory = 20,
         color = "p.adjust", #pvalue, p.adjust, qvalue
         label_format = 10, fontsize = 4,
         #添加颜色框
         hilight.params = list(hilight = T, #是否添加颜色框
                               align = "both"),# 颜色框的对齐方式’none’,’left’,’right’,’both’
         #控制各个元素之间的距离
         offset.params = list(bar_tree = rel(4.0),tiplab = rel(5.0),
                              extend = 0.4, hexpand = 0.1),
         
         #聚类相关参数，见前面的解释
         cluster.params = list(method = "ward.D", n  = 5,
                               color =  NULL, #控制颜色框的颜色，
                               label_words_n = 2,
                               label_format = 10))+
  theme(legend.position = "right",
        legend.justification = c(1, 0.5),  # Center the legend
        legend.box.spacing = unit(0.3, "cm"),  # Adjust the distance between the legend and plot
        legend.key.width = unit(0.5, "cm"),  # Adjust the width of the legend key
        legend.title = element_text(size = 12, color = "black", face = "bold"))

ggsave("GO_result/GO_treeplot.pdf", width = 8, height = 7)
ggsave("GO_result/GO_treeplot.png", width = 8, height = 7, dpi = 300)

# 计算通路间的相似性矩阵
go_sim_matrix <- pairwise_termsim(go_down)
# 可以查看通路之间的相似性矩阵
head(go_sim_matrix@termsim)
# 提取相似性矩阵
sim_matrix <- go_sim_matrix@termsim
library(pheatmap)
# 可视化相似性矩阵
pheatmap(sim_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         display_numbers = FALSE,
         color = colorRampPalette(c("white", "steelblue"))(100),
         main = "GO Term Similarity Matrix")
# 提取BP top20相似性矩阵
top_terms <- go_down@result$Description[1:20]
sim_top <- sim_matrix[top_terms, top_terms]
write.csv(as.data.frame(sim_top), file = "GO_result/go_sim_top.csv")

# 统一聚类方式（例如 Euclidean + complete linkage）
dist_top <- dist(sim_top, method = "euclidean")
hc_top <- hclust(dist_top, method = "complete")

# pheatmap
pheatmap(sim_top,
         cluster_rows = hc_top,
         cluster_cols = hc_top,
         show_colnames = FALSE,
         fontsize = 12,
         color = colorRampPalette(c("white", "steelblue"))(100),
         treeheight_row = 20,   # ← 控制行树高度，默认50
         treeheight_col = 20)   # ← 控制列树高度，默认50

# ComplexHeatmap 
library(ComplexHeatmap)
library(circlize)
pdf("GO_result/go_top_similarity.pdf", width = 9, height = 5)
Heatmap(sim_top,
        cluster_rows = as.dendrogram(hc_top),
        cluster_columns = as.dendrogram(hc_top),
        show_column_names = FALSE,
        col = colorRamp2(c(0, 1), c("white", "darkblue")),
        row_names_gp = gpar(fontsize = 12, fontface = "bold"),
        row_dend_width = unit(0.5, "cm"),
        column_dend_height = unit(0.5, "cm"),
        row_names_max_width = unit(11.2, "cm"),
        border = TRUE,
        rect_gp = gpar(col = "grey80", lwd = 0.5),
        border_gp = gpar(col = "grey80", lwd =0.5),
        heatmap_legend_param = list(
          at = seq(0, 1, by = 0.2),
          title = "Similarity",            # 或者 title = ""
          legend_direction = "vertical",
          legend_width = unit(0, "cm"),          # 色条宽度，垂直时一般不宜太宽
          legend_height = unit(5, "cm"),        # 调整色条长度
          title_gp = gpar(fontsize = 12, fontface = "bold"),
          labels_gp = gpar(fontsize = 12)
        )
)
dev.off()

png("GO_result/go_top_similarity.png", width = 9, height = 5, units = "in", res = 300)
Heatmap(sim_top,
        cluster_rows = as.dendrogram(hc_top),
        cluster_columns = as.dendrogram(hc_top),
        show_column_names = FALSE,
        col = colorRamp2(c(0, 1), c("white", "darkblue")),
        row_names_gp = gpar(fontsize = 12, fontface = "bold"),
        row_dend_width = unit(0.5, "cm"),
        column_dend_height = unit(0.5, "cm"),
        row_names_max_width = unit(11.2, "cm"),
        border = TRUE,
        rect_gp = gpar(col = "grey80", lwd = 0.5),
        border_gp = gpar(col = "grey80", lwd =0.5),
        heatmap_legend_param = list(
          at = seq(0, 1, by = 0.2),
          title = "Similarity",            # 或者 title = ""
          legend_direction = "vertical",
          legend_width = unit(0, "cm"),          # 色条宽度，垂直时一般不宜太宽
          legend_height = unit(5, "cm"),        # 调整色条长度
          title_gp = gpar(fontsize = 12, fontface = "bold"),
          labels_gp = gpar(fontsize = 12)
        )
)
dev.off()
dev.new()

# 差异蛋白热图
# 创建带有样本分组的列注释
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(sig_matrix)  

ann_colors <- list(Group = c(Y = "skyblue", O = "salmon"))

pdf("ribosomal_deg_pheatmap.pdf", width = 6, height = 10)
pheatmap(ribo_matrix,
         scale = "row",
         name = "log2FC",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize = 14,
         fontface = "bold",
         border_color = NA,
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

png("ribosomal_deg_pheatmap.png", width = 6, height = 10, units = "in", res = 300)
pheatmap(ribo_matrix,
         scale = "row",
         name = "log2FC",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize = 14,
         fontface = "bold",  # 加粗行名
         border_color = NA,
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()
dev.new()

library(clusterProfiler)
library(org.Hs.eg.db)

# mitochondrial matrix 的 GO ID
go_id <- "GO:0005759"

# 映射该GO term涉及的所有基因
mito_matrix <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = go_id,
                                     keytype = "GOALL",
                                     columns = c("SYMBOL", "ENTREZID")) %>%
  dplyr::filter(ONTOLOGYALL == "CC") %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

deg_mito_matrix <- diff_result %>% filter(Gene %in% mito_matrix$SYMBOL)
deg_mito_matrix <- deg_mito_matrix %>% filter(regulation != "NotSig")
deg_mito_matrix_down <- deg_mito_matrix%>% filter(regulation == "Down")
# 保存靶向蛋白列表
write.csv(deg_mito_matrix, "deg_mito_proteins.csv", row.names = FALSE)
write.csv(deg_mito_matrix_down, "deg_mito_down_proteins.csv", row.names = FALSE)

# mitochondrial gene expression 的 GO ID
go_id <- "GO:0140053"

# 映射该GO term涉及的所有基因
mito_gene_expr <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = go_id,
                                     keytype = "GOALL",
                                     columns = c("SYMBOL", "ENTREZID")) %>%
  dplyr::filter(ONTOLOGYALL == "BP") %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

deg_mito_gene_expr <- diff_result %>% filter(Gene %in% mito_gene_expr$SYMBOL)
deg_mito_gene_expr <- deg_mito_gene_expr %>% filter(regulation != "NotSig")

#  mitochondrial translation 的 GO ID
go_id <- "GO:0032543"

# 映射该GO term涉及的所有基因
mito_trans <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = go_id,
                                        keytype = "GOALL",
                                        columns = c("SYMBOL", "ENTREZID")) %>%
  dplyr::filter(ONTOLOGYALL == "BP") %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

deg_mito_trans <- diff_result %>% filter(Gene %in% mito_trans$SYMBOL)
deg_mito_trans <- deg_mito_trans %>% filter(regulation != "NotSig")

# translation 的 GO ID
go_id <- "GO:0006412"

# 映射该GO term涉及的所有基因
trans <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = go_id,
                                    keytype = "GOALL",
                                    columns = c("SYMBOL", "ENTREZID")) %>%
  dplyr::filter(ONTOLOGYALL == "BP") %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

deg_trans <- diff_result %>% filter(Gene %in% trans$SYMBOL)
deg_trans <- deg_trans %>% filter(regulation != "NotSig")

# peptide biosynthetic process的 GO ID
go_id <- "GO:0043043"

# 映射该GO term涉及的所有基因
pep <- AnnotationDbi::select(org.Hs.eg.db,
                               keys = go_id,
                               keytype = "GOALL",
                               columns = c("SYMBOL", "ENTREZID")) %>%
  dplyr::filter(ONTOLOGYALL == "BP") %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

deg_pep <- diff_result %>% filter(Gene %in% pep$SYMBOL)
deg_pep <- deg_pep %>% filter(regulation != "NotSig")

#  organellar ribosome的 GO ID
go_id <- "GO:0000313"

# 映射该GO term涉及的所有基因
org_ribo <- AnnotationDbi::select(org.Hs.eg.db,
                                       keys = go_id,
                                       keytype = "GOALL",
                                       columns = c("SYMBOL", "ENTREZID")) %>%
  dplyr::filter(ONTOLOGYALL == "CC") %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

deg_org_ribo <- diff_result %>% filter(Gene %in% org_ribo$SYMBOL)
deg_org_ribo <- deg_org_ribo %>% filter(regulation != "NotSig")

# mitochondrial ribosome的 GO ID
go_id <- "GO:0005761"

# 映射该GO term涉及的所有基因
mito_ribo <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = go_id,
                                  keytype = "GOALL",
                                  columns = c("SYMBOL", "ENTREZID")) %>%
  dplyr::filter(ONTOLOGYALL == "CC") %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

deg_mito_ribo <- diff_result %>% filter(Gene %in% mito_ribo$SYMBOL)
deg_mito_ribo<- deg_mito_ribo %>% filter(regulation != "NotSig")
 
# ribosomal subunit的 GO ID
go_id <- "GO:0044391"

# 映射该GO term涉及的所有基因
ribo_sub <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = go_id,
                                   keytype = "GOALL",
                                   columns = c("SYMBOL", "ENTREZID")) %>%
  dplyr::filter(ONTOLOGYALL == "CC") %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

deg_ribo_sub <- diff_result %>% filter(Gene %in% ribo_sub$SYMBOL)
deg_ribo_sub<- deg_ribo_sub %>% filter(regulation != "NotSig")

# ribosome 的 GO ID
go_id <- "GO:0005840"

# 映射该GO term涉及的所有基因
ribo <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = go_id,
                                  keytype = "GOALL",
                                  columns = c("SYMBOL", "ENTREZID")) %>%
  dplyr::filter(ONTOLOGYALL == "CC") %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

deg_ribo <- diff_result %>% filter(Gene %in% ribo$SYMBOL)
deg_ribo <- deg_ribo %>% filter(regulation != "NotSig")
 
# large ribosomal subunit 的 GO ID
go_id <- "GO:0015934"

# 映射该GO term涉及的所有基因
large_ribo <- AnnotationDbi::select(org.Hs.eg.db,
                              keys = go_id,
                              keytype = "GOALL",
                              columns = c("SYMBOL", "ENTREZID")) %>%
  dplyr::filter(ONTOLOGYALL == "CC") %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

deg_large_ribo <- diff_result %>% filter(Gene %in% large_ribo$SYMBOL)
deg_large_ribo <- deg_large_ribo %>% filter(regulation != "NotSig")

# 提取每组的蛋白 symbol 向量
large_ribo_genes <- deg_large_ribo$Gene
mito_gene <- deg_mito_gene_expr$Gene
matrix_gene <- deg_mito_matrix$Gene
mito_ribo_gene <- deg_mito_ribo$Gene
mito_trans_gene <- deg_mito_trans$Gene
org_ribo_gene <- deg_mito_ribo$Gene
pep_gene <- deg_pep$Gene
ribo_gene <- deg_ribo$Gene
ribo_sub_gene <- deg_ribo_sub$Gene
trans_gene <- deg_trans$Gene

library(ComplexUpset)
library(ggplot2)

# 构建 presence/absence 矩阵（每行一个基因，每列一个group）
all_genes <- unique(unlist(list(
  large_ribo_genes, mito_gene, matrix_gene,
  mito_ribo_gene, mito_trans_gene, org_ribo_gene,
  pep_gene, ribo_gene, ribo_sub_gene, trans_gene
)))

binary_matrix <- data.frame(
  Gene = all_genes,
  large_ribosomal_subunit = all_genes %in% large_ribo_genes,
  mitochondrial_gene_expression = all_genes %in% mito_gene,
  mitochondrial_matrix = all_genes %in% matrix_gene,
  mitochondrial_ribosome = all_genes %in% mito_ribo_gene,
  mitochondrial_translation = all_genes %in% mito_trans_gene,
  organellar_ribosome = all_genes %in% org_ribo_gene,
  peptide_biosynthetic_process = all_genes %in% pep_gene,
  ribosome = all_genes %in% ribo_gene,
  ribosomal_subunit = all_genes %in% ribo_sub_gene,
  translation = all_genes %in% trans_gene
)

# 绘图
ComplexUpset::upset(binary_matrix,
                    intersect = c("large_ribosomal_subunit", "mitochondrial_gene_expression", 
                                  "mitochondrial_matrix", "mitochondrial_ribosome",
                                  "mitochondrial_translation", "organellar_ribosome", 
                                  "peptide_biosynthetic_process", "ribosome", "ribosomal_subunit", 
                                  "translation"),
                    name = "Gene Intersections",
                    base_annotations = list(
                      'Intersection size' = intersection_size(
                        text = list(size = 5, angle = 0),
                        fill = "#1f30b4"  # 更改柱状图填充颜色
                      )
                    ),
                    set_sizes = FALSE)+
  theme(
   # axis.text.x = element_text(size = 10, face = "bold"),   # x轴字体加粗
    axis.text.y = element_text(size = 14, face = "bold", color = "black"),   # y轴字体加粗
    text = element_text(size = 14, face = "bold", color = "black"),        # 全局字体设置
    strip.text = element_text(size = 12, face = "bold", color = "black"),    # 条带标题（如图例）加粗
    plot.title = element_text(size = 14, face = "bold", color = "black")     # 图标题加粗
  ) + scale_color_manual(values = c("TRUE" = "#1f30b4", "FALSE" = "lightgrey"))

ggsave("GO_result/go_pathway_intersection_protein.pdf", width = 10, height = 7)
ggsave("GO_result/go_pathway_intersection_protein.png", width = 10, height = 7, dpi = 300)

# 去掉第一列（Gene），统计每行 TRUE 的数量
binary_matrix$Count <- rowSums(binary_matrix[ , -1])

# 提取至少在 5 个通路中出现的蛋白
overlap_genes_5plus <- binary_matrix %>%
  filter(Count >= 5) %>%
  pull(Gene)

# 查看或保存结果
print(overlap_genes_5plus)

# 从原始差异表达表中提取交集蛋白的信息
overlap_info <- diff_result %>%
  filter(Gene %in% overlap_genes_5plus)

write.csv(overlap_info, "overlap_mito_proteins.csv", row.names = FALSE)

overlap_info_all <- df_filtered %>%
  filter(`Gene name` %in% overlap_genes_5plus)
write.csv(overlap_info_all, "overlap_mito_proteins_all.csv", row.names = FALSE)

# 分类蛋白种类与调控方向
# 标注类别
library(dplyr)

# 创建一个注释表
overlap_info <- overlap_info %>%
  mutate(
    Type = case_when(
      grepl("^MRP", Gene, ignore.case = TRUE) ~ "MitoRibo",
      grepl("^RP", Gene, ignore.case = TRUE)  ~ "CytoRibo",
      TRUE ~ "Other Mitochondrial"
    )
  )

# 调控方向本身已有（Up/Down），直接组合标签
overlap_info$Category <- paste0(overlap_info$Type, "_", overlap_info$regulation)

# 按顺序设置 factor 确保热图时不会打乱顺序
overlap_info$Category <- factor(overlap_info$Category, levels = c(
  "MitoRibo_Up", "MitoRibo_Down", "CytoRibo_Up", "CytoRibo_Down", "Other Mitochondrial_Up", "Other Mitochondrial_Down"))

# 提取并排序表达矩阵
# 确保 gene symbol 是行名
overlap_genes <- log2_expr[rownames(log2_expr) %in% overlap_genes_5plus, ]
# 按你指定的类别顺序排序
overlap_genes_sorted <- overlap_info %>% arrange(Category)
overlap_matrix <- overlap_genes[overlap_genes_sorted$Gene, ]

# 差异蛋白热图
# 创建带有样本分组的列注释
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(overlap_matrix)  

ann_colors <- list(Group = c(Y = "steelblue", O ="firebrick"))
# 定义颜色
col_fun <- colorRamp2(c(-2, 0, 2),
                      c("blue", "white", "red"))
library(ComplexHeatmap)
library(circlize)  # 用于 colorRamp2
Heatmap(overlap_matrix,
        name = "log2FC",
        col = col_fun,
        cluster_rows = FALSE,          # 不聚类行
        cluster_columns = FALSE,       # 不聚类列
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 12, fontface = "bold"),
        column_names_gp = gpar(fontsize = 12, fontface = "bold"),
        border = NA,
        top_annotation = if(!is.null(annotation_col)) HeatmapAnnotation(df = annotation_col, 
                                                                        col = ann_colors) else NULL,
        heatmap_legend_param = list(
          title = "log2FC",
          title_gp = gpar(fontsize = 12, fontface = "bold"),
          labels_gp = gpar(fontsize = 12)
        )
)

pdf("overlap_gene_pheatmap.pdf", width = 4.5, height = 8)
pheatmap(overlap_matrix,
         scale = "row",
         name = "log2FC",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize = 14,
         fontface = "bold",
         border_color = NA,
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

png("overlap_gene_pheatmap.png", width = 4.5, height = 8, units = "in", res = 300)
pheatmap(overlap_matrix,
         scale = "row",
         name = "log2FC",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize = 14,
         fontface = "bold",  # 加粗行名
         border_color = NA,
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()
dev.new()

