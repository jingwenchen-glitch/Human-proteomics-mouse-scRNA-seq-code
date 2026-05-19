# 清除系统环境变量，加载R包：
rm(list=ls())
options(stringsAsFactors = F) 
setwd("/Users/jingwenchen/Desktop/Ph.D/Jin Lab/scRNA-seq/Results/")

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)

# 读取文件
tcell_clean <- readRDS("tcell_clean_celltype.rds")

# 挖掘子群的DEG基因
Idents(tcell_clean) <- "orig.ident"
tcell_DEG <- FindMarkers(tcell_clean, ident.1 = "OE",ident.2 = "NC", logfc.threshold = 0.25,  min_pct = 0.1)
write.csv(tcell_DEG, "tcell_DEG/OE_vs_NC_DEG.csv", row.names = TRUE) 

saveRDS(tcell_DEG, file = "tcell_DEG.rds")
tcell_DEG <- readRDS(file = "tcell_DEG.rds")

# 加载必要的包
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(Seurat)

# 创建排序基因列表（按logFC降序）
gene_rank <- tcell_DEG %>%
  arrange(desc(avg_log2FC)) %>%
  mutate(gene = rownames(.)) %>%
  dplyr::select(gene, avg_log2FC) %>%
  {setNames(.$avg_log2FC, .$gene)}

# 获取MSigDB基因集（小鼠）
msigdb_sets <- list(
  Hallmark = msigdbr(species = "Mus musculus", category = "H"),
  KEGG = msigdbr(species = "Mus musculus", category = "C2") %>% 
    filter(grepl("KEGG", gs_subcat)),
  GO_BP = msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP"),
  GO_MF = msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:MF"),
  GO_CC = msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:CC"),
  Immunologic = msigdbr(species = "Mus musculus", category = "C7"),
  Oncogenic = msigdbr(species = "Mus musculus", category = "C6")
)

# 转换为fgsea需要的格式
pathways <- lapply(msigdb_sets, function(x) split(x$gene_symbol, x$gs_name))

# 执行GSEA分析
gsea_results <- list()
output_dir <- "tcell_DEG"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (set_name in names(pathways)) {
  cat("Running GSEA for:", set_name, "\n")
  
  # 运行fgsea
  fgsea_res <- fgseaMultilevel(
    pathways = pathways[[set_name]],
    stats = gene_rank,
    minSize = 10,
    maxSize = 1000,
    eps = 0,
    nproc = 8
  )
  
  # 添加富集方向列
  fgsea_res$Enrichment <- ifelse(fgsea_res$NES > 0, "Up", "Down")
  
  # 创建可保存的副本（转换所有列表列为字符串）
  fgsea_res_save <- fgsea_res
  if ("leadingEdge" %in% colnames(fgsea_res_save)) {
    fgsea_res_save$leadingEdge <- sapply(
      fgsea_res_save$leadingEdge, 
      function(x) paste(x, collapse = ";")
    )
  }
  
  # 筛选显著结果 (FDR < 0.25)
  sig_res <- fgsea_res_save[fgsea_res_save$padj < 0.25, ]
  
  if (nrow(sig_res) > 0) {
    # 保存完整结果
    write.csv(
      fgsea_res_save, 
      file.path(output_dir, paste0("Full_", set_name, "_GSEA.csv")), 
      row.names = FALSE
    )
    
    # 保存显著结果
    write.csv(
      sig_res, 
      file.path(output_dir, paste0("Significant_", set_name, "_GSEA.csv")), 
      row.names = FALSE
    )
  }
  
  gsea_results[[set_name]] <- fgsea_res
}

# 保存所有结果
saveRDS(gsea_results, file.path(output_dir, "All_GSEA_Results.rds"))

# NC vs OE 富集分析
library(clusterProfiler)
library(org.Mm.eg.db)  
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# 上调和下调 DEG 筛选
up_genes <- tcell_DEG %>%
  filter(avg_log2FC > 0.25, p_val <= 0.01) %>%
  rownames() %>% unique()
down_genes <- tcell_DEG %>%
  filter(avg_log2FC < -0.25, p_val <= 0.01) %>%
  rownames() %>% unique()
deg_genes <- tcell_DEG %>%
  filter(abs(avg_log2FC) > 0.25, p_val <= 0.01) %>%
  rownames() %>% unique()
# SYMBOL 转 ENTREZID
up_gene_df <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
down_gene_df <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
deg_gene_df <- bitr(deg_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
# 富集分析
# GO
go_up <- enrichGO(gene = up_gene_df$ENTREZID, OrgDb = org.Mm.eg.db, ont = "ALL",
                  pAdjustMethod = "BH", pvalueCutoff = 0.3, readable = TRUE)
go_down <- enrichGO(gene = down_gene_df$ENTREZID, OrgDb = org.Mm.eg.db, ont = "ALL",
                    pAdjustMethod = "BH", pvalueCutoff = 0.3, readable = TRUE)
go_deg <- enrichGO(gene = deg_gene_df$ENTREZID, OrgDb = org.Mm.eg.db, ont = "ALL",
                   pAdjustMethod = "BH", pvalueCutoff = 0.3, readable = TRUE)
write.csv(as.data.frame(go_up), "tcell_DEG/GO_Up.csv", row.names = FALSE)
write.csv(as.data.frame(go_down), "tcell_DEG/GO_Down.csv", row.names = FALSE)
write.csv(as.data.frame(go_deg), "tcell_DEG/GO_Deg.csv", row.names = FALSE)

# KEGG
kegg_up <- enrichKEGG(gene = up_gene_df$ENTREZID, organism = "mmu", pvalueCutoff = 0.3)
kegg_up <- setReadable(kegg_up, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
kegg_down <- enrichKEGG(gene = down_gene_df$ENTREZID, organism = "mmu", pvalueCutoff = 0.3)
kegg_down <- setReadable(kegg_down, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
kegg_deg <- enrichKEGG(gene = deg_gene_df$ENTREZID, organism = "mmu", pvalueCutoff = 0.3)
kegg_deg <- setReadable(kegg_deg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
write.csv(as.data.frame(kegg_up), "tcell_DEG/KEGG_Up.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_down), "tcell_DEG/KEGG_Down.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_deg), "tcell_DEG/KEGG_Deg.csv", row.names = FALSE)

# GSEA
library(clusterProfiler)
library(org.Mm.eg.db)  # 小鼠注释数据库
library(dplyr)
library(enrichplot)    # 可视化
library(ggplot2)
# 构建有序基因向量
gene_list <- tcell_DEG %>%
  filter(!is.na(avg_log2FC)) %>%
  arrange(desc(avg_log2FC)) %>%  # 可改为升序/降序看具体需求
  pull(avg_log2FC)

names(gene_list) <- rownames(tcell_DEG)[match(gene_list, tcell_DEG$avg_log2FC)]
head(gene_list)
# SYMBOL 转 ENTREZID
gene_df <- bitr(names(gene_list),
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)

# 重新映射 gene_list 到 ENTREZID
gene_list_named <- gene_list[gene_df$SYMBOL]
names(gene_list_named) <- gene_df$ENTREZID

# GSEA GO 分析（所有 ontology）
gsea_go <- gseGO(
  geneList = gene_list_named,
  OrgDb = org.Mm.eg.db,
  ont = "ALL",         # "BP", "MF", "CC", 或 "ALL"
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.3,
  verbose = FALSE
)
gsea_go <- setReadable(gsea_go, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
write.csv(as.data.frame(gsea_go), file = "tcell_DEG/GSEA_GO.csv", row.names = FALSE)

# GSEA KEGG 分析
gsea_kegg <- gseKEGG(
  geneList = gene_list_named,
  organism = "mmu",     # 人类为 "hsa"
  keyType = "ncbi-geneid",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.3,
  verbose = FALSE
)
gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
write.csv(as.data.frame(gsea_kegg), file = "tcell_DEG/GSEA_KEGG.csv", row.names = FALSE)

# 可视化 GSEA 结果
library(enrichplot)
library(clusterProfiler)

# GSEA条形图（可选 topN）
dotplot(gsea_go, showCategory = 5, split = "ONTOLOGY") +
  ggplot2::facet_grid(~ONTOLOGY)

barplot(gsea_go, showCategory = 15, split = "ONTOLOGY") + 
  ggplot2::facet_grid(~ONTOLOGY)
barplot(gsea_kegg, showCategory = 15)
# ridge plot
ridgeplot(gsea_go) + labs(title = "GO GSEA Ridgeplot")
# enrichment plot for a specific pathway
gseaplot2(gsea_go, geneSetID = "GO:0031577", title = gsea_go$Description[1])
# 取某个 cluster 画 GSEA 富集图
term_id <- "GO:0031577"
desc <- gsea_go@result[term_id, "Description"]

gseaplot2(gsea_go,
          geneSetID = term_id,
          title = desc)


library(Seurat)
library(msigdbr)
library(dplyr)
# 确认 cluster 列
Idents(tcell_clean) <- "celltype"

# 只保留 cluster 0
sub <- subset(tcell_clean, idents = c("TCM"))

# 获取 KEGG 通路基因集
msig_kegg <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "KEGG_LEGACY")

# 取相关通路
gly_kegg  <- msig_kegg %>% filter(gs_name == "KEGG_GLYCOLYSIS_GLUCONEOGENESIS") %>% pull(gene_symbol)
oxphos_kegg <- msig_kegg %>% filter(gs_name == "KEGG_OXIDATIVE_PHOSPHORYLATION") %>% pull(gene_symbol)
ribo_kegg <- msig_kegg %>% filter(gs_name == "KEGG_RIBOSOME") %>% pull(gene_symbol)


# 计算单细胞模块分数
sub <- AddModuleScore(sub, features = list(gly_kegg), name = "KEGG_GLY_")
sub <- AddModuleScore(sub, features = list(oxphos_kegg), name = "KEGG_OXPHOS_")
sub <- AddModuleScore(sub, features = list(ribo_kegg), name = "KEGG_RIBOSOME_")


# 可视化（按 cluster 分组、按 condition 分面）
VlnPlot(sub, features = c("KEGG_GLY_1"), group.by = "celltype", split.by = "orig.ident", pt.size = 0.1) +
  labs(title="KEGG Glycolysis Activity")
VlnPlot(sub, features = c("KEGG_OXPHOS_1"), group.by = "celltype", split.by = "orig.ident", pt.size = 0.1) +
  labs(title="KEGG OxPhos Activity")
VlnPlot(sub, features = c("KEGG_RIBOSOME_1"), group.by = "celltype", split.by = "orig.ident", pt.size = 0.1) +
  labs(title="KEGG Ribosome Activity")


# 获取 GOCC 通路基因集
msig_gocc <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "CC")

# 取相关通路
ribo_gocc  <- msig_gocc %>% filter(gs_name == "GOCC_RIBOSOME") %>% pull(gene_symbol)
res_gocc <- msig_gocc %>% filter(gs_name == "GOCC_RESPIRATORY_CHAIN_COMPLEX") %>% pull(gene_symbol)
org_ribo_gocc <- msig_gocc %>% filter(gs_name == "GOCC_ORGANELLAR_RIBOSOME") %>% pull(gene_symbol)

# 计算单细胞模块分数
sub <- AddModuleScore(sub, features = list(ribo_gocc), name = "GOCC_RIBO_")
sub <- AddModuleScore(sub, features = list(res_gocc), name = "GOCC_RESPI_")
sub <- AddModuleScore(sub, features = list(org_ribo_gocc), name = "GOCC_ORGAN_RIBO_")

# 可视化（按 cluster 分组、按 condition 分面）
VlnPlot(sub, features = c("GOCC_RIBO_1"), group.by = "celltype", split.by = "orig.ident", pt.size = 0.1) +
  labs(title="GOCC Ribosome Activity")
VlnPlot(sub, features = c("GOCC_RESPI_1"), group.by = "celltype", split.by = "orig.ident", pt.size = 0.1) +
  labs(title="GOCC Respiratory chain Activity")
VlnPlot(sub, features = c("GOCC_ORGAN_RIBO_1"), group.by = "celltype", split.by = "orig.ident", pt.size = 0.1) +
  labs(title="GOCC Organellar Ribosome Activity") 


# 获取 GOBP 通路基因集
msig_gobp <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")

# 取相关通路
mitogene_gobp  <- msig_gobp %>% filter(gs_name == "GOBP_MITOCHONDRIAL_GENE_EXPRESSION") %>% pull(gene_symbol)
mitotrans_gobp <- msig_gobp %>% filter(gs_name == "GOBP_MITOCHONDRIAL_TRANSLATION") %>% pull(gene_symbol)
oxphos_gobp <- msig_gobp %>% filter(gs_name == "GOBP_OXIDATIVE_PHOSPHORYLATION") %>% pull(gene_symbol)
fao_gobp <- msig_gobp %>% filter(gs_name == "GOBP_FATTY_ACID_BETA_OXIDATION") %>% pull(gene_symbol)

# 计算单细胞模块分数
sub <- AddModuleScore(sub, features = list(ribo_gocc), name = "GOBP_MITOGENE_")
sub <- AddModuleScore(sub, features = list(res_gocc), name = "GOBP_MITOTRANS_")
sub <- AddModuleScore(sub, features = list(org_ribo_gocc), name = "GOBP_OXPHOS_")
sub <- AddModuleScore(sub, features = list(fao_gobp), name = "GOBP_FAO_")

# 可视化（按 cluster 分组、按 condition 分面）
VlnPlot(sub, features = c("GOBP_MITOGENE_1"), group.by = "celltype", split.by = "orig.ident", pt.size = 0.1) +
  labs(title="GOBP Mitochondrial gene expression Activity")
VlnPlot(sub, features = c("GOBP_MITOTRANS_1"), group.by = "celltype", split.by = "orig.ident", pt.size = 0.1) +
  labs(title="GOBP Mitochondrial translation Activity")
VlnPlot(sub, features = c("GOBP_OXPHOS_1"), group.by = "celltype", split.by = "orig.ident", pt.size = 0.1) +
  labs(title="GOBP OXPHOS Activity") 
VlnPlot(sub, features = c("GOBP_FAO_1"), group.by = "celltype", split.by = "orig.ident", pt.size = 0.1) +
  labs(title="GOBP Fatty Acid Oxidation Activity")

library(ggpubr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(scales)
library(tidyr)

# 多通路
score_cols <- c( "GOBP_MITOTRANS_1", "GOCC_RIBO_1", "GOCC_ORGAN_RIBO_1", "GOBP_MITOGENE_1", "GOCC_RESPI_1", "KEGG_OXPHOS_1", "KEGG_GLY_1")
meta_multi <- sub@meta.data[, c("celltype", "orig.ident", score_cols)]
meta_multi$celltype <- as.factor(meta_multi$celltype)
meta_multi$orig.ident <- as.factor(meta_multi$orig.ident)
meta_melt <- meta_multi %>%
  pivot_longer(
    cols = score_cols,         # 要展开的列
    names_to = "Pathway",      # 新列名（通路名）
    values_to = "Score"        # 数值列名（打分）
  )

# 生成比较组
my_comparisons <- list(c("NC", "OE"))

# 每个分面的 y 位置
pos_df <- meta_melt %>%
  group_by(Pathway, celltype) %>%
  summarise(y_pos = max(Score, na.rm = TRUE) * 0.7, .groups = "drop")

# 定义映射表
pathway_names <- c(
  "GOBP_MITOTRANS_1"  = "Mitochondrial\ntranslation",
  "GOCC_RIBO_1"       = "Ribosome",
  "GOCC_ORGAN_RIBO_1" = "Organellar\nribosome",
  "GOBP_MITOGENE_1"   = "Mitochondrial\ngene expression",
  "GOCC_RESPI_1"      = "Respiratory\nchain complex",
  "KEGG_OXPHOS_1"     = "Oxidative\nphosphorylation",
  "KEGG_GLY_1"        = "Glycolysis"
)
# 按给定顺序设置 factor
meta_melt$Pathway <- factor(
  recode(meta_melt$Pathway, !!!pathway_names),
  levels = pathway_names  # 按原始定义的顺序排列
)

ggplot(meta_melt, aes(x = orig.ident, y = Score, fill = orig.ident)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  facet_wrap(~ Pathway, scales = "fixed", nrow = 2) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.35,
                     label.y.npc = "top",
                     vjust = 0.95,
                     size = 7,
                     fontface = "bold") +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  scale_fill_manual(
    values = c("NC" = "#EA3323", "OE" = "#0000F5")) + # 红色
  theme_bw(base_size = 12) +
  labs(y = "Module Score", x = "") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 0.8),
    strip.text.x = element_text(size = 14, face = "bold", color = "black", angle = 0, hjust = 0.5),
    strip.text.y = element_text(size = 14, face = "bold", color = "black"),
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none",
    panel.spacing = unit(0.2, "lines")
  )

ggsave("enrichment_results/interested_pathway_violin_plot_cluster0.pdf", width = 7.5, height = 6)
ggsave("enrichment_results/interested_pathway_violin_plot_cluster0.png", width = 7.5, height = 6, dpi = 300)