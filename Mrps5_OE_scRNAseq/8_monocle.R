setwd("/Users/jingwenchen/Desktop/Ph.D/Jin Lab/scRNA-seq/Results/")
rm(list= ls())

library(Seurat)
library(monocle3)
library(ggplot2)
library(patchwork)

sce <- readRDS("tcell_clean_celltype.rds")
p1 <- DimPlot(sce,pt.size = 0.8,group.by = "celltype",label=T)

Idents(sce) <- sce$celltype
levels(Idents(sce))

expression_matrix<- GetAssayData(sce,assay='RNA',layer='counts')
cell_metadata <-  sce@meta.data
gene_annotation <- data.frame(gene_short_name=rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata=cell_metadata,
                         gene_metadata=gene_annotation)

cds <- preprocess_cds(cds, num_dim = 30)
plot_pc_variance_explained(cds)

reducedDim(cds,"UMAP") <- Embeddings(sce,"umap")
cds <- cluster_cells(cds)

cds<- learn_graph(cds,verbose=T,
                   use_partition=T)

plot_cells(cds,color_cells_by = "celltype",label_groups_by_cluster=FALSE,label_leaves=FALSE,label_branch_points=FALSE)



start_celltypes <- c("TCM", "TEM") # 替换为你的两个起点
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds), ]) # 细胞→主图顶点映射
root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name # 所有主图顶点名称

root_vertices <- c()
for (start_type in start_celltypes) {
  start_cells <- colnames(cds)[colData(cds)$celltype == start_type]
  vertex_freq <- table(closest_vertex[start_cells, ])
  top_vertex <- as.numeric(names(which.max(vertex_freq)))
  root_vertices <- c(root_vertices, root_pr_nodes[top_vertex])
}
cds = order_cells(cds, root_pr_nodes=root_vertices)

#start = c("TEM")
#closest_vertex = cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
#closest_vertex = as.matrix(closest_vertex[colnames(cds), ])
#root_pr_nodes = igraph::V(principal_graph(cds)[["UMAP"]])$name
#flag = closest_vertex[as.character(colData(cds)$celltype) %in% start,]
#flag = as.numeric(names(which.max(table(flag))))
#root_pr_nodes = root_pr_nodes[flag]
#cds = order_cells(cds, root_pr_nodes=root_pr_nodes)

# 可视化
p2 <- plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=F,
           label_groups_by_cluster=F,
           label_roots=F,
           label_leaves=F,
           label_branch_points=F,
           cell_size=0.8,
           group_label_size=4,
           trajectory_graph_color='black',
           trajectory_graph_segment_size=0.5
           )

p1 + p2
ggsave("./monocle3_twotree.pdf",width = 14, height = 7)

