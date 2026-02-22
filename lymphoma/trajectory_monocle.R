### monocle 2 analysis ###

library(SeuratData)
library(SeuratDisk)
library(Seurat)
library(ggplot2)
library(dplyr)
Convert("C5C1_subset.h5ad", dest = "h5seurat", overwrite = TRUE)
novaseq <- LoadH5Seurat("C5C1_subset.h5seurat", meta.data = FALSE, misc = FALSE)

#now start analysis
library(Seurat)
library(ggplot2)
library(dplyr)
library(monocle)
metadata <- read.csv("signac_output/EpiTrace_SignacPeak_Convergence_out.csv")

filtered_metadata <- metadata[metadata$cell %in% Cells(novaseq), ]
ordered_metadata <- filtered_metadata[match(Cells(novaseq), filtered_metadata$cell), ]

novaseq$celltype <- ordered_metadata$celltype
novaseq$EpiTraceAge_iterative <- ordered_metadata$EpiTraceAge_iterative

metadata2 <- read.csv("snapatac_out/MALT_meta.csv")
filtered_metadata2 <- metadata2[metadata2$X %in% Cells(novaseq), ]
ordered_metadata2 <- filtered_metadata2[match(Cells(novaseq), filtered_metadata2$X), ]
novaseq$n_fragment <- ordered_metadata2$n_fragment
novaseq$frac_mito <- ordered_metadata2$frac_mito

novaseq<-subset(novaseq, subset = n_fragment > 3000)

data <- as(as.matrix(novaseq@assays$RNA@data), 'sparseMatrix')
pData <- data.frame(cell_name = colnames(data), row.names = colnames(data),
                    cluster_use=novaseq@meta.data$celltype, EpiTraceAge_iterative=novaseq@meta.data$EpiTraceAge_iterative,
                    n_fragment=novaseq@meta.data$n_fragment, frac_mito=novaseq@meta.data$frac_mito)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
pd <- new("AnnotatedDataFrame", data = pData)
fd <- new("AnnotatedDataFrame", data = fData)
cds <- newCellDataSet(data,phenoData = pd,featureData = fd)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))
length(expressed_genes)

fData(cds)$use_for_ordering <-fData(cds)$num_cells_expressed > 0.1 * ncol(cds)
#table(fData(cds)$use_for_ordering)

clustering_DEG_genes <-
  differentialGeneTest(cds[expressed_genes,],
                       fullModelFormulaStr = '~cluster_use',
                       cores = 5)

HSMM_ordering_genes <-
  row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
HSMM_myo <-
  setOrderingFilter(cds,
                    ordering_genes = HSMM_ordering_genes)
HSMM_myo <-
  reduceDimension(HSMM_myo, max_components = 3, method = 'DDRTree')

plot_cell_trajectory(HSMM_myo, color_by = "cluster_use",
                     cell_size = 1, cell_link_size = 1, show_branch_points=FALSE)+
labs(color = "Cluster")+
theme(legend.position = "right",legend.text=element_text(size=10),
        legend.key = element_rect(fill = "white"))+
guides(color = guide_legend(override.aes = list(size = 4)))+
  scale_color_manual(breaks = c('C1','C5'), values=c('darkorchid','firebrick'))

plot_cell_trajectory(HSMM_myo, color_by = "EpiTraceAge_iterative",
                     cell_size = 1, cell_link_size = 1, show_branch_points=FALSE)+
labs(color = "EpiTrace Age")+
theme(legend.position = "right",legend.text=element_text(size=10),
        legend.key = element_rect(fill = "white"))+
scale_color_gradientn(colors = c("red", "lightgrey", "green"))


HSMM_myo <- clusterCells(HSMM_myo, num_clusters = 8, verbose = F)
plot_cell_trajectory(HSMM_myo, color_by = "Cluster",
                     cell_size = 1, cell_link_size = 1, show_branch_points=FALSE)+
scale_color_brewer(palette = "Set1")+labs(color = "Monocle Cluster")+
theme(legend.position = "right",legend.text=element_text(size=11),
        legend.key = element_rect(fill = "white"))+guides(color = guide_legend(override.aes = list(size = 4)))


df <- data.frame(classes = HSMM_myo$Cluster, values = HSMM_myo$EpiTraceAge_iterative)
summary_df <- df %>%
  group_by(classes) %>%
  summarise(mean = mean(values),sem = sd(values) / sqrt(n()))
# Plot
ggplot(summary_df, aes(x = classes, y = mean, fill = classes)) +
  geom_bar(stat = "identity") +
  #scale_fill_manual(values = c("1" = "blue", "2" = "green", "3" = "red", "4" = "red", "5" = "red", "6" = "red", "7" = "red")) +
scale_fill_brewer(palette = "Set1") +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
  labs(x = "Monocle Cluster", y = "EpiTrace Age\n") +
  coord_cartesian(ylim = c(0.5, 0.82)) +  # Set y-axis limits without removing data
  theme_minimal() + 
theme(legend.position = "none",
      axis.title.x = element_text(size = 12.5), axis.text = element_text(size = 12),
      axis.title.y = element_text(size = 14))



p=plot_cell_trajectory(HSMM_myo, color_by = "EpiTraceAge_iterative",
                     cell_size = 1, cell_link_size = 1, show_branch_points=FALSE)+
labs(color = "EpiTrace Age")+
theme(legend.position = "right",legend.text=element_text(size=10),
        legend.key = element_rect(fill = "white"))+
scale_color_gradientn(colors = c("#FF0000", "#D3D3D3", "#008000"))
ggsave(filename = 'signac_output/Monocle2/C5C1_Monocle2_age.png', plot = p, width = 7, height = 4.5, units = "in", dpi = 300, bg = "white")


p=plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)')+
labs(color = "Monocle Cluster")+theme(legend.position = "right",legend.text=element_text(size=12),
        legend.key = element_rect(fill = "white"))+guides(color = guide_legend(override.aes = list(size = 3)))+
scale_color_brewer(palette = "Set1")
ggsave(filename = 'signac_output/Monocle2/C5C1_Monocle2_tsne.png', plot = p, width = 8, height = 6, dpi = 300, bg = "transparent")




### monocle 3 analysis ###

library(Seurat)
library(monocle3)
library(ggplot2)
library(dplyr)

load('C5C1.Rdata')

metadata2 <- read.csv("MALT_snapatac_meta.csv")
filtered_metadata2 <- metadata2[metadata2$X %in% Cells(novaseq), ]
ordered_metadata2 <- filtered_metadata2[match(Cells(novaseq), filtered_metadata2$X), ]
novaseq$n_fragment <- ordered_metadata2$n_fragment
novaseq$frac_mito <- ordered_metadata2$frac_mito
novaseq<-subset(novaseq, subset = n_fragment > 3000)

data <- as(as.matrix(novaseq@assays$RNA@data), 'sparseMatrix')
pData <- data.frame(cell_name = colnames(data), row.names = colnames(data),
                    cluster_use=novaseq@meta.data$celltype, EpiTraceAge_iterative=novaseq@meta.data$EpiTraceAge_iterative,
                    n_fragment=novaseq@meta.data$n_fragment, frac_mito=novaseq@meta.data$frac_mito)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

cds <- new_cell_data_set(data,
                         cell_metadata = pData,
                         gene_metadata = fData)

cds <- preprocess_cds(cds, num_dim = 6)
plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds, reduction_method='UMAP', preprocess_method='PCA', umap.min_dist = 0.4)

plot_cells(cds, cell_size = 2, color_cells_by = "cluster_use", label_cell_groups=FALSE)+
theme(legend.position = "right",legend.text=element_text(size=10),
        legend.key = element_rect(fill = "white"))+
guides(color = guide_legend(override.aes = list(size = 4)))+
  scale_color_manual(breaks = c('C1','C5'),
                     values=c('darkorchid','firebrick'))

plot_cells(cds, cell_size = 2, color_cells_by = "EpiTraceAge_iterative", label_cell_groups=FALSE)+
theme(legend.position = "right",legend.text=element_text(size=10),
        legend.key = element_rect(fill = "white"))+
scale_color_gradientn(colors = c("#FF0000", "#D3D3D3", "#008000"))

colData(cds)$clusters <- cds@clusters@listData$UMAP$clusters
cds <- learn_graph(cds, close_loop=FALSE)

plot_cells(cds,
           color_cells_by = "cluster_use",cell_size = 2, alpha=1,
           trajectory_graph_segment_size=2,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, label_roots = FALSE)+
  theme(legend.position = "right",legend.text=element_text(size=15),legend.title=element_blank(),
        legend.key = element_rect(fill = "white"))+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  theme(axis.title = element_text(size = 16))+
  scale_color_manual(breaks = c('C1','C5'),
                     values=c('darkorchid','firebrick'))

plot_cells(cds,
           color_cells_by = "EpiTraceAge_iterative",cell_size = 2, alpha=1,
           trajectory_graph_segment_size=3,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, label_roots = FALSE)+
  theme(legend.position = "right",legend.text=element_text(size=15),legend.title=element_blank(),
        legend.key = element_rect(fill = "white"))+
  theme(axis.title = element_text(size = 16))+
  scale_color_gradientn(colors = c("#FF0000", "#D3D3D3", "#008000"))