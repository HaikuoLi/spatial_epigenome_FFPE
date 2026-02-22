options(warn=-1)

library(SeuratData)
library(SeuratDisk)
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)
library(SeuratWrappers)
library(NICHES)
library(viridis)


novaseq <- LoadH5Seurat("H3K27me3_data.h5seurat", meta.data = FALSE, misc = FALSE)
metadata <- read.csv("./snapatac_out/H3K27me3_snapatac_meta.csv")
filtered_metadata <- metadata[metadata$X %in% Cells(novaseq), ]
ordered_metadata <- filtered_metadata[match(Cells(novaseq), filtered_metadata$X), ]
novaseq$n_fragment <- ordered_metadata$n_fragment
novaseq$frac_dup <- ordered_metadata$frac_dup
novaseq$frac_mito <- ordered_metadata$frac_mito
novaseq$coord <- ordered_metadata$coord
novaseq$x_coord <- ordered_metadata$x_coord
novaseq$y_coord <- ordered_metadata$y_coord
novaseq$tsse <- ordered_metadata$tsse
novaseq$cluster_use <- ordered_metadata$cluster_use

pos_sample_df_cor <- data.frame(
  y_coord = novaseq@meta.data$y_coord,
  x_coord = novaseq@meta.data$x_coord
)
rownames(pos_sample_df_cor)<-rownames(novaseq@meta.data)
novaseq[['whole']] <- new(Class = 'SlideSeq', assay = "RNA", coordinates = pos_sample_df_cor)

SpatialDimPlot(novaseq, group.by = "cluster_use", pt.size.factor = 5) +
  scale_fill_manual(values = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b'))

Idents(novaseq)<-novaseq$cluster_use
novaseq@meta.data$x <- novaseq@meta.data$x_coord
novaseq@meta.data$y <- novaseq@meta.data$y_coord

novaseq <- NormalizeData(novaseq)

novaseq <- SeuratWrappers::RunALRA(novaseq)

NICHES_output <- RunNICHES(object = novaseq,
                           LR.database = "fantom5",
                           species = "human",
                           assay = "alra",
                           position.x = 'x',
                           position.y = 'y',
                           k = 4, 
                           cell_types = "cluster_use",
                           min.cells.per.ident = 0,
                           min.cells.per.gene = NULL,
                           meta.data.to.map = c('cluster_use'),
                           CellToCell = F,CellToSystem = F,SystemToCell = F,
                           CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = T)


niche <- NICHES_output[['NeighborhoodToCell']]

# Scale and visualize
niche <- ScaleData(niche)
niche <- FindVariableFeatures(niche,selection.method = "disp")
niche <- RunPCA(niche)
ElbowPlot(niche,ndims = 50)

niche <- RunUMAP(niche,dims = 1:10)

levels(niche)<-c('C0','C1','C2','C3','C4','C5')
DimPlot(niche,reduction = 'umap',pt.size = 2) +ggtitle('Cellular Microenvironment') +
  scale_color_manual(values = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b'))

# Find markers
mark <- FindAllMarkers(niche,min.pct = 0.25,only.pos = T,test.use = "roc")
GOI_niche <- mark %>% group_by(cluster) %>% top_n(5,myAUC)
DoHeatmap(niche,features = unique(GOI_niche$gene))+ 
  scale_fill_gradientn(colors = c("grey","white", "blue"))

# Add Niches output as an assay
niches.data <- GetAssayData(object =  niche[['NeighborhoodToCell']], slot = 'data')
colnames(niches.data) <- niche[['ReceivingCell']]$ReceivingCell

novaseq[["NeighborhoodToCell"]] <- CreateAssayObject(data = niches.data )
DefaultAssay(novaseq) <- "NeighborhoodToCell"
novaseq <- ScaleData(novaseq)

p <- SpatialFeaturePlot(
  novaseq,
  features = "IL23A—IL12RB1", pt.size.factor = 5,
  slot = "scale.data",
  min.cutoff = "q1",
  max.cutoff = "q95"
) + theme(legend.position = "right") + labs(fill = "IL23A-IL12RB1")

# Save the plot as a PDF with 300 DPI
ggsave("niches/20250702_CT1_0605_niches_IL23A—IL12RB1.pdf", plot = p, dpi = 300, width = 8, height = 6)

write.csv(mark,"niches/H3K27me3_mark.csv")