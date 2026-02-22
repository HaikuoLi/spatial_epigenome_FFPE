.libPaths(c("/gpfs/gibbs/project/fan/hl737/R/4.3",'/vast/palmer/home.mccleary/hl737/.conda/envs/ds_notebook/lib/R/library',.libPaths()))
load('signac_output/mouse_brain_atac_20240828.Rdata') #you can generate this seurat object with the previous processing scripts (we named it as 'novaseq')

library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

#https://github.com/stuart-lab/signac/issues/486
main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
keep.peaks <- as.logical(seqnames(granges(novaseq)) %in% main.chroms)
novaseq

novaseq <- novaseq[keep.peaks, ]
novaseq


# add motif information
novaseq <- AddMotifs(
  object = novaseq,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

novaseq <- RunChromVAR(
  object = novaseq,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

save(novaseq,keep.peaks,main.chroms,pfm,file='signac_output/mouse_brain_atac_20240905_chromvar.Rdata')

DefaultAssay(novaseq) <- 'chromvar'

da_motifs <- FindAllMarkers(
  object = novaseq,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.1,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  logfc.threshold = 0 #find all cluster-specific motifs
)
Sys.time()

# Add motif name
DefaultAssay(novaseq) <- "peaks"
motif.names <- da_motifs$gene
da_motifs$gene <- ConvertMotifID(novaseq, id = motif.names)
da_motifs$motif <- motif.names

colnames(da_motifs) <- c("motif.p_val","motif.avg_diff.","motif.pct.1","motif.pct.2","motif.p_val_adj","cluster","gene","motif")
write.csv(da_motifs,"signac_output/chromvar_DEmotif_minpct0.1.csv")

motif_gene <- ConvertMotifID(novaseq, id = rownames(novaseq@assays$chromvar@meta.features))
write.csv(data.frame(motif = rownames(novaseq@assays$chromvar@meta.features), gene = motif_gene),
          "signac_output/chromvar_motif_to_gene.csv")

library(SeuratData)
library(SeuratDisk)
SaveH5Seurat(novaseq, filename = "signac_output/chromvar.h5Seurat")

Convert("signac_output/chromvar.h5Seurat", dest = "h5ad", assay='chromvar')