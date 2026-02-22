#After finishing processing

celltype <- read.csv("snapatac_out/mouse_brain_atac_cluster_use.csv")
#this is our meta output from python, including spatial coord and cluster ID

colnames(celltype)[1] <- 'cellid'

matches <- match(rownames(novaseq@meta.data), celltype$cellid)
novaseq@meta.data$cluster_use <- celltype$cluster_use[matches]
# Replace any NAs in the new column with "NA"
novaseq@meta.data$cluster_use[is.na(novaseq@meta.data$cluster_use)] <- "NA"

matches <- match(rownames(novaseq@meta.data), celltype$cellid)
novaseq@meta.data$coord <- celltype$coord[matches]
novaseq@meta.data$coord[is.na(novaseq@meta.data$coord)] <- "NA"

novaseq@meta.data$x_coord <- celltype$x_coord[matches]
novaseq@meta.data$x_coord[is.na(novaseq@meta.data$x_coord)] <- "NA"
novaseq@meta.data$y_coord <- celltype$y_coord[matches]
novaseq@meta.data$y_coord[is.na(novaseq@meta.data$y_coord)] <- "NA"

novaseq <- subset(novaseq, subset = cluster_use != 'NA')

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style
seqlevelsStyle(annotations) <- 'UCSC'

#seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(novaseq) <- annotations

# compute nucleosome signal score per cell
novaseq <- NucleosomeSignal(object = novaseq)

# compute TSS enrichment score per cell
novaseq <- TSSEnrichment(object = novaseq, fast = FALSE)

#Fig1 TSS plot
p=TSSPlot(novaseq)+theme(axis.title = element_text(size = 16),
                      axis.text = element_text(size = 15))+labs(title=NULL)

ggsave(filename = 'signac_output/R_mouse_brain_atac_TSSplot.png', plot = p, width = 7, height = 7, units = "in", dpi = 300)


total_fragments <- CountFragments('/gpfs/gibbs/pi/fan/hl737/Bo_FFPE_ATAC/COATAC/fragments.tsv.gz')
rownames(total_fragments)<-total_fragments$CB
novaseq$fragments <- total_fragments[colnames(novaseq), "frequency_count"]
novaseq$reads_count <- total_fragments[colnames(novaseq), "reads_count"]
novaseq$mononucleosomal <- total_fragments[colnames(novaseq), "mononucleosomal"]
novaseq$nucleosome_free <- total_fragments[colnames(novaseq), "nucleosome_free"]

#Calculate FRiP and blacklist ratio

novaseq <- FRiP(
  object = novaseq,
  assay = 'peaks',
  total.fragments = 'fragments'
)

novaseq$blacklist_fraction <- FractionCountsInRegion(
  object = novaseq, 
  assay = 'peaks',
  regions = blacklist_mm10
)

#Violin plots in supp (across clusters)
p=VlnPlot(object = novaseq,features = c('FRiP'),pt.size = 0.02, group.by='cluster_use',cols=cmap)+ 
  theme(legend.position = "none")+labs(title='FRiP')
ggsave(filename = 'signac_output/R_mouse_brain_atac_FRiP_cluster.png', plot = p, width = 8, height = 4, units = "in", dpi = 300)


p=VlnPlot(object = novaseq,features = c('fragments'),pt.size = 0.02, group.by='cluster_use',cols=cmap)+ 
  theme(legend.position = "none")+labs(title='Fragment count')
ggsave(filename = 'signac_output/R_mouse_brain_atac_fragments_cluster.png', plot = p, width = 8, height = 4, units = "in", dpi = 300)

Idents(novaseq)<-novaseq$cluster_use

levels(novaseq)<-c('0','1','2','3','4','5','6','7','8','9')

cmap<-c('0'='steelblue',
'1'= 'orange',
'2'= '#279e68',
'3'= '#d62728',
'4'= 'darkviolet',
'5'= '#8c564b',
'6'= 'navy',
'7'= '#b5bd61',
'8'= '#17becf',
'9'= '#aec7e8')

#An example of coverage plot
options(repr.plot.width = 8, repr.plot.height = 6, repr.plot.res = 200)
CoveragePlot(
  object = novaseq, region = "Apc",
  extend.upstream = 5000, extend.downstream = -80000) & scale_fill_manual(values = cmap)