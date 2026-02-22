##in linux:
#tail -n +2 MALT_clustering_meta.csv | cut -d',' -f1 > 8444cells.txt
#zgrep -Ff 8444cells.txt MALT100.chromap.sorted.bed.gz | gzip > 8444subset.bed.gz
#gunzip -c 8444subset.bed.gz > 8444subset.bed


#R-bundle-Bioconductor/3.19-foss-2022b-R-4.4.1
library(ggdendro)
library(devtools)
library(BiocManager)
library(GenomicAlignments)
library(SummarizedExperiment)
library(plyranges)
library(Rsamtools)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(Biostrings)
library(BiocGenerics)
library(S4Vectors)
library(GenomicFeatures)
#so we have all these packages

install_github("colomemaria/epiAneufinder")

library(epiAneufinder)
#epiAneufinder_1.0.3

epiAneufinder(input="8444subset.bed",
              outdir="epiAneufinder_8444subset_500k_sbatch",
              blacklist="hg38-blacklist.v2.bed",
              windowSize=5e5, 
              genome="BSgenome.Hsapiens.UCSC.hg38",
              exclude=c('chrX','chrY','chrM'), 
              reuse.existing=TRUE,
              title_karyo="Karyogram of sample data", 
              ncores=20,
              minFrags=4000,
              minsizeCNV=0,
              k=4,
              plotKaryo=TRUE)


#Load result table
res_table<-read.table("epiAneufinder_results/results_table.tsv")
#Need to reformat column names as dash are replaced by points in R
colnames(res_table)<-gsub("\\.","-",colnames(res_table))


subclones<-split_subclones(res_table,tree_depth=8,
                           plot_tree=TRUE,
                           plot_path="epiAneufinder_results/8subclones.pdf",
                           plot_width=9,
                           plot_height=4)


annot_dt<-data.frame(cell=subclones$cell,
                     annot=paste0("Clone",subclones$subclone))


#Reformat somy dataframe
library(data.table)
res_table<-as.data.table(res_table)
somies.dt <- res_table[,-c("start","end")]
somies.dt$rn <- 1:nrow(somies.dt)
somies_melted <- melt(somies.dt, id.vars=c('rn','seq'))
somies_melted$value <- as.factor(paste0(somies_melted$value,'-somy'))

#Sort the karyogram chromosomes correctly (assuming that the order in res_table is correct)
somies_melted$seq<-factor(somies_melted$seq,levels=unique(res_table$seq))

library(ggplot2)
library(ggdendro)
library(cowplot)
#Create a hierarchical clustering and plot clustering
counts_t <- t(somies.dt[ ,.SD, .SDcols=patterns('cell-')])
dist_matrix <- dist(counts_t)
dist_matrix[is.na(dist_matrix)] <- 0
hc_counts <- hclust(dist_matrix, method = "ward.D")
ord <- hc_counts$order
dhc <- stats::as.dendrogram(hc_counts)
ddata <- ggdendro::dendro_data(dhc)
ggdndr <- ggplot(ddata$segments) + 
  geom_segment(aes(x = x, xend = xend, y = y, yend = yend)) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_reverse()
ggdndr <- ggdndr + coord_flip() + labs(title = "", subtitle = "")
ggdndr <- ggdndr + theme(panel.background=element_blank(), 
                         axis.ticks=element_blank(), 
                         axis.text=element_blank(), 
                         axis.line=element_blank(), 
                         axis.title=element_blank(),
                         panel.grid.minor=element_blank(),
                         plot.title = element_text(size=18),
                         plot.subtitle=element_text(size=12))

#Make sure the karyogram is ordered according to the dendogram
somies_melted$variable <- factor(somies_melted$variable,
                                 levels = rownames(counts_t)[ord])

#Also order the cell type annotations accordingly
annot_dt$cell<-factor(as.character(annot_dt$cell),
                      levels=rownames(counts_t)[ord])


somycolours <- c(`0-somy` = "royalblue4",
                 `1-somy` = "gray90",
                 `2-somy` = "red3")

ggsomy <- ggplot(somies_melted, aes(x=rn, y=variable, fill=value)) + 
  geom_tile() +
  facet_grid(cols=vars(seq), scales = 'free_x', space = 'free') +
  labs(x="", fill=NULL) +  # Remove legend title
  scale_fill_manual(values=somycolours, labels = c("loss", "normal", "gain")) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size=18),
        strip.text.x = element_text(size = 15, angle = 90),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())

library(RColorBrewer)
annot_dt$type<-"Clone"
ggannot<-ggplot(annot_dt, aes(x=cell, y=1, fill=annot)) +
  geom_tile() +
  coord_flip() +
  facet_grid(~type, scales = 'free', space = 'free') +
  labs(title = "", fill="Karyotype") +
  theme(legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        plot.title = element_text(size=18)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_fill_manual(values = c('#d6b09c', '#e41a1c', '#f781bf', '#984ea3', 
                               '#ff7f00', '#ffff33', '#377eb8', '#4daf4a'), levels(annot_dt$annot))

gglegend<-get_legend(ggannot)

ggannot<-ggannot+
  ylab("")+
  theme(axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, angle = 90),
        legend.position="none")

library(patchwork)


##get a new column for spatial clusters
annot_dt2<-annot_dt
annot_dt2$celltype<-"Clusters"
annot_dt2$type<-NULL
annot_dt2$cell_id <- gsub("cell-", "", annot_dt2$cell)

metadata <- read.csv("MALT_clustering_meta.csv")

annot_dt2 <- merge(annot_dt2, metadata[, c("X", "cluster_use")], by.x = "cell_id", by.y = "X", all.x = TRUE)
annot_dt2$cell_id<-NULL

ggannot2<-ggplot(annot_dt2, aes(x=cell, y=1, fill=annot)) +
  geom_tile() +
  coord_flip() +
  facet_grid(~cluster_use, scales = 'free', space = 'free') +
  labs(title = "", fill="Karyotype") +
  theme(legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        plot.title = element_text(size=18)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_fill_manual(values = c('#d6b09c', '#e41a1c', '#f781bf', '#984ea3', 
                               '#ff7f00', '#ffff33', '#377eb8', '#4daf4a'), levels(annot_dt$annot))
gglegend2<-get_legend(ggannot2)
ggannot2<-ggannot2+
  ylab("")+
  theme(axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, angle = 90),
        legend.position="none")



#this aligns better
combiplot <- ggdndr + ggsomy + ggannot + ggannot2 + 
  plot_layout(ncol = 4, widths = c(0.08, 1, 0.02, 0.15))

combiplot <- cowplot::plot_grid(combiplot, gglegend, ncol = 1, rel_heights = c(1, 0.05))+ theme(
  panel.background = element_rect(fill = "white", color = "white"),
  plot.background = element_rect(fill = "white", color = "white")
)

ggsave('epiAneufinder_results/karyo_annotated_8clones_clusters_use.png', combiplot, width = 30, height=16, units = "in")


#annot_dt2$cell[annot_dt2$cluster_use=='C1']
#outdir<-'20241105_epiAneufinder_8444subset_bin1M_sbatch/epiAneufinder_results/plot_single_cell_profile/'
#plot_single_cell_profile('epiAneufinder_results/',
#                         threshold_blacklist_bins=0.85,
#                         selected_cells=annot_dt2$cell[annot_dt2$cluster_use=='C1'][1:3],
#                         plot_path=paste0(outdir,"C1_selected.pdf"),
#                         plot_width=25,
#                         plot_height=15)

count_summary <- readRDS('epiAneufinder_results/count_summary.rds')
cnv_calls <- readRDS('epiAneufinder_results/cnv_calls.rds')
counts_gc_corrected <- readRDS('epiAneufinder_results/counts_gc_corrected.rds')
results_gc_corrected <- readRDS('epiAneufinder_results/results_gc_corrected.rds')

write.csv(annot_dt, "8444_bin500k_8clone_annot_dt.csv", row.names = FALSE)
write.csv(counts_gc_corrected, "8444_bin500k_8clone_counts_gc_corrected.csv", row.names = FALSE)


