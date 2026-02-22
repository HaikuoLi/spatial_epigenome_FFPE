options(bitmapType='cairo')
library(infercnv)

# a gene order file was generated from the gene matrix
gene_order <- read.table("gene_matrix_gene_order_file.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(gene_order) <- c("gene", "chromosome", "start", "end")

# Define the chromosome ordering
chromosome_order <- c(paste0("chr", 1:22), "chrX", "chrY")

# Sort the data by the chromosome order
gene_order_sorted <- gene_order[order(factor(gene_order$chromosome, levels = chromosome_order)), ]

# Save the sorted file
write.table(gene_order_sorted, file = "gene_matrix_filtered_sorted_gene_order_file.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

filtered_genes <- read.table("gene_matrix_gene_order_file.txt", header = FALSE, sep = "\t")
gene_list <- filtered_genes[[1]]
subset_matrix <- gene_matrix[rownames(gene_matrix) %in% gene_list, ]


#gene_matrix_annotations_file.txt has two columns - cell ID and leiden ID
infercnv_obj = CreateInfercnvObject(delim = '\t',
                                    chr_exclude = c("chrX", "chrY", "chrM"),
                                    raw_counts_matrix = subset_matrix,
                                    annotations_file = 'gene_matrix_annotations_file.txt',
                                    gene_order_file = 'gene_matrix_filtered_sorted_gene_order_file.txt',
                                    ref_group_names = NULL
)
infercnv_obj = infercnv::run(infercnv_obj,
                             num_threads = 16, cutoff=0.1,
                             out_dir='result_final/',
                             cluster_by_groups=TRUE, denoise=TRUE,
                             HMM=TRUE)

aa = readRDS('result_final/run.final.infercnv_obj')
plot_cnv(infercnv_obj=aa, out_dir='result_final/',output_filename='used_this_run.final',
         k_obs_groups=1, cluster_by_groups=TRUE, title='', output_format= "pdf",
         obs_title = "", contig_cex = 2,
         cluster_references=TRUE, plot_chr_scale=FALSE, chr_lengths=FALSE,
         x.center=1, x.range=c(0.9,1.1), png_res=300, useRaster=TRUE)

