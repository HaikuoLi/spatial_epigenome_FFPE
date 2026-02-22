load(system.file(package = "circlize", "extdata", "DMR.RData"))
library(readr)

##check Python_4_visualization for how to split the bulk fragment bed into bed for each cluster.


file_path_0 <- "0.bed.gz"
bed_data <- read_delim(file_path_0, delim = "\t", col_names = FALSE)
# Select the first three columns and assign appropriate column names
bed_df_0 <- bed_data[, 1:3]
colnames(bed_df_0) <- c("chr", "start", "end")
bed_df_0 <- bed_df_0[bed_df_0$chr != "chrY", ]
bed_df_0 <- bed_df_0[sample(nrow(bed_df_0), round(length(rownames(bed_df_0))/30)), ]

file_path_1 <- "1.bed.gz"
bed_data <- read_delim(file_path_1, delim = "\t", col_names = FALSE)
# Select the first three columns and assign appropriate column names
bed_df_1 <- bed_data[, 1:3]
colnames(bed_df_1) <- c("chr", "start", "end")
bed_df_1 <- bed_df_1[bed_df_1$chr != "chrY", ]
bed_df_1 <- bed_df_1[sample(nrow(bed_df_1), round(length(rownames(bed_df_1))/30)), ]

file_path_2 <- "2.bed.gz"
bed_data <- read_delim(file_path_2, delim = "\t", col_names = FALSE)
# Select the first three columns and assign appropriate column names
bed_df_2 <- bed_data[, 1:3]
colnames(bed_df_2) <- c("chr", "start", "end")
bed_df_2 <- bed_df_2[bed_df_2$chr != "chrY", ]
bed_df_2 <- bed_df_2[sample(nrow(bed_df_2), round(length(rownames(bed_df_2))/30)), ]

file_path_3 <- "3.bed.gz"
bed_data <- read_delim(file_path_3, delim = "\t", col_names = FALSE)
# Select the first three columns and assign appropriate column names
bed_df_3 <- bed_data[, 1:3]
colnames(bed_df_3) <- c("chr", "start", "end")
bed_df_3 <- bed_df_3[bed_df_3$chr != "chrY", ]
bed_df_3 <- bed_df_3[sample(nrow(bed_df_3), round(length(rownames(bed_df_3))/30)), ]

bed_data <- read_delim("4.bed.gz", delim = "\t", col_names = FALSE)
bed_df_4 <- bed_data[, 1:3]
colnames(bed_df_4) <- c("chr", "start", "end")
bed_df_4 <- bed_df_4[bed_df_4$chr != "chrY", ]
bed_df_4 <- bed_df_4[sample(nrow(bed_df_4), round(length(rownames(bed_df_4))/30)), ]

bed_data <- read_delim("5.bed.gz", delim = "\t", col_names = FALSE)
bed_df_5 <- bed_data[, 1:3]
colnames(bed_df_5) <- c("chr", "start", "end")
bed_df_5 <- bed_df_5[bed_df_5$chr != "chrY", ]
bed_df_5 <- bed_df_5[sample(nrow(bed_df_5), round(length(rownames(bed_df_5))/30)), ]

bed_data <- read_delim("6.bed.gz", delim = "\t", col_names = FALSE)
bed_df_6 <- bed_data[, 1:3]
colnames(bed_df_6) <- c("chr", "start", "end")
bed_df_6 <- bed_df_6[bed_df_6$chr != "chrY", ]
bed_df_6 <- bed_df_6[sample(nrow(bed_df_6), round(length(rownames(bed_df_6))/30)), ]


bed_list = list(bed_df_0,bed_df_1,bed_df_2,bed_df_3,bed_df_4, bed_df_5,bed_df_6)


pdf(file = "fragment_subset_C2.pdf", width = 6, height = 6)
circos.initializeWithIdeogram(chromosome.index = 'chr2')
circos.genomicDensity(bed_list[[1]], col='orange', track.height = 0.08)
circos.genomicDensity(bed_list[[2]], col='#279e68', track.height = 0.08)
circos.genomicDensity(bed_list[[3]], col='steelblue', track.height = 0.08)
circos.genomicDensity(bed_list[[4]], col='darkviolet', track.height = 0.08)
circos.genomicDensity(bed_list[[5]], col='#d62728', track.height = 0.08)
circos.genomicDensity(bed_list[[6]], col='#e377c2', track.height = 0.08)
circos.genomicDensity(bed_list[[7]], col='#8c564b', track.height = 0.08)
circos.clear()
dev.off()

## Pipeline described in: Li, H., Li, D. & Humphreys, B.D. Chromatin conformation and histone modification profiling across human kidney anatomic regions. Sci Data 11, 797 (2024). https://doi.org/10.1038/s41597-024-03648-8