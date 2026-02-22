library(ChIPseeker)
library(EpiTrace)

##load a Seurat object called novaseq, containing 102293 peaks (following the Signac processing)
load('signac_output/MALT_8444x102293.Rdata')
novaseq_initiated_peaks <- Init_Peakset(novaseq@assays$peaks@ranges)
novaseq_initiated_peaks_df <- as.data.frame(novaseq_initiated_peaks,row.names = NULL)

paste0(novaseq_initiated_peaks_df$seqnames,':',novaseq_initiated_peaks_df$start,'-',novaseq_initiated_peaks_df$end) -> novaseq_initiated_peaks_df$peakName

as(novaseq@assays$peaks@data, "sparseMatrix") -> novaseq_mtx

initiated_mm <- Init_Matrix(cellname = colnames(novaseq_mtx),
                            peakname = novaseq_initiated_peaks_df$peakName,
                            matrix = novaseq_mtx)

library(dplyr)
library(tidyr)
library(readr)
library(GenomicRanges)
library(reshape2)
library(openxlsx)
library(ggplot2)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(ggtree)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#library(ArchR)

initiated_peaks <- NULL
clock_gr <- NULL
original_clk_peakset <- NULL



EpiTraceAge_Convergence <- function (peakSet, matrix, celltype = NULL, min.cutoff = 50,lsi_dim = 2:50, fn.k.param = 21, ref_genome = "hg38", sep_string = c(":","-"), clock_gr = plyranges::reduce_ranges(c(clock_gr_list[[1]],clock_gr_list[[2]])), non_standard_clock = F, qualnum = 10,Z_cutoff = 3, mean_error_limit = 0.01, ncore_lim = 12, parallel = T,iterative_time = 2,remove_peaks_number=10,normalization_method='randomized',select_minimal_percentage = 0.05, select_size_of_dispersion = 3000){
  norm_meth = normalization_method
  original_clk_peakset <- clock_gr
  if (ref_genome %in% "hg38") {
    original_clk_peakset <- easyLift::easyLiftOver(original_clk_peakset,
                                                   "hg19_hg38")
  }
  if (ref_genome != "hg19" & ref_genome != "hg38") {
    message("please make double sure your ref genome, peak set and cells are similar.")
  }
  iterative_count = 1
  na_vector_current <- ncol(matrix)
  na_vector_previous <- ncol(matrix)
  error <- rep(1, ncol(matrix))
  mean_error <- mean(error)
  iterative_GR_list <- list(iterative = original_clk_peakset)
  overlap_with_clk <- findOverlaps(peakSet, original_clk_peakset)@from %>%
    unique()
  initial_matrix_clk <- matrix[overlap_with_clk, ]
  initial_peakSet_clk <- peakSet[overlap_with_clk, ]
  message("Preparing obj...")
  epitrace_obj <- EpiTrace_prepare_object(initial_peakSet_clk,
                                          initial_matrix_clk, celltype, ref_genome = "hg38", non_standard_clock = T,
                                          clock_gr_list = iterative_GR_list, sep_string = sep_string,
                                          fn.k.param = fn.k.param, lsi_dim = lsi_dim, qualnum = qualnum,
                                          min.cutoff = min.cutoff, run_reduction = T,remove_peaks_number=remove_peaks_number)
  message("Finished prepare obj")
  epitrace_obj_original_metadata <- epitrace_obj@meta.data
  message("Estimating age for initialization...")
  if (parallel == F) {
    epitrace_obj_age_estimated <- RunEpiTraceAge(epitrace_obj,normalization_method=norm_meth,select_minimal_percentage=select_minimal_percentage,select_size_of_dispersion=select_size_of_dispersion) %>%
      suppressMessages() %>% suppressWarnings()
  } else {
    epitrace_obj_age_estimated <- RunEpiTraceAge(epitrace_obj,
                                                 ncores = ncore_lim, parallel = T,normalization_method=norm_meth,select_minimal_percentage=select_minimal_percentage,select_size_of_dispersion=select_size_of_dispersion) %>% suppressMessages() %>%
      suppressWarnings()
  }
  message("Finished age estimation")
  age_current <- epitrace_obj_age_estimated@meta.data$EpiTraceAge_iterative
  names(age_current) <- epitrace_obj_age_estimated@meta.data$cell
  na_vector_current <- is.na(age_current)
  next_peakset <- iterative_GR_list$iterative
  cell_current <- epitrace_obj_age_estimated@meta.data$cell[!is.na(epitrace_obj_age_estimated@meta.data$EpiTraceAge_iterative)]
  age_current <- age_current[cell_current]
  to_be_associated_mtx <- matrix[,cell_current]
  while ((sum(na_vector_current) <= sum(na_vector_previous)) &
         (mean_error >= mean_error_limit) & (iterative_count <=
                                             iterative_time)) {
    message("Iterating ", iterative_count)
    age_previous <- age_current
    na_vector_previous <- na_vector_current
    iterative_count = iterative_count + 1
    cell_current <- cell_current[!na_vector_current] %>% na.omit() 
    age_current <- age_current[cell_current]
    to_be_associated_mtx <- matrix[,cell_current]
    message("Performing Corr ")
    if (parallel == T) {
      res_list_cor <- parallel::mclapply(c(1:ceiling(dim(to_be_associated_mtx)[1]/1000)),
                                         mc.cores = ncore_lim, function(x) {
                                           target_mtx <- t(to_be_associated_mtx[((1000 *
                                                                                    (x - 1)) + 1):min(dim(to_be_associated_mtx)[1],
                                                                                                      1000 * x), ]) %>% as.matrix()
                                           WGCNA::cor(x = target_mtx, y = age_current)
                                         })
      cor_res_PT <- res_list_cor %>% unlist
    }
    else {
      res_list_cor <- lapply(c(1:ceiling(dim(to_be_associated_mtx)[1]/1000)),
                             function(x) {
                               target_mtx <- t(to_be_associated_mtx[((1000 *
                                                                        (x - 1)) + 1):min(dim(to_be_associated_mtx)[1],
                                                                                          1000 * x), ]) %>% as.matrix()
                               WGCNA::cor(x = target_mtx, y = age_current)
                             })
      cor_res_PT <- res_list_cor %>% unlist
    }
    cor_res_PT[is.na(cor_res_PT)] <- 0
    scaled_cor_res_PT <- scale(cor_res_PT)
    message("Finished Corr ")
    updated_peakset <- peakSet
    updated_peakset$correlation_of_EpiTraceAge <- cor_res_PT
    updated_peakset$scaled_correlation_of_EpiTraceAge <- scaled_cor_res_PT
    peaks_overlap_with_clock <- findOverlaps(updated_peakset,
                                             original_clk_peakset)@from %>% unique
    updated_peakset$locus_type <- "Others"
    updated_peakset$locus_type[peaks_overlap_with_clock] <- "Chronology"
    iterative_clock_gr_list <- list(iterative = updated_peakset[(updated_peakset$scaled_correlation_of_EpiTraceAge >
                                                                   Z_cutoff) | (updated_peakset$locus_type %ni% "Others") %>%
                                                                  as.vector() | (updated_peakset$peakId %in% next_peakset$peakId),
    ])
    overlap_with_clk <- findOverlaps(peakSet, iterative_clock_gr_list[[1]])@from %>%
      unique()
    initial_matrix_clk <- matrix[overlap_with_clk, ]
    initial_peakSet_clk <- peakSet[overlap_with_clk, ]
    message("Calculate iterative age")
    if (parallel == T) {
      message("Parallel run")
      tryCatch({
        nrow(initial_matrix_clk) %>% as.numeric() -> row_num 
        ncol(initial_matrix_clk) %>% as.numeric() -> col_num
        if ((row_num*col_num) >
            (50000 * 2000)) {
          batch_size = 400
        }
        else {
          batch_size = 1000
        }
      },error=function(e){batch_size=400})
      res1 <- EpiTraceAge(initial_matrix_clk, parallel = parallel,
                          ncores = ncore_lim, size = batch_size, subsamplesize = batch_size,normalization_method=normalization_method,select_minimal_percentage=select_minimal_percentage,select_size_of_dispersion=select_size_of_dispersion) %>%
        suppressMessages() %>% suppressWarnings()
    }
    else {
      message("Single thread run")
      res1 <- EpiTraceAge(initial_matrix_clk, parallel = F,
                          ncores = 1,normalization_method=normalization_method,select_minimal_percentage=select_minimal_percentage,select_size_of_dispersion=select_size_of_dispersion) %>% suppressMessages() %>% suppressWarnings()
    }
    colnames(res1)[2:4] <- c("EpiTraceAge_iterative", "Accessibility_iterative",
                             "AccessibilitySmooth_iterative")
    message("Update dataframe")
    epitrace_obj_original_metadata_update <- left_join(epitrace_obj_original_metadata,
                                                       res1)
    epitrace_obj_iterative_age_estimated <- epitrace_obj_age_estimated
    rownames(epitrace_obj_original_metadata_update) <- rownames(epitrace_obj_original_metadata)
    epitrace_obj_iterative_age_estimated@meta.data <- epitrace_obj_original_metadata_update
    age_current <- epitrace_obj_iterative_age_estimated$EpiTraceAge_iterative
    cell_current <- epitrace_obj_iterative_age_estimated$cell[!is.na(epitrace_obj_iterative_age_estimated$EpiTraceAge_iterative)]
    names(age_current) <- cell_current
    error <- (age_current - age_previous[cell_current])
    tryCatch({
      error[is.infinite(error)] <- 0
    }, error = function(e) {
      message("")
    })
    tryCatch({
      error[is.na(error)] <- 1
    }, error = function(e) {
      message("")
    })
    mean_error = mean(abs(error), na.rm = T)
    na_vector_current <- is.na(age_current)
    age_previous <- age_current
    epitrace_obj_iterative_age_estimated@misc$iterative_count <- iterative_count
    epitrace_obj_iterative_age_estimated@misc$mean_error <- mean_error
    epitrace_obj_iterative_age_estimated@misc$next_ref_peak_size <- length(iterative_clock_gr_list$iterative)
    message("mean_error = ", mean_error)
    next_peakset <- iterative_clock_gr_list$iterative
  }
  epitrace_obj_iterative_age_estimated$EpiTraceAge_Clock_initial <- epitrace_obj_age_estimated@meta.data$EpiTraceAge_iterative
  epitrace_obj_iterative_age_estimated$Accessibility_initial <- epitrace_obj_age_estimated@meta.data$Accessibility_iterative
  epitrace_obj_iterative_age_estimated$AccessibilitySmooth_initial <- epitrace_obj_age_estimated@meta.data$AccessibilitySmooth_iterative
  return(epitrace_obj_iterative_age_estimated)
}
#in EpiTrace_prepare_object change to hg38
#non_standard_clock must be T in EpiTrace_prepare_object to let it run

epitrace_obj_age_estimated <- EpiTraceAge_Convergence(novaseq_initiated_peaks, initiated_mm, celltype = novaseq$cluster_use, 
                                                      Z_cutoff = 2.5, iterative_time = 10, parallel = F,
                                                      qualnum = 5, remove_peaks_number = 5, 
                                                      ref_genome = 'hg38', non_standard_clock = F,
                                                      sep_string = c(":","-"))


data_clean <- epitrace_obj_age_estimated@meta.data[!is.na(epitrace_obj_age_estimated@meta.data$EpiTraceAge_iterative), ]
data_clean %>% group_by(celltype) %>% summarise(avg_EpiTraceAge = mean(EpiTraceAge_iterative, na.rm = TRUE))

data_clean <- epitrace_obj_age_estimated@meta.data[!is.na(epitrace_obj_age_estimated@meta.data$EpiTraceAge_Clock_initial ), ]
data_clean %>% group_by(celltype) %>% summarise(avg_EpiTraceAge = mean(EpiTraceAge_Clock_initial , na.rm = TRUE))

data_clean <- epitrace_obj_age_estimated@meta.data[!is.na(epitrace_obj_age_estimated@meta.data$EpiTraceAge_iterative), ]
write.csv(data_clean, "signac_output/EpiTrace_SignacPeak_Convergence_out.csv", row.names = FALSE)




#Compute association between regulatory region accessibility and cell age

associated_res_MALT <- AssociationOfPeaksToAge(subset(epitrace_obj_age_estimated,celltype %in% c('C1','C5')),epitrace_age_name = "EpiTraceAge_iterative")
#in R 4.3: BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
associated_res_MALT <- separate(associated_res_MALT,col='locus',into=c('chr','start','end'),remove=F,convert=T)
associated_res_MALT_gr <- makeGRangesFromDataFrame(associated_res_MALT)
findOverlaps(associated_res_MALT_gr,plyranges::reduce_ranges(c(clock_gr_list[[1]],clock_gr_list[[2]])))@from %>% unique -> peaks_overlap_with_clock
associated_res_MALT_gr$locus_type <- 'Others'
associated_res_MALT_gr$locus_type[peaks_overlap_with_clock] <- 'Clock-like DML'
save(epitrace_obj_age_estimated,associated_res_MALT,associated_res_MALT_gr,txdb,
     file='signac_output/EpiTrace_Convergence.RData')
write.csv(associated_res_MALT, "signac_output/EpiTrace_SignacPeak_Convergence_correlation.csv", row.names = FALSE)
