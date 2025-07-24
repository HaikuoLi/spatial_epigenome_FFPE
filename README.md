# Spatial epigenomics in FFPE tissues
This repository documents scripts used for data analysis and visualization in an unpublished manuscript. Please access all compressed files with the reviewer access token. The GitHub page will be completely publicly accessible upon manuscript's acceptance.

![schematic](https://github.com/HaikuoLi/spatial_epigenome_FFPE/blob/main/workflow.jpeg)

## general - repository of scripts for general data processing and analysis
1. Fragment file processing and spatial barcode mapping in Python (SnapATAC2)
2. Fragment data QC and clustering analysis in Python
3. Gene activity computation and differential analysis in Python
4. General data visualization scripts
5. Fragment file processing and peak calling in R (Signac)
6. Fragment data QC and visualization in R
7. Motif analysis with chromvar in R
8. Visuazation of SPOTlight deconvolution outputs
9. Cross-sample comparison related to Figure 1

## lymphoma - repository of scripts for human lymphoma data analysis
1. Tissue architecture super-resolution inference using iStar, including input file preparation in Python, iStar bash scripts and a mask tissue image
2. Copy number variation inference using epiAneufinder, including epiAneufinder analysis in R and downstream analysis in Python
3. Cell cycle scoring and gene module scoring analysis
4. Trajectory inference with Monocle2 and Monocle3
5. Mitotic age inference with Epitrace, including Epitrace analysis in R and downstream analysis in Python
6. Analysis of Patho-DBiT (spatial-RNA-seq) data, including dimension reduction/clustering in Scanpy and inferCNV analysis in R
7. Visualization of fragment data with circos
