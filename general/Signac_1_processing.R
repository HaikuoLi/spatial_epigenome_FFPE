library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)

fragments<-CreateFragmentObject('fragments.tsv.gz')

features<-CallPeaks(fragments,format='BED',
                    outdir = "signac_output/MACS2/",
                    name='bulk', effective.genome.size = 2.5e+09,
                    macs2.path ='/gpfs/gibbs/pi/fan/hl737/conda/py27/bin/macs2',cleanup=FALSE)

matrix<-FeatureMatrix(fragments,features,process_n = 10000)

chrom_assay <- CreateChromatinAssay(
  counts = matrix,
  genome = 'mm10',
  fragments = fragments,
    min.cells = 1, min.features = 1, ranges=features)

novaseq <- CreateSeuratObject(
  counts = chrom_assay,
  assay = 'peaks')
