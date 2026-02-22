module purge
module load miniconda
conda activate istar
module load NCCL/2.16.2-GCCcore-12.2.0-CUDA-12.0.0
module load cuDNN/8.8.0.121-CUDA-12.0.0
module load UCX-CUDA/1.13.1-GCCcore-12.2.0-CUDA-12.0.0
module load magma/2.7.1-foss-2022b-CUDA-12.0.0
module list

prefix="MALT100"
set -e
device="cuda"  # "cuda" or "cpu"
pixel_size=0.5  # desired pixel size for the whole analysis
n_genes=20000  # number of most variable genes to impute
echo $pixel_size > ${prefix}pixel-size.txt


#Unless specified, .py scripts were adapted from the original iSTAR package
python rescale.py ${prefix} --image

python preprocess.py ${prefix} --image

python extract_features.py ${prefix} --device=${device} --random-weights

#we provide the png file of the mask
python select_genes.py --n-top=${n_genes} "${prefix}cnts.tsv" "${prefix}gene-names.txt"

python rescale.py ${prefix} --locs --radius

python impute.py ${prefix} --epochs=400 --device=${device}  # train model from scratch

#Note: we changed the original visual.py file
python cluster.py --filter-size=12 --min-cluster-size=20 --n-clusters=20 --mask=${prefix}mask-small.png ${prefix}embeddings-gene.pickle ${prefix}clusters-gene/

echo "~~~ visualize imputed gene expression ~~~"
python plot_imputed.py ${prefix}

# differential analysis by clusters
echo "~~~ differential analysis by clusters ~~~"
python aggregate_imputed.py ${prefix}

python reorganize_imputed.py ${prefix}

python differential.py ${prefix}


##Note: we modified the original plotting function
python plot_imputed-with-color-bar.py ${prefix}


##finally, let's convert feature names from genome bins to genes using Homer:
#python code to generate a homer input:

#input_file = "MALT100/gene-names.txt"
#output_file = "homer_input.txt"

#with open(input_file, "r") as infile, open(output_file, "w") as outfile:
#    for i, line in enumerate(infile):
#        line = line.strip()
#        # Parse chromosome, start, and end positions
#        chrom, pos = line.split(":")
#        start, end = pos.split("-")
#
#        # Format the output line
#        peak_id = f"peak_{i+1}"       # Unique Peak ID
#        strand = "0"                  # Set strand to "+" (0) by default
#
#        # Write to output file
#        outfile.write(f"{peak_id}\t{chrom}\t{start}\t{end}\t{strand}\n")

annotatePeaks.pl homer_input.txt hg38 > homer_output.txt

sort -k1.6n homer_output.txt -o homer_output_sort.txt

# Skip the header in homer_input_sorted.txt and add gene-names.txt as the last column
tail -n +2 homer_output_sort.txt | paste - MALT100/gene-names.txt > homer_input_final.txt

# Add the header back to the new file
head -n 1 homer_output_sort.txt | awk '{print $0 "\tOriginal_Genomic_Position"}' > temp_header.txt
cat temp_header.txt homer_input_final.txt > homer_output_sort_use.txt
rm temp_header.txt