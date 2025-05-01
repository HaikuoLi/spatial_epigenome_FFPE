#!/bin/bash
#SBATCH --partition=day 
#SBATCH --job-name=yourjobname
#SBATCH --ntasks=1 --cpus-per-task=12
#SBATCH --mem=36g
#SBATCH --time=3:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haikuo.li@yale.edu

#change directories and parameters below
raw_folder=/gpfs/gibbs/pi/fan/userid/FFPE_epi/epiPatho/raw/01.RawData/MB/
fastq_intput_1=$raw_folder/data_1.fq.gz
fastq_intput_2=$raw_folder/data_2.fq.gz

output_folder=/gpfs/gibbs/pi/fan/userid/FFPE_epi/epiPatho/process/MB/
out_file=/gpfs/gibbs/pi/fan/userid/FFPE_epi/epiPatho/process/MB/MB.bed
chromap_input_2=/gpfs/gibbs/pi/fan/userid/FFPE_epi/epiPatho/process/MB/bbduk/bbduk_L2_R1.fastq.gz

BC_process_script=/home/userid/project/20240503_ATAC_practice/process/BC_process_BCB_wUMI.py
index_file=/gpfs/gibbs/pi/fan/userid/chromap_index/GRCm38/GRCm38.primary_assembly.genome.index
fa_file=/gpfs/gibbs/pi/fan/userid/chromap_index/GRCm38/GRCm38.primary_assembly.genome.fa
bc_file=/home/userid/palmer_pi_fan/20240614_2500barcode.txt


###########-----------------###########
#usually no need to change below
date
module purge
module load miniconda
conda activate spatial-atac
module load Java

###1. Linker filter: Filter L1 and then filter L2
mkdir -p $output_folder
bbduk_folder=$output_folder/bbduk
mkdir -p $bbduk_folder

/gpfs/gibbs/pi/fan/userid/Downloads/bbmap/bbduk.sh in1=$fastq_intput_1 \
in2=$fastq_intput_2 \
outm1=$bbduk_folder/bbduk_L1_R1.fastq.gz \
outm2=$bbduk_folder/bbduk_L1_R2.fastq.gz \
k=30 mm=f rcomp=f restrictleft=108 skipr1=t hdist=3 \
stats=$bbduk_folder/bbduk_stats_L1.txt \
threads=$core literal=GTGGCCGATGTTTCGCATCGGCGTACGACT

/gpfs/gibbs/pi/fan/userid/Downloads/bbmap/bbduk.sh in1=$bbduk_folder/bbduk_L1_R1.fastq.gz \
in2=$bbduk_folder/bbduk_L1_R2.fastq.gz \
outm1=$bbduk_folder/bbduk_L2_R1.fastq.gz \
outm2=$bbduk_folder/bbduk_L2_R2.fastq.gz \
k=30 mm=f rcomp=f restrictleft=70 skipr1=t hdist=3 \
stats=$bbduk_folder/bbduk_stats_L2.txt \
threads=$core literal=ATCCACGTGCTTGAGAGGCCAGAGCATTCG

rm -r $bbduk_folder/bbduk_L1_R1.fastq.gz $bbduk_folder/bbduk_L1_R2.fastq.gz


###2. BC_process.
date
chromap_input_folder=$output_folder/chromap_input/
mkdir -p $chromap_input_folder
python $BC_process_script \
--input $bbduk_folder/bbduk_L2_R2.fastq.gz \
--output_R1 $chromap_input_folder/sample_S1_L001_R1_001.fastq \
--output_R2 $chromap_input_folder/sample_S1_L001_R2_001.fastq

gzip $chromap_input_folder/sample_S1_L001_R1_001.fastq
#117:
gzip $chromap_input_folder/sample_S1_L001_R2_001.fastq
#16bp bc
rm $bbduk_folder/bbduk_L2_R2.fastq.gz


###3. Run chromap
date
cd $output_folder

chromap --preset atac -x $index_file -r $fa_file -1 $chromap_input_folder/sample_S1_L001_R1_001.fastq.gz -2 $chromap_input_2 \
-b $chromap_input_folder/sample_S1_L001_R2_001.fastq.gz --barcode-whitelist $bc_file -t 12 -o $out_file
