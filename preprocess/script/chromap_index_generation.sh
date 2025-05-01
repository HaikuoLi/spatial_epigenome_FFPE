#By Haikuo Li
https://www.gencodegenes.org/mouse/release_M23.html
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
#gunzip

conda activate spatial-atac
chromap -i -r /gpfs/gibbs/pi/fan/userid/chromap_index/GRCm38/GRCm38.primary_assembly.genome.fa -o /gpfs/gibbs/pi/fan/userid/chromap_index/GRCm38/GRCm38.primary_assembly.genome.index