import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import scanpy as sc
sc.settings.verbosity = 3
sc.logging.print_header()
from matplotlib import rcParams

#https://github.com/scverse/scanpy_usage/tree/master/180209_cell_cycle
cell_cycle_genes = [x.strip() for x in open('regev_lab_cell_cycle_genes.txt')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]

cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)


# Initialize an empty list to collect data
data = []

# Loop through clusters C0 to C8
for cluster in [f'C{i}' for i in range(9)]:
    # Get counts of phases (G1, G2M, S) for each cluster
    counts = pd.value_counts(adata.obs[adata.obs['cluster_use'] == cluster]['phase'])
    
    # Add the counts to the list, making sure to handle any missing phases
    for phase in ['G1', 'G2M', 'S']:
        count = counts.get(phase, 0)  # Use 0 if the phase is not present
        data.append({'cluster': cluster, 'phase': phase, 'count': count})

# Create a DataFrame from the collected data
df = pd.DataFrame(data)

# Group by 'cluster' and calculate the sum of counts for each cluster
cluster_sums = df.groupby('cluster')['count'].transform('sum')

# Divide each count by the sum for its corresponding cluster to get proportions
df['proportion'] = df['count'] / cluster_sums

# Pivot the DataFrame to create a better structure for plotting
df_pivot = df.pivot(index='cluster', columns='phase', values='proportion').fillna(0)

# Create a stacked bar plot
df_pivot.plot(kind='bar', stacked=True, figsize=(5,4), color=['#1f77b4', '#ff7f0e', '#2ca02c'])

# Add labels and title
plt.ylabel('Fraction', fontsize=11)
plt.xlabel('Cluster', fontsize=11)
plt.title('Cell Cycle Phase Distribution', fontsize=12)
plt.legend(loc='right',bbox_to_anchor=(1.3, 0.5), fontsize=10)
plt.xticks(rotation=0)

# Show the plot
plt.tight_layout()
plt.savefig('snapatac_out/cellcycle_distribution.png', dpi=300, bbox_inches='tight')
plt.show()




###for epithelial cell score:
#Epithelial all
gene_list=['GPA33','ACE','ALCAM','RNPEP','ANPEP','ANPEP','AMN','ICOSLG','CD276','MUC16','MUC1','CDH1','CD1A','CD1D','CD1D','CD46','CD74','CEACAM1','CEACAM3','CEACAM4','CEACAM5','CEACAM6','CEACAM7','COL1A1','COL1A2','C1QTNF5','CUBN','DDR1','DDR1','DDR2','DEFB4A','DEFB103A','DEFA1','DEFA5','HSPG2','EPCAM','FASLG','GKN1','SCGB3A1','HAS1','HAS2','HAS3','CADM4','ITGA4','ITGA4','ITGB1','ITGA4','ITGB7','F11R','JAM2','JAM3','L1CAM','LAMA1','LAMB1','LAMC1','MFGE8','MST1R','MUC1','MUC19','MUC4','PVRL1','PVRL2','PVRL3','PVRL4','NID1','OCLN','CD274','PLET1','PGF','PRSS8','SLURP2','TFRC','SCGB3A2',
    'CLDN1','CLDN3','CLDN4','CLDN6','CLDN10','CLDN12','CLDN17','CLDN19','CRNN','KRT8','KRT14','KRT18','KRT19','FOXJ1','FOXN1','KLF4','KLF5','KLF10','LRRC1','TCF7L1','TCF3']
sc.tl.score_genes(adata, gene_list, ctrl_size=len(gene_list), score_name='epithelial_all_score')


###for B cell activation score:
#https://maayanlab.cloud/Harmonizome/gene_set/B+cell+activation/GO+Biological+Process+Annotations+2023
gene_list=[
    'ADAM17','ADGRG3','AICDA','AKAP17A','AQP8','BANK1','BCL2','BLNK','BTK','CASP8','CD22','CD40','CD40LG','CD70','CD79A','CD79B','CD86','CEBPG','CHRNA4','CHRNB2','CLCF1','CR2','CTPS1','DCAF1','DNAJB9','EPHB2','FCRLA','FLT3','FNIP1','GPR183','GPS2','HDAC4','HDAC5','HDAC9','HHEX','HSPD1','IFNA1','IFNA10','IFNA14','IFNA16','IFNA17','IFNA2','IFNA21','IFNA4','IFNA5','IFNA6','IFNA7','IFNA8','IFNB1','IFNE','IFNK','IFNW1','IL10','IL11','IL4','ITGA4','ITGB1','JAK3','KIT','KLF6','LAT2','LAX1','LYL1','MALT1','MEF2C','MS4A1','MSH2','NFAM1','NHEJ1','PAX5','PHB','PHB2','PIK3CD','PIK3R1','PIK3R2','PIK3R3','PLCG2','PRKCB','PTPN2','PTPRC','PTPRJ','RABL3','RAG1','RAG2','RASGRP1','SLC25A5','SLC39A7','TCF3','TOP2B','TPD52','VCAM1','ZAP70','ZBTB7A'
]
sc.tl.score_genes(adata, gene_list, ctrl_size=len(gene_list), score_name='Bcell_activation_score')
