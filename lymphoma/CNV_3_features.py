import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import scanpy as sc
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
from matplotlib import rcParams


df=pd.read_csv('epiAneufinder_results/results_table.tsv', sep=' ')

index=[]
for i in range(len(df)):
    index.append(list(df['seq'])[i]+":"+str(list(df['start'])[i])+"-"+str(list(df['end'])[i]))
df2=df.copy()
df2.index=index
del df2['seq']; del df2['start']; del df2['end']
columns=[]
for i in df2.columns:
    columns.append(i.split('-')[1])
df2.columns=columns

output_file = "CNV_homer_input.txt"
i=1
with open(output_file, "w") as outfile:
    for line in df2.index:
        line = line.strip()
        # Parse chromosome, start, and end positions
        chrom, pos = line.split(":")
        start, end = pos.split("-")

        # Format the output line
        peak_id = f"peak_{i+1}"       # Unique Peak ID
        strand = "0"                  # Set strand to "+" (0) by default

        # Write to output file
        outfile.write(f"{peak_id}\t{chrom}\t{start}\t{end}\t{strand}\n")
        i+=1

print(f"Converted file saved as '{output_file}'")


#Using Homer in Linux
#annotatePeaks.pl CNV_homer_input.txt hg38 > CNV_homer_output.txt
!sort -k1.6n CNV_homer_output.txt -o CNV_homer_output_sort.txt

homer_df=pd.read_csv('CNV_homer_output_sort.txt', sep='\t')

position=[]
for i in range(len(homer_df)):
    position.append(homer_df['Chr'][i]+':'+str(homer_df['Start'][i])+'-'+str(homer_df['End'][i]))
homer_df['position']=position

adata=sc.read('MALT_clustering_8444x250k.h5ad')

clone_df=pd.read_csv("20241110_8444_bin500k_8clone_annot_dt.csv")
clone_df.index=clone_df['cell']
cnv_clone=[]
for i in adata.obs.index:
    id='cell-'+str(i)
    if id in clone_df['cell']:
        cnv_clone.append(clone_df.loc[id]['annot'])
    else:
        cnv_clone.append('NA')
adata.obs['cnv_clone']=cnv_clone

adata2=adata[adata.obs["cnv_clone"] != 'NA']
df_cluster=adata2.obs[['cluster_use','cnv_clone']]

new_index=[]
for i in range(len(homer_df)):
    new_index.append(str(homer_df['Gene Name'][i])+' ('+homer_df['position'][i]+')')

df2.index=new_index

df2_gain=df2.applymap(lambda x: 1 if x == 2 else 0)

df_merged = df2_gain.T.merge(df_cluster[['cluster_use']], left_index=True, right_index=True)
df_summed = df_merged.groupby('cluster_use').sum()
df_summed = df_summed.T

# Sum up columns 'C1' and 'C5' and divide by 2
df_summed['C1_C5_avg'] = df_summed[['C1', 'C5']].sum(axis=1) / 2

# Sum up all other columns and divide by 7
other_columns = df_summed.columns.difference(['C1', 'C5'])
df_summed['other_avg'] = df_summed[other_columns].sum(axis=1) / 7

# Optionally, keep only the new averaged columns
df_result = df_summed[['C1_C5_avg', 'other_avg']]

# Calculate the ratio of C1_C5_avg to other_avg
df_result['ratio'] = df_result['C1_C5_avg'] / df_result['other_avg']

# Sort the DataFrame by the 'ratio' column in descending order
df_result_sorted = df_result.sort_values(by='ratio', ascending=False)

##Find top 10 gain CNVs
for bin_id in df_result_sorted.index[:10]:
    df_tmp=df2.loc[bin_id]
    bin_class=[]
    for i in adata.obs.index:
        if i in df2.columns:
            if df_tmp[i]==2:
                bin_class.append('Gain')
            elif df_tmp[i]==1:
                bin_class.append('NA')
            elif df_tmp[i]==0:
                bin_class.append('Loss')
        else:
            bin_class.append('NA')
    
    category_order=['Gain','Loss','NA']
    adata.obs['bin_class']=bin_class
    adata.obs['bin_class'] = pd.Categorical(
        adata.obs['bin_class'], 
        categories=category_order, 
        ordered=True)
    
    out_id='_gain_'+str(bin_id)+'.png'
    sc.pl.spatial(adata, color=["bin_class"],img_key=None, spot_size=1, title=bin_id,
                 palette=['#cd0000','#27408b','#e5e5e5'], save=out_id, show=False)



df2_loss=df2.applymap(lambda x: 1 if x == 0 else 0)
df_merged = df2_loss.T.merge(df_cluster[['cluster_use']], left_index=True, right_index=True)
df_summed = df_merged.groupby('cluster_use').sum()
df_summed = df_summed.T

# Sum up columns 'C1' and 'C5' and divide by 2
df_summed['C1_C5_avg'] = df_summed[['C1', 'C5']].sum(axis=1) / 2
# Sum up all other columns and divide by 7
other_columns = df_summed.columns.difference(['C1', 'C5'])
df_summed['other_avg'] = df_summed[other_columns].sum(axis=1) / 7
# Optionally, keep only the new averaged columns
df_result = df_summed[['C1_C5_avg', 'other_avg']]
# Calculate the ratio of C1_C5_avg to other_avg
df_result['ratio'] = df_result['C1_C5_avg'] / df_result['other_avg']
# Sort the DataFrame by the 'ratio' column in descending order
df_result_sorted = df_result.sort_values(by='ratio', ascending=False)

##find top 10 loss CNVs
for bin_id in df_result_sorted.index[:10]:
    out_id='_loss_'+str(bin_id)+'.png'
    df_tmp=df2.loc[bin_id]
    bin_class=[]
    for i in adata.obs.index:
        if i in df2.columns:
            if df_tmp[i]==2:
                bin_class.append('Gain')
            elif df_tmp[i]==1:
                bin_class.append('NA')
            elif df_tmp[i]==0:
                bin_class.append('Loss')
        else:
            bin_class.append('NA')
    
    category_order=['Gain','Loss','NA']
    adata.obs['bin_class']=bin_class
    adata.obs['bin_class'] = pd.Categorical(
        adata.obs['bin_class'], 
        categories=category_order, 
        ordered=True)
    
    sc.pl.spatial(adata, color=["bin_class"],img_key=None, spot_size=1, title=bin_id,
                 palette=['#cd0000','#27408b','#e5e5e5'], save=out_id, show=False)