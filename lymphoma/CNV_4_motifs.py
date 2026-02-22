import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import scanpy as sc
sc.settings.verbosity = 3
sc.logging.print_header()
from matplotlib import rcParams
from scipy.stats import ttest_ind

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

adata=sc.read('MALT_clustering_8444x250k.h5ad')

df2_gain=df2.applymap(lambda x: 1 if x == 2 else 0)

clone_df=pd.read_csv("8444_bin500k_8clone_annot_dt.csv")
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

##### Check OTX2 (Score 500) (C1C5/others ratio > 3) ######
#using UCSC-bigbedtobed

BIGBEDTOBED = "/gpfs/gibbs/pi/fan/hl737/Downloads/bigBedToBed"
BIGBED_URL = "http://hgdownload.soe.ucsc.edu/gbdb/hg38/jaspar/JASPAR2020.bb"

df = pd.DataFrame(index=df_result_sorted[df_result_sorted['ratio']>3].index)

# Parse index into a DataFrame for sorting and merging
intervals = []
for idx in df.index:
    chrom, pos = idx.split(":")
    start, end = map(int, pos.split("-"))
    intervals.append((chrom, start, end))

# Convert to DataFrame and sort
intervals_df = pd.DataFrame(intervals, columns=["chrom", "start", "end"]).sort_values(["chrom", "start"])

# Merge consecutive intervals
merged_intervals = []
current_chrom, current_start, current_end = intervals_df.iloc[0]

for _, row in intervals_df.iloc[1:].iterrows():
    chrom, start, end = row
    if chrom == current_chrom and start == current_end + 1:  # If consecutive
        current_end = end  # Extend current interval
    else:
        merged_intervals.append(f"{current_chrom}:{current_start}-{current_end}")
        current_chrom, current_start, current_end = chrom, start, end

# Append the last interval
merged_intervals.append(f"{current_chrom}:{current_start}-{current_end}")
print(len(merged_intervals))

counts_dict = {}  # Initialize dictionary

for interval in merged_intervals:
    chrom, pos = interval.split(":")
    start, end = pos.split("-")

    command = f"""
    {BIGBEDTOBED} {BIGBED_URL} -chrom={chrom} -start={start} -end={end} stdout | \
    awk '$4 == "OTX2" {{ if ($5 > 500) print }}' | wc -l
    """
    
    # Run the command in Jupyter
    #print(f"Running command for {interval} ...")
    output = !{command}  # Run the command and capture output
    count = int(output[0])  # Convert output to integer
    counts_dict[interval] = count  # Store in dictionary

print(np.sum(list(counts_dict.values())))
print(np.sum(list(counts_dict.values()))/len(df))
df_OTX2_gain = pd.DataFrame(list(counts_dict.items()), columns=['Interval', 'Count'])


df = pd.DataFrame(index=list(df_result_sorted.index)[-227:])

intervals = []
for idx in df.index:
    chrom, pos = idx.split(":")
    start, end = map(int, pos.split("-"))
    intervals.append((chrom, start, end))

# Convert to DataFrame and sort
intervals_df = pd.DataFrame(intervals, columns=["chrom", "start", "end"]).sort_values(["chrom", "start"])

# Merge consecutive intervals
merged_intervals = []
current_chrom, current_start, current_end = intervals_df.iloc[0]

for _, row in intervals_df.iloc[1:].iterrows():
    chrom, start, end = row
    if chrom == current_chrom and start == current_end + 1:  # If consecutive
        current_end = end  # Extend current interval
    else:
        merged_intervals.append(f"{current_chrom}:{current_start}-{current_end}")
        current_chrom, current_start, current_end = chrom, start, end

# Append the last interval
merged_intervals.append(f"{current_chrom}:{current_start}-{current_end}")
print(len(merged_intervals))

counts_dict = {}  # Initialize dictionary

for interval in merged_intervals:
    chrom, pos = interval.split(":")
    start, end = pos.split("-")

    command = f"""
    {BIGBEDTOBED} {BIGBED_URL} -chrom={chrom} -start={start} -end={end} stdout | \
    awk '$4 == "OTX2" {{ if ($5 > 500) print }}' | wc -l
    """
    
    # Run the command in Jupyter
    #print(f"Running command for {interval} ...")
    output = !{command}  # Run the command and capture output
    count = int(output[0])  # Convert output to integer
    counts_dict[interval] = count  # Store in dictionary

print(np.sum(list(counts_dict.values())))
print(np.sum(list(counts_dict.values()))/len(df))
df_OTX2_nogain = pd.DataFrame(list(counts_dict.items()), columns=['Interval', 'Count'])


t_stat, p_value = ttest_ind(df_OTX2_gain["Count"], df_OTX2_nogain["Count"], equal_var=False)
print(f"T-test results: t-statistic = {t_stat:.4f}, p-value = {p_value:.4f}")

# Calculate means and SEM
mean_df1 = np.mean(df_OTX2_gain["Count"])
mean_df2 = np.mean(df_OTX2_nogain["Count"])
sem_df1 = np.std(df_OTX2_gain["Count"], ddof=1) / np.sqrt(len(df_OTX2_gain["Count"]))
sem_df2 = np.std(df_OTX2_nogain["Count"], ddof=1) / np.sqrt(len(df_OTX2_nogain["Count"]))

# Create a plot
fig, ax = plt.subplots(figsize=(3, 4))
x = np.arange(2)

ax.bar(x[0], mean_df1, yerr=sem_df1, capsize=10, color="skyblue")
ax.bar(x[1], mean_df2, yerr=sem_df2, capsize=10, color="lightgreen")

jitter_strength = 0.1
ax.scatter(np.random.normal(x[0], jitter_strength, len(df_OTX2_gain)), df_OTX2_gain["Count"], 
           color="blue", alpha=0.2)
ax.scatter(np.random.normal(x[1], jitter_strength, len(df_OTX2_nogain)), df_OTX2_nogain["Count"], 
           color="green", alpha=0.2)

ax.set_xticks(x)
ax.set_xticklabels(["w/ gained copies", "w/o gained copies"], rotation=30)
ax.set_ylabel("Motif number in regions with tumor CNV")
ax.set_title("Frequency of OTX2 motifs")
plt.ylim((0,800))
plt.tight_layout()
#plt.show()
plt.savefig("./CNV_motif/tumor_OTX2.png", dpi=300)