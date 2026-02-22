import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import scanpy as sc
sc.settings.verbosity = 3
sc.logging.print_header()
from matplotlib import rcParams
from scipy.stats import sem
from scipy.stats import ttest_ind


df=pd.read_csv('signac_output/EpiTrace_SignacPeak_Convergence_out.csv')
df.index=df.cell
EpiTraceAge_iterative=[]

for i in adata.obs.index:
    if i not in df.index:
        EpiTraceAge_iterative.append(-1)
    else:
        EpiTraceAge_iterative.append(df.loc[i]['EpiTraceAge_iterative'])
adata.obs['EpiTraceAge_iterative']=EpiTraceAge_iterative

adata2=adata[adata.obs['EpiTraceAge_iterative'] != -1, :]

rcParams['axes.grid'] = False
summary_df = adata2.obs.groupby('cluster_use')['EpiTraceAge_iterative'].agg(['mean', sem]).reset_index()

summary_df.columns = ['cluster_use', 'mean', 'sem']

# Assuming your DataFrame is named 'df' and has 'cluster_use' and 'EpiTraceAge_MitosisClock' columns
# Define the two groups
group1 = adata2.obs[adata2.obs['cluster_use'] == 'C1']['EpiTraceAge_iterative']
group2 = adata2.obs[adata2.obs['cluster_use'] == 'C5']['EpiTraceAge_iterative']

# Perform an independent t-test
t_stat, p_value = ttest_ind(group1, group2, equal_var=False)  # Use equal_var=False for Welch's t-test if variances are unequal

# Display the results
print(f"T-statistic: {t_stat}")
print(f"P-value: {p_value}")


summary_df = adata2.obs.groupby('cluster_use')['EpiTraceAge_iterative'].agg(['mean', sem]).reset_index()

# Rename columns for clarity
summary_df.columns = ['cluster_use', 'mean', 'sem']
plt.figure(figsize=(4.5,4))
ax = sns.barplot(
    data=summary_df,
    x='cluster_use',
    y='mean',
    palette=adata.uns['cluster_use_colors'],
    ci=None  # Disable seaborn's confidence intervals
)

# Add error bars for SEM
plt.errorbar(
    x=range(len(summary_df)),
    y=summary_df['mean'],
    yerr=summary_df['sem'],
    fmt='none',  # No line connecting markers
    ecolor='black',
    capsize=5  # Length of error bar caps
)

# Plot p-value bar between the two bars
x1, x2 = 1, 5  # Coordinates of the bars to compare (adjust based on position of 'Group1' and 'Group2')
y, h, col = summary_df['mean'].max() + summary_df['sem'].max()+0.015, 0.03, 'black'  # Set y position above bars and height

plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1, color=col)
plt.text((x1 + x2) / 2, y + h, "****", ha='center', va='bottom', color=col, fontsize=12)

# Add labels
plt.xlabel('Cluster')
plt.ylabel('EpiTrace Age', fontsize=16)
plt.ylim((0.3,0.72))
plt.title('')
plt.xticks(ticks=range(len(summary_df)), labels=summary_df['cluster_use'], rotation=0)  # Rotate x-ticks if needed
plt.tight_layout()
plt.show()