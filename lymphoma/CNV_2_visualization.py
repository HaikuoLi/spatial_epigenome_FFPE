import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import scanpy as sc
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
from matplotlib import rcParams

adata=sc.read('MALT_clustering_8444x250k.h5ad')

df=pd.read_csv('epiAneufinder_results/results_table.tsv', sep=' ')

cnv_gain=[]
cnv_loss=[]
cnv_normal=[]
cnv_varied=[]
for i in adata.obs.index:
    id='cell-'+str(i)
    if id in df.columns:
        cnv_gain.append(list(df[id]).count(2))
        cnv_normal.append(list(df[id]).count(1))
        cnv_loss.append(list(df[id]).count(0))
        cnv_varied.append(list(df[id]).count(0)+list(df[id]).count(2))
    else:
        cnv_loss.append(0)
        cnv_gain.append(0)
        cnv_normal.append(0)
        cnv_varied.append(0)

adata.obs['cnv_gain']=cnv_gain
adata.obs['cnv_loss']=cnv_loss
adata.obs['cnv_normal']=cnv_normal
adata.obs['cnv_varied']=cnv_varied

adata.obs["cnv_loss"] = adata.obs["cnv_loss"].astype("uint64")
adata.obs["cnv_gain"] = adata.obs["cnv_gain"].astype("uint64")
adata.obs["cnv_normal"] = adata.obs["cnv_normal"].astype("uint64")
adata.obs["cnv_varied"] = adata.obs["cnv_varied"].astype("uint64")

clone_df=pd.read_csv("8444_bin500k_8clone_annot_dt.csv")
#data frame written from R

clone_df.index=clone_df['cell']
cnv_clone=[]
for i in adata.obs.index:
    id='cell-'+str(i)
    if id in clone_df['cell']:
        cnv_clone.append(clone_df.loc[id]['annot'])
    else:
        cnv_clone.append('NA')
adata.obs['cnv_clone']=cnv_clone


print(sns.color_palette("Set2", 8).as_hex())

sc.settings.set_figure_params(facecolor=(0, 0, 0, 0), dpi=300, fontsize=12)
sc.pl.spatial(adata, color=["cnv_clone"],img_key=None, spot_size=1, title='Karyotypes',
              palette=[ '#d6b09c',
'#e41a1c',
 '#f781bf',
 '#984ea3',
 '#ff7f00',
 '#ffff33',
 '#377eb8',
 '#4daf4a',
 'gainsboro'])

sc.settings.set_figure_params(dpi=300, facecolor=(0, 0, 0, 0),fontsize=12)

palette_cmap=['gainsboro']*9
original_cmap=['#d6b09c',
'#e41a1c',
 '#f781bf',
 '#984ea3',
 '#ff7f00',
 '#ffff33',
 '#377eb8',
 '#4daf4a',
 'gainsboro']

for i in range(8):
    palette_cmap[i]=original_cmap[i]
    sc.pl.spatial(adata, color=["cnv_clone"],img_key=None, spot_size=0.92, title='Karyotype Clone'+str(i+1),
             palette=palette_cmap, legend_loc=None)
    palette_cmap=['gainsboro']*9

adata2=adata[adata.obs["cnv_normal"] != 0]

# Count occurrences of each combination of class1 and class2
count_df = adata2.obs[['cluster_use','cnv_clone']].groupby(['cluster_use', 'cnv_clone']).size().unstack(fill_value=0)
# Plot with reversed stacking order by reordering columns only for the plot
ax = count_df[count_df.columns[::-1]].plot(kind='bar', stacked=True, color=['#d6b09c',
'#e41a1c',
 '#f781bf',
 '#984ea3',
 '#ff7f00',
 '#ffff33',
 '#377eb8',
 '#4daf4a'][::-1], figsize=(8,5))
plt.xlabel("Spatial clusters",fontsize=16)
plt.xticks(rotation=0,fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel("Karyotype clone count",fontsize=16)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], title="Karyotypes", bbox_to_anchor=(1.03, 0.82), loc='upper left')
ax.grid(False)
plt.tight_layout()
plt.savefig("clone_per_cluster.png", dpi=300)
plt.show()


count_df = adata2.obs[['cluster_use','cnv_clone']].groupby(['cnv_clone','cluster_use']).size().unstack(fill_value=0)

count_df_normalized = count_df.div(count_df.sum(axis=1), axis=0)
ax = count_df_normalized[count_df_normalized.columns[::-1]].plot(kind='bar', stacked=True, color=adata.uns['cluster_use_colors'][::-1], figsize=(8,5))
plt.xlabel("Karyotypes",fontsize=16)
plt.xticks(rotation=45,fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel("Spatial cluster proportion",fontsize=16)
plt.ylim(-0.01,1.01)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], title="Clusters", bbox_to_anchor=(1.03, 0.82), loc='upper left')
ax.grid(False)
plt.tight_layout()
plt.savefig("cluster_per_clone_ratio.png", dpi=300)
plt.show()

