import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import scanpy as sc
sc.settings.verbosity = 3
sc.logging.print_header()
from matplotlib import rcParams
#use scipy 1.13.1

adata=sc.read('MALT_clustering_8444x250k.h5ad')

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
del adata2.obsm['X_umap']

x_coord=[]
y_coord=[]
for i in adata2.obs['coord']:
    i=i.split('x')
    x_coord.append(int(i[0]))
    y_coord.append(int(i[1]))
adata2.obs['x_coord']=x_coord
adata2.obs['y_coord']=y_coord

adata2.uns['spatial']={}
adata2.uns['spatial']['whole']={

                           }
adata2.obsm['spatial']=pd.DataFrame.to_numpy(adata2.obs[['x_coord','y_coord']].astype('int64'))

adata2.obsm['X_umap']=adata2.obsm['spatial']


import sys
import matplotlib.pyplot as plt
import scvelo as scv
import scanpy as sc
import cellrank as cr
from typing import Any
from copy import copy
from anndata import AnnData
import numpy as np
import scipy.sparse as sp
import petsc4py

class AgeKernel(cr.kernels.Kernel):
    def __init__(self, adata: AnnData, obs_key: str = "EpiTraceAge_iterative", **kwargs: Any):
        super().__init__(adata=adata, obs_key=obs_key, **kwargs)

    def _read_from_adata(self, obs_key: str, **kwargs: Any) -> None:
        super()._read_from_adata(**kwargs)
        print(f"Reading `adata.obs[{obs_key!r}]`")
        self.pseudotime = 1 - self.adata.obs[obs_key].values

    def compute_transition_matrix(self, some_parameter: float = 0.5) -> "AgeKernel":
        if self._reuse_cache({"some_parameter": some_parameter}):
            print("Using cached values for parameters:", self.params)
            return self
        # Manually create a transition matrix as an example
        transition_matrix = sp.diags((some_parameter,) * self.adata.shape[0], dtype=np.float64)
        # Store the transition matrix as an attribute for later use
        self.transition_matrix = transition_matrix
        return self

    def copy(self) -> "AgeKernel":
        return copy.copy(self)
    
    def backward(self) -> None:
        # Implement the required abstract method. This example is minimal.
        print("Executing backward method (this is a placeholder).")


adata2.layers['spliced'] = adata2.X

scv.pp.filter_and_normalize(adata2)

scv.pp.moments(adata2)
#now with obsp['connectivities']

sc.pp.pca(adata2)
sc.pp.neighbors(adata2, random_state=0)
adata2.obsm['X_umap'] = adata2.obsm['X_umap'].astype(np.float64)
ak = AgeKernel(adata2).compute_transition_matrix()

from matplotlib.colors import LinearSegmentedColormap
cmap=LinearSegmentedColormap.from_list('name',["red","lightgrey","green"], N=256)
sc.settings.set_figure_params(facecolor=(0, 0, 0, 0), dpi=300, fontsize=12)
scv.pl.velocity_embedding_stream(adata2, basis='umap', color='EpiTraceAge_iterative', vkey="T_bwd",
                                 legend_loc='right', color_map=cmap, colorbar=False,
                                 size=50, alpha=0.6, linewidth=0.5, title='EpiTrace Age')
scv.pl.velocity_embedding_stream(adata2, basis='umap', color='cluster_use', vkey="T_bwd",
                                 size=50, alpha=0.6, linewidth=0.5)





#below is an example how to visualize mitotic age vs cholesterol genes

adata2=adata[adata.obs['EpiTraceAge_iterative'] != -1, :]
# Extract EpiTraceAge_iterative
epi_trace_age = adata2.obs['EpiTraceAge_iterative']

# Calculate tertile thresholds (33rd and 66th percentiles)
low_threshold, high_threshold = np.percentile(epi_trace_age, [33.33, 66.67])
print(f"Low threshold (33rd percentile): {low_threshold:.3f}")
print(f"High threshold (66th percentile): {high_threshold:.3f}")

# Create group labels: low (< 33rd), medium (33rd to 66th), high (> 66th)
group_labels = np.select(
    [epi_trace_age < low_threshold, 
     (epi_trace_age >= low_threshold) & (epi_trace_age <= high_threshold), 
     epi_trace_age > high_threshold],
    ['Low mitotic age', 'Medium mitotic age', 'High mitotic age']
)

# Add the group labels to adata2.obs
adata2.obs['EpiTraceAge_group'] = group_labels

adata2.obs['EpiTraceAge_group'] = adata2.obs['EpiTraceAge_group'].cat.reorder_categories(['High mitotic age','Medium mitotic age','Low mitotic age'])

sc.settings.set_figure_params(facecolor=(0, 0, 0, 0), dpi=300, fontsize=12)
markers = ['SREBF1', 'SREBF2', 'LDLR', 'HMGCS1', 'HMGCR', 'INSIG1','NPC1','NPC2']
sc.pl.matrixplot(adata2, markers, groupby='EpiTraceAge_group', use_raw=False, standard_scale='var',
                colorbar_title='Mean activity\nin group', cmap='Reds')

plt.rcParams['axes.grid'] = False
sc.pl.violin(adata2,  markers, groupby='EpiTraceAge_group', use_raw=False,
            stripplot=False,  inner="box")