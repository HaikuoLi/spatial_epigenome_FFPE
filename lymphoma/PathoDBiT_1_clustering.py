import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import scanpy as sc
sc.settings.verbosity = 3
sc.logging.print_header()
from matplotlib import rcParams
sc.settings.set_figure_params(dpi=100, facecolor="white")

adata=sc.read('../dataset/LM19RNA/LM19_stdata_names.tsv',delimiter='\t')

x_coord=[]
y_coord=[]
for i in adata.obs.index:
    i=i.split('x')
    x_coord.append(i[0])
    y_coord.append(i[1])
adata.obs['x_coord']=x_coord
adata.obs['y_coord']=y_coord

adata.uns['spatial']={}
adata.uns['spatial']['whole']={

                           }
adata.obsm['spatial']=pd.DataFrame.to_numpy(adata.obs[['x_coord','y_coord']].astype('int64'))
adata.obsm['spatial'].shape

sc.pp.filter_cells(adata, min_counts=1)
print(adata)
sc.pp.filter_cells(adata, min_genes=1)
print(adata)

sc.pp.filter_genes(adata, min_cells=20)

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

sc.pp.normalize_total(adata, target_sum=None)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=15000)
sc.pl.highly_variable_genes(adata)

adata.raw = adata
adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ["n_counts",'n_genes','total_counts'], n_jobs=4)

sc.pp.scale(adata, max_value=None)
sc.tl.pca(adata, svd_solver="arpack")


sc.pl.pca(adata, color=['n_counts','pct_counts_mt'])
sc.pl.pca_variance_ratio(adata, log=True)

sc.pp.neighbors(adata, n_neighbors=20, n_pcs=20)

sc.tl.umap(adata, min_dist=0.1)

sc.tl.leiden(adata, resolution=0.85)
sc.pl.umap(adata, color=["leiden"])
sc.pl.spatial(adata, color=['leiden'],img_key=None, spot_size=1)

marker_genes=['PAPOLG','REL','PUS10','PEX13','KIAA1841','USP34','XPO1']
sc.pl.matrixplot(adata, marker_genes ,groupby = 'leiden',standard_scale='var',
                 cmap='Oranges', swap_axes=True)