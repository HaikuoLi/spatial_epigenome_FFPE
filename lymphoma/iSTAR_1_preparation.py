import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import scanpy as sc
sc.settings.verbosity = 3
sc.logging.print_header()
from matplotlib import rcParams
import snapatac2 as snap

adata=sc.read('MALT_raw_8444.h5ad')

adata=adata[adata.obs['n_fragment']<15000, :]

snap.pp.add_tile_matrix(adata, bin_size=10000)

snap.pp.select_features(adata, n_features=50000, inplace = True)

adata = adata[:,adata.var['selected']==True]

actual_position=[]
for i in range(len(adata.obs)):
    actual_position.append(str(adata.obs['x_coord'][i])+'x'+str(adata.obs['y_coord'][i]))

df = pd.DataFrame(adata.X.toarray(), index=actual_position, columns=adata.var_names)
df.index.name='spot'

#this is one of the iSTAR input files
df.to_csv('cnts.tsv', sep='\t')


# Generate row names
row_names = [f"{i}x{j}" for i in range(1, 101) for j in range(1, 101)]

# Create a DataFrame with zeros
locs = pd.DataFrame(0, index=row_names, columns=['x-coord', 'y-coord'])

# tissue HE coord for the x axis
x_start = 4364
x_end = 12923
x_increment = (x_end - x_start) / 100
y_value = 2845

for i in range(1, 101):
    for j in range(1, 101):
        row_name = f"{i}x{j}"
        locs.at[row_name, 'x-coord'] = x_start + x_increment/2 + x_increment * (i - 1)
        locs.at[row_name, 'y-coord'] = y_value


# tissue HE coord for the y axis
y_start = 2845
y_end = 11153
y_increment = (y_end - y_start) / 100

for i in range(1, 101):
    for j in range(1, 101):
        row_name = f"{i}x{j}"
        locs.at[row_name, 'y-coord'] = y_start + y_increment/2 + y_increment * (j - 1)

in_tissue=[]
for i in locs.index:
    if i in df.index:
        in_tissue.append(1)
    else:
        in_tissue.append(0)
locs['in_tissue']=in_tissue
locs=locs[locs['in_tissue'] == 1]

locs=locs.reindex(df.index)
del locs['in_tissue']

#this is one of the iSTAR input files
locs[['x-coord','y-coord']].to_csv('locs-raw.tsv', sep='\t',header=['x', 'y'],index_label='spot')

#pixel size is 0.4744
#radius is 42.1675