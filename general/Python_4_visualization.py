# followed by gene activity computation

#dot plot
from matplotlib.colors import LinearSegmentedColormap

marker_genes=['Satb2','Pdzrn3','Mobp','Sox10','Ankub1','Tshz2','Plekhg1','Stat3','Snx21','Pde7b','Meis2',
             'Ttc23','Atr','Ccdc33','Slc8a2','Impa2','Art3','Szt2','Gdi1','Apc','Nwd2'
 ]
cmap=LinearSegmentedColormap.from_list('name',["w","lightblue","orangered"], N=256)

sc.pl.DotPlot(gene_matrix, marker_genes, groupby='cluster_use',
             standard_scale='var',expression_cutoff=0.85)\
.style(cmap=cmap,largest_dot=350,dot_edge_lw=0,size_exponent=1.4)\
.legend(width=1.65,colorbar_title='Mean activity',size_title='Fraction of spots\n').swap_axes().show()

#feature plot
sc.pl.spatial(gene_matrix, color=["Sox10"],img_key=None, spot_size=1, use_raw=False, vmax=0.16, vmin=0.092, cmap='coolwarm')

##to visualize chromvar output with scanpy
adata=sc.read('signac_output/chromvar.h5ad')
motif_gene=pd.read_csv('signac_output/chromvar_motif_to_gene.csv', index_col=0)

adata.var['motif_gene']=list(motif_gene['gene'])

motif_and_gene=[]
for i in range(len(adata.var)):
    motif_and_gene.append(adata.var['features'][i]+'_'+adata.var['motif_gene'][i])
adata.var['motif_and_gene']=motif_and_gene
adata.var.index=adata.var['motif_and_gene']

markers=['MA0476.1_FOS','MA0490.2_JUNB',
    
'MA0442.2_SOX10',
'MA0867.2_SOX4',
    
'MA0677.1_Nr2f6',
'MA0512.2_Rxra',

'MA0466.2_CEBPB',
'MA0746.2_SP3',
    
'MA0162.4_EGR1',
'MA0052.4_MEF2A',
         
'MA1575.1_THRB(var.2)',
'MA1637.1_EBF3',

'MA1109.1_NEUROD1',
'MA0827.1_OLIG3',
         
'MA1546.1_PAX3(var.2)','MA0711.1_OTX1',

'MA1514.1_KLF17',
'MA0502.2_NFYB',

'MA0768.1_LEF1',
'MA1474.1_CREB3L4',

       ]


sc.pl.matrixplot(adata, markers, groupby = 'cluster_use', standard_scale='var',
                 cmap='coolwarm',colorbar_title='chromVAR Activity')

sc.settings.set_figure_params(dpi=300, facecolor='white',fontsize=12)
sc.pl.spatial(adata, color=['MA0442.2'], img_key=None, spot_size=1, title='SOX10 (MA0442.2)', vmin=-0.5, vmax=2.5)


#Split the bulk fragment bed into bed for each cluster. Good for Browser visualization

celltype_csv=pd.read_csv('snapatac_out/meta.csv',index_col=0)

celltype_dict={}
for i in set(celltype_csv['cluster_use']):
    celltype_dict[i]=celltype_csv[celltype_csv['cluster_use']==i].index.to_list()

!gunzip -c fragments.tsv.gz > fragments.tsv

fragment_file=pd.read_csv('fragments.tsv',sep='\t', header=None, comment='#')

for i,j in celltype_dict.items():
    subset_df=fragment_file[fragment_file[3].isin(j)]
    filename='snapatac_out/cluster_use_fragment/'+str(i)+'.bed'
    subset_df.to_csv(filename, sep='\t', header=False, index=False)

#module load tabix
#for file in *.bed; do bgzip "$file"; done
#for file in *.bed.gz; do tabix -p bed "$file"; done


##ODC.MOL and ODC.MFOL fragments can be downloaded from: http://catlas.org/mousebrain/#!/downloads

#zcat NonN.OGC.MFOL.bedpe.gz | awk '{print $1"\t"$2"\t"$3}' | sort -k1,1 -k2,2n -S 50G | bgzip > NonN.OGC.MFOL.sort.bed.gz
#tabix -p bed NonN.OGC.MFOL.sort.bed.gz

#zcat NonN.OGC.MOL.bedpe.gz | awk '{print $1"\t"$2"\t"$3}' | sort -k1,1 -k2,2n -S 50G | bgzip > NonN.OGC.MOL.sort.bed.gz
#tabix -p bed NonN.OGC.MOL.sort.bed.gz

## An example of correlation plot (fragment level)
adata=sc.read('MALT_raw_8444.h5ad')
snap.pp.add_tile_matrix(adata, bin_size=10000)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.settings.set_figure_params(dpi=300, facecolor=(0, 0, 0, 0),fontsize=12)
rcParams['axes.grid'] = False
sc.pl.correlation_matrix(adata, groupby='cluster_use',dendrogram=True, show_correlation_numbers=True)