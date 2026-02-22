##followed by clustering. Run in Jupyter

%%time
gene_matrix = snap.pp.make_gene_matrix(data2, snap.genome.mm10)

sc.pp.filter_genes(gene_matrix, min_cells = 20)

sc.pp.normalize_total(gene_matrix)
sc.pp.log1p(gene_matrix)
sc.pp.highly_variable_genes(gene_matrix, min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.settings.set_figure_params(dpi=100, facecolor='white',fontsize=12)
sc.pl.highly_variable_genes(gene_matrix)

sc.external.pp.magic(gene_matrix, solver="approximate")

# Copy over UMAP embedding
gene_matrix.obsm["X_umap"] = data2.obsm["X_umap"]

sc.tl.rank_genes_groups(gene_matrix, groupby='cluster_use',
                        method='wilcoxon',key_added='rank_genes_groups_cluster_use')

writer = pd.ExcelWriter('snapatac_out/rank_genes_groups_cluster_use.xlsx', engine='xlsxwriter')

# Write each dataframe to a different worksheet.
pd.DataFrame(gene_matrix.uns['rank_genes_groups_cluster_use']['names']).to_excel(writer, sheet_name='names')
pd.DataFrame(gene_matrix.uns['rank_genes_groups_cluster_use']['scores']).to_excel(writer, sheet_name='scores')
pd.DataFrame(gene_matrix.uns['rank_genes_groups_cluster_use']['pvals']).to_excel(writer, sheet_name='pvals')
pd.DataFrame(gene_matrix.uns['rank_genes_groups_cluster_use']['pvals_adj']).to_excel(writer, sheet_name='pvals_adj')
pd.DataFrame(gene_matrix.uns['rank_genes_groups_cluster_use']['logfoldchanges']).to_excel(writer, sheet_name='logfoldchanges')
writer.close()

