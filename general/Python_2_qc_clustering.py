##followed by processing. Run in Jupyter

fig = snap.pl.frag_size_distr(adata, show=False)
fig.update_layout(
    xaxis_title="Fragment size (bp)", yaxis_title="Count",
    width=400,
    height=450
)
fig.update_layout(
    plot_bgcolor='rgba(0,0,0,0)',  # Remove background color (transparent)
    paper_bgcolor='rgba(0,0,0,0)',  # Remove the background of the entire figure (transparent)
    xaxis=dict(
        showline=True,  # Show frame on x-axis
        showgrid=False,
        linewidth=1,  # Thickness of the frame line
        linecolor='black',  # Color of the frame line
        ticks='outside',  # Place ticks outside the axis line
        tickwidth=2,  # Thickness of the ticks
        tickcolor='black',  # Color of the ticks
        ticklen=10,  # Length of the ticks
        tickmode='auto',  # Can be 'auto' or 'array'
        # Optional: Define specific tick positions
        # tickvals=[0, 1, 2, 3, 4],
        # Optional: Define custom tick labels
        # ticktext=['A', 'B', 'C', 'D', 'E'],
        mirror=True,
        title_font=dict(
            family="Arial",  # Font family
            size=18,  # Font size
            color="black"  # Font color
        ),tickfont=dict(
            family="Arial",  # Font family
            size=16,  # Font size
            color="black"  # Font color
        )
    ),
    yaxis=dict(
        showline=True,  # Show frame on y-axis
        showgrid=False,
        linewidth=1,  # Thickness of the frame line
        linecolor='black',  # Color of the frame line
        ticks='outside',  # Place ticks outside the axis line
        tickwidth=2,  # Thickness of the ticks
        tickcolor='black',  # Color of the ticks
        ticklen=10,  # Length of the ticks
        tickmode='auto',  # Can be 'auto' or 'array'
        mirror=True,title_font=dict(
            family="Arial",  # Font family
            size=18,  # Font size
            color="black"  # Font color
        ),tickfont=dict(
            family="Arial",  # Font family
            size=16,  # Font size
            color="black"  # Font color
        )
    ),
)
fig.show()
fig.write_image('fragment_dis.png', scale=4)

%%time
snap.metrics.tsse(data, snap.genome.mm10)


fig = snap.pl.tsse(data, interactive=False, show=False)
fig.update_layout(
    width=600,
    height=500
)
fig.update_layout(
    plot_bgcolor='rgba(0,0,0,0)',  # Remove background color (transparent)
    paper_bgcolor='rgba(0,0,0,0)',  # Remove the background of the entire figure (transparent)
    xaxis=dict(
        showline=True,  # Show frame on x-axis
        showgrid=False,
        linewidth=1,  # Thickness of the frame line
        linecolor='black',  # Color of the frame line
        ticks='outside',  # Place ticks outside the axis line
        tickwidth=2,  # Thickness of the ticks
        tickcolor='black',  # Color of the ticks
        ticklen=10,  # Length of the ticks
        tickmode='auto',  # Can be 'auto' or 'array'
        # Optional: Define specific tick positions
        # tickvals=[0, 1, 2, 3, 4],
        # Optional: Define custom tick labels
        # ticktext=['A', 'B', 'C', 'D', 'E'],
        mirror=True,
        title_font=dict(
            family="Arial",  # Font family
            size=18,  # Font size
            color="black"  # Font color
        ),tickfont=dict(
            family="Arial",  # Font family
            size=16,  # Font size
            color="black"  # Font color
        )
    ),
    yaxis=dict(
        showline=True,  # Show frame on y-axis
        showgrid=False,
        linewidth=1,  # Thickness of the frame line
        linecolor='black',  # Color of the frame line
        ticks='outside',  # Place ticks outside the axis line
        tickwidth=2,  # Thickness of the ticks
        tickcolor='black',  # Color of the ticks
        ticklen=10,  # Length of the ticks
        tickmode='auto',  # Can be 'auto' or 'array'
        mirror=True,title_font=dict(
            family="Arial",  # Font family
            size=18,  # Font size
            color="black"  # Font color
        ),tickfont=dict(
            family="Arial",  # Font family
            size=16,  # Font size
            color="black"  # Font color
        )
    ),
)
fig.show()
fig.write_image('tsse.png', scale=4)

%%time
snap.pp.add_tile_matrix(data)

snap.pp.select_features(data, n_features=500000, inplace = True)

data2 = data[:,data.var['selected']==True]

snap.tl.spectral(data2, n_comps=20)

snap.tl.umap(data2, min_dist=0.02)

sc.settings.set_figure_params(dpi=300, facecolor='white',fontsize=12)
sc.pl.umap(data2, color=["n_fragment",'frac_dup','tsse'], use_raw=False, size=25, wspace=0.25)

%%time
snap.pp.knn(data2)
data2

snap.tl.leiden(data2, resolution=1.4, key_added='leiden')
sc.pl.umap(data2,color='leiden', size=25)


## Generate a spatial obsm

data2.uns['spatial']={}
data2.uns['spatial']['whole']={

                           }
data2.obsm['spatial']=pd.DataFrame.to_numpy(data2.obs[['x_coord','y_coord']].astype('int64'))

sc.pl.spatial(data2, color=["leiden"],img_key=None, spot_size=1)

sc.pl.spatial(data2, color=["n_fragment"],img_key=None, spot_size=1, vmax=10000
             )
sc.pl.spatial(data2, color=["tsse"],img_key=None, spot_size=1
             )



palette_cmap=['gainsboro']*10
original_cmap=['steelblue', 'orange', '#279e68', '#d62728', 'darkviolet', '#8c564b', 'navy', '#b5bd61', '#17becf', '#aec7e8']

for i in range(10):
    palette_cmap[i]=original_cmap[i]
    sc.pl.spatial(data2, color=["leiden"],img_key=None, spot_size=0.92, title='Cluster #'+str(i),
             palette=palette_cmap, legend_loc=None)
    palette_cmap=['gainsboro']*10