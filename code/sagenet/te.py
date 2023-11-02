import squidpy as sq                                            
import scanpy as sc
import pandas as pd
import anndata 


bc_atlas_1 = sc.read('../../xe_bc_public/xe_bc_analysis/data_raw/bc_atlas/BC_atlas_xe.h5ad')
bc_atlas_1.obs['dataset'] = 'atlas'

# sc.pp.filter_cells(sc_adata, min_genes=200)
# sc.pp.filter_genes(sc_adata, min_cells=100)

adata_frp2 = sc.read_10x_h5('../../xe_bc_public/xe_bc_analysis/data_raw/frp/count_sample_filtered_feature_bc_matrix.h5')
adata_frp.var = adata_frp.var.reset_index()
adata_frp = adata_frp[:, adata_frp.var['index'].isin(sc_adata.var_names)] 
adata_frp.var_names = adata_frp.var['index']
# adata_frp = adata_frp[:, sc_adata.var_names]
obs = pd.read_csv('../../xe_bc_public/xe_bc_analysis/data_raw/frp/count_sample_filtered_barcodes.csv', header=None)
adata_frp.obs = obs

adata_frp.X = adata_frp.X.toarray()
adata_frp.obs['dataset'] = 'frp'
sc_adata = sc_adata[:, adata_frp.var_names]

sc.pp.filter_cells(adata_frp, min_genes=10)
sc.pp.normalize_total(adata_frp, target_sum=1e4)
sc.pp.log1p(adata_frp)


sc.pp.pca(sc_adata)
sc.pp.neighbors(sc_adata)
sc.tl.umap(sc_adata)
# sc.pl.umap(sc_adata, color='celltype_major')
sc.tl.ingest(adata_frp, sc_adata, obs='celltype_major')
# adata_frp.obs['celltype_minor'] = 'NA'

ad = {}
for t in adata_frp.obs['celltype_major'].unique():
    print(t)
    ad3 = adata_frp[adata_frp.obs['celltype_major']==t]
    sc2 = sc_adata[sc_adata.obs['celltype_major']==t]
    sc.pp.pca(sc2)
    sc.pp.neighbors(sc2)
    sc.tl.umap(sc2)
    # sc.pl.umap(sc_adata[sc_adata.obs['celltype_major']==t], color='celltype_major')
    sc.tl.ingest(ad3, sc2, obs='celltype_minor', inplace=True)
    ad[t] = ad3
	# adata_frp.obs['celltype_minor'] = 'NA'

adata_frp = anndata.concat(ad)
adata_frp.write('data_tidy/adata_frp.h5ad')
frp = adata_frp

ad_xe_1 = sc.read('data_tidy/xe_human_bc/xe_rep1.h5ad')
ad_xe_2 = sc.read('data_tidy/xe_human_bc/xe_rep2.h5ad')

sc.pp.subsample(ad_xe_1, 0.2)
sc.pp.subsample(ad_xe_2, 0.2)
sc.pp.subsample(bc_atlas_1, 0.2)
sc.pp.subsample(frp, 0.2)

frp.var_names_make_unique()
genes = list(set(frp.var_names) & set(ad_xe_1.var_names))
ad_xe_1    = ad_xe_1[:, genes]
ad_xe_2    = ad_xe_2[:, genes]
bc_atlas_1 = bc_atlas_1[:, genes]
frp = frp[:, genes]

sc.pp.normalize_total(ad_xe_1 , target_sum=1)
sc.pp.normalize_total(bc_atlas_1 , target_sum=1)
sc.pp.normalize_total(ad_xe_2 , target_sum=1)
sc.pp.normalize_total(frp, target_sum=1)



