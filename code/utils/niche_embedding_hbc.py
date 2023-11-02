import squidpy as sq                                            
import scanpy as sc
import pandas as pd
import numpy as np
import anndata 
import os, sys
from pathlib import Path
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import OneHotEncoder   
from sklearn.metrics import pairwise_distances  




models_f = 'models/human_breast_cancer'

type_dists = np.empty((0, 9))
barcodes   = []
samples    = []
adatas     = []
for f in os.listdir(models_f):
	print(f)
	if f != 'atlas_subset' and f != 'CytAssist_Fresh_Frozen_Human_Breast_Cancer':
		adata_f = Path(models_f, f, 'sp.h5ad')
		adata   = sc.read(adata_f)
		sq.gr.spatial_neighbors(adata)            
		ohe = OneHotEncoder()                                                                           
		type_dist = np.dot(adata.obsp['spatial_connectivities'].toarray(), adata.obsm['q95_cell_abundance_w_sf']) + adata.obsm['q95_cell_abundance_w_sf']
		# type_dist /= type_dist.sum(axis=1)      
		type_dists = np.append(type_dists, type_dist, axis = 0) 
		barcodes.append(adata.obs_names)
		samples.append([f] * adata.shape[0])
		adatas.append(adata)

	# adata.obsm['type_dist'] = np.array(type_dist)
	# # sc.pp.neighbors(adata, use_rep='type_dist', metric='cosine')
	# A = kneighbors_graph(np.array(type_dist), 100, mode='connectivity', include_self=True, metric='cosine')
	# adata.obsp['niche_connectivities_type'] = A 
	# sc.tl.leiden(adata, adjacency=adata.obsp['niche_connectivities_type'], key_added='type_niche', resolution=0.5) 
	# # adata.write('data_tidy/codex_niched_v2.h5ad')
	# sc.pl.spatial(adata, spot_size=150, color = 'type_niche', save='test_visium_niche.pdf')

barcodes = [item for sublist in barcodes for item in sublist]
samples  = [item for sublist in samples for item in sublist]

meta = pd.DataFrame({'barcode': barcodes, 'sample': samples})
embedding = pd.DataFrame(type_dists)

df = pd.concat([meta, embedding], axis=1)
df.to_csv('data_tidy/human_breast_cancer/visium/niche_embedding.csv')



adata_comb = anndata.concat(adatas)

adata_comb.obs['sample'] = samples

adata_comb.obsm['type_dist'] = type_dists
sc.pp.neighbors(adata_comb, use_rep='type_dist', metric='cosine')
A = kneighbors_graph(np.array(type_dists), 100, mode='connectivity', include_self=True, metric='cosine')
adata_comb.obsp['niche_connectivities_type'] = A 
sc.tl.leiden(adata_comb, adjacency=adata_comb.obsp['niche_connectivities_type'], key_added='type_niche', resolution=0.5) 



for f in os.listdir(models_f):
	print(f)
	if f != 'atlas_subset' and f != 'CytAssist_Fresh_Frozen_Human_Breast_Cancer':
		adata_f = Path(models_f, f, 'sp.h5ad')
		adata   = sc.read(adata_f)
		adata.obs = adata.obs.merge(adata_comb[adata_comb.obs['sample'] == f].obs)
		sc.pl.spatial(adata, spot_size=150, color = 'type_niche', save= f +'_niched.pdf')






# sc.pl.spatial(adata, spot_size=150, save='test_visium_vis.pdf') 