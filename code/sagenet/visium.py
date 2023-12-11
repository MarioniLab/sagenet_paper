import sagenet as sg
import scanpy as sc
import squidpy as sq
import anndata as ad
import random
import torch
from sagenet.utils import glasso
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import os


visium = sq.read.visium('../../xe_bc_public/xe_bc_analysis/data_raw/visium')
visium.var_names_make_unique()
visium.X = visium.X.toarray()
out_spots = np.asarray(np.where(np.isnan(visium.obsm['spatial'])))[0][0]
visium = visium[~np.isnan(visium.obsm['spatial'])[:,0],:]
visium = visium[:,~(visium.X.var(0) == 0)]

xe1    = sc.read('data_tidy/xe_human_bc/xe_rep1.h5ad')
xe2    = sc.read('data_tidy/xe_human_bc/xe_rep2.h5ad')
atlas  = sc.read('data_tidy/xe_human_bc/bc_atlas_raw.h5ad')

genes  = list(set(visium.var_names) & set(xe1.var_names) & set(atlas.var_names))
visium = visium[:, genes]
xe1    = xe1[:, genes]
xe2    = xe2[:, genes]
atlas  = atlas[:, genes]

sc.pp.normalize_total(xe1 , target_sum=1)
sc.pp.normalize_total(xe2 , target_sum=1)
sc.pp.normalize_total(atlas , target_sum=1)



# visium.varm['adj'] = xe1.varm['adj']
# del visium.uns
sq.gr.spatial_neighbors(visium, coord_type="grid", n_rings=1)
sc.tl.leiden(visium, resolution=.01, random_state=0, key_added='leiden_0.01', adjacency=visium.obsp["spatial_connectivities"])
sc.tl.leiden(visium, resolution=.05, random_state=0, key_added='leiden_0.05', adjacency=visium.obsp["spatial_connectivities"])
sc.tl.leiden(visium, resolution=.1, random_state=0, key_added='leiden_0.1', adjacency=visium.obsp["spatial_connectivities"])
sc.tl.leiden(visium, resolution=.5, random_state=0, key_added='leiden_0.5', adjacency=visium.obsp["spatial_connectivities"])
sc.tl.leiden(visium, resolution=1, random_state=0, key_added='leiden_1', adjacency=visium.obsp["spatial_connectivities"])
sc.tl.leiden(visium, resolution=2, random_state=0, key_added='leiden_2', adjacency=visium.obsp["spatial_connectivities"])
sc.tl.leiden(visium, resolution=5, random_state=0, key_added='leiden_5', adjacency=visium.obsp["spatial_connectivities"])
sc.tl.leiden(visium, resolution=10, random_state=0, key_added='leiden_10', adjacency=visium.obsp["spatial_connectivities"])
sc.tl.leiden(visium, resolution=10, random_state=0, key_added='leiden_20', adjacency=visium.obsp["spatial_connectivities"])
sc.tl.leiden(visium, resolution=10, random_state=0, key_added='leiden_15', adjacency=visium.obsp["spatial_connectivities"])
sc.pp.filter_cells(visium, np.quantile(visium.X.sum(1), .25))
sc.pp.normalize_total(visium, target_sum=1)
glasso(visium)

glasso(xe1)
xe1.obsm['spatial'] = np.array(xe1.obs[['x_centroid', 'y_centroid']])
sq.gr.spatial_neighbors(xe1, coord_type="generic")
sc.tl.leiden(xe1, resolution=.01, random_state=0, key_added='leiden_0.01', adjacency=xe1.obsp["spatial_connectivities"])
sc.tl.leiden(xe1, resolution=.05, random_state=0, key_added='leiden_0.05', adjacency=xe1.obsp["spatial_connectivities"])
sc.tl.leiden(xe1, resolution=.1, random_state=0, key_added='leiden_0.1', adjacency=xe1.obsp["spatial_connectivities"])
sc.tl.leiden(xe1, resolution=.5, random_state=0, key_added='leiden_0.5', adjacency=xe1.obsp["spatial_connectivities"])
sc.tl.leiden(xe1, resolution=1, random_state=0, key_added='leiden_1', adjacency=xe1.obsp["spatial_connectivities"])
sc.tl.leiden(xe1, resolution=1, random_state=0, key_added='leiden_2', adjacency=xe1.obsp["spatial_connectivities"])


glasso(xe2)
xe2.obsm['spatial'] = np.array(xe2.obs[['x_centroid', 'y_centroid']])
sq.gr.spatial_neighbors(xe2, coord_type="generic")
sc.tl.leiden(xe2, resolution=.01, random_state=0, key_added='leiden_0.01', adjacency=xe2.obsp["spatial_connectivities"])
sc.tl.leiden(xe2, resolution=.05, random_state=0, key_added='leiden_0.05', adjacency=xe2.obsp["spatial_connectivities"])
sc.tl.leiden(xe2, resolution=.1, random_state=0, key_added='leiden_0.1', adjacency=xe2.obsp["spatial_connectivities"])
sc.tl.leiden(xe2, resolution=.5, random_state=0, key_added='leiden_0.5', adjacency=xe2.obsp["spatial_connectivities"])
sc.tl.leiden(xe2, resolution=1, random_state=0, key_added='leiden_1', adjacency=xe2.obsp["spatial_connectivities"])
sc.tl.leiden(xe2, resolution=1, random_state=0, key_added='leiden_2', adjacency=xe2.obsp["spatial_connectivities"])


# sc.pp.subsample(xe1, 0.2)
# sc.pp.subsample(xe2, 0.2)
# sc.pp.subsample(atlas, 0.2)

sc.pl.spatial(
    visium,
    color='leiden_2',
    # color=['sink', 'transcript_counts', 'leiden_0.01', 'leiden_0.1', 'leiden_1'],
    # palette=celltype_colours, # Color cells based on 'cell_type'
    # color_map=cell_type_color_map,  # Use the custom color map
    # library_id='r1_mapping',  # Use 'r1_mapping' coordinates
    # title='Spatial Plot with Cell Type Coloring',
    save='_test_visium.pdf',
    alpha=0.5,
    spot_size=200
    # ,
    # spot_size=30
)

if torch.cuda.is_available():  
  dev = "cuda:0"  
else:  
  dev = "cpu"  
  

device = torch.device(dev)
print(device)

sg_obj = sg.sage.sage(device=device)
sg_obj.add_ref(visium, comm_columns=['leiden_1', 'leiden_2', 'leiden_5'], tag='visium', epochs=1000, verbose = True, classifier='GraphSAGE')
sg_obj.add_ref(visium, comm_columns=['leiden_10', 'leiden_15', 'leiden_20'], tag='visium2', epochs=1000, verbose = True, classifier='GraphSAGE')
sg_obj.add_ref(xe1, comm_columns=['leiden_0.1', 'leiden_1', 'leiden_1'], tag='xe1', epochs=30, verbose = True, classifier='GraphSAGE')
sg_obj.add_ref(xe2, comm_columns=['leiden_0.1', 'leiden_1', 'leiden_1'], tag='visium', epochs=30, verbose = True, classifier='GraphSAGE')
os.makedirs('models/bc_ens_all_model')
sg_obj.save_model_as_folder('models/bc_ens_all_model')
sg_obj.load_model_as_folder('models/bc_ens_all_model')

sg_obj.map_query(visium, save_prob=True)
sg_obj.map_query(xe1, save_prob=True)
sg_obj.map_query(xe2, save_prob=True)
atlas.X = atlas.X.toarray()
sg_obj.map_query(atlas, save_prob=True)

import anndata
from functools import reduce
import numpy as np

def prob_con(adata):
    # Get a list of obsm matrices with names starting with "prob"
    prob_matrices = [matrix_name for matrix_name in adata.obsm.keys() if matrix_name.startswith("prob")]
    # Define a function to concatenate two matrices
    def concatenate_matrices(matrix1, matrix2):
        return np.concatenate((matrix1, matrix2), axis=1)
    # Use functools.reduce to concatenate all matrices in prob_matrices
    if prob_matrices:
        concatenated_matrix = reduce(concatenate_matrices, [adata.obsm[matrix] for matrix in prob_matrices])
        adata.obsm["prob_concatenated"] = concatenated_matrix
    else:
        print("No 'prob' matrices found in the AnnData object.")
    return adata

visium = prob_con(visium)
xe1 = prob_con(xe1)
xe2 = prob_con(xe2)
atlas = prob_con(atlas)
# frp = prob_con(frp)



from torch.nn import Softmax
    
reference_embeddings = visium.obsm['prob_concatenated']
# reference_embeddings= reference_embeddings/reference_embeddings.sum(axis=1)[:,None]
kdtree_r1 = cKDTree(reference_embeddings)
target_embeddings = atlas.obsm['prob_concatenated']
# target_embeddings= target_embeddings/target_embeddings.sum(axis=1)[:,None]
k_neighbors = 10  # You can adjust this as needed
distances, indices = kdtree_r1.query(target_embeddings, k=k_neighbors)


m = Softmax(dim=1)
probs = m(-torch.tensor(distances))
dist=torch.distributions.categorical.Categorical(probs=probs)
idx = dist.sample().numpy()
indices = indices[np.arange(len(indices)), idx]
# x_coordinates = visium.obsm['spatial'][indices]
# y_coordinates = visium.obs['y_centroid'][indices]
# coordinates = np.column_stack((x_coordinates, y_coordinates))
atlas.obsm['spatial'] = visium.obsm['spatial'][indices] 
atlas.obsm['spatial'] += np.random.uniform(-50, 50, atlas.obsm['spatial'].shape)
visium.obs['sink'] = 0
ind, counts = np.unique(indices, return_counts=True)
visium.obs['sink'][ind] = np.log10(counts)
# swap.obs['sink']



major_colors = {
        'B-cells' : '#d8f55e',
        'CAFs' : '#532C8A',
        'Cancer Epithelial' : '#C72228',
        'Endothelial' : '#9e6762',
        'Myeloid' : '#ffe012',
        'T-cells' : '#3cb44b',
        'Normal Epithelial' : '#0F4A9C',
        'PVL': '#c09d9a',
        'Plasmablasts' : '#000075'
}


sc.pl.spatial(
    atlas,
    color='celltype_major',
    # color=['sink', 'transcript_co-unts'],
    palette=major_colors, # Color cells based on 'cell_type'
    # color_map=cell_type_color_map,  # Use the custom color map
    # library_id='r1_mapping',  # Use 'r1_mapping' coordinates
    title='xe1 -> visium',
    save='_xe1_visium5.pdf',
    alpha=0.5,
    spot_size=100
)


import tangram as tg
tg.pp_adatas(atlas, visium, genes=None)
ad_map = tg.map_cells_to_space(atlas, visium, device=device)
dist_np = ad_map.X
mapped_spots = np.argmax(dist_np, axis=1)
# x_coordinates = visium.obs['x_centroid'][mapped_spots]
# y_coordinates = visium.obs['y_centroid'][mapped_spots]
# coordinates = np.column_stack((x_coordinates, y_coordinates))
# coordinates = np.column_stack((x_coordinates, y_coordinates))
atlas.obsm['spatial'] = visium.obsm['spatial'][mapped_spots]

visium.obs['leiden_2'] = visium.obs['leiden_2'].astype(str).astype('category')
sc.pl.spatial(
    visium,
    color='leiden_2',
    # color=['sink', 'transcript_counts', 'leiden_0.01', 'leiden_0.1', 'leiden_1'],
    # palette=celltype_colours, # Color cells based on 'cell_type'
    # color_map=cell_type_color_map,  # Use the custom color map
    # library_id='r1_mapping',  # Use 'r1_mapping' coordinates
    title='Spatial Plot with Cell Type Coloring',
    save='_atlas_visium_tangram.pdf',
    spot_size=100
)







reference_embeddings = atlas.obsm['prob_concatenated']
# reference_embeddings= reference_embeddings/reference_embeddings.sum(axis=1)[:,None]
kdtree_r1 = cKDTree(reference_embeddings)
target_embeddings = visium.obsm['prob_concatenated']
# target_embeddings= target_embeddings/target_embeddings.sum(axis=1)[:,None]
k_neighbors = 5  # You can adjust this as needed
distances, indices = kdtree_r1.query(target_embeddings, k=k_neighbors)


visium.obs['type0'] = atlas.obs['celltype_major'][indices[:,0]].values
visium.obs['type1'] = atlas.obs['celltype_major'][indices[:,1]].values
visium.obs['type2'] = atlas.obs['celltype_major'][indices[:,2]].values
visium.obs['type3'] = atlas.obs['celltype_major'][indices[:,3]].values
visium.obs['type4'] = atlas.obs['celltype_major'][indices[:,4]].values

sc.pl.spatial(
    visium,
    color=['type0', 'type1', 'type2', 'type3', 'type4'],
    # color=['sink', 'transcript_counts', 'leiden_0.01', 'leiden_0.1', 'leiden_1'],
    palette=major_colors, # Color cells based on 'cell_type'
    # color_map=cell_type_color_map,  # Use the custom color map
    # library_id='r1_mapping',  # Use 'r1_mapping' coordinates
    # title='Spatial Plot with Cell Type Coloring',
    save='_visium_atlas_types.pdf',
    spot_size=200
)