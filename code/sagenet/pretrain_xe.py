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


ad_xe_1 = sc.read('data_tidy/xe_human_bc/xe_rep1.h5ad')
ad_xe_2 = sc.read('data_tidy/xe_human_bc/xe_rep2.h5ad')
bc_atlas_1 = sc.read('data_tidy/xe_human_bc/bc_atlas_raw.h5ad')
frp = sc.read('data_tidy/adata_frp.h5ad')
frp = adata_frp
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



glasso(ad_xe_1, [0.25, 0.5])
ad_xe_1.obsm['spatial'] = np.array(ad_xe_1.obs[['x_centroid', 'y_centroid']])
sq.gr.spatial_neighbors(ad_xe_1, coord_type="generic")
sc.tl.leiden(ad_xe_1, resolution=.01, random_state=0, key_added='leiden_0.01', adjacency=ad_xe_1.obsp["spatial_connectivities"])
sc.tl.leiden(ad_xe_1, resolution=.05, random_state=0, key_added='leiden_0.05', adjacency=ad_xe_1.obsp["spatial_connectivities"])
sc.tl.leiden(ad_xe_1, resolution=.1, random_state=0, key_added='leiden_0.1', adjacency=ad_xe_1.obsp["spatial_connectivities"])
sc.tl.leiden(ad_xe_1, resolution=.5, random_state=0, key_added='leiden_0.5', adjacency=ad_xe_1.obsp["spatial_connectivities"])
sc.tl.leiden(ad_xe_1, resolution=1, random_state=0, key_added='leiden_1', adjacency=ad_xe_1.obsp["spatial_connectivities"])
sc.tl.leiden(ad_xe_1, resolution=1, random_state=0, key_added='leiden_2', adjacency=ad_xe_1.obsp["spatial_connectivities"])


glasso(ad_xe_2, [0.25, 0.5])
ad_xe_2.obsm['spatial'] = np.array(ad_xe_2.obs[['x_centroid', 'y_centroid']])
sq.gr.spatial_neighbors(ad_xe_2, coord_type="generic")
sc.tl.leiden(ad_xe_2, resolution=.01, random_state=0, key_added='leiden_0.01', adjacency=ad_xe_2.obsp["spatial_connectivities"])
sc.tl.leiden(ad_xe_2, resolution=.05, random_state=0, key_added='leiden_0.05', adjacency=ad_xe_2.obsp["spatial_connectivities"])
sc.tl.leiden(ad_xe_2, resolution=.1, random_state=0, key_added='leiden_0.1', adjacency=ad_xe_2.obsp["spatial_connectivities"])
sc.tl.leiden(ad_xe_2, resolution=.5, random_state=0, key_added='leiden_0.5', adjacency=ad_xe_2.obsp["spatial_connectivities"])
sc.tl.leiden(ad_xe_2, resolution=1, random_state=0, key_added='leiden_1', adjacency=ad_xe_2.obsp["spatial_connectivities"])
sc.tl.leiden(ad_xe_2, resolution=1, random_state=0, key_added='leiden_2', adjacency=ad_xe_2.obsp["spatial_connectivities"])


if torch.cuda.is_available():  
  dev = "cuda:0"  
else:  
  dev = "cpu"  
  

device = torch.device(dev)
print(device)

sg_obj = sg.sage.sage(device=device)
sg_obj.add_ref(ad_xe_1, comm_columns=['leiden_0.01', 'leiden_0.1', 'leiden_1'], tag='xe1', epochs=20, verbose = True, classifier='GraphSAGE')
sg_obj.add_ref(ad_xe_2, comm_columns=['leiden_0.01', 'leiden_0.1', 'leiden_1'], tag='xe2', epochs=20, verbose = True, classifier='GraphSAGE')
os.makedirs('models/xe_model')
sg_obj.save_model_as_folder('models/xe_model')
sg_obj.load_model_as_folder('models/xe_model')

sg_obj.map_query(ad_xe_1, save_prob=True)
frp.X = frp.X.toarray()
sg_obj.map_query(frp, save_prob=True)
sg_obj.map_query(ad_xe_2, save_prob=True)
bc_atlas_1.X = bc_atlas_1.X.toarray()
sg_obj.map_query(bc_atlas_1, save_prob=True)

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


ad_xe_1 = prob_con(ad_xe_1)
ad_xe_2 = prob_con(ad_xe_2)
bc_atlas_1 = prob_con(bc_atlas_1)
frp = prob_con(frp)



from torch.nn import Softmax
    
reference_embeddings = ad_xe_2.obsm['prob_concatenated']
# reference_embeddings= reference_embeddings/reference_embeddings.sum(axis=1)[:,None]
kdtree_r1 = cKDTree(reference_embeddings)
target_embeddings = frp.obsm['prob_concatenated']
# target_embeddings= target_embeddings/target_embeddings.sum(axis=1)[:,None]
k_neighbors = 1  # You can adjust this as needed
distances, indices = kdtree_r1.query(target_embeddings, k=k_neighbors)


# m = Softmax(dim=1)
# probs = m(-torch.tensor(distances))
# dist=torch.distributions.categorical.Categorical(probs=probs)
# idx = dist.sample().numpy()
# indices = indices[np.arange(len(indices)), idx]
x_coordinates = ad_xe_2.obs['x_centroid'][indices]
y_coordinates = ad_xe_2.obs['y_centroid'][indices]
coordinates = np.column_stack((x_coordinates, y_coordinates))
frp.obsm['spatial'] = coordinates
ad_xe_2.obs['sink'] = 0
ind, counts = np.unique(indices, return_counts=True)
ad_xe_2.obs['sink'][ind] = np.log10(counts)
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
    ad_xe_2,
    color='celltype_major',
    # color=['sink', 'transcript_counts'],
    # palette=major_colors, # Color cells based on 'cell_type'
    # color_map=cell_type_color_map,  # Use the custom color map
    # library_id='r1_mapping',  # Use 'r1_mapping' coordinates
    title='Xenium Replicate 2',
    save='_xe2_true.pdf',
    spot_size=30
)