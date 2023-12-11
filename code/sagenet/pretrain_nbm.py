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
from scipy.sparse.csgraph import connected_components


coldt = pd.read_csv('/omics/groups/OE0606/internal/hluo/data/Neuroblastoma/internal/Neuroblastoma2/new_annots.csv')
coldt.index = coldt.iloc[:,0]

ad_xe_827 = sc.read('data_tidy/NBS827/nbx827.h5ad')
ad_xe_827.obsm['spatial'] = np.array(ad_xe_827.obs[['cell_centroid_x', 'cell_centroid_y']])
ad_sc_827 = sc.read('data_tidy/NBS827/scrnaseq.h5ad')
ad_sc_827.obs['clusters_0'] = coldt.loc[ad_sc_827.obs_names,'clusters']
ad_sc_827.obs['clusters_1'] = coldt.loc[ad_sc_827.obs_names,'seurat_clusters']
ad_sc_827.obs['Annotation_new'] = [decoding[item] for item in ad_sc_827.obs['clusters_1'].values.astype('str')]
with open('/omics/groups/OE0606/internal/hluo/data/Neuroblastoma/internal/Neuroblastoma2/feature_names.txt', 'r') as file:
    # Read all the lines of the file into a list
    lines = file.read().splitlines() 
    
    
    
ad_sc_827.var_names = lines
gg = (set(ad_sc_827.var_names) & set(ad_xe_827.var_names))

ad_xe_827 = ad_xe_827[:, list(gg)]
ad_sc_827 = ad_sc_827[:, list(gg)]
sc.pp.normalize_total(ad_xe_827 , target_sum=1)
sc.pp.normalize_total(ad_sc_827 , target_sum=1)
    
    
sq.gr.spatial_neighbors(ad_xe_827, radius=200, coord_type = 'generic')
ad_xe_827.obs['core'] = connected_components(ad_xe_827.obsp['spatial_connectivities'])[1].astype('str')
core_counts = ad_xe_827.obs['core'].value_counts()
# Filter rows with core counts >= 100
cores_to_keep = core_counts[core_counts >= 10000].index
# Filter AnnData object to keep only rows with cores >= 100
ad_xe_827 = ad_xe_827[ad_xe_827.obs['core'].isin(cores_to_keep)]
# ad_xe_827.obs['core'] = ad_xe_827.obs[['sample', 'core']].agg('_core_'.join, axis=1)






if torch.cuda.is_available():  
  dev = "cuda:0"  
else:  
  dev = "cpu"  
  

device = torch.device(dev)
print(device)

sg_obj = sg.sage.sage(device=device)

for core in ad_xe_827.obs['core'].unique()[1:]:
    print('core')
    ad_xe = ad_xe_827[ad_xe_827.obs['core'] == core]
    glasso(ad_xe)
    ad_xe.obsm['spatial'] = np.array(ad_xe.obs[['cell_centroid_x', 'cell_centroid_y']])
    sq.gr.spatial_neighbors(ad_xe, coord_type="generic")
    sc.tl.leiden(ad_xe, resolution=.01, random_state=0, key_added='leiden_0.01', adjacency=ad_xe.obsp["spatial_connectivities"])
    sc.tl.leiden(ad_xe, resolution=.1, random_state=0, key_added='leiden_0.1', adjacency=ad_xe.obsp["spatial_connectivities"])
    sc.tl.leiden(ad_xe, resolution=1, random_state=0, key_added='leiden_1', adjacency=ad_xe.obsp["spatial_connectivities"])
    sg_obj.add_ref(ad_xe, comm_columns=['leiden_0.01', 'leiden_0.1', 'leiden_1'], tag=core, epochs=50, verbose = True, classifier='GraphSAGE')
    # sg_obj.add_ref(ad_xe_827, comm_columns=['leiden_0.1', 'leiden_0.5', 'leiden_1'], tag='xe1', epochs=20, verbose = True, classifier='GraphSAGE')
    # sg_obj.add_ref(ad_xe_827, comm_columns=['leiden_0.1', 'leiden_0.5', 'leiden_1'], tag='xe1', epochs=20, verbose = True, classifier='GraphSAGE')






os.makedirs('models/nbx827_model_2023_11_28')
sg_obj.save_model_as_folder('models/nbx827_model_2023_11_28')
sg_obj.load_model_as_folder('models/nbx827_model_subset')



ad_xe = {}
for core in ad_xe_827.obs['core'].unique():
    print(core)
    ad_xe = ad_xe_827[ad_xe_827.obs['core'] == core]
    ad_xe.X = ad_xe.X.toarray()
    # sc.pp.subsample(ad_xe, .25)
    sg_obj.map_query(ad_xe, save_prob=True, save_dist=False)
    ad_xe.write('data_tidy/NBS827/xe_827_%s.h5ad' % core)

ad_sc_827.X = ad_sc_827.X.toarray()
sg_obj.map_query(ad_sc_827, save_prob=True)
ad_sc_827.write('data_tidy/NBS827/sc_827.h5ad')

ad_sc_827 = sc.read('data_tidy/NBS827/sc_827.h5ad')
ad_xe_0 = sc.read('data_tidy/NBS827/xe_827_0.h5ad')
ad_xe_1 = sc.read('data_tidy/NBS827/xe_827_1.h5ad')
ad_xe_2 = sc.read('data_tidy/NBS827/xe_827_2.h5ad')

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


ad_xe_0 = prob_con(ad_xe_0)
ad_sc_827 = prob_con(ad_sc_827)



from torch.nn import Softmax
    
reference_embeddings = ad_sc_827.obsm['prob_concatenated']
# reference_embeddings= reference_embeddings/reference_embeddings.sum(axis=1)[:,None]
kdtree_r1 = cKDTree(reference_embeddings)
target_embeddings = ad_xe_0.obsm['prob_concatenated']
# target_embeddings= target_embeddings/target_embeddings.sum(axis=1)[:,None]
k_neighbors = 1  # You can adjust this as needed
distances, indices = kdtree_r1.query(target_embeddings, k=k_neighbors)


# m = Softmax(dim=1)
# probs = m(-torch.tensor(distances))
# dist=torch.distributions.categorical.Categorical(probs=probs)
# idx = dist.sample().numpy()
# indices = indices[np.arange(len(indices)), idx]
# x_coordinates = ad_xe_0.obs['cell_centroid_x'][indices]
# y_coordinates = ad_xe_0.obs['cell_centroid_y'][indices]
# coordinates = np.column_stack((x_coordinates, y_coordinates))
# ad_xe_0.obsm['spatial'] = coordinates

ad_xe_0.obs['clusters'] = ad_sc_827.obs['clusters_0'][indices].values.astype('str')
# ad_sc_827.obs['sink'] = 0
# ind, counts = np.unique(indices, return_counts=True)
# ad_sc_827.obs['sink'][ind] = np.log10(counts)
# swap.obs['sink']





sc.pl.spatial(
    ad_xe_0,
    color='clusters',
    # color=['sink', 'transcript_counts'],
    # palette=major_colors, # Color cells based on 'cell_type'
    # color_map=cell_type_color_map,  # Use the custom color map
    # library_id='r1_mapping',  # Use 'r1_mapping' coordinates
    title='Xenium Locs',
    save='_xe0_sagenet_clusters.pdf',
    spot_size=10
)

ad_xe_0.obsm['X_spatial'] = ad_xe_0.obsm['spatial']
sc.tl.embedding_density(ad_xe_0, basis='spatial', groupby='clusters')
sc.pl.embedding_density(ad_xe_0, basis='spatial', key='spatial_density_clusters', save='xe0_sagenet.pdf')

decoding = {
    '0': 'Cycling',
    '1': 'Tumor',
    '2': 'Tumor',
    '3': 'Tumor',
    '4': 'Tumor',
    '5': 'Tumor',
    '6': 'Tumor',
    '7': 'Tumor',
    '8': 'DBH+',
    '9': 'Tumor',
    '10': 'Tumor',
    '11': 'Cycling',
    '12': 'Tumor',
    '13': 'Tumor',
    '14': 'Cycling',
    '15': 'Tumor',
    '16': 'Tumor',
    '17': 'Tumor',
    '18': 'Tumor',
    '19': 'Tumor',
    '20': 'DBH+ DLK1+',
    '21': 'Tumor',
    '22': 'Tumor',
    '23': 'Cycling',
    '24': 'M@',
    '25': 'Tumor',
    '26': 'Endothelial',
    '27': 'Lymphocyte',
    '28': 'Mesenchymal',
    '29': 'Cycling',
    '30': 'Cycling 2'
    
}


major_colors = {
        'M@' : '#d8f55e',
        'Mesenchymal' : '#532C8A',
        'DBH+ DLK1+' : '#C72228',
        'Endothelial' : '#9e6762',
        'Myeloid' : '#ffe012',
        'Lymphocyte' : '#3cb44b',
        'Cycling' : '#8DB5CE',
        'Cycling 2' : '#C9EBFB',
        'DBH+': '#ff891c',
        'Tumor' : '#EF5A9D'
}


  "Anterior Primitive Streak" = "#c19f70",
                     "Notochord" = "#0F4A9C",
                     "Def. endoderm" = "#F397C0",
                     "Definitive endoderm" = "#F397C0",
                     "Gut" = "#EF5A9D",
                     "Gut tube" = "#EF5A9D",
                     
                     "Nascent mesoderm" = "#C594BF",
                     "Mixed mesoderm" = "#DFCDE4",
                     "Intermediate mesoderm" = "#139992",
                     "Caudal Mesoderm" = "#3F84AA",
                     "Paraxial mesoderm" = "#8DB5CE",
                     "Somitic mesoderm" = "#005579",
                     "Pharyngeal mesoderm" = "#C9EBFB",
                     "Splanchnic mesoderm" = "#C9EBFB",
                     "Cardiomyocytes" = "#B51D8D",
                     "Allantois" = "#532C8A",
                     "ExE mesoderm" = "#8870ad",
                     "Lateral plate mesoderm" = "#8870ad",
                     "Mesenchyme" = "#cc7818",
                     "Mixed mesenchymal mesoderm" = "#cc7818",
                     
                     "Haematoendothelial progenitors" = "#FBBE92",
                     "Endothelium" = "#ff891c",
                     "Blood progenitors 1" = "#f9decf",
                     "Blood progenitors 2" = "#c9a997",
                     
                     "Erythroid1" = "#C72228",
                     "Erythroid2" = "#f79083",
                     "Erythroid3" = "#EF4E22",
                     
                     "Erythroid" = "#f79083",
                     "Blood progenitors" = "#f9decf",
                     
                     "NMP" = "#8EC792",
                     
                     "Rostral neurectoderm" = "#65A83E",
                     "Caudal neurectoderm" = "#354E23",
                     "Neural crest" = "#C3C388",
                     "Forebrain/Midbrain/Hindbrain" = "#647a4f",
                     "Spinal cord" = "#CDE088",
                     
                     "Surface ectoderm" = "#f7f79e",
                     