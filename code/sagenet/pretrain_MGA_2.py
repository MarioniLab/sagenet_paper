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

random.seed(1996)


adata_r1 = sg.MGA_data.seqFISH1_1()
adata_r2 = sg.MGA_data.seqFISH2_1()
adata_r3 = sg.MGA_data.seqFISH3_1()
adata_q  = sg.MGA_data.scRNAseq()


glasso(adata_r1, [0.25, 0.5])
adata_r1.obsm['spatial'] = np.array(adata_r1.obs[['x','y']])

sq.gr.spatial_neighbors(adata_r1, coord_type="generic")
sc.tl.leiden(adata_r1, resolution=.01, random_state=0, key_added='leiden_0.01', adjacency=adata_r1.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r1, resolution=.05, random_state=0, key_added='leiden_0.05', adjacency=adata_r1.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r1, resolution=.1, random_state=0, key_added='leiden_0.1', adjacency=adata_r1.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r1, resolution=.5, random_state=0, key_added='leiden_0.5', adjacency=adata_r1.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r1, resolution=1, random_state=0, key_added='leiden_1', adjacency=adata_r1.obsp["spatial_connectivities"])
# adata_r1.obs.to_csv('int_data/seqfish_mouse_embryo/col_dt_embryo1_2.txt', sep='\t')

glasso(adata_r2, [0.25, 0.5])
adata_r2.obsm['spatial'] = np.array(adata_r2.obs[['x','y']])
sq.gr.spatial_neighbors(adata_r2, coord_type="generic")
sc.tl.leiden(adata_r2, resolution=.01, random_state=0, key_added='leiden_0.01', adjacency=adata_r2.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r2, resolution=.05, random_state=0, key_added='leiden_0.05', adjacency=adata_r2.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r2, resolution=.1, random_state=0, key_added='leiden_0.1', adjacency=adata_r2.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r2, resolution=.5, random_state=0, key_added='leiden_0.5', adjacency=adata_r2.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r2, resolution=1, random_state=0, key_added='leiden_1', adjacency=adata_r2.obsp["spatial_connectivities"])
# adata_r2.obs.to_csv('int_data/seqfish_mouse_embryo/col_dt_embryo2_2.txt', sep='\t')

glasso(adata_r3, [0.25, 0.5])
adata_r3.obsm['spatial'] = np.array(adata_r3.obs[['x','y']])
sq.gr.spatial_neighbors(adata_r3, coord_type="generic")
sc.tl.leiden(adata_r3, resolution=.01, random_state=0, key_added='leiden_0.01', adjacency=adata_r3.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r3, resolution=.05, random_state=0, key_added='leiden_0.05', adjacency=adata_r3.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r3, resolution=.1, random_state=0, key_added='leiden_0.1', adjacency=adata_r3.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r3, resolution=.5, random_state=0, key_added='leiden_0.5', adjacency=adata_r3.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r3, resolution=1, random_state=0, key_added='leiden_1', adjacency=adata_r3.obsp["spatial_connectivities"])
# adata_r3.obs.to_csv('int_data/seqfish_mouse_embryo/col_dt_embryo3_2.txt', sep='\t')

if torch.cuda.is_available():  
  dev = "cuda:0"  
else:  
  dev = "cpu"  


device = torch.device(dev)
print(device)

sg_obj = sg.sage.sage(device=device)
sg_obj.add_ref(adata_r1, comm_columns=['leiden_0.01', 'leiden_0.05', 'leiden_0.1','leiden_0.5', 'leiden_1'], tag='embryo1_2', epochs=20, verbose = True, classifier='GraphSAGE')

# # adata_r1.var.to_csv('int_data/seqfish_mouse_embryo/gene_dt_embryo1_2.txt', sep='\t')
# # sg_obj.save_model_as_folder('models/seqfish_mouse_embryo/embryo1_2')
# # sg_obj = sg.sage.sage(device=device)
# sg_obj.add_ref(adata_r2, comm_columns=['leiden_0.01', 'leiden_0.1', 'leiden_1'], tag='embryo2_2', epochs=20, verbose = False, classifier='GraphSAGE')
# # sg_obj.save_model_as_folder('models/seqfish_mouse_embryo/embryo2_2')
# # adata_r2.var.to_csv('int_data/seqfish_mouse_embryo/gene_dt_embryo2_2.txt', sep='\t')
# # sg_obj = sg.sage.sage(device=device)
# sg_obj.add_ref(adata_r3, comm_columns=['leiden_0.01', 'leiden_0.1', 'leiden_1'], tag='embryo3_2', epochs=20, verbose = False, classifier='GraphSAGE')
# # sg_obj.save_model_as_folder('models/seqfish_mouse_embryo/embryo3_2')
# # adata_r3.var.to_csv('int_data/seqfish_mouse_embryo/gene_dt_embryo3_2.txt', sep='\t')


sg_obj.map_query(adata_r1, save_prob=True)
sg_obj.map_query(adata_r2, save_prob=True)
sg_obj.map_query(adata_r3, save_prob=True)
sg_obj.map_query(adata_q, save_prob=True)

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


adata_r1 = prob_con(adata_r1)
adata_r2 = prob_con(adata_r2)
adata_r3 = prob_con(adata_r3)
adata_q = prob_con(adata_q)



# sg_obj = sg.sage.sage(device=device)
# sg_obj.add_ref(adata_r3, comm_columns=['leiden_0.01', 'leiden_0.1', 'leiden_1'], tag='embryo3_2', epochs=20, verbose = False, classifier='GraphSAGE')
# sg_obj.add_ref(adata_r3, comm_columns=['leiden_0.01', 'leiden_0.1', 'leiden_1'], tag='embryo3_2', epochs=20, verbose = False, classifier='GraphSAGE')
# sg_obj.add_ref(adata_r3, comm_columns=['leiden_0.01', 'leiden_0.1', 'leiden_1'], tag='embryo3_2', epochs=20, verbose = False, classifier='GraphSAGE')
# os.makedirs('models/seqfish_mouse_embryo/embryo_all')
# sg_obj.save_model_as_folder('models/seqfish_mouse_embryo/embryo_all')
# # adata_r3.var.to_csv('int_data/seqfish_mouse_embryo/gene_dt_embryo3_2.txt', sep='\t')
# # bsub -o .logs/sagenet/seqfish_mouse_embryo  -q gpu -gpu "num=1:gmem=20000" "python3 code/experiments/run_sagenet.py"



from torch.nn import Softmax
    
reference_embeddings = adata_r1.obsm['prob_concatenated']
kdtree_r1 = cKDTree(reference_embeddings)
target_embeddings = adata_q.obsm['prob_concatenated']
k_neighbors = 10  # You can adjust this as needed
distances, indices = kdtree_r1.query(target_embeddings, k=k_neighbors)


m = Softmax(dim=1)
probs = m(-torch.tensor(distances))
dist=torch.distributions.categorical.Categorical(probs=probs)
idx = dist.sample().numpy()
indices = indices[np.arange(len(indices)), idx]
x_coordinates = adata_r1.obs['x'][indices]
y_coordinates = adata_r1.obs['y'][indices]
coordinates = np.column_stack((x_coordinates, y_coordinates))
adata_q.obsm['spatial'] = coordinates



sc.pl.spatial(
    adata_q,
    color='cell_type',
    palette=celltype_colours, # Color cells based on 'cell_type'
    # color_map=cell_type_color_map,  # Use the custom color map
    # library_id='r1_mapping',  # Use 'r1_mapping' coordinates
    title='Spatial Plot with Cell Type Coloring',
    save='_sagenet_q3.pdf',
    spot_size=.1
)

sc.pl.spatial(
    adata_r1,
    color='cell_type',
    palette=celltype_colours, # Color cells based on 'cell_type'
    # color_map=cell_type_color_map,  # Use the custom color map
    # library_id='r1_mapping',  # Use 'r1_mapping' coordinates
    title='Spatial Plot with Cell Type Coloring',
    save='_adata_r1_true.pdf',
    spot_size=.1
)





celltype_colours = {
  "Epiblast" : "#635547",
  "Primitive Streak" : "#DABE99",
  "Caudal epiblast" : "#9e6762",
  "PGC" : "#FACB12",
  "Anterior Primitive Streak" : "#c19f70",
  "Notochord" : "#0F4A9C",
  "Def. endoderm" : "#F397C0",
  "Definitive endoderm" : "#F397C0",
  "Gut" : "#EF5A9D",
  "Gut tube" : "#EF5A9D",
  "Nascent mesoderm" : "#C594BF",
  "Mixed mesoderm" : "#DFCDE4",
  "Intermediate mesoderm" : "#139992",
  "Caudal Mesoderm" : "#3F84AA",
  "Paraxial mesoderm" : "#8DB5CE",
  "Somitic mesoderm" : "#005579",
  "Pharyngeal mesoderm" : "#C9EBFB",
  "Splanchnic mesoderm" : "#C9EBFB",
  "Cardiomyocytes" : "#B51D8D",
  "Allantois" : "#532C8A",
  "ExE mesoderm" : "#8870ad",
  "Lateral plate mesoderm" : "#8870ad",
  "Mesenchyme" : "#cc7818",
  "Mixed mesenchymal mesoderm" : "#cc7818",
  "Haematoendothelial progenitors" : "#FBBE92",
  "Endothelium" : "#ff891c",
  "Blood progenitors 1" : "#f9decf",
  "Blood progenitors 2" : "#c9a997",
  "Erythroid1" : "#C72228",
  "Erythroid2" : "#f79083",
  "Erythroid3" : "#EF4E22",
  "Erythroid" : "#f79083",
  "Blood progenitors" : "#f9decf",
  "NMP" : "#8EC792",
  "Rostral neurectoderm" : "#65A83E",
  "Caudal neurectoderm" : "#354E23",
  "Neural crest" : "#C3C388",
  "Forebrain/Midbrain/Hindbrain" : "#647a4f",
  "Spinal cord" : "#CDE088",
  "Surface ectoderm" : "#f7f79e",
  "Visceral endoderm" : "#F6BFCB",
  "ExE endoderm" : "#7F6874",
  "ExE ectoderm" : "#989898",
  "Parietal endoderm" : "#1A1A1A",
  "Unknown" : "#FFFFFF",
  "Low quality" : "#e6e6e6",
  # somitic and paraxial types
  # colour from T chimera paper Guibentif et al Developmental Cell 2021
  "Cranial mesoderm" : "#77441B",
  "Anterior somitic tissues" : "#F90026",
  "Sclerotome" : "#A10037",
  "Dermomyotome" : "#DA5921",
  "Posterior somitic tissues" : "#E1C239",
  "Presomitic mesoderm" : "#9DD84A"
}


def plot_confusion_matrix(confusion_matrix, labels, output_file):
    row_normalized_cm = confusion_matrix.astype('float') / confusion_matrix.sum(axis=1)[:, np.newaxis]
    plt.figure(figsize=(15, 12))
    sns.heatmap(confusion_matrix, annot=False, cmap="Blues", xticklabels=labels, yticklabels=labels)
    plt.xlabel('Predicted Labels')
    plt.ylabel('True Labels')
    plt.title('Confusion Matrix')
    plt.savefig(output_file, format="pdf")  # Save the plot to a PDF file

import tangram as tg
tg.pp_adatas(adata_q, adata_r1, genes=None)
ad_map = tg.map_cells_to_space(adata_q, adata_r1, device=device)
dist_np = ad_map.X
mapped_spots = np.argmax(dist_np, axis=1)
x_coordinates = adata_r1.obs['x'][mapped_spots]
y_coordinates = adata_r1.obs['y'][mapped_spots]
coordinates = np.column_stack((x_coordinates, y_coordinates))
coordinates = np.column_stack((x_coordinates, y_coordinates))
adata_q.obsm['spatial'] = coordinates


sc.pl.spatial(
    adata_r2,
    color='cell_type',
    palette=celltype_colours,# Color cells based on 'cell_type'
    # color_map=cell_type_color_map,  # Use the custom color map
    # library_id='r1_mapping',  # Use 'r1_mapping' coordinates
    title='Spatial Plot with Cell Type Coloring',
    save='_adata_r2.pdf',
    spot_size=.1
)


import novosparc

tissue = novosparc.cm.Tissue(dataset=adata_q, locations=adata_r1.obs[['x', 'y']])
tissue.setup_reconstruction(atlas_matrix=adata_r1.X)
tissue.reconstruct(alpha_linear=0.5, epsilon=5e-3)

gw = tissue.gw
ngw = (gw.T / gw.sum(1)).T

dist_np = ad_map.X
mapped_spots = ngw.argmax(1)
x_coordinates = adata_r1.obs['x'][mapped_spots]
y_coordinates = adata_r1.obs['y'][mapped_spots]
coordinates = np.column_stack((x_coordinates, y_coordinates))
coordinates = np.column_stack((x_coordinates, y_coordinates))
adata_q.obsm['spatial'] = coordinates


sc.pl.spatial(
    adata_q,
    color='cell_type',
    palette=celltype_colours,# Color cells based on 'cell_type'
    # color_map=cell_type_color_map,  # Use the custom color map
    # library_id='r1_mapping',  # Use 'r1_mapping' coordinates
    title='Spatial Plot with Cell Type Coloring',
    save='_novosparc_q.pdf',
    spot_size=.1
)