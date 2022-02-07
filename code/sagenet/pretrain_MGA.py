import sagenet as sg
import scanpy as sc
import squidpy as sq
import anndata as ad
import random
import torch
from sagenet.utils import glasso
import numpy as np
import pandas as pd

random.seed(1996)


adata_r1 = sg.MGA_data.seqFISH1_1()
adata_r2 = sg.MGA_data.seqFISH2_1()
adata_r3 = sg.MGA_data.seqFISH3_1()


glasso(adata_r1, [0.25, 0.5])
adata_r1.obsm['spatial'] = np.array(adata_r1.obs[['x','y']])
sq.gr.spatial_neighbors(adata_r1, coord_type="generic")
sc.tl.leiden(adata_r1, resolution=.01, random_state=0, key_added='leiden_0.01', adjacency=adata_r1.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r1, resolution=.05, random_state=0, key_added='leiden_0.05', adjacency=adata_r1.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r1, resolution=.1, random_state=0, key_added='leiden_0.1', adjacency=adata_r1.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r1, resolution=.5, random_state=0, key_added='leiden_0.5', adjacency=adata_r1.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r1, resolution=1, random_state=0, key_added='leiden_1', adjacency=adata_r1.obsp["spatial_connectivities"])
adata_r1.obs.to_csv('int_data/seqfish_mouse_embryo/col_dt_embryo1_2.txt', sep='\t')

glasso(adata_r2, [0.25, 0.5])
adata_r2.obsm['spatial'] = np.array(adata_r2.obs[['x','y']])
sq.gr.spatial_neighbors(adata_r2, coord_type="generic")
sc.tl.leiden(adata_r2, resolution=.01, random_state=0, key_added='leiden_0.01', adjacency=adata_r2.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r2, resolution=.05, random_state=0, key_added='leiden_0.05', adjacency=adata_r2.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r2, resolution=.1, random_state=0, key_added='leiden_0.1', adjacency=adata_r2.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r2, resolution=.5, random_state=0, key_added='leiden_0.5', adjacency=adata_r2.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r2, resolution=1, random_state=0, key_added='leiden_1', adjacency=adata_r2.obsp["spatial_connectivities"])
adata_r2.obs.to_csv('int_data/seqfish_mouse_embryo/col_dt_embryo2_2.txt', sep='\t')

glasso(adata_r3, [0.25, 0.5])
adata_r3.obsm['spatial'] = np.array(adata_r3.obs[['x','y']])
sq.gr.spatial_neighbors(adata_r3, coord_type="generic")
sc.tl.leiden(adata_r3, resolution=.01, random_state=0, key_added='leiden_0.01', adjacency=adata_r3.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r3, resolution=.05, random_state=0, key_added='leiden_0.05', adjacency=adata_r3.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r3, resolution=.1, random_state=0, key_added='leiden_0.1', adjacency=adata_r3.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r3, resolution=.5, random_state=0, key_added='leiden_0.5', adjacency=adata_r3.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r3, resolution=1, random_state=0, key_added='leiden_1', adjacency=adata_r3.obsp["spatial_connectivities"])
adata_r3.obs.to_csv('int_data/seqfish_mouse_embryo/col_dt_embryo3_2.txt', sep='\t')

if torch.cuda.is_available():  
  dev = "cuda:0"  
else:  
  dev = "cpu"  
device = torch.device(dev)
print(device)

sg_obj = sg.sage.sage(device=device)
sg_obj.add_ref(adata_r1, comm_columns=['leiden_0.01', 'leiden_0.05', 'leiden_0.1', 'leiden_0.5', 'leiden_1'], tag='embryo1_2', epochs=20, verbose = True)

adata_r1.var.to_csv('int_data/seqfish_mouse_embryo/gene_dt_embryo1_2.txt', sep='\t')
sg_obj.save_model_as_folder('models/seqfish_mouse_embryo/embryo1_2')
sg_obj = sg.sage.sage(device=device)
sg_obj.add_ref(adata_r2, comm_columns=['leiden_0.01', 'leiden_0.05', 'leiden_0.1', 'leiden_0.5', 'leiden_1'], tag='embryo2_2', epochs=20, verbose = False)
sg_obj.save_model_as_folder('models/seqfish_mouse_embryo/embryo2_2')
adata_r2.var.to_csv('int_data/seqfish_mouse_embryo/gene_dt_embryo2_2.txt', sep='\t')
sg_obj = sg.sage.sage(device=device)
sg_obj.add_ref(adata_r3, comm_columns=['leiden_0.01', 'leiden_0.05', 'leiden_0.1', 'leiden_0.5', 'leiden_1'], tag='embryo3_2', epochs=20, verbose = False)
sg_obj.save_model_as_folder('models/seqfish_mouse_embryo/embryo3_2')
adata_r3.var.to_csv('int_data/seqfish_mouse_embryo/gene_dt_embryo3_2.txt', sep='\t')
# bsub -o .logs/sagenet/seqfish_mouse_embryo  -q gpu -gpu "num=1:gmem=20000" "python3 code/experiments/run_sagenet.py"