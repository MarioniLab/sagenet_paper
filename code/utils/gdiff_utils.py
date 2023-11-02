import glob
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from copy import copy
from numpy import newaxis
import numpy as np




DATA_DIR = "/Users/b450-admin/Desktop/spverse/txsim_out" # ADJUST
DATASET = "10x_Xenium_Janesick_breast-cancer_2023_spverse_rep1"
SUBSETS = {"2k":"_crop_2000px", "4k":"_crop_4000px", "8k":"_crop_8000px"}

spatial_h5ads = [
    "counts_baysor-0_area-0.h5ad",
    "counts_binning-0_basic-0_area-0.h5ad",
    "counts_cellpose-0_basic-0_area-0.h5ad",
    "counts_custom-0_basic-0_area-0.h5ad",
    "counts_custom-1_basic-0_area-0.h5ad",
    "counts_watershed-0_basic-0_area-0.h5ad",
]
names = [n.split("_",1)[1].rsplit("_",1)[0] for n in spatial_h5ads]
names = [n.split("_",1)[0] if "basic" in n else n for n in names]

adata_sc = sc.read(Path(DATA_DIR, DATASET + SUBSETS["8k"], "sc_normalized.h5ad"))
adata_sp = sc.read(Path(DATA_DIR, DATASET + SUBSETS["8k"], "counts_custom-0_basic-0_area-0.h5ad"))
adata_sc = adata_sc[:, adata_sp.var_names]
adatas_sp = {name: sc.read(Path(DATA_DIR, DATASET + SUBSETS["8k"], h5ad)) for name, h5ad in zip(names, spatial_h5ads)}









sc_corr = np.corrcoef(adata_sc.layers['raw'].toarray().T)


# temp = sc.pp.normalize_total(
#     adata_sp, target_sum=1e4, exclude_highly_expressed=True,
#     max_fraction=0.2, inplace=False
# )['X']
temp = copy(adata_sp.layers['raw'])

temp_1 = temp[:, :, newaxis]
temp_2 = temp[:, newaxis, :]

gg = np.matmul(temp_1, temp_2)
# del temp, temp_1, temp_2
res = np.multiply(gg, sc_corr).sum((1,2))




adata_sp.obs['sum'] = adata_sp.X.sum(1)

adata_sp.obs['cell_corr_qual'] = np.log(res) 


spots_df = pd.DataFrame(adata_sp.uns['spots']).dropna()
spots_df['cell_id'] = spots_df['cell'].astype('int')
spots_df = spots_df.groupby(['cell_id']).mean(['x', 'y'])
adata_sp.obs = spots_df.merge(adata_sp.obs, on='cell_id')  

adata_sp.obsm['spatial'] = np.array(adata_sp.obs[['x','y']])
# adata_sp.uns['spatial'] = np.array(adata_sp.obs[['x','y']])


sc.pl.spatial(adata_sp, color=['cell_corr_qual'], spot_size=50)





