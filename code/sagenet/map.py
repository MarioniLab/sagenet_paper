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


sg_obj = sg.sage.sage(device=device)

sg_obj.load_model_as_folder('models/seqfish_mouse_embryo/embryo1_2')
adata_q = sc.read_h5ad('data_tidy/seqfish_mouse_embryo/query_data.h5ad')
adata_q.obs['class_'] = 0
sg_obj.map_query(adata_q)
adata_q.write('output/seqfish_mosue_embryo/sagenet/preds_embryo1_2_query.h5ad')