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


adata_q  = sg.MGA_data.scRNAseq()
adata_q.obs['x'] = adata_q.obs['y'] = np.nan
adata_q.obs['embryo'] = 'scRNAseq'
adata_q = ad.concat([adata_q1, adata_q2, adata_q3, adata_q4, adata_q5, adata_q6, adata_q], join="inner")
adata_q.obs = adata_q.obs.drop('class_', axis=1)  
adata_q.write('data_tidy/seqfish_mouse_embryo/query_scRNAseq.h5ad')