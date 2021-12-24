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


adata_q1 = sg.MGA_data.seqFISH1_1()
adata_q2 = sg.MGA_data.seqFISH2_1()
adata_q3 = sg.MGA_data.seqFISH3_1()
adata_q4 = sg.MGA_data.seqFISH1_2()
adata_q5 = sg.MGA_data.seqFISH2_2()
adata_q6 = sg.MGA_data.seqFISH3_2()
adata_q  = sg.MGA_data.scRNAseq()
sc.pp.subsample(adata_q1, fraction=0.2)
sc.pp.subsample(adata_q2, fraction=0.2)
sc.pp.subsample(adata_q3, fraction=0.2)
sc.pp.subsample(adata_q4, fraction=0.2)
sc.pp.subsample(adata_q5, fraction=0.2)
sc.pp.subsample(adata_q6, fraction=0.2)
sc.pp.subsample(adata_q, fraction=0.25)
adata_q.obs['x'] = adata_q.obs['y'] = np.nan
adata_q.obs['embryo'] = 'scRNAseq'
adata_q = ad.concat([adata_q1, adata_q2, adata_q3, adata_q4, adata_q5, adata_q6, adata_q], join="inner")
adata_q.obs = adata_q.obs.drop('class_', axis=1)  
adata_q.write('data_tidy/seqfish_mouse_embryo/query_data.h5ad')