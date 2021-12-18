import sagenet as sg
import scanpy as sc
import squidpy as sq
import anndata as ad
import random
import torch
import re

random.seed(1996)

if torch.cuda.is_available():  
  dev = "cuda:0" 
else:  
  dev = "cpu"  
device = torch.device(dev)
print(device)

sg_obj = sg.sage.sage()

sg_obj.load_model_as_folder('models/seqFISH')
adata_q1 = sg.datasets.seqFISH1()
adata_q2 = sg.datasets.seqFISH2()
adata_q3 = sg.datasets.seqFISH3()
adata_q4 = sg.datasets.MGA_scRNAseq()
sc.pp.subsample(adata_q1, fraction=0.25)
sc.pp.subsample(adata_q2, fraction=0.01)
sc.pp.subsample(adata_q3, fraction=0.25)
adata_q = ad.concat([adata_q1, adata_q2, adata_q3], join="inner")
adata_q.obs['class_'] = 0
sg_obj.map_query(adata_q)
adata_q.write('int_data/test_query.h5ad')