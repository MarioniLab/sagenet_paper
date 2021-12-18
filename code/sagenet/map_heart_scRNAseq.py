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

sg_obj.load_model_as_folder('models/visium')
adata_q = sg.datasets.DHH_scRNAseq()
adata_q.obs['class_'] = 0
sg_obj.map_query(adata_q)
adata_q.write('int_data/mapped_DHH_scRNAseq.h5ad')