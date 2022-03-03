import sagenet as sg
import scanpy as sc
import squidpy as sq
import anndata as ad
import random
import numpy as np

rna_ad = sc.read('data_tidy/ST_human_heart/scRNAseq.h5ad')
rna_ad.obs = rna_ad.obs[['cell_id']].assign(class__ = rna_ad.obs.celltype, sample = rna_ad.obs.experiment, dataset='scRNAseq')

st_ad  = sc.read('data_tidy/ST_human_heart/ST.h5ad')
st_ad.obs = st_ad.obs[['cell_id', 'class__', 'sample']].assign(dataset='ST', class__ = str(st_ad.obs.class__))

query_ad = st_ad.concatenate(rna_ad)

query_ad.write('data_tidy/ST_human_heart/query.h5ad')
