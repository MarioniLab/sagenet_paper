import argparse
import json
import os, sys
import sys
import torch
import tangram as tg
import scanpy as sc
import pandas as pd
import sagenet as sg
import scanpy as sc
import squidpy as sq
import anndata as ad
import random
import torch
import re
import numpy as np
from sklearn.metrics import pairwise_distances

currentdir = os.path.dirname(os.path.abspath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 



parser = argparse.ArgumentParser()
parser.add_argument('--tag', type=str, default='Syn', help='Input dataset. Can either be the name of the synthetic dataset generation method, \
                                              or the path to a real dataset.')
parser.add_argument('--tag_ref', type=str, default='train', help='Input dataset. Can either be the name of the synthetic dataset generation method, \
                                              or the path to a real dataset.')
parser.add_argument('--tag_query', type=str, default='test', help='Input dataset. Can either be the name of the synthetic dataset generation method, \
                                              or the path to a real dataset.')
parser.add_argument('-i', type=str,default='../../data_tidy', help='Input dataset. Can either be the name of the synthetic dataset generation method, \
                                              or the path to a real dataset.')
parser.add_argument('--oo', type=str,default='../../output', help='Output directory.')

args = parser.parse_args()
print(args)


if not os.path.exists(args.oo):
  os.makedirs(args.oo)
args.oo = os.path.join(args.oo, args.tag)
if not os.path.exists(args.oo):
  os.makedirs(args.oo)
args.oo = os.path.join(args.oo, 'sagenet')
if not os.path.exists(args.oo):
  os.makedirs(args.oo)

if torch.cuda.is_available():  
  dev = "cuda:0" 
else:  
  dev = "cpu"  


device = torch.device(dev)
print(device)


# device = torch.device('cpu')

print('1')
models_path = os.path.join(os.path.join('models', args.tag), args.tag_ref)

sg_obj = sg.sage.sage(device)

models_path = 'models/seqfish_mouse_embryo/embryo1_2'
sg_obj.load_model_as_folder(models_path)
adata_q_1 = sc.read('data_tidy/seqfish_mouse_embryo/embryo1_2.h5ad')
adata_q_2 = sc.read('data_tidy/seqfish_mouse_embryo/query_data.h5ad')

print('2')
path_q = os.path.join(args.i, args.tag, args.tag_query) + '.h5ad'
adata_q = sc.read_h5ad(path_q)
adata_q.obs['class_'] = 0
print('3')
sg_obj.map_query(adata_q_1, save_prob=True)
sg_obj.map_query(adata_q_2, save_prob=True)
print('4')
preds_f = "_".join(['preds', args.tag_ref, args.tag_query]) + ".h5ad"
preds_f = os.path.join(args.oo, preds_f)
print('5')
adata_q.write(filename=preds_f)



sg_obj.add_ref(embryo1_2, tag='embryo1_2', comm_columns=['leiden_1', 'leiden_0.05', 'leiden_0.1', 'leiden_0.5', 'leiden_1'], classifier='GraphSAGE')


sg_obj.save_model_as_folder('models/seqfish_mouse_embryo/embryo1_2')


adata_q.obsm['prob_combined'] = np.column_stack((
    np.nan_to_num(adata_q.obsm['prob_embryo1_2_leiden_1'],0), 
    np.nan_to_num(adata_q.obsm['prob_embryo1_2_leiden_0.1'],0),
    np.nan_to_num(adata_q.obsm['prob_embryo1_2_leiden_0.01'],0),
    np.nan_to_num(adata_q.obsm['prob_embryo1_2_leiden_0.05'], 0),
    np.nan_to_num(adata_q.obsm['prob_embryo1_2_leiden_0.5'], 0)

  )) 


base_dists = np.zeros((adata_q.obsm['prob_combined'].shape[0], adata_q.obsm['prob_combined'].shape[0]))
prob_list = ['prob_embryo1_2_leiden_1', 'prob_embryo1_2_leiden_0.1', 'prob_embryo1_2_leiden_0.05', 'prob_embryo1_2_leiden_0.5']
for prob in prob_list:
  print(prob)
  pd = pairwise_distances(adata_q.obsm[prob])
  del adata_q.obsm[prob]
  pd /= np.linalg.norm(pd, 'fro')
  base_dists += pd

adata_q.obsp['sagenet_dist'] = base_dists





    # def add_ref(self, 
    #     adata, 
    #     tag = None,
    #     comm_columns = 'class_',
    #     classifier   = 'TransformerConv',
    #     num_workers  = 0,
    #     batch_size   = 32,
    #     epochs       = 10,
    #     n_genes      = 10,
    #     verbose      = False):




