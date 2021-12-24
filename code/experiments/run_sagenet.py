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
sg_obj.load_model_as_folder(models_path)

print('2')
path_q = os.path.join(args.i, args.tag, args.tag_query) + '.h5ad'
adata_q = sc.read_h5ad(path_q)
adata_q.obs['class_'] = 0
print('3')
sg_obj.map_query(adata_q)
print('4')
preds_f = "_".join(['preds', args.tag_ref, args.tag_query]) + ".h5ad"
preds_f = os.path.join(args.oo, preds_f)
print('5')
adata_q.write(filename=preds_f)


