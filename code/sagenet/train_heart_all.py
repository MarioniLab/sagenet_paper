import sagenet as sg
import scanpy as sc
import squidpy as sq
import anndata as ad
import random
import numpy as np
random.seed(10)


def grid_array(df, n = 2):
	min_x = np.min(df['x'])
	max_x = np.max(df['x'])
	to_ret_x = np.ones(len(df))
	for i in range(1, n):
		qr = min_x + i*((max_x - min_x)/n)
		to_ret_x[np.where(df['x'].values > qr)] = i + 1
	min_y = np.min(df['y'])
	max_y = np.max(df['y'])
	to_ret_y = np.ones(len(df))
	for i in range(1, n):
		qr = min_y + i*((max_y - min_y)/n)
		to_ret_y[np.where(df['y'].values > qr)] = i + 1
	temp = to_ret_x * (to_ret_y**3)
	df['grid_' + str(n)] = temp
	return df

def grid_adata(adata, n=[2, 3, 4]):
	df  = adata.obs
	for i in n:
		tag = 'grid_' + str(i)
		df_temp = df.groupby('sample').apply(grid_array, n = i)
		temp = np.searchsorted(np.unique(df_temp[tag].values), df_temp[tag].values)
		sg.utils.save_adata(adata, attr='obs', key=tag, data=temp)

adata_r = sg.datasets.DHH_visium_()
grid_adata(adata_r, [2, 3, 4])
print('1')
sg.utils.glasso(adata_r, [0.5, 0.75, 1])
print('2')
import torch
if torch.cuda.is_available():  
	dev = "cuda:0" 
else:  
	dev = "cpu"  
device = torch.device(dev)
print(device)

sg_obj = sg.sage.sage(device=device)
sg_obj.add_ref(adata_r, comm_columns=['grid_2', 'grid_3', 'grid_4'], tag='visium_ref', epochs=20, verbose = False)
print('3')
sg_obj.save_model_as_folder('models/visium')
print('4')
# bsub -o .logs/visium  -q production -n 2 -M 8000 -R rusage[mem=8000] "python3 scPotter/train_heart_all.py"