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
	df = adata.obs
	for i in n:
		tag = 'grid_' + str(i)
		df_temp = df.groupby('sample').apply(grid_array, n = i)
		temp = np.searchsorted(np.unique(df_temp[tag].values), df_temp[tag].values)
		sg.utils.save_adata(adata, attr='obs', key=tag, data=temp)

adata_r = sg.DHH_data.ST()
grid_adata(adata_r, [2, 3, 4])
sg.utils.glasso(adata_r, [0.5, 0.75])
adata_r.obs['class_'] = 0
import torch
if torch.cuda.is_available():  
	dev = "cuda:0" 
else:  
	dev = "cpu"  
device = torch.device(dev)
print(device)

sg_obj = sg.sage.sage(device=device)
sg_obj.add_ref(adata_r, comm_columns=['grid2', 'grid_3', 'grid_4'], tag='ST_all', epochs=20, verbose = False)
adata_r.var.to_csv('int_data/ST_human_heart/ST_gene_dt.txt', sep='\t')
adata_r.obs.to_csv('int_data/ST_human_heart/col_dt_ST.txt', sep='\t')
sg_obj.save_model_as_folder('models/ST_human_heart')


for i in ['CN73_C2', 'CN73_D2', 'CN73_E1', 'CN73_E2', 'CN74_C1', 'CN74_D1', 'CN74_D2', 'CN74_E1', 'CN74_E2']:
	print(i)
	adata_rx = adata_r[adata_r.obs['sample'] != i]
	sg.utils.glasso(adata_rx, [0.5, 0.75])
	print(i)
	sg_obj = sg.sage.sage(device=device)
	sg_obj.add_ref(adata_rx, comm_columns = ['grid_2', 'grid_3', 'grid_4'], tag = i, epochs=10, verbose = False)
	adata_q = adata_r[adata_r.obs['sample'] == i].copy()
	sg_obj.map_query(adata_q)
	preds_f = "_".join(['preds', str(i + '_lo'), i]) + ".h5ad"
	preds_f = os.path.join('output/ST_human_heart/sagenet', preds_f)
	adata_q.write(filename=preds_f)
