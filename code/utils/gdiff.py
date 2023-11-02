import os
# import sagenet as sg
import scanpy as sc
import squidpy as sq
import anndata as ad
import random
import anndata as ad 
import re 
random.seed(10)
from scipy import sparse
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import pairwise_distances

from gglasso.helper.data_generation import generate_precision_matrix, group_power_network, sample_covariance_matrix
from gglasso.problem import glasso_problem
from gglasso.helper.basic_linalg import adjacency_matrix
import numpy as np


from sklearn.covariance import empirical_covariance
from sklearn.metrics import *
from sklearn.covariance import GraphicalLassoCV, graphical_lasso, GraphicalLasso
from sklearn.preprocessing import StandardScaler
from scipy import sparse
from sklearn.metrics import pairwise_distances as pd

from squidpy.datasets._utils import AMetadata
# from sagenet.utils import save_adata


from numpy import zeros, newaxis

import leidenalg as la

from sklearn.metrics import jaccard_score
import igraph as ig


from copy import copy 
from sklearn.covariance import GraphicalLasso


sc_adata = AMetadata(
		name="scRNAseq",
		doc_header="",
		# shape=(270876, 43),
		url="https://figshare.com/ndownloader/files/31767704",
).download()
sc.pp.subsample(sc_adata, 0.5)
sp_adata_1 = AMetadata(
    name="seqFISH1_1",
    doc_header="",
    # shape=(270876, 43),
    url="https://figshare.com/ndownloader/files/31716029",
).download()
sp_adata_2 = AMetadata(
    name="seqFISH3_1",
    doc_header="",
    # shape=(270876, 43),
    url="https://figshare.com/ndownloader/files/31716089",
).download()





def glasso(adata, alphas):
	N = adata.shape[1]
	scaler = StandardScaler()
	data = scaler.fit_transform(adata.X)
	S    = empirical_covariance(data)
	P    = glasso_problem(S, N, latent = False, do_scaling = True)
	# lambda1_range = np.logspace(-0.1, -1, 10)
	lambda1_range = np.logspace(-1000,-10,10)
	mu1_range = np.logspace(-1000,-10,10)
	modelselect_params = {'lambda1_range': lambda1_range, 'mu1_range': mu1_range}
	P.model_selection(modelselect_params = modelselect_params, method = 'eBIC', gamma = 0.5, tol=1e-7)
	sol = P.solution.precision_
	P.solution.calc_adjacency(t = 1e-6)
	# save_adata(adata, attr='varm', key='adj', data=sparse.csr_matrix(P.solution.precision_))
	return P





sp_adj_1 = GraphicalLasso(max_iter=1000, alpha=0.0001).fit(sp_adata_1.X).precision_
np.fill_diagonal(sp_adj_1, 0)
sp_adj_2 = GraphicalLasso(max_iter=1000, alpha=0.0001).fit(sp_adata_2.X).precision_
np.fill_diagonal(sp_adj_2, 0)
sc_adj = GraphicalLasso(max_iter=1000, alpha=0.01).fit(sc_adata.X).precision_
np.fill_diagonal(sc_adj, 0)
# sp_P_1 = glasso(sp_adata_1, [1, 2, 5])
# sp_P_2 = glasso(sp_adata_2, [0.001, 0.1, 0.01])
# sc_P = glasso(sc_adata, [0.001, 0.1, 0.01])

# sp_adj_1 = sp_P_1.solution.precision_
# sp_adj_2 = sp_P_2.solution.precision_
# sc_adj   = sc_P.solution.precision_





def gdiff_adj(adj_1, adj_2, ord='fro'):
	return np.linalg.norm((adj_1 - adj_2), ord)




import scanpy as sc
def partition_graph(adj):
	g = sc._utils.get_igraph_from_adjacency(adj)
	partition = la.find_partition(g, la.ModularityVertexPartition)
	return partition.membership


def plot_graph(adj, target='output/test_partitions.pdf'):
	g = sc._utils.get_igraph_from_adjacency(adj)
	partition = la.find_partition(g, la.ModularityVertexPartition)
	ig.plot(partition, target = target)
# plot_graph(sc_adj)

# plot_graph(sp_adj_2, target='output/test_partitions.pdf')



def gdiff_jaccard(adj_1, adj_2):
	part_1 = partition_graph(adj_1)
	part_2 = partition_graph(adj_2)
	return 1 - jaccard_score(part_1, part_2, average='micro')




gdiff_adj(sp_adj_1, sp_adj_2)
gdiff_jaccard(sp_adj_2, sp_adj_1)

gene_metric = pd(sp_adj_1, sp_adj_2).diagonal()


a = np.dot(sp_adata_1.X[0][:, None], sp_adata_1.X[0][None, :]) * sp_adj_1




xx = copy(sp_adata_1.X)

sc_corr = np.corrcoef(sc_adata.X.T)

xx_1 = xx[:, :, newaxis]
xx_2 = xx[:, newaxis, :]
gg = np.matmul(xx_1, xx_2)
res = np.multiply(gg, sc_corr).sum((1,2))


sp_adata_1.obs['cell_corr_qual'] = res

sp_adata_1.obsm['spatial'] = np.array(sp_adata_1.obs[['x','y']])
sp_adata_1.uns['spatial'] = np.array(sp_adata_1.obs[['x','y']])
sc.pl.spatial(sp_adata_1, color='cell_corr_qual', spot_size=0.08)
import torch
if torch.cuda.is_available():  
  dev = "cuda:0" 
else:  
  dev = "cpu"  

device = torch.device(dev)
print(device)

sg_obj = sg.sage.sage(device=device)



adata_dic = {}


# filter = df['Event Name'].str.contains(patternDel)

for adata_f in os.listdir(data_path):
	tag = re.sub('\\.h5ad', '', adata_f) 
	print(tag)
	adata_dic[tag] = sc.read(os.path.join(data_path, adata_f))
	normal_pat = '^Ctr|^Normal|^unassigned'
	normal_spots = adata_dic[tag].obs['nodule'].str.contains(normal_pat)
	adata_dic[tag] = adata_dic[tag][~normal_spots,:]
	print(adata_dic[tag].shape)
	adata_dic[tag].obs['section'] = tag
	adata_dic[tag] = adata_dic[tag][adata_dic[tag].X.sum(1) != 0, :]
	sc.pp.normalize_total(adata_dic[tag], target_sum=1e4)
	sc.pp.log1p(adata_dic[tag])
	adata_dic[tag].varm['adj'] = sparse.csr_matrix(adata_dic[tag].uns['adj'])
	le = LabelEncoder()
	adata_dic[tag].obs['nodule_encoded'] = le.fit_transform(adata_dic[tag].obs['nodule'])
	sg_obj.add_ref(adata_dic[tag], comm_columns=['nodule_encoded'], tag=tag, epochs=100, verbose = False)


