from gglasso.helper.data_generation import generate_precision_matrix, group_power_network, sample_covariance_matrix
from gglasso.problem import glasso_problem
from gglasso.helper.basic_linalg import adjacency_matrix


from sklearn.covariance import empirical_covariance
from sklearn.metrics import *
from sklearn.covariance import GraphicalLassoCV, graphical_lasso, GraphicalLasso
from sklearn.preprocessing import StandardScaler
from scipy import sparse

from sagenet.utils import save_adata


embryo1_2.obsm['spatial'] = embryo1_2.obs[['x', 'y']]
sq.gr.spatial_neighbors(embryo1_2)

sc.tl.leiden(embryo1_2, resolution=.01, random_state=0, key_added='leiden_0.01', adjacency=embryo1_2.obsp["spatial_connectivities"])
sc.tl.leiden(embryo1_2, resolution=.05, random_state=0, key_added='leiden_0.05', adjacency=embryo1_2.obsp["spatial_connectivities"])
sc.tl.leiden(embryo1_2, resolution=.1, random_state=0, key_added='leiden_0.1', adjacency=embryo1_2.obsp["spatial_connectivities"])
sc.tl.leiden(embryo1_2, resolution=.5, random_state=0, key_added='leiden_0.5', adjacency=embryo1_2.obsp["spatial_connectivities"])
sc.tl.leiden(embryo1_2, resolution=1, random_state=0, key_added='leiden_1', adjacency=embryo1_2.obsp["spatial_connectivities"])



embryo1_2.write('data_tidy/seqfish_mouse_embryo/embryo1_2.h5ad')



def glasso(adata):
    N = visium.shape[1]
    scaler = StandardScaler()
    data = scaler.fit_transform(visium.X)
    S    = empirical_covariance(data)
    P    = glasso_problem(S, N, latent = False, do_scaling = True)
    # lambda1_range = np.logspace(-0.1, -1, 10)
    lambda1_range = np.logspace(-10, -1,10)
    modelselect_params = {'lambda1_range': lambda1_range}
    P.model_selection(modelselect_params = modelselect_params, method = 'eBIC', gamma = 0.1, tol=1e-7)
    sol = P.solution.precision_
    P.solution.calc_adjacency(t = 1e-4)
    save_adata(adata, attr='varm', key='adj', data=sparse.csr_matrix(P.solution.precision_))
