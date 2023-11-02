import sagenet as sg
import scanpy as sc
import squidpy as sq
import anndata as ad
import random
random.seed(10)




from copy import copy
adata_r1 = sg.MGA_data.seqFISH1_1()
adata_r2 = sg.MGA_data.seqFISH2_1()
adata_r3 = sg.MGA_data.seqFISH3_1()
adata_q1 = copy(adata_r1)
adata_q1.obs['section'] = 'embryo1'
adata_q2 = copy(adata_r2)
adata_q2.obs['section'] = 'embryo2'
adata_q3 = copy(adata_r3)
adata_q3.obs['section'] = 'embryo3'
adata_q4 = sg.MGA_data.scRNAseq()
adata_q4.obs['section'] = 'scRNAseq'
sc.pp.subsample(adata_q1, fraction=0.1)
sc.pp.subsample(adata_q2, fraction=0.1)
sc.pp.subsample(adata_q3, fraction=0.1)
sc.pp.subsample(adata_q4, fraction=0.25)
adata_q = ad.concat([adata_q1, adata_q2, adata_q3, adata_q4], join="inner")
del adata_q1 
del adata_q2 
del adata_q3 
del adata_q4


from sagenet.utils import glasso
import numpy as np



glasso(adata_r1, [0.5, 0.75, 1])
adata_r1.obsm['spatial'] = np.array(adata_r1.obs[['x','y']])
sq.gr.spatial_neighbors(adata_r1, coord_type="generic")
sc.tl.leiden(adata_r1, resolution=.01, random_state=0, key_added='leiden_0.01', adjacency=adata_r1.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r1, resolution=.05, random_state=0, key_added='leiden_0.05', adjacency=adata_r1.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r1, resolution=.1, random_state=0, key_added='leiden_0.1', adjacency=adata_r1.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r1, resolution=.5, random_state=0, key_added='leiden_0.5', adjacency=adata_r1.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r1, resolution=1, random_state=0, key_added='leiden_1', adjacency=adata_r1.obsp["spatial_connectivities"])
glasso(adata_r2, [0.5, 0.75, 1])
adata_r2.obsm['spatial'] = np.array(adata_r2.obs[['x','y']])
sq.gr.spatial_neighbors(adata_r2, coord_type="generic")
sc.tl.leiden(adata_r2, resolution=.01, random_state=0, key_added='leiden_0.01', adjacency=adata_r2.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r2, resolution=.05, random_state=0, key_added='leiden_0.05', adjacency=adata_r2.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r2, resolution=.1, random_state=0, key_added='leiden_0.1', adjacency=adata_r2.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r2, resolution=.5, random_state=0, key_added='leiden_0.5', adjacency=adata_r2.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r2, resolution=1, random_state=0, key_added='leiden_1', adjacency=adata_r2.obsp["spatial_connectivities"])
glasso(adata_r3, [0.5, 0.75, 1])
adata_r3.obsm['spatial'] = np.array(adata_r3.obs[['x','y']])
sq.gr.spatial_neighbors(adata_r3, coord_type="generic")
sc.tl.leiden(adata_r3, resolution=.01, random_state=0, key_added='leiden_0.01', adjacency=adata_r3.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r3, resolution=.05, random_state=0, key_added='leiden_0.05', adjacency=adata_r3.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r3, resolution=.1, random_state=0, key_added='leiden_0.1', adjacency=adata_r3.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r3, resolution=.5, random_state=0, key_added='leiden_0.5', adjacency=adata_r3.obsp["spatial_connectivities"])
sc.tl.leiden(adata_r3, resolution=1, random_state=0, key_added='leiden_1', adjacency=adata_r3.obsp["spatial_connectivities"])


import torch
if torch.cuda.is_available():  
  dev = "cuda:0" 
else:  
  dev = "cpu"  

device = torch.device(dev)
print(device)



sg_obj = sg.sage.sage(device=device)

sg_obj.add_ref(adata_r1, comm_columns=['leiden_0.01', 'leiden_0.05', 'leiden_0.1', 'leiden_0.5'], tag='seqFISH_ref1', epochs=15, verbose = False)
# sg_obj.add_ref(adata_r2, comm_columns=['leiden_0.01', 'leiden_0.05', 'leiden_0.1', 'leiden_0.5', 'leiden_1'], tag='seqFISH_ref2', epochs=15, verbose = False)
# sg_obj.add_ref(adata_r3, comm_columns=['leiden_0.01', 'leiden_0.05', 'leiden_0.1', 'leiden_0.5', 'leiden_1'], tag='seqFISH_ref3', epochs=15, verbose = False)


ind = np.argsort(-(adata_r1.var['seqFISH_ref1_entropy']+ adata_r2.var['seqFISH_ref2_entropy'] + adata_r3.var['seqFISH_ref3_entropy']))[0:12]
with rc_context({'figure.figsize': (4, 4)}):
  sc.pl.spatial(adata_r1, color=list(adata_r1.var_names[ind]), ncols=4, spot_size=0.03, legend_loc=None)



sg_obj.map_query(adata_q, save_prob=True)
dist_map = np.zeros((adata_q.shape[0], adata_q.shape[0]))

from sklearn.metrics import pairwise_distances

for prob_mat_tag in adata_q.obsm:
	prob_mat = adata_q.obsm[prob_mat_tag]
	dist_map += pairwise_distances(prob_mat)



dist_adata = ad.AnnData(dist_map, obs = adata_q.obs)
dist_adata = dist_adata[:, dist_adata.obs['section'] == 'embryo1']
dist_np = dist_adata.X
mapped_spots = np.argmin(dist_np, axis=1)
mapped_spots = adata_q.obs_names[mapped_spots]
adata_q.obs['x'] = np.random.normal(adata_r1[mapped_spots,:].obs['x'].values, 0.000001)
adata_q.obs['y'] = np.random.normal(adata_r1[mapped_spots,:].obs['y'].values, 0.000001)
adata_q.obsm['spatial'] = np.array(adata_q.obs[['x', 'y']])
sc.pl.spatial(adata_q, color='cell_type', palette=celltype_colours, frameon=False, spot_size=.1, save='_sagenet_spot_map.pdf')
sc.pl.spatial(adata_q, color='section', frameon=False, spot_size=.1, save='_spot_map_sections.pdf')



celltype_colours = {
  "Epiblast" : "#635547",
  "Primitive Streak" : "#DABE99",
  "Caudal epiblast" : "#9e6762",
  "PGC" : "#FACB12",
  "Anterior Primitive Streak" : "#c19f70",
  "Notochord" : "#0F4A9C",
  "Def. endoderm" : "#F397C0",
  "Definitive endoderm" : "#F397C0",
  "Gut" : "#EF5A9D",
  "Gut tube" : "#EF5A9D",
  "Nascent mesoderm" : "#C594BF",
  "Mixed mesoderm" : "#DFCDE4",
  "Intermediate mesoderm" : "#139992",
  "Caudal Mesoderm" : "#3F84AA",
  "Paraxial mesoderm" : "#8DB5CE",
  "Somitic mesoderm" : "#005579",
  "Pharyngeal mesoderm" : "#C9EBFB",
  "Splanchnic mesoderm" : "#C9EBFB",
  "Cardiomyocytes" : "#B51D8D",
  "Allantois" : "#532C8A",
  "ExE mesoderm" : "#8870ad",
  "Lateral plate mesoderm" : "#8870ad",
  "Mesenchyme" : "#cc7818",
  "Mixed mesenchymal mesoderm" : "#cc7818",
  "Haematoendothelial progenitors" : "#FBBE92",
  "Endothelium" : "#ff891c",
  "Blood progenitors 1" : "#f9decf",
  "Blood progenitors 2" : "#c9a997",
  "Erythroid1" : "#C72228",
  "Erythroid2" : "#f79083",
  "Erythroid3" : "#EF4E22",
  "Erythroid" : "#f79083",
  "Blood progenitors" : "#f9decf",
  "NMP" : "#8EC792",
  "Rostral neurectoderm" : "#65A83E",
  "Caudal neurectoderm" : "#354E23",
  "Neural crest" : "#C3C388",
  "Forebrain/Midbrain/Hindbrain" : "#647a4f",
  "Spinal cord" : "#CDE088",
  "Surface ectoderm" : "#f7f79e",
  "Visceral endoderm" : "#F6BFCB",
  "ExE endoderm" : "#7F6874",
  "ExE ectoderm" : "#989898",
  "Parietal endoderm" : "#1A1A1A",
  "Unknown" : "#FFFFFF",
  "Low quality" : "#e6e6e6",
  # somitic and paraxial types
  # colour from T chimera paper Guibentif et al Developmental Cell 2021
  "Cranial mesoderm" : "#77441B",
  "Anterior somitic tissues" : "#F90026",
  "Sclerotome" : "#A10037",
  "Dermomyotome" : "#DA5921",
  "Posterior somitic tissues" : "#E1C239",
  "Presomitic mesoderm" : "#9DD84A"
}




import argparse
import json
import os, sys
import sys
import torch
import tangram as tg
import scanpy as sc
import pandas as pd

tg.pp_adatas(adata_q1, adata_q, genes=None)
ad_out = tg.map_cells_to_space(
	adata_q1,
	adata_q,
	device=device
)


ad_q = ad_out.T
mapped_spots = np.argmax(ad_q.X, axis=1)
mapped_spots = ad_q.var_names[mapped_spots]
ad_q.obs['x'] = np.random.normal(adata_r1[mapped_spots,:].obs['x'].values, 0.000001)
ad_q.obs['y'] = np.random.normal(adata_r1[mapped_spots,:].obs['y'].values, 0.000001)
ad_q.obsm['spatial'] = np.array(ad_q.obs[['x', 'y']])
sc.pl.spatial(ad_q, color='cell_type', palette=celltype_colours, frameon=False, spot_size=.1, save='_tangram_spot_map.pdf')
sc.pl.spatial(ad_q, color='section', frameon=False, spot_size=.1, save='_tangram_spot_map_sections.pdf')





def calc_bins(preds):
  # Assign each prediction to a bin
  num_bins = 10
  bins = np.linspace(0.1, 1, num_bins)
  binned = np.digitize(preds, bins)

  # Save the accuracy, confidence and size of each bin
  bin_accs = np.zeros(num_bins)
  bin_confs = np.zeros(num_bins)
  bin_sizes = np.zeros(num_bins)

  for bin in range(num_bins):
    bin_sizes[bin] = len(preds[binned == bin])
    if bin_sizes[bin] > 0:
      bin_accs[bin] = (labels_oneh[binned==bin]).sum() / bin_sizes[bin]
      bin_confs[bin] = (preds[binned==bin]).sum() / bin_sizes[bin]

  return bins, binned, bin_accs, bin_confs, bin_sizes



import torch
from torch import nn, optim
from torch.nn import functional as F


class ModelWithTemperature(nn.Module):
    """
    A thin decorator, which wraps a model with temperature scaling
    model (nn.Module):
        A classification neural network
        NB: Output of the neural network should be the classification logits,
            NOT the softmax (or log softmax)!
    """
    def __init__(self, model):
        super(ModelWithTemperature, self).__init__()
        self.model = model
        self.temperature = nn.Parameter(torch.ones(1) * 1.5)

    def forward(self, input):
        logits = self.model(input)
        return self.temperature_scale(logits)

    def temperature_scale(self, logits):
        """
        Perform temperature scaling on logits
        """
        # Expand temperature to match the size of logits
        temperature = self.temperature.unsqueeze(1).expand(logits.size(0), logits.size(1))
        return logits / temperature

    # This function probably should live outside of this class, but whatever
    def set_temperature(self, valid_loader):
        """
        Tune the tempearature of the model (using the validation set).
        We're going to set it to optimize NLL.
        valid_loader (DataLoader): validation set loader
        """
        self.cuda()
        nll_criterion = nn.CrossEntropyLoss().cuda()
        ece_criterion = _ECELoss().cuda()
        # First: collect all the logits and labels for the validation set
        logits_list = []
        labels_list = []
        with torch.no_grad():
            for input, label in valid_loader:
                input = input.cuda()
                logits = self.model(input)
                logits_list.append(logits)
                labels_list.append(label)
            logits = torch.cat(logits_list).cuda()
            labels = torch.cat(labels_list).cuda()
        # Calculate NLL and ECE before temperature scaling
        before_temperature_nll = nll_criterion(logits, labels).item()
        before_temperature_ece = ece_criterion(logits, labels).item()
        print('Before temperature - NLL: %.3f, ECE: %.3f' % (before_temperature_nll, before_temperature_ece))
        # Next: optimize the temperature w.r.t. NLL
        optimizer = optim.LBFGS([self.temperature], lr=0.01, max_iter=50)
        def eval():
            optimizer.zero_grad()
            loss = nll_criterion(self.temperature_scale(logits), labels)
            loss.backward()
            return loss
        optimizer.step(eval)
        # Calculate NLL and ECE after temperature scaling
        after_temperature_nll = nll_criterion(self.temperature_scale(logits), labels).item()
        after_temperature_ece = ece_criterion(self.temperature_scale(logits), labels).item()
        print('Optimal temperature: %.3f' % self.temperature.item())
        print('After temperature - NLL: %.3f, ECE: %.3f' % (after_temperature_nll, after_temperature_ece))
        return self


class _ECELoss(nn.Module):
    """
    Calculates the Expected Calibration Error of a model.
    (This isn't necessary for temperature scaling, just a cool metric).

    The input to this loss is the logits of a model, NOT the softmax scores.

    This divides the confidence outputs into equally-sized interval bins.
    In each bin, we compute the confidence gap:

    bin_gap = | avg_confidence_in_bin - accuracy_in_bin |

    We then return a weighted average of the gaps, based on the number
    of samples in each bin

    See: Naeini, Mahdi Pakdaman, Gregory F. Cooper, and Milos Hauskrecht.
    "Obtaining Well Calibrated Probabilities Using Bayesian Binning." AAAI.
    2015.
    """
    def __init__(self, n_bins=15):
        """
        n_bins (int): number of confidence interval bins
        """
        super(_ECELoss, self).__init__()
        bin_boundaries = torch.linspace(0, 1, n_bins + 1)
        self.bin_lowers = bin_boundaries[:-1]
        self.bin_uppers = bin_boundaries[1:]

    def forward(self, logits, labels):
        softmaxes = F.softmax(logits, dim=1)
        confidences, predictions = torch.max(softmaxes, 1)
        accuracies = predictions.eq(labels)
        ece = torch.zeros(1, device=logits.device)
        for bin_lower, bin_upper in zip(self.bin_lowers, self.bin_uppers):
            # Calculated |confidence - accuracy| in each bin
            in_bin = confidences.gt(bin_lower.item()) * confidences.le(bin_upper.item())
            prop_in_bin = in_bin.float().mean()
            if prop_in_bin.item() > 0:
                accuracy_in_bin = accuracies[in_bin].float().mean()
                avg_confidence_in_bin = confidences[in_bin].mean()
                ece += torch.abs(avg_confidence_in_bin - accuracy_in_bin) * prop_in_bin
        return ece







