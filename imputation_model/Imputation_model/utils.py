## from scGPT (https://github.com/bowang-lab/scGPT/blob/main/scgpt/model/model.py)
import torch
import torch.nn.functional as F
from torch import Tensor, nn
import scipy.sparse as sp
from scipy.sparse import csr_matrix
from typing import Dict, Mapping, Optional, Tuple, Any, Union
import numpy as np
import scanpy as sc

from .normalization import ref_query_norm, RNA_scale



def load_preprocess(ref_PATH, query_PATH, ref_subsample_num, log1p, scale_type):
    
    if ref_PATH.split('.')[-1] == 'h5ad':
        sc_RNA_ref = sc.read_h5ad(ref_PATH)
    elif ref_PATH.split('.')[-1] == 'csv':
        sc_RNA_ref = sc.read_csv(ref_PATH, delimiter=',', first_column_names=True, dtype='float32')

    if query_PATH.split('.')[-1] == 'h5ad':
        sc_RNA_query = sc.read_h5ad(query_PATH)
    elif query_PATH.split('.')[-1] == 'csv':
        sc_RNA_query = sc.read_csv(query_PATH, delimiter=',', first_column_names=True, dtype='float32')
        
    if ref_subsample_num != None:
        sc.pp.subsample(sc_RNA_ref, n_obs=ref_subsample_num, random_state=0, copy=False) 

    sc_RNA_ref.obs['total_counts'] = sc_RNA_ref.X.sum(axis=1)
    sc_RNA_query.obs['total_counts'] = sc_RNA_query.X.sum(axis=1)
    sc_RNA_ref.obs['dataset'] = 'ref'
    sc_RNA_query.obs['dataset'] = 'query'

    sc_RNA_ref.X = sp.csr_matrix(sc_RNA_ref.X)
    sc_RNA_query.X = sp.csr_matrix(sc_RNA_query.X)
    
    sc_RNA, sc_RNA_ref_normalized, sc_RNA_query_normalized, common_genes, pred_genes = ref_query_norm(sc_RNA_ref, sc_RNA_query)
    
    if log1p == True:
        sc.pp.log1p(sc_RNA)

    if scale_type != False:
        sc_RNA = RNA_scale(sc_RNA, common_genes, sc_RNA_ref_normalized.shape[0], scale_type)
        
    return sc_RNA, sc_RNA_ref_normalized, sc_RNA_query_normalized, common_genes, pred_genes



## adjacency matrix normalization, from VGAE paper(https://github.com/DaehanKim/vgae_pytorch/tree/master)
## extract coor, values and shape of a sparse matrix
def sparse_to_tuple(sparse_mx, clip_min=None, clip_max=None):
    if not sp.isspmatrix_coo(sparse_mx):
        sparse_mx = sparse_mx.tocoo()
    coords = np.vstack((sparse_mx.row, sparse_mx.col)).transpose()
    values = sparse_mx.data
    if clip_min != None or clip_max != None:
        values = np.clip(values, a_min=clip_min, a_max=clip_max)
    shape = sparse_mx.shape
    return coords, values, shape



class weighted_MSELoss(nn.Module):
    def __init__(self, reduction='mean'):
        super().__init__()
        self.reduction = reduction
    def forward(self, inputs, targets, weights):
        if self.reduction == 'mean':
            return (((inputs - targets)**2 ) * weights).sum() / weights.sum()
        else:
            return (((inputs - targets)**2 ) * weights).sum()
        


def get_weight_matrix(features, power=1):
    
    n_nodes, feat_RNA_dim = features.shape
    
    features = sp.csr_matrix(features)
    features = sparse_to_tuple(features)
    features = torch.sparse.FloatTensor(torch.LongTensor(features[0].T), 
                                torch.FloatTensor(features[1]), 
                                torch.Size(features[2]))

    #loss_norm = features.shape[0] * features.shape[0] / float((features.shape[0] * features.shape[0] - torch.sparse.sum(features)) * 2)
    pos_weight = float(features.shape[0] * features.shape[1] - len(features.coalesce().values())) / len(features.coalesce().values())
    pos_weight = pos_weight * power
    print("Weight of non-zero values:", pos_weight)
    if pos_weight > 1:
        weight_mask = features.to_dense().view(-1) != 0
        weight_tensor = torch.ones(weight_mask.size(0)) 
        weight_tensor[weight_mask] = pos_weight ## non-zero positions marked as pos_weight, zero positions marked as 1
    else:
        weight_tensor = torch.ones_like(features)
    weight_tensor = weight_tensor.view(n_nodes, feat_RNA_dim)
    
    return weight_tensor


