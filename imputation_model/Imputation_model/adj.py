import torch
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse as sp
#from sklearn.metrics import DistanceMetric
from sklearn.neighbors import KDTree
from sklearn.metrics import DistanceMetric
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.neighbors import NearestNeighbors


def find_mutual_nn(data1, 
                   data2, 
                   dist_method, 
                   k1, 
                   k2, 
                  ):
    if dist_method == 'cosine':
        cos_sim1 = cosine_similarity(data1, data2)
        cos_sim2 = cosine_similarity(data2, data1)
        k_index_1 = torch.topk(torch.tensor(cos_sim2), k=k2, dim=1)[1]
        k_index_2 = torch.topk(torch.tensor(cos_sim1), k=k1, dim=1)[1]
    else:
        dist = DistanceMetric.get_metric(dist_method)
        k_index_1 = KDTree(data1, metric=dist).query(data2, k=k2, return_distance=False)
        k_index_2 = KDTree(data2, metric=dist).query(data1, k=k1, return_distance=False)
    mutual_1 = []
    mutual_2 = []
    mutual = []
    for index_2 in range(data2.shape[0]):
        for index_1 in k_index_1[index_2]:
            if index_2 in k_index_2[index_1]: 
                mutual_1.append(index_1)
                mutual_2.append(index_2)
                mutual.append([index_1, index_2])
    return mutual



def intra_exp_adj(input_feature, 
                  find_neighbor_method='KNN', 
                  dist_method='euclidean', 
                  corr_dist_neighbors=20, 
                  ):
    
    n_samples = input_feature.shape[0]
    row_indices = []
    col_indices = []
        
    if find_neighbor_method == 'KNN':
        if dist_method == 'cosine':
            cos_sim = cosine_similarity(input_feature, input_feature)
            k_index = torch.topk(torch.tensor(cos_sim), k=corr_dist_neighbors, dim=1)[1]
        else:
            dist = DistanceMetric.get_metric(dist_method)
            k_index = KDTree(input_feature, metric=dist).query(input_feature, k=corr_dist_neighbors, return_distance=False)
        for i in range(k_index.shape[0]):
            for j in k_index[i]:
                if i != j:
                    row_indices.append(i)
                    col_indices.append(j)
                    row_indices.append(j) 
                    col_indices.append(i)
    elif find_neighbor_method == 'MNN':
        mut = find_mutual_nn(input_feature, input_feature, dist_method=dist_method, k1=corr_dist_neighbors, k2=corr_dist_neighbors)
        mut = pd.DataFrame(mut, columns=['data1', 'data2'])
        for i in mut.index:
            data1 = mut.loc[i, 'data1']
            data2 = mut.loc[i, 'data2']
            if data1 != data2:
                row_indices.append(data1)
                col_indices.append(data2)
                row_indices.append(data2)
                col_indices.append(data1)

    data = np.ones(len(row_indices))
    A_exp = sp.csr_matrix((data, (row_indices, col_indices)), shape=(n_samples, n_samples))
    
    return A_exp



def inter_adj(data1,
              data2,
              find_neighbor_method='MNN',
              dist_method='euclidean',
              corr_dist_neighbors=20, 
             ):
    data1_num = data1.shape[0]
    data2_num = data2.shape[0]
    row_indices = []
    col_indices = []
    
    if find_neighbor_method == 'KNN':
        if dist_method == 'cosine':
            cos_sim = cosine_similarity(data1, data2)
            k_index = torch.topk(torch.tensor(cos_sim), k=corr_dist_neighbors, dim=1)[1]
        else:
            dist = DistanceMetric.get_metric(dist_method)
            k_index = KDTree(data2, metric=dist).query(data1, k=corr_dist_neighbors, return_distance=False)
        A_exp = np.zeros((data1_num+data2_num, data1_num+data2_num), dtype=float)
        for i in range(k_index.shape[0]):
            for j in k_index[i]:
                row_indices.append(i)
                col_indices.append(data1_num + j)  # Map to data2 indices
                row_indices.append(data1_num + j)  # Ensure symmetry
                col_indices.append(i)
    elif find_neighbor_method == 'MNN':
        mut = find_mutual_nn(data2, data1, dist_method=dist_method, k1=corr_dist_neighbors, k2=corr_dist_neighbors)
        mut = pd.DataFrame(mut, columns=['data2', 'data1'])
        for i in mut.index:
            row_indices.append(mut.loc[i, 'data1'])
            col_indices.append(data1_num + mut.loc[i, 'data2'])
            row_indices.append(data1_num + mut.loc[i, 'data2'])
            col_indices.append(mut.loc[i, 'data1'])

    data = np.ones(len(row_indices))
    A_exp = sp.csr_matrix((data, (row_indices, col_indices)), shape=(data1_num + data2_num, data1_num + data2_num))
    
    return A_exp



## adjacency matrix normalization, from STdGCN
def adj_normalize(mx, symmetry=True):

    mx = sp.csr_matrix(mx)
    rowsum = np.array(mx.sum(1))
    r_inv = np.power(rowsum, -1).flatten() # flatten(): dimension [m,n] into one dimension [m*n].
    r_inv[np.isinf(r_inv)] = 0. 
    if symmetry == True:
        r_mat_inv = sp.diags(np.sqrt(r_inv)) # generate a sparse matrix, the diag is 'sqrt(r_inv)'
        mx = r_mat_inv.dot(mx).dot(r_mat_inv)
    else:
        r_mat_inv = sp.diags(r_inv) 
        mx = r_mat_inv.dot(mx)
    
    return mx


