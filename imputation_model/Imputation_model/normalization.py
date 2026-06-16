import scipy.sparse as sp
import scanpy as sc
from scipy.sparse import csr_matrix
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import issparse



def RNA_scale(sc_RNA, common_genes, ref_cell_num, scale_type):

    sc_RNA.obs['cell_max_common_genes'] = sc_RNA[:, common_genes].X.max(axis=1).A.squeeze()
    sc_RNA.var['gene_max_ref'] = sc_RNA[:ref_cell_num].X.max(axis=0).A.squeeze()

    if scale_type == 'gene':
        sc_RNA.X = sp.csr_matrix(sc_RNA.X / sc_RNA.var['gene_max_ref'].values)
    elif scale_type == 'cell':
        sc_RNA.obs['cell_max_common_genes'] = sc_RNA.X.max(axis=1).A.squeeze()
        sc_RNA.X = sp.csr_matrix((sc_RNA.X.T / sc_RNA.obs['cell_max_common_genes'].values).T)
        
    return sc_RNA



def ref_query_norm(sc_RNA_ref, sc_RNA_query):
    ## Extract common and unique genes
    common_genes = sc_RNA_query.var_names.intersection(sc_RNA_ref.var_names)
    unique_genes_sc_RNA = sc_RNA_ref.var_names.difference(common_genes)

    ## Subset both datasets to the common genes
    sc_RNA_query_common = sc_RNA_query[:, common_genes].copy()
    sc_RNA_common = sc_RNA_ref[:, common_genes].copy()

    median_value = np.median(np.array(sc_RNA_query_common.X.sum(axis=1)))
    common_target_sum = median_value

    ## Calculate scaling factors based on the shared genes in sc_RNA
    sc_RNA_common_total_counts = np.array(sc_RNA_common.X.sum(axis=1)).squeeze()
    scaling_factors = (sc_RNA_ref.obs["total_counts"] - sc_RNA_common_total_counts) / sc_RNA_common_total_counts

    ## Normalize both datasets to the same total counts (shared genes)
    sc.pp.normalize_total(sc_RNA_query_common, target_sum=common_target_sum)
    sc.pp.normalize_total(sc_RNA_common, target_sum=common_target_sum)

    ## Normalize the other 9000 genes in sc_RNA using the scaling factors
    sc_RNA_unique = sc_RNA_ref[:, unique_genes_sc_RNA].copy()
    sc.pp.normalize_total(sc_RNA_unique, target_sum = common_target_sum)
    sc_RNA_unique.X = sp.csr_matrix((sc_RNA_unique.X.A.T * scaling_factors.values).T)

    ## Combine normalized shared and unique genes back into a single dataset for sc_RNA
    sc_RNA_ref_normalized = sc.AnnData(
        X=sp.hstack([sc_RNA_common.X, sc_RNA_unique.X]),
        obs=sc_RNA_ref.obs,
        var=pd.concat([sc_RNA_common.var, sc_RNA_unique.var])
    )

    ## Verify that the genes are in the correct order
    #sc_RNA_normalized.var_names_make_unique()
    
    padding_columns = sc_RNA_ref_normalized.shape[1] - sc_RNA_query_common.shape[1]
    zero_padding = sp.csr_matrix((sc_RNA_query_common.shape[0], padding_columns))
    sc_RNA_query_padded = sp.hstack([sc_RNA_query_common.X, zero_padding])
    
    sc_RNA = sp.vstack([sc_RNA_ref_normalized.X, sc_RNA_query_padded])
    
    sc_RNA = sc.AnnData(
        X = sc_RNA,
        obs = pd.concat([sc_RNA_ref_normalized.obs, sc_RNA_query_common.obs]),
        var = sc_RNA_ref_normalized.var
        )
    
    return sc_RNA, sc_RNA_ref_normalized, sc_RNA_query_common, common_genes, unique_genes_sc_RNA