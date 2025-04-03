import scanpy as sc
import pandas as pd
import numpy as np

def select_highly_variable_genes(adata, n_top_genes=5000):

    """
    Selects the top highly variable genes from an AnnData object.

    This function computes the variance of each gene from the expression data in the AnnData object,
    and then returns a list of the gene names corresponding to the top n_top_genes with the highest variance.

    Parameters:
        adata (anndata.AnnData):
            An AnnData object containing gene expression data. The function converts the expression matrix
            to a pandas DataFrame to compute variances.
        
        n_top_genes (int, optional):
            The number of top highly variable genes to select. Defaults to 5000.

    Returns:
        list:
            A list of gene names (strings) corresponding to the top n_top_genes with the highest variance.
    """

    gene_variances = adata.to_df().var()
    highly_variable_genes = gene_variances.nlargest(n_top_genes).index
    return list(highly_variable_genes)

def extract_labels_from_scdata(adata, celltype_col='celltype_major', sample_col='Patient'):

    """
    Extract cell type labels, assign numeric clusters, barcodes, and renamed sample IDs from an AnnData object.

    Parameters:
    - adata: AnnData object with `.obs` containing cell type and patient/sample information.
    - celltype_col: Name of the column in `adata.obs` indicating cell types.
    - sample_col: Name of the column in `adata.obs` indicating patient/sample IDs.

    Returns:
    - sc_labels: DataFrame with columns ['cell_type', 'cluster', 'barcode', 'sample']
    """
    
    sc_labels = adata.obs[celltype_col].to_frame()
    sc_labels = sc_labels.rename(columns={celltype_col: 'cell_type'})
    
    sample_id_old = list(adata.obs[sample_col])
    id_unqs = adata.obs[sample_col].unique()
    
    sample_id_new = [f'sample{i+1}' for i in range(len(id_unqs))]
    new_labels = dict(zip(id_unqs, sample_id_new))
    
    sample_id_label = [new_labels[x] for x in sample_id_old]
    
    celltype = list(set(sc_labels['cell_type']))
    sc_labels['cluster'] = [celltype.index(x) + 1 for x in sc_labels['cell_type']]
    sc_labels['barcode'] = sc_labels.index
    sc_labels['sample'] = sample_id_label
    
    return sc_labels

