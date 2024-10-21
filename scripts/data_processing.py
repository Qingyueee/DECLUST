import scanpy as sc

# Selected highly variable genes from scRNA-seq data and spatial data overlapped genes
st_adata = sc.read_h5ad('./data/st_adata.h5ad')
st_adata.var_names_make_unique()

sc_adata = sc.read_h5ad('./data/sc_adata.h5ad')
sc_adata.var_names_make_unique()

st_gene_variances = st_adata.to_df().var()
st_highly_variable_genes = st_gene_variances.nlargest(5000) 

sc_gene_variances = sc_adata.to_df().var()
sc_highly_variable_genes = sc_gene_variances.nlargest(5000) 

high_variable_common_genes = st_highly_variable_genes.index.intersection(sc_highly_variable_genes.index)

sc_data_high_variable_overlapped = sc_adata[:, sc_adata.var_names.isin(high_variable_common_genes)].copy()
sc_data_high_variable_overlapped.to_df().to_csv('./data/sc_adata_high_variable_genes.csv')

# Save the labels of the scRNA-seq data for marker gene selection
sc_labels = sc_data_high_variable_overlapped.obs['celltype_major'].to_frame()
sc_labels = sc_labels.rename(columns={'celltype_major': 'cell_type'})
sample_id_old = list(sc_data_high_variable_overlapped.obs['Patient'])
id_unqs = sc_data_high_variable_overlapped.obs['Patient'].unique()

sample_id_new = []
for i, ids in enumerate (id_unqs):
    sample_id_new.append(f'sample{i+1}')
    
new_labels = dict(zip(id_unqs, sample_id_new))
sample_id_label = []
for x in sample_id_old:
    sample_id_label.append(new_labels[x])

celltype = list(set(sc_labels.cell_type))
sc_labels['cluster'] = [celltype.index(x)+1 for x in sc_labels.cell_type]
sc_labels['barcode'] = sc_labels.index
sc_labels['sample'] = sample_id_label

sc_labels.to_csv('./data/sc_labels.csv')