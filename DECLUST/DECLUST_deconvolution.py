import numpy as np
import scanpy as sc
import pandas as pd


def multiple_linear_regression_OLS(Y_star, X1_star):
    beta_hat = np.linalg.lstsq(Y_star, X1_star, rcond=None)[0]
    return beta_hat

def nonnegative_constraint(beta_hat):
    beta_hat[beta_hat < 0] = 0
    beta_normalized = beta_hat / np.sum(beta_hat)
    return beta_normalized

def deconvolution(st_adata, sc_marker_gene, sc_adata, sc_adata_marker_h5ad, label_df):
    unique_labels = sorted(label_df['label'].unique())
    index_lists = {}
    for label in unique_labels:
        index_lists[f'cluster_index_list_{label}'] = label_df[label_df['label'] == label].index.tolist()

    pseudo_bulk_df = pd.DataFrame()
    for label in unique_labels:
        cluster_index_list = index_lists.get(f'cluster_index_list_{label}', [])
        selected_rows = st_adata.to_df().loc[cluster_index_list]
        sum_values = selected_rows.sum(axis=0)
        pseudo_bulk_df = pd.concat([pseudo_bulk_df, sum_values.to_frame().T])
    pseudo_bulk_df.index = unique_labels

    unique_cell_types = sc_adata_marker_h5ad.obs['celltype_major'].unique()
    celltype_indexes = {cell_type: cell_type for cell_type in unique_cell_types}
    mean_cell_type_df = pd.DataFrame(index=celltype_indexes.keys(), columns=sc_adata_marker_h5ad.var_names)

    for celltype, label in celltype_indexes.items():
        indexes = sc_adata_marker_h5ad.obs['celltype_major'][sc_adata_marker_h5ad.obs['celltype_major'] == label].index
        selected_rows = sc_adata_marker_h5ad[indexes, :]
        mean_cell_type_df.loc[celltype] = np.mean(selected_rows.X.toarray(), axis=0)

    common_columns = mean_cell_type_df.columns.intersection(pseudo_bulk_df.columns)
    mean_cell_type_df = mean_cell_type_df.loc[:, common_columns]
    pseudo_bulk_df = pseudo_bulk_df.loc[:, common_columns]
    mean_cell_type_df = mean_cell_type_df.astype(float)
    pseudo_bulk_df = pseudo_bulk_df.astype(float)

    normalized_coefficients_df = pd.DataFrame(index=pseudo_bulk_df.index, columns=mean_cell_type_df.index)
    for index, row in pseudo_bulk_df.iterrows():
        beta_hat = multiple_linear_regression_OLS(mean_cell_type_df.T, row)
        beta_normalized = nonnegative_constraint(beta_hat)
        normalized_coefficients_df.loc[index] = beta_normalized

    label_dataframes = {}
    for label in label_df['label'].unique():
        indexes_with_label = label_df.index[label_df['label'] == label].tolist()
        new_df_label = pd.DataFrame(index=indexes_with_label, columns=normalized_coefficients_df.columns)
        values_for_new_df_label = normalized_coefficients_df.loc[label]

        for index in new_df_label.index:
            new_df_label.loc[index] = values_for_new_df_label

        label_dataframes[label] = new_df_label
    DECLUST_df = pd.concat(label_dataframes.values())

    DECLUST_df = DECLUST_df.reindex(st_adata.obs_names)
    return DECLUST_df


# # Loading data
# st_adata = sc.read_h5ad('./data/st_adata.h5ad')
# st_adata.var_names_make_unique()
# sc.pp.normalize_total(st_adata, target_sum=1e4)
# sc.pp.log1p(st_adata)
# sc_marker_gene = pd.read_csv('./data/marker_genes.csv', index_col='Gene')
# sc_adata= sc.read_h5ad('./data/sc_adata.h5ad')
# sc.pp.normalize_total(sc_adata, target_sum=1e4)
# sc.pp.log1p(sc_adata)
# sc_adata_marker_h5ad = sc_adata[:, sc_adata.var_names.isin(sc_marker_gene.index)].copy()
# label_df = pd.read_csv('./srg_results.csv', index_col=0)

# unique_labels = sorted(label_df['label'].unique())
# index_lists = {}
# for label in unique_labels:
#     index_lists[f'cluster_index_list_{label}'] = label_df[label_df['label'] == label].index.tolist()

# pseudo_bulk_df = pd.DataFrame()
# for label in unique_labels:
#     cluster_index_list = index_lists.get(f'cluster_index_list_{label}', [])
#     selected_rows = st_adata.to_df().loc[cluster_index_list]
#     sum_values = selected_rows.sum(axis=0)
#     pseudo_bulk_df = pd.concat([pseudo_bulk_df, sum_values.to_frame().T])
# pseudo_bulk_df.index = unique_labels

# unique_cell_types = sc_adata_marker_h5ad.obs['celltype_major'].unique()
# celltype_indexes = {cell_type: cell_type for cell_type in unique_cell_types}
# mean_cell_type_df = pd.DataFrame(index=celltype_indexes.keys(), columns=sc_adata_marker_h5ad.var_names)

# for celltype, label in celltype_indexes.items():
#     indexes = sc_adata_marker_h5ad.obs['celltype_major'][sc_adata_marker_h5ad.obs['celltype_major'] == label].index
#     selected_rows = sc_adata_marker_h5ad[indexes, :]
#     mean_cell_type_df.loc[celltype] = np.mean(selected_rows.X.toarray(), axis=0)

# common_columns = mean_cell_type_df.columns.intersection(pseudo_bulk_df.columns)
# mean_cell_type_df = mean_cell_type_df.loc[:, common_columns]
# pseudo_bulk_df = pseudo_bulk_df.loc[:, common_columns]
# mean_cell_type_df = mean_cell_type_df.astype(float)
# pseudo_bulk_df = pseudo_bulk_df.astype(float)

# normalized_coefficients_df = pd.DataFrame(index=pseudo_bulk_df.index, columns=mean_cell_type_df.index)
# for index, row in pseudo_bulk_df.iterrows():
#     beta_hat = multiple_linear_regression_OLS(mean_cell_type_df.T, row)
#     beta_normalized = nonnegative_constraint(beta_hat)
#     normalized_coefficients_df.loc[index] = beta_normalized

# label_dataframes = {}
# for label in label_df['label'].unique():
#     indexes_with_label = label_df.index[label_df['label'] == label].tolist()
#     new_df_label = pd.DataFrame(index=indexes_with_label, columns=normalized_coefficients_df.columns)
#     values_for_new_df_label = normalized_coefficients_df.loc[label]

#     for index in new_df_label.index:
#         new_df_label.loc[index] = values_for_new_df_label

#     label_dataframes[label] = new_df_label
# DECLUST_df = pd.concat(label_dataframes.values())

# DECLUST_df = DECLUST_df.reindex(st_adata.obs_names)
# print(DECLUST_df.head(5))