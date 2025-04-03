import numpy as np
import pandas as pd
import os

def generate_pseudo_bulk(st_adata, label_df, save_csv=False, output_path="pseudo_bulk.csv"):
    """
    Generate a pseudo-bulk expression profile by summing expression across spots in each cluster.

    Parameters:
    - st_adata: AnnData object containing spatial transcriptomics data.
    - label_df: DataFrame with cluster labels. Must have a column named 'label'.
    - save_csv: Boolean, whether to save the result as a CSV file. Default is False.
    - output_path: Path to save the CSV file if save_csv is True. Default is 'pseudo_bulk.csv'.

    Returns:
    - pseudo_bulk_df: DataFrame of summed expression values per cluster.
    """
    unique_labels = sorted(label_df['label'].unique())
    index_lists = {
        f'cluster_index_list_{label}': label_df[label_df['label'] == label].index.tolist()
        for label in unique_labels
    }

    pseudo_bulk_df = pd.DataFrame()
    for label in unique_labels:
        cluster_index_list = index_lists[f'cluster_index_list_{label}']
        selected_rows = st_adata.to_df().loc[cluster_index_list]
        sum_values = selected_rows.sum(axis=0)
        pseudo_bulk_df = pd.concat([pseudo_bulk_df, sum_values.to_frame().T])

    pseudo_bulk_df.index = unique_labels

    if save_csv:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        pseudo_bulk_df.to_csv(output_path)
        print(f"Pseudo-bulk expression data saved to: {output_path}")

    return pseudo_bulk_df

def multiple_linear_regression_OLS(Y_star, X1_star):
    beta_hat = np.linalg.lstsq(Y_star, X1_star, rcond=None)[0]
    return beta_hat

def nonnegative_constraint(beta_hat):
    beta_hat[beta_hat < 0] = 0
    beta_normalized = beta_hat / np.sum(beta_hat)
    return beta_normalized

def ols(st_adata, sc_adata_marker_h5ad, label_df):
    """
    Perform deconvolution using ordinary least squares (OLS) regression to estimate 
    cell-type proportions in spatial transcriptomics data.

    Parameters:
        st_adata (anndata.AnnData):
            Spatial transcriptomics data in an AnnData object. The .obs attribute should 
            contain the spatial spot names (obs_names) used for reindexing the final output.
        
        sc_adata_marker_h5ad (anndata.AnnData):
            Single-cell RNA-seq data with marker information in an AnnData object. Its .obs 
            attribute must include a 'celltype_major' column indicating the major cell type for 
            each cell. The .var_names represent the gene names.
        
        label_df (pandas.DataFrame):
            A DataFrame containing labels for the spatial spots. It should have a 'label' column 
            and its index corresponds to the spatial spots in st_adata.

    Returns:
        pandas.DataFrame:
            A DataFrame (DECLUST_df) with estimated deconvolution coefficients for each spatial 
            spot. The rows are reindexed to match st_adata.obs_names and the columns represent 
            cell types.
    """   
    print("Running deconvolution...")

    pseudo_bulk_df = generate_pseudo_bulk(st_adata, label_df)

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
