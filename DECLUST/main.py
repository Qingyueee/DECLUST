import argparse
import scanpy as sc
import pandas as pd

from DECLUST.Hierarchical_clustering import hierarchical_clustering, selected_number_of_clusters
from DECLUST.DBSCAN_clustering import dbscan_clustering
from DECLUST.SRG import srg_process
from DECLUST.DECLUST_deconvolution import deconvolution

def load_and_preprocess_data(st_path, sc_path, marker_genes_path):
    """Load and preprocess spatial transcriptomics and single-cell RNA-seq data."""
    # Loading spatial transcriptomics data
    print(f"Loading spatial transcriptomics data from {st_path}...")
    st_adata = sc.read_h5ad(st_path)
    st_adata.var_names_make_unique()
    sc.pp.normalize_total(st_adata, target_sum=1e4)
    sc.pp.log1p(st_adata)
    st_highly_variable_genes_df = st_adata[:, st_adata.to_df().var().nlargest(500).index].to_df()
    coords = st_adata.obs[['array_row', 'array_col']].rename(columns={'array_row': 'x', 'array_col': 'y'}).copy()

    # Loading single-cell RNA-seq data
    print(f"Loading single-cell RNA-seq data from {sc_path}...")
    sc_adata = sc.read_h5ad(sc_path)
    sc_adata.var_names_make_unique()
    sc.pp.normalize_total(sc_adata, target_sum=1e4)
    sc.pp.log1p(sc_adata)
    
    # Loading marker genes
    print(f"Loading marker genes from {marker_genes_path}...")
    sc_marker_gene = pd.read_csv(marker_genes_path, index_col='Gene')
    sc_adata_marker_h5ad = sc_adata[:, sc_adata.var_names.isin(sc_marker_gene.index)].copy()

    return st_adata, st_highly_variable_genes_df, coords, sc_marker_gene, sc_adata, sc_adata_marker_h5ad

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Run DECLUST on spatial transcriptomics and single-cell RNA-seq data.")
    parser.add_argument('--st_adata', required=True, help='Path to spatial transcriptomics data (h5ad file).')
    parser.add_argument('--sc_adata', required=True, help='Path to single-cell RNA-seq data (h5ad file).')
    parser.add_argument('--marker_genes', required=True, help='Path to marker genes CSV file.')
    parser.add_argument('--visualize', action='store_true', help='Set this flag to visualize DBSCAN clustering results.')
    args = parser.parse_args()

    # Load and preprocess data
    print("Loading and preprocessing data...")
    st_adata, st_highly_variable_genes_df, coords, sc_marker_gene, sc_adata, sc_adata_marker_h5ad = load_and_preprocess_data(
        args.st_adata, args.sc_adata, args.marker_genes
    )

    # Step 1: Hierarchical Clustering
    print("Running hierarchical clustering...")
    selected_number_of_clusters(st_highly_variable_genes_df)
    while True:
        try:
            num_clusters = int(input("Please enter the number of clusters based on the elbow method: "))
            if num_clusters > 1:
                break
            else:
                print("The number of clusters must be greater than 1. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a valid integer for the number of clusters.")

    hierarchical_results = hierarchical_clustering(st_highly_variable_genes_df, num_clusters, coords)

    # Step 2: DBSCAN Clustering
    print("Running DBSCAN clustering...")
    dbscan_centers_df = dbscan_clustering(hierarchical_results, coords, visualize=args.visualize)

    # Step 3: Seeded Region Growing
    print("Running seeded region growing...")
    neighbours_rules = [(-2, 0), (2, 0), (0, -2), (0, 2), (-1, -1), (1, -1), (-1, 1), (1, 1)]
    label_df = srg_process(dbscan_centers_df, coords, st_highly_variable_genes_df, neighbours_rules)

    # Step 4: Deconvolution
    print("Running deconvolution...")
    DECLUST_df = deconvolution(st_adata, sc_marker_gene, sc_adata, sc_adata_marker_h5ad, label_df)

    print("Deconvolution results:")
    print(DECLUST_df.head())

    # Optionally, save the results to a CSV file
    output_file = './DECLUST_results.csv'
    print(f"Saving results to {output_file}...")
    DECLUST_df.to_csv(output_file, index=False)
    print("Results saved.")

if __name__ == '__main__':
    main()
