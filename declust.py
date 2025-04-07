import os
import argparse
import scanpy as sc
import pandas as pd
import declust as dc

# -----------------------------------------
# Utility: create folders
# -----------------------------------------
os.makedirs('data', exist_ok=True)
os.makedirs('results', exist_ok=True)

# -----------------------------------------
# Argument parser
# -----------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--module', type=str, choices=['marker', 'cluster', 'pseudo_bulk', 'deconv', 'visualize'], required=True)
parser.add_argument('--celltype_col', type=str, help='Column name in scRNA-seq data for cell type')
parser.add_argument('--sample_col', type=str, help='Column name in scRNA-seq data for sample')

parser.add_argument('--num_clusters', type=int, default=None,
                    help='Manually specify number of clusters. If not set, auto-selected by elbow method.')
parser.add_argument('--visualize_selection', action='store_true', help='Show elbow and silhouette plots for cluster selection')
parser.add_argument('--visualize_hierarchical', action='store_true', help='Show final clustering result')
parser.add_argument('--visualize_dbscan', action='store_true', help='Show DBSCAN clustering result')
parser.add_argument('--visualize_srg', action='store_true', help='Show final SRG result plot')

args = parser.parse_args()

# -----------------------------------------
# Load data
# -----------------------------------------
print("\U0001F4E5 Loading scRNA-seq and ST data...")
sc_adata = sc.read_h5ad('data/sc_adata.h5ad')
st_adata = sc.read_h5ad('data/st_adata.h5ad')
st_adata.var_names_make_unique()
sc_adata.var_names_make_unique()
coords_df = st_adata.obs.rename(columns={'array_row': 'x', 'array_col': 'y'})[['x', 'y']]

sc.pp.normalize_total(st_adata, target_sum=1e4)
sc.pp.log1p(st_adata)
sc.pp.normalize_total(sc_adata, target_sum=1e4)
sc.pp.log1p(sc_adata)

st_top_500_genes_df = st_adata[:, dc.preprocessing.select_highly_variable_genes(st_adata, n_top_genes=500)].to_df()

# -----------------------------------------
# Step 1: Marker selection
# -----------------------------------------
def prepare_marker_genes():
    if not args.celltype_col or not args.sample_col:
        print("\n❌ Error: You must specify both --celltype_col and --sample_col when running marker gene selection.")
        print("📌 Available options in scRNA-seq data:", sc_adata.obs.columns.tolist())
        exit(1)

    if not os.path.exists('data/sc_adata_overlapped.csv') or not os.path.exists('data/sc_labels.csv'):
        print("⚙️ Generating marker gene input files...")
        common_genes = set(dc.preprocessing.select_highly_variable_genes(st_adata, n_top_genes=5000)).intersection(sc_adata.var_names)
        sc_adata_overlapped = sc_adata[:, list(common_genes)].copy()
        sc_labels = dc.preprocessing.extract_labels_from_scdata(
            sc_adata_overlapped,
            celltype_col=args.celltype_col,
            sample_col=args.sample_col
        )
        sc_adata_overlapped.to_df().to_csv('data/sc_adata_overlapped.csv')
        sc_labels.to_csv('data/sc_labels.csv')

    print("🔍 Selecting marker genes...")
    return dc.marker_selection.generate_marker_genes(
        'data/sc_adata_overlapped.csv',
        'data/sc_labels.csv',
        'results/marker_genes.csv'
    )


# -----------------------------------------
# Step 2: Clustering
# -----------------------------------------
def run_clustering():
    if not os.path.exists('results/marker_genes.csv'):
        print("⚠️ Marker gene file not exit. Auto-running marker step...")
        prepare_marker_genes()
    print("📦 Running clustering...")
    sc_marker_gene_df = pd.read_csv('results/marker_genes.csv', index_col='Gene')
    sc_adata_marker = sc_adata[:, sc_adata.var_names.isin(sc_marker_gene_df.index)].copy()
    hier_df = dc.hierarchical.clustering(
        st_top_500_genes_df,
        coords_df,
        num_clusters=args.num_clusters,
        show_plot=args.visualize_hierarchical,
        save_path='results/hierarchical' if not args.visualize_hierarchical else None,
        show_selection_plot=args.visualize_selection
    )
    dbscan_df = dc.dbscan.clustering(
        hier_df,
        coords_df,
        visualize=args.visualize_dbscan,
        plot_save_dir=None if args.visualize_dbscan else 'results/dbscan' 
    )

    # srg_df = dc.srg.clustering(dbscan_df, coords_df, st_top_500_genes_df)
    srg_df = dc.srg.clustering(
        dbscan_df,
        coords_df,
        st_top_500_genes_df,
        show_plot=args.visualize_srg,
        save_path=None if args.visualize_srg else 'results/srg'
    )
    srg_df.to_csv('results/srg_df.csv')
    return srg_df, sc_adata_marker

# -----------------------------------------
# Step 3: Pseudo bulk
# -----------------------------------------
def generate_pseudo_bulk():
    if not os.path.exists('results/srg_df.csv'):
        print("⚠️ Clustering result not exit. Auto-running clustering step...")
        run_clustering()
    srg_df = pd.read_csv('results/srg_df.csv', index_col=0)
    print("🧪 Generating pseudo bulk expression...")
    dc.deconvolution.generate_pseudo_bulk(st_adata, srg_df, save_csv=True, output_path="results/pseudo_bulk.csv")

# -----------------------------------------
# Step 4: Deconvolution
# -----------------------------------------
def run_deconvolution():
    if not os.path.exists('results/marker_genes.csv'):
        print("⚠️ Marker gene file not exit. Auto-running marker step...")
        prepare_marker_genes()
    if not os.path.exists('results/srg_df.csv'):
        print("⚠️ Clustering result not exit. Auto-running clustering step...")
        run_clustering()
    srg_df = pd.read_csv('results/srg_df.csv', index_col=0)
    sc_marker_gene_df = pd.read_csv('results/marker_genes.csv', index_col='Gene')
    sc_adata_marker = sc_adata[:, sc_adata.var_names.isin(sc_marker_gene_df.index)].copy()
    print("🧩 Performing deconvolution with DECLUST...")
    DECLUST_df = dc.deconvolution.ols(st_adata, sc_adata_marker, srg_df)
    DECLUST_df.to_csv('results/DECLUST_result.csv')

# -----------------------------------------
# Step 5: Visualization
# -----------------------------------------
def visualize():
    marker_path = 'results/marker_genes.csv'
    declust_path = 'results/DECLUST_result.csv'

    if not os.path.exists(marker_path) or not os.path.exists(declust_path):
        print("❌ Required files for visualization not found. Please run deconv step first.")
        return

    sc_marker_gene_df = pd.read_csv(marker_path, index_col='Gene')
    DECLUST_df = pd.read_csv(declust_path, index_col=0)

    print("\n📊 Visualization Options:")
    print("1. Marker Gene Boxplot")
    print("2. DECLUST Deconvolution Results")
    while True:
        choice = input("Select option (1 or 2): ").strip()
        if choice in ['1', '2']:
            break

    if choice == '1':
        if 'tstat' in sc_marker_gene_df.columns:
            sorted_genes = sc_marker_gene_df.sort_values('tstat', ascending=False)
            print("\n🔝 Top 5 marker genes (by t-statistic):", sorted_genes.index[:5].tolist())
            see_all = input("See all marker genes? (yes/no): ").strip().lower()
            if see_all in ['yes', 'y', 'all']:
                print(sorted_genes.index.tolist())
        while True:
            gene_name = input("Enter gene name for boxplot: ").strip()
            if gene_name in sc_adata.var_names:
                break
        dc.visualize.declust_marker_boxplot(sc_adata, sc_marker_gene_df, gene_name)

    elif choice == '2':
        print("\n🔎 Available cell types:", DECLUST_df.columns.tolist())
        while True:
            cell_type = input("Enter cell type: ").strip()
            if cell_type in DECLUST_df.columns:
                break
        print("\nMethods: \n1. Single gene\n2. Sum of markers\n3. Mean of markers")
        method = input("Choose method (1/2/3): ").strip()
        if method == '1':
            markers = sc_marker_gene_df[sc_marker_gene_df['maxgroup'] == cell_type]
            if 'tstat' in markers.columns:
                markers = markers.sort_values('tstat', ascending=False)
            print(f"\nTop 5 markers for {cell_type}: {markers.index[:5].tolist()}")
            see_all = input("See all markers? (yes/no): ").strip().lower()
            if see_all in ['yes', 'y', 'all']:
                print(markers.index.tolist())
            while True:
                gene_name = input("Enter gene name to visualize: ").strip()
                if gene_name in markers.index:
                    break
            dc.visualize.declust_results_visualize(st_adata, sc_marker_gene_df, DECLUST_df, coords_df, gene_name=gene_name, cell_type=cell_type)
        elif method == '2':
            dc.visualize.declust_results_visualize(st_adata, sc_marker_gene_df, DECLUST_df, coords_df, cell_type=cell_type, agg_method='sum')
        elif method == '3':
            dc.visualize.declust_results_visualize(st_adata, sc_marker_gene_df, DECLUST_df, coords_df, cell_type=cell_type, agg_method='mean')

# -----------------------------------------
# Stage dispatching
# -----------------------------------------
if args.module == 'marker':
    prepare_marker_genes()
    print("✅ Marker Gene Selection Completed")

elif args.module == 'cluster':
    run_clustering()
    print("✅ Clustering Completed")

elif args.module == 'pseudo_bulk':
    generate_pseudo_bulk()
    print("✅ Pseudo Bulk Gene Profile Saving Completed")

elif args.module == 'deconv':
    run_deconvolution()
    print("✅ Deconvolution Completed")

elif args.module == 'visualize':
    visualize()
    # print("✅ Stage completed: visualize")