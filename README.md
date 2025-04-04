#  <img src="./logo.png" align="left" height="150" /></a>

<strong>DECLUST</strong> is a Python package developed to identify spatially coherent clusters of spots by integrating gene expression profiles with spatial coordinates in spatial transcriptomics data. It also enables accurate estimation of cell-type compositions within each cluster. The recommended starting point for using DECLUST is to follow the provided <a href="https://github.com/Qingyueee/DECLUST/blob/main/tutorial.ipynb" target="_blank">tutorial</a>。


<br> 

## 🌟 Features

 **Spatially-aware clustering:** Combines gene expression and spatial coordinates.

 **Robust deconvolution:** Aggregates signals over clusters to enhance cell type detection.

 **Easy to install:** Available via pip.

 **Visualization:** Includes modules for visualizing clustering and marker gene expression.

## ⏬ Installation

We recommend using a separate conda environment:

- Create a conda environment and install the DECLUST package

```bash
   conda create -n declust_env python=3.9
   conda activate declust_env

   pip install declust
```

- Additional dependencies can be installed using the `install_dependencies.sh` script:

```bash
   sh install_dependencies.sh
```

## 📊 Data Requirement

DECLUST requires AnnData files for processing and analysis. 

|             | Description                                                                                 |
|----------------------|---------------------------------------------------------------------------------------------|
| **sc_adata**    | Single-cell RNA sequencing (scRNA-seq) data, including gene expression profiles and cell type annotations.|
| **st_adata**    | Spatial transcriptomics (ST) data, including spatial gene expression profiles and spot coordinates. |

- Example Data Download  

   Download the `data.zip` file from [Google Drive](https://drive.google.com/uc?export=download&id=1LrSQYf1_IqQzxx7GeJrbBsEyuLLHHERC).


## 🔜 How to run DECLUST

The main workflow consists of:

- Data preprocessing
- Construction of Reference Matrix from Annotated Single-Cell Transcriptomic Data
- Identification of spatial clusters of spots from ST data
- Deconvolution
- Visualization

To start using DECLUST, import `declust` in your jupyter notebooks or/and scripts

```python
   import declust as dc
```

#### 1. Data Preprocessing

Load your spatial data and your single cell data (which should be in AnnData format)

```python
   sc_adata = sc.read_h5ad(path)
   sc_adata.var_names_make_unique()

   st_adata = sc.read_h5ad(path)
   st_adata.var_names_make_unique()

   coords_df = st_adata.obs.rename(
      columns={'array_row': 'x', 'array_col': 'y'}
   )[['x', 'y']]
```

Normalization

```python
   sc.pp.normalize_total(st_adata, target_sum=1e4)
   sc.pp.log1p(st_adata)

   sc.pp.normalize_total(sc_adata, target_sum=1e4)
   sc.pp.log1p(sc_adata)
```

And pre-process data using `dc.preprocessing.select_highly_variable_genes` and `dc.preprocessing.extract_labels_from_scdata`:

```python
   st_top_500_genes_df = st_adata[
      :, dc.preprocessing.select_highly_variable_genes(st_adata, n_top_genes=500)
   ].to_df()

   high_variable_common_genes = set(
      dc.preprocessing.select_highly_variable_genes(st_adata, n_top_genes=5000)
   ).intersection(sc_adata.var_names)

   sc_adata_overlapped = sc_adata[:, list(high_variable_common_genes)].copy()
   sc_labels = dc.preprocessing.extract_labels_from_scdata(
      sc_adata_overlapped, celltype_col, sample_col
   )
```

#### 2. Built Reference Matrix

Build a cell type reference matrix from annotated scRNA seq data.

```python
   sc_marker_gene = dc.marker_selection.generate_marker_genes(
      sc_adata_overlappe_path, sc_labels_path, output_path
   )

   sc_adata_marker = sc_adata[:, sc_adata.var_names.isin(sc_marker_gene_df.index)].copy()
```

#### 3. Clustering

DECLUST provides a three-step clustering process to identify spatial clusters. The process involves:

-  **Hierarchical Clustering**  

   This function clusters spatial transcriptomics spots based on the similarity of their highly variable gene expression profiles.
   
```python
   hierarchical_df = dc.hierarchical.clustering(
      st_top_500_genes_df, coords_df
   )
```

- **DBSCAN Clustering**  

   This function refines the initial clustering by applying DBSCAN to the results of the hierarchical clustering. It also selects the initial seeds for the SRG. You can choose to visualize the DBSCAN clustering results.
   
```python
   dbscan_centers_df = dc.dbscan.clustering(
      hierarchical_df, coords_df, visualize=True
   )
```

- **Seeded Region Growing (SRG)**  
   Finally, the SRG function further refines the clusters.
   
```python
   srg_df = dc.srg.clustering(
      dbscan_centers_df, coords_df, st_top_500_genes_df
   )
```

#### 4. Deconvolution

The final stage of the DECLUST pipeline is deconvolution, where cell type proportions are estimated from the clustered spatial transcriptomics data. By default, we apply an ordinary least squares (OLS) regression model to infer these proportions.

Alternatively, users may export the pseudo-bulk gene expression profiles aggregated by spatial clusters by `dc.deconvolution.generate_pseudo_bulk` function:

```python
   dc.deconvolution.generate_pseudo_bulk(
      st_adata, srg_df, save_csv=True, output_path
   )
```
This exported data can be used with any preferred deconvolution method.

To perform OLS-based deconvolution using DECLUST, the following function is provided:

```python
   decon_result_df = dc.deconvolution.ols(st_adata, sc_marker_adata, seg_df)
```

#### 5. Visualization

The DECLUST package provides two visualization functions to assist in interpreting the deconvolution results.

- **Cell-Type-Specific Marker Gene Expression**

   The `dc.visualize.declust_marker_boxplot` function displays the expression levels of top cell-type-specific marker genes across different cell types in the scRNA-seq reference data. This helps validate the specificity of selected markers.

```python
   dc.visualize.declust_marker_boxplot(sc_adata, sc_marker_adata, gene_idx)
```

- **Spatial Visualization of Deconvolution Results**

   The `dc.visualize.declust_results_visualize` function displays the spatial distribution of deconvolved cell type proportions, along with expression patterns of marker genes in the spatial transcriptomics data.

```python
   dc.visualize.declust_results_visualize(
      st_adata, sc_marker_adata, decon_result_df, coords_df, gene_idx, cell_type_idx
   )
```

## License  

GNU General Public License v3.0
