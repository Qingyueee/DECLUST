# DECLUST

**DECLUST** is a cluster-based cell-type deconvolution method designed for spatial transcriptomic data. By integrating gene expression profiles with spatial coordinates, DECLUST preserves the tissue's spatial structure and overcomes the low expression levels observed in individual spots. The method first identifies spatial clusters and then applies deconvolution on the aggregated gene expression within each cluster, leading to improved cell type identification.

![workflow](./overview.jpg)

## Features

🌟 **Spatially-aware clustering:** Combines gene expression and spatial coordinates.

🌟 **Robust deconvolution:** Aggregates signals over clusters to enhance cell type detection.

🌟 **Easy to install:** Available via pip.

🌟 **Visualization:** Includes modules for visualizing clustering and marker gene expression.

## Installation

We suggest using a separate conda environment for installing DECLUST.

- Create conda environment and install DECLUST package

```bash
   conda create -n declust_env python=3.9
   conda activate declust_env

   pip install declust
```



Some dependencies need to be installed separately for optimal compatibility:


- Scanpy

   Install via pip:
```python
   pip install scanpy
```

- rpy2

   It is recommended to install rpy2 via conda-forge to avoid build issues:
   
```python
   conda install -c conda-forge rpy2
```

- R Packages

   DECLUST requires the dplyr R package. After installing rpy2, you can install dplyr via Python:
   
```python
   from rpy2.robjects.packages import importr
   utils = importr('utils')
   utils.install_packages('dplyr')
```

## Usage
DECLUST provides a Python API. Below are quick start instructions.

#### 1. Preprocess Data 

The DECLUST package provides two utility functions for preprocessing your data.

**Selecting Highly Variable Genes**

The `declust.preprocessing.select_highly_variable_genes` function computes the variance of each gene from the expression data in an AnnData object and returns a list of gene names corresponding to the top variable genes.

```python
   declust.preprocessing.select_highly_variable_genes(anndata, n_top_genes=5000)
```

**Extracting Labels from Single-Cell Data**

The `declust.preprocessing.extract_labels_from_scdata` function extracts cell type labels, assigns numeric clusters, retains barcodes, and renames sample IDs from an AnnData object. It returns a DataFrame with the following columns: cell_type, cluster, barcode, and sample.

```python
   declust.preprocessing.extract_labels_from_scdata(anndata, celltype_col='celltype_major',
                                                   sample_col='Patient')
```

#### 2. Generating Marker Genes

Once you have preprocessed your single-cell data and extracted the necessary labels, you can generate marker genes by running the following function:

```python
   declust.marker_selection.generate_marker_genes(sc_data_path, sc_labels_path, output_path)
```

#### 3. Clustering

DECLUST provides a three-step clustering process to identify spatial clusters and refine the deconvolution results. The process involves:

a. **Hierarchical Clustering**  
   This function clustering spots based on the spatial transcriptome highly variable gene expression profiles.
   
   ```python
      hierarchical_df = declust.hierarchical.clustering(st_highly_variable_genes, coords)
   ```

b. **DBSCAN Clustering**  
   This function refines the initial clustering by applying DBSCAN to the results of the hierarchical clustering. It also selects the initial seeds for the SRG. You can choose to visualize the DBSCAN clustering results.
   
   ```python
      dbscan_centers_df = declust.dbscan.clustering(hierarchical_df, coords, visualize=True)
   ```

c. **Seeded Region Growing (SRG)**  
   Finally, the SRG function further refines the clusters.
   
   ```python
      srg_df = declust.srg.clustering(dbscan_centers_df, coords, st_highly_variable_genes)
   ```

Each step builds on the previous one:

- Hierarchical Clustering: Spots are initially grouped based on their overall similarity.

- DBSCAN Clustering: Within these hierarchical groups, DBSCAN identifies denser, local clusters, and selects initial seeds.

- Seeded Region Growing: Finally, the SRG approach refines the cluster labels and completes the clusters by leveraging spatial information.

#### 4. Deconvolution

The final stage of the DECLUST pipeline is deconvolution. Here, we use a simple ordinary least squares (OLS) model to estimate cell type proportions from the clustered spatial transcriptomics data. Alternatively, if you prefer another deconvolution method, you can export the pseudo bulk data generated during clustering and use it as input for your own model.

```python
   declust.deconvolution.generate_pseudo_bulk(st_anndata, label_df, save_csv=True, 
                                             output_path="results/pseudo_bulk.csv")

   declust.deconvolution.ols(st_anndata, sc_marker_anndata, label_df)
```
The function `generate_pseudo_bulk` aggregates the gene expression data according to the clusters and optionally saves the result as a CSV file. The `ols` function then uses this aggregated (pseudo bulk) data to perform deconvolution via an OLS model, producing the estimated cell type proportions.

#### 5. Visualization

The DECLUST package provides two visualization functions to help you explore the deconvolution results.

a. **Marker Gene Expression Boxplot**

Use the `declust.visualize.declust_marker_boxplot` function to display a boxplot of a marker gene's expression across different cell types.

```python
   declust.visualize.declust_marker_boxplot(sc_anndata, sc_marker_anndata, gene_idx)
```

b. **Spatial Visualization of Deconvolution Results**

Use the `declust.visualize.declust_results_visualize` function to generate a deconvolution result visualization:

```python
   declust.visualize.declust_results_visualize(st_anndata, sc_marker_anndata, decon_result_df,
                                             coords, gene_idx, cell_type_idx)
```

## Input Files  

DECLUST requires several input files for processing and analysis. Ensure these files are located in the `data` folder as specified:

| File Name            | Description                                                                                 |
|----------------------|---------------------------------------------------------------------------------------------|
| **sc_adata.h5ad**    | This file contains single-cell RNA sequencing (scRNA-seq) data, including gene expression profiles and cell type annotations.|
| **st_adata.h5ad**    | This file contains spatial transcriptomics (ST) data, including spatial gene expression profiles and spot coordinates. |

## Example Data Download  

To access the example data:

1. Download the `data.zip` file from [Google Drive](https://drive.google.com/uc?export=download&id=1LrSQYf1_IqQzxx7GeJrbBsEyuLLHHERC).


## License  

GNU General Public License v3.0