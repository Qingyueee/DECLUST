#  <img src="./logo.png" align="left" height="150" /></a>

<strong>DECLUST</strong> is a Python package developed to identify spatially coherent clusters of spots by integrating gene expression profiles with spatial coordinates in spatial transcriptomics data. It also enables accurate estimation of cell-type compositions within each cluster. The recommended starting point for using DECLUST is to follow the provided <a href="https://github.com/Qingyueee/DECLUST/blob/main/tutorial.ipynb" target="_blank">**tutorial**</a>.


<br> 

## 🌟 Features

 **Spatially-aware clustering:** Combines gene expression and spatial coordinates.

 **Robust deconvolution:** Aggregates signals over clusters to enhance cell type detection.

 **Easy to install:** Available via pip.

 **Visualization:** Includes modules for visualizing clustering and marker gene expression.

## ⏬ Installation

We recommend using a separate conda environment in a Linux environment:

- Create a conda environment and install the DECLUST package

```bash
   conda create -n declust_env python=3.9
   conda activate declust_env

   pip install declust
```
- Following dependencies are required to installed in advanace: scanpy, rpy2, and R version >= 4.3 with dplyr R-packages. These dependencies can be installed using the [`install_dependencies.sh` script](https://github.com/Qingyueee/DECLUST/blob/main/install_dependencies.sh):

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


## ⚙️ DECLUST Pipeline

The pipeline consists of 5 modular stages. You can flexibly run any stage. Before start running the pipeline, place your input files(.h5ad) in the data/ folder first.
The main workflow consists of:


- Construction of Reference Matrix from Annotated Single-Cell Transcriptomic Data
  
```bash
   python main.py --stage marker
```

- Identification of spatial clusters of spots from ST data

```bash
   python main.py --stage cluster --visualize_dbscan True
```

- Pseudo-Bulk Generation
  
```bash
   python main.py --stage pseudo_bulk
```
- Deconvolution
  
```bash
   python main.py --stage deconv
```

- Visualization
  
```bash
   python main.py --stage visualize
```

## License  

GNU General Public License v3.0
