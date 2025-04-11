#  <img src="./logo.png" align="left" height="150" /></a>

<strong>DECLUST</strong> is a Python package developed to identify spatially coherent clusters of spots by integrating gene expression profiles with spatial coordinates in spatial transcriptomics data. It also enables accurate estimation of cell-type compositions within each cluster.


<br> 

## 🌟 Features

 **Spatially-aware clustering:** Combines gene expression and spatial coordinates.

 **Robust deconvolution:** Aggregates signals over clusters to enhance cell type detection.

 **Easy to install:** Available via pip.

 **Visualization:** Includes modules for visualizing clustering and marker gene expression.

## ⏬ Installation

We recommend using a separate Conda environment. Information about Conda and how to install it can be found in the [anaconda webpage](https://www.anaconda.com/docs/getting-started/miniconda/main).

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

The DECLUST package has been installed successfully on Operating systems:

- macOS Sequoia 15.3.2
- SUSE Linux Enterprise Server 15 SP5 (Dardel HPC system)

## 📊 Data Input

DECLUST uses `.h5ad` files, which are [AnnData](https://anndata.readthedocs.io/en/latest/) objects commonly used for storing annotated data matrices in single-cell and spatial transcriptomics analysis.

Each `.h5ad` file includes:

#### **sc_adata.h5ad** (Single-cell RNA-seq data)
- `.X`: Gene expression matrix (cells × genes)
- `.obs`: Cell type annotation of single cells

#### **st_adata.h5ad** (Spatial transcriptomics data)
- `.X`: Spatial gene expression matrix (spots × genes)
- `.obs`: Spots coordinates

> 💡 Both datasets should **originate from the same tissue** and have **overlapping gene sets** to ensure proper implementation of DECLUST.

## 🔗 Example Data Download  

- Download the [Real Data Example](https://drive.google.com/uc?export=download&id=1LrSQYf1_IqQzxx7GeJrbBsEyuLLHHERC). 
   
- Download the [Simulation Data Example](https://drive.google.com/uc?export=download&id=1VY_vIuZalCBe2IhNCNBSQwo5m5Da8aFw).


## ⚙️ Usage
DECLUST can be embedded into python scripts or used independently as a tool.
A guide of how to use it in python scripts is provided in this <a href="https://github.com/Qingyueee/DECLUST/blob/main/tutorial.ipynb" target="_blank">**tutorial**</a>. In this section, we introduce how to use it as a bioinformatics pipeline.

Run the pipeline using the following command:

```bash
python declust.py --module <module_name> [other options]
```

- Available Modules

| Module       | Description                                                                    |
|--------------|--------------------------------------------------------------------------------|
| `marker`     | Construction of Reference Matrix from Annotated Single-Cell Transcriptomic Data|
| `cluster`    | Identification of spatial clusters of spots from ST data                       |
| `pseudo_bulk`| Generate pseudo-bulk ST profiles per cluster                                   |
| `deconv`     | Run deconvolution by Ordinary Least Squares                                    |
| `visualize`  | Visualize markers or deconvolution results                                     |

Type `python declust.py --help` in the terminal to see a list of available commands.

## 🧬 DECLUST pipeline

1. Download DECLUST:

```bash
   wget https://github.com/Qingyueee/DECLUST/archive/refs/tags/0.1.1.tar.gz
   tar -xvf 0.1.1.tar.gz
```

2. Unpack data:

```bash
   cd DECLUST-0.1.1
   unzip data.zip
```
3. Marker gene selection:

```bash
   python declust.py --module marker \
   --celltype_col \
   --sample_col
```

Outputs:

- `sc_data_overlapped.csv` and `sc_label.csv` in the `data/` folder

- `marker_genes.csv` in the `results/` folder

4. Clustering:

```bash
   python declust.py --module cluster
```

Performs Hierarchical Clustering → DBSCAN → Seeded Region Growing (SRG). Saves:

- `srg_df.csv` and clustering plots in `results/`

5. Deconvolution:

```bash
   python declust.py --module deconv
```

Performs OLS-based deconvolution and outputs:

- `DECLUST_result.csv` in `results/`

You can run each step individually or execute the entire pipeline by running the deconvolution script.

To export pseudo-bulk profiles for external methods:

```bash
   python declust.py --module pseduo_bulk
```

- Generates `pseudo_bulk.csv` in the `results/` folder.

#### 💡 Custom Marker Genes

Users can provide their own marker gene list in one of two formats:

- **CSV file** containing two columns:
   - `Gene`: gene names  
   - `maxgroup`: corresponding cell type annotations

```bash
   --custom_marker_genes file_path
```

- **Comma-separated gene list**, along with a corresponding **comma-separated list of cell types**:

```bash
   --custom_marker_genes "DCN, LUM, C1S, AGR2, PPDPF, ..."
   --custom_marker_celltype "CAFs, CAFs, CAFs, Cancer Epithelial, Cancer Epithelial, ..."
```
> ⚠️ The provided marker genes and cell type annotations must exist in the single-cell dataset.

## 📬 Quick example to run DECLUST on a simulated data

```bash
# 1. Download DECLUST
   wget https://github.com/Qingyueee/DECLUST/archive/refs/tags/0.1.1.tar.gz
   tar -xvf 0.1.1.tar.gz
   cd DECLUST-0.1.1

# 2. Configuring environment and install dependencies
   conda create -n declust_env python=3.9
   conda activate declust_env
   pip install declust
   sh install_dependencies.sh

# 3. Download and unpack simulated data
   wget "https://drive.usercontent.google.com/download?id=1VY_vIuZalCBe2IhNCNBSQwo5m5Da8aFw&export=download&authuser=0&confirm=t&uuid=93730baf-2a12-49d7-b475-ab715a3644c3&at=APcmpow759exSs6opQk4zSMVbjXf%3A1744370330609" -O simulation_data.zip
   unzip simulation_data.zip

# 4. Run pipeline - it may take about 2 minutes to complete on a personal computer
   python declust.py --module deconv \
      --data_dir simulation_data \
      --results_dir simulation_results \
      --sc_file sc_adata_200_per_celltype.h5ad \
      --st_file st_simu_adata.h5ad \
      --celltype_col celltype_major \
      --sample_col Patient

# 5. Results visulization
   python declust.py --module visualize \
      --data_dir simulation_data \
      --results_dir simulation_results \
      --sc_file sc_adata_200_per_celltype.h5ad \
      --st_file st_simu_adata.h5ad
```

## 📁 Output Structure

```bash
   project/
   │
   ├── data/
   │   ├── sc_adata_overlapped.h5ad
   │   ├── sc_labels.csv
   │   └── ...
   │
   ├── results/
   │   ├── marker_genes.csv
   │   ├── srg_df.csv
   │   ├── pseudo_bulk.csv
   │   ├── DECLUST_result.csv
   │   └── [visualization plots]
```

## License  

GNU General Public License v3.0
