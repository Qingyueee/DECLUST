import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def declust_marker_boxplot(sc_adata, sc_marker_gene, gene_index):

    """
    Create a boxplot showing the expression of a specific marker gene across different cell types.

    This function extracts the expression values for a specified gene from single-cell data (sc_adata)
    and plots a boxplot of its expression across the major cell types defined in the dataset. The gene
    to plot can be specified by its index (integer) or its name (string).

    Parameters:
        sc_adata (anndata.AnnData):
            Single-cell RNA-seq data in an AnnData object. The object should have:
                - `obs['celltype_major']`: A column indicating the major cell type for each cell.
                - Gene expression data accessible via `sc_adata.to_df()`, where columns are gene names.
                
        sc_marker_gene (pandas.DataFrame):
            A DataFrame of marker genes, where the index contains gene names.
            
        gene_index (int or str):
            The identifier for the gene to plot. If an integer, it is used as the index position
            in sc_marker_gene's index. If a string, it is interpreted as the gene name. If the string
            is not found in the DataFrame index, the function will attempt to locate it in the DataFrame values.
    """

    box_color = '#92B5CA'

    if isinstance(gene_index, int):
        gene_name = sc_marker_gene.index[gene_index]
    elif isinstance(gene_index, str):
        if gene_index in sc_marker_gene.index:
            gene_name = gene_index
        elif gene_index in sc_marker_gene.values:
            gene_name = sc_marker_gene[sc_marker_gene == gene_index].index[0]
        else:
            raise ValueError(f"Gene name '{gene_index}' not found in marker_genes.")
    else:
        raise TypeError("gene_index must be either an integer or a string.")

    sc_marker_df = pd.DataFrame({
        'CellType': sc_adata.obs['celltype_major'],
        'gene_Expression': sc_adata.to_df()[gene_name]
    })

    plt.figure(figsize=(8, 6))
    sns.boxplot(
        x='CellType', 
        y='gene_Expression', 
        data=sc_marker_df,
        showfliers=False,
        hue='CellType',
        palette={cell_type: box_color for cell_type in sc_marker_df['CellType'].unique()}
    )

    plt.xlabel('')
    plt.title(gene_name)
    plt.xticks(rotation=90)

    for spine in plt.gca().spines.values():
        spine.set_linewidth(0.5)

    plt.tight_layout()
    plt.show()



def declust_results_visualize(st_adata, sc_marker_gene, DECLUST_df, coords,
                              idx=None, gene_name=None, cell_type=None,
                              agg_method='sum', save=False, save_path=None):
    
    """
    Visualize deconvolution results by plotting marker gene expression and estimated cell type proportion 
    on spatial coordinates.

    This function generates a two-panel figure:
      - The left panel displays the expression of one or a set of marker genes (aggregated using the 
        specified method) across spatial coordinates.
      - The right panel shows the estimated proportion for a specified cell type as obtained from the 
        deconvolution results (DECLUST_df).
    
    The function determines the gene(s) and cell type based on the provided arguments:
      - If 'idx' is provided (and gene_name and cell_type are None), then the gene name is taken from 
        sc_marker_gene at that index, and the corresponding cell type is taken from sc_marker_gene['maxgroup'].
      - If gene_name is not provided but cell_type is provided, all marker genes corresponding to that cell 
        type are used.
      - If gene_name is provided without cell_type, an error is raised.
    
    Parameters:
        st_adata (anndata.AnnData):
            Spatial transcriptomics data with gene expression values. Expression data is accessed via 
            st_adata.to_df(), with gene names as columns.
            
        sc_marker_gene (pandas.DataFrame):
            DataFrame of marker genes. The index contains gene names, and it must include a column 
            'maxgroup' that indicates the associated cell type for each marker gene.
            
        DECLUST_df (pandas.DataFrame):
            DataFrame containing deconvolution results. It should have a column for each cell type 
            representing the estimated proportion at each spatial spot.
            
        coords (pandas.DataFrame):
            DataFrame containing spatial coordinates for each spot. It must include at least 'x' and 'y' columns.
            
        idx (int, optional):
            Index into sc_marker_gene to select a gene and its corresponding cell type if gene_name and 
            cell_type are not provided.
            
        gene_name (str or list of str, optional):
            The name of the gene to visualize, or a list of gene names. If a list is provided, the 
            expression values will be aggregated according to agg_method.
            
        cell_type (str, optional):
            The cell type associated with the marker gene(s). Must be provided if gene_name is given.
            
        agg_method (str, optional):
            Aggregation method for multiple genes. Options are 'sum' or 'mean'. Defaults to 'sum'.
            
        save (bool, optional):
            If True, the figure is saved to file using save_path.
            
        save_path (str, optional):
            Path and filename to save the figure. If not provided and save is True, the file will be 
            saved as "{cell_type}.png" (spaces replaced with underscores).
    """

    if gene_name is None and cell_type is None and idx is not None:
        gene_name = sc_marker_gene.index[idx]
        cell_type = sc_marker_gene.iloc[idx]['maxgroup']
    elif gene_name is None and cell_type is not None:
        gene_name = sc_marker_gene[sc_marker_gene['maxgroup'] == cell_type].index.tolist()
        if len(gene_name) == 0:
            raise ValueError(f"No marker genes found for cell type '{cell_type}'.")
    elif gene_name is not None and cell_type is None:
        raise ValueError("If gene_name is given, you must also provide cell_type.")

    if isinstance(gene_name, list):
        gene_expr_df = st_adata.to_df()[gene_name]
        if agg_method == 'sum':
            gene_expr = gene_expr_df.sum(axis=1)
        elif agg_method == 'mean':
            gene_expr = gene_expr_df.mean(axis=1)
        else:
            raise ValueError("agg_method must be 'sum' or 'mean'")

        marker_genes_all = sc_marker_gene[sc_marker_gene['maxgroup'] == cell_type].index.tolist()
        if set(gene_name) == set(marker_genes_all):
            title_expr = "All markers" if agg_method == 'sum' else "mean(All markers)"
        else:
            title_expr = ' + '.join(gene_name) if agg_method == 'sum' else f'mean({", ".join(gene_name)})'
    else:
        gene_expr = st_adata.to_df()[gene_name]
        title_expr = gene_name

    fig, axs = plt.subplots(1, 2, figsize=(12, 5))

    scatter1 = axs[0].scatter(coords['x'], coords['y'], c=gene_expr, cmap='viridis', s=10,
                              vmin=gene_expr.min(), vmax=gene_expr.max())
    fig.colorbar(scatter1, ax=axs[0])
    axs[0].set_title(title_expr)

    cell_prop = DECLUST_df[cell_type]
    scatter2 = axs[1].scatter(coords['x'], coords['y'], c=cell_prop, cmap='viridis', s=10,
                              vmin=cell_prop.min(), vmax=cell_prop.max() * 1.1)
    fig.colorbar(scatter2, ax=axs[1])
    axs[1].set_title(f'{cell_type} proportion')

    for ax in axs:
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xlabel('')
        ax.set_ylabel('')

    plt.tight_layout()

    if save:
        if save_path is None:
            filename = f"{cell_type}.png".replace(" ", "_")
        else:
            filename = save_path
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {filename}")

    plt.show()
