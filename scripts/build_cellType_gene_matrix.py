import pandas as pd

def build_ct_gene_mat(adata, convert_to_gene_symbols=True):
    """
    Computes the average gene expression for each cell type in an AnnData object
    and returns the result as a DataFrame.

    Parameters:
        adata (AnnData): AnnData object containing the single-cell data.
        convert_to_gene_symbols (bool): Whether to convert gene codes to gene symbols.
    Returns:
        pd.DataFrame: DataFrame containing the average gene expression for each cell type.
    """
      
    # Create an empty dataframe
    celltype_gene_matrix = pd.DataFrame(index=adata.var.index, columns=adata.obs['cell_type'].unique())
    
    for celltype in adata.obs['cell_type'].unique():
        # Subset data to only include the cells of the current cell type
        cdata = adata[adata.obs['cell_type'] == celltype]
        
        # Calculate the average expression of each gene
        avg_expression = cdata.X.mean(axis=0).A1
        
        # Store in the celltype-gene matrix
        celltype_gene_matrix[celltype] = avg_expression
    
    # Convert gene codes to gene symbols
    if not convert_to_gene_symbols:
        return celltype_gene_matrix

    from mygene import MyGeneInfo
    mg = MyGeneInfo()
    genes = list(celltype_gene_matrix.index)
    mapping = mg.querymany(genes, scopes="ensembl.gene", fields="symbol", species="human")

    # Convert to dictionary
    gene_code_to_symbol = {g["query"]: g.get("symbol", g["query"]) for g in mapping}
    
    celltype_gene_matrix.rename(index=gene_code_to_symbol, inplace=True)

    # Remove rows without gene symbols (assumed to start with 'ENSG00000')
    celltype_gene_matrix = celltype_gene_matrix[~celltype_gene_matrix.index.astype(str).str.startswith("ENSG00000")]


    return celltype_gene_matrix
