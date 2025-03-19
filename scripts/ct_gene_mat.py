def build_celltype_gene_matrix(input_file: str, output_file: str, celltype_column: str = 'cell_type', feature_name_column: str = 'feature_name'):
    """
    Computes the average gene expression for each cell type in an AnnData object
    and saves the result as a CSV file.

    Parameters:
        input_file (str): Path to the input .h5ad file.
        output_file (str): Path to save the output CSV file.
        celltype_column (str): Column name in .obs that contains cell type labels.
        feature_name_column (str): Column name in .var that contains gene names.
    """
    # Load the AnnData object
    adata = sc.read(input_file)
    
    # Get unique cell types
    celltypes = adata.obs[celltype_column].unique()
    
    # Create an empty dataframe
    celltype_gene_matrix = pd.DataFrame(index=adata.var[feature_name_column], columns=celltypes)
    
    for celltype in celltypes:
        # Subset data to only include the cells of the current cell type
        cdata = adata[adata.obs[celltype_column] == celltype]
        
        # Calculate the average expression of each gene
        avg_expression = cdata.X.mean(axis=0).A1
        
        # Store in the matrix
        celltype_gene_matrix[celltype] = avg_expression
    
    # Remove rows without gene symbols (assumed to start with 'ENSG00000')
    celltype_gene_matrix = celltype_gene_matrix[~celltype_gene_matrix.index.astype(str).str.startswith("ENSG00000")]
    
    # save the matrix
    celltype_gene_matrix.to_csv(output_file, index=True)

    return celltype_gene_matrix