import pandas as pd


def build_ct_fun_mat(celltype_gene_matrix, gene_functions):
    """
    Constructs a celltype-function matrix from a celltype-gene matrix and a gene-function JSON file.
    
    Parameters:
    - celltype_gene_matrix (pd.DataFrame): The celltype-gene matrix.
    - gene_functions (dict): The gene-function mapping.
    
    Returns:
    - pd.DataFrame: The computed celltype-function matrix.
    """
     
    # Initialize the celltype-function matrix
    celltype_function_matrix = pd.DataFrame(
        index=gene_functions.keys(),
        columns=celltype_gene_matrix.keys(),
        dtype=int
    )
    
    # Precompute expressed genes as sets for each cell type
    celltype_expressed_sets = {
        celltype: set(expressed_genes.index[expressed_genes == 1.0])
        for celltype, expressed_genes in celltype_gene_matrix.items()
    }
    
    ## COMPUTE FUNCTION-CELLTYPE ASSOCIATIONS
    # Count the genes in each function that are expressed in each cell type

    # Get the total number of functions
    total_functions = len(gene_functions)
    # Loop through each function
    for i, (function, function_data) in enumerate(gene_functions.items(), start=1):
        if i % 500 == 0 or i == total_functions:
            print(f"Processing function {i}/{total_functions}...")
        genes = set(function_data['geneSymbols'])  # Convert genes to a set
        # Compute intersection with each cell type's expressed genes
        for celltype, expressed_gene_set in celltype_expressed_sets.items():
            celltype_function_matrix.loc[function, celltype] = len(genes & expressed_gene_set)

    ## Normalize the celltype-function matrix by dividing by the number of genes in each function
    ## calculate the number of genes in each function using gene_functions
    num_genes_in_function = {function: len(function_data['geneSymbols']) for function, function_data in gene_functions.items()}
    
    # Divide each celltype-function matrix entry by the number of genes in the function
    celltype_function_matrix = celltype_function_matrix.div(num_genes_in_function.values(), axis=0)
    
    return celltype_function_matrix
