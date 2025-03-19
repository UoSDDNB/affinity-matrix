def plot_gene_expression_histogram(celltype_gene_matrix, output_path):
    """
    Plots a histogram of gene expression values from a cell type gene matrix.
    
    Parameters:
    celltype_gene_matrix (pd.DataFrame): A DataFrame containing gene expression values.
    output_path (str): Path to save the generated histogram plot.
    """
    plt.hist(celltype_gene_matrix.values.flatten(), bins=100)
    plt.xlim(0, 30)
    plt.ylim(0, 2000)
    plt.title("Histogram of gene expression")
    plt.xlabel("Gene expression")
    plt.ylabel("Frequency")
    plt.savefig(output_path)
    plt.clf()