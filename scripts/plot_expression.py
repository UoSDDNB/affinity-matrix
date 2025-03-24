import matplotlib.pyplot as plt

def plt_gene_exp_hist(celltype_gene_matrix, output_path, xlim=30, ylim=2000, bins=100):
    """
    Plots a histogram of gene expression values from a cell type gene matrix.
    
    Parameters:
    celltype_gene_matrix (pd.DataFrame): A DataFrame containing gene expression values.
    output_path (str): Path to save the generated histogram plot.
    """
    plt.hist(celltype_gene_matrix.values.flatten(), bins=bins)
    plt.xlim(0, xlim)
    plt.ylim(0, ylim)
    plt.title("Histogram of gene expression")
    plt.xlabel("Gene expression")
    plt.ylabel("Frequency")
    plt.savefig(output_path)
    plt.clf()