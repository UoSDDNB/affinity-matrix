import pandas as pd
from matplotlib import pyplot as plt

def query_celltype_gene_matrix_binarized(ct_gene_bin):
    """
    Query the gene_ct_matrix to return the number of cell types that express each gene
    and the number of genes expressed by each cell type, generate and save histograms
    for both distributions.
    
    Args:
        gene_ct_matrix (pd.DataFrame): Matrix with genes as rows and cell types as columns.
        output_dir (str): Directory to save the output plots.
        
    Returns:
        None
    """
    # Number of cell types that express each gene
    num_ct_expressing_gene = ct_gene_bin.sum(axis=1)
    
    # Number of genes expressed by each cell type
    num_genes_expressed_by_ct = ct_gene_bin.sum(axis=0)

    # Plot the number of cell types that express each gene
    plt.close()
    num_ct_expressing_gene.plot(kind='hist', bins=44, color='green', alpha=0.7)
    plt.title('Number of CellTypes that express each gene')
    plt.xlabel('Number of CellTypes')
    plt.ylabel('Number of Genes')
    plt.tight_layout()
    plt.savefig('/scratch/mtn1n22/affinity-matrix/output/figures/query/celltype_gene_matrix_binarized/num_cell_types_expressing_gene.png')
    plt.close()

    # Plot the number of genes expressed by each cell type in a bar plot
    num_genes_expressed_by_ct.sort_values(ascending=False).plot(kind='barh', color='blue', alpha=0.7, figsize=(12, 12))
    plt.title('Number of genes expressed by each CellType')
    plt.xlabel('Number of Genes')
    plt.ylabel('CellType')
    plt.tight_layout()
    plt.savefig('/scratch/mtn1n22/affinity-matrix/output/figures/query/celltype_gene_matrix_binarized/num_genes_expressed_by_cell_type.png')
    plt.close()

    # plot a histogram of the variance of each gene (row)
    ct_gene_bin.var(axis=1).sort_values(ascending=True).plot(kind='hist', bins=100, color='purple', alpha=0.7)
    plt.title('Variance of each gene (Binarized)')
    plt.xlabel('Variance')
    plt.ylabel('Frequency')
    plt.tight_layout()
    # draw a line at the mean
    mean = ct_gene_bin.var(axis=1).mean()
    plt.axvline(mean, color='red', linestyle='dashed', linewidth=1)
    plt.savefig('/scratch/mtn1n22/affinity-matrix/output/figures/query/celltype_gene_matrix_binarized/variance_of_each_gene.png')
    plt.close()


    