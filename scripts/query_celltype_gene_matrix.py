import pandas as pd
from matplotlib import pyplot as plt

def query_celltype_gene_matrix(ct_gene):
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
    # Sum of expression across each gene
    sum_exp_per_gene = ct_gene.sum(axis=1)
    
    # Sum of expression across each cell type
    sum_exp_per_ct = ct_gene.sum(axis=0)

    # Plot the sum of expression per gene
    sum_exp_per_gene.plot(kind='hist', bins=100, color='orange', alpha=0.9)
    plt.title('Sum of expression per gene')
    plt.xlabel('Sum of expression')
    plt.ylabel('Frequency')
    plt.tight_layout()
    # draw a line at the mean
    mean = sum_exp_per_gene.mean()
    plt.axvline(mean, color='red', linestyle='dashed', linewidth=1)
    plt.savefig('/scratch/mtn1n22/affinity-matrix/output/figures/query/celltype_gene_matrix/sum_exp_per_gene.png')
    plt.close()

    # Plot the Sum of expression by CellType
    sum_exp_per_ct.sort_values(ascending=False).plot(kind='barh', color='purple', alpha=0.7, figsize=(12, 12))
    plt.title('Sum of expression by CellType')
    plt.xlabel('Sum of expression')
    plt.ylabel('CellType')
    plt.tight_layout()
    plt.savefig('/scratch/mtn1n22/affinity-matrix/output/figures/query/celltype_gene_matrix/sum_exp_per_cell_type.png')
    plt.close()

    # Plot a histogram of the variance of each gene (row)
    ct_gene.var(axis=1).sort_values(ascending=True).plot(kind='hist', bins=100, color='purple', alpha=0.7)
    plt.title('Variance of each gene')
    plt.xlabel('Variance')
    plt.ylabel('Frequency')
    plt.tight_layout()
    # draw a line at the mean
    mean = ct_gene.var(axis=1).mean()
    plt.axvline(mean, color='red', linestyle='dashed', linewidth=1)
    plt.savefig('/scratch/mtn1n22/affinity-matrix/output/figures/query/celltype_gene_matrix/variance_of_each_gene.png')
    plt.close()

    # Plot a histogram of all values in the DataFrame
    pd.Series(ct_gene.values.flatten()).plot(kind='hist', bins=100, color='darkgreen', alpha=0.7)
    plt.title('Histogram Gene Expression Values')
    plt.xlabel('Expression value')
    plt.ylabel('Frequency')
    plt.tight_layout()
    # draw a line at the mean
    mean = ct_gene.values.flatten().mean()
    plt.axvline(mean, color='red', linestyle='dashed', linewidth=1)
    plt.savefig('/scratch/mtn1n22/affinity-matrix/output/figures/query/celltype_gene_matrix/Gene_EXP_hist.png')
    plt.close()

# # load binarized celltype-gene matrix
# ct_gene_mat_bin = pd.read_csv('output/celltype_gene_matrix_binarized.csv', index_col=0)
# num_ct_expressing_gene = ct_gene_mat_bin.sum(axis=1)
# num_genes_expressed_by_ct.sort_values(ascending=False)

# # Plot the Sum of expression by CellType again,
# # this time ordered to match num_genes_expressed_by_ct.sort_values(ascending=False)
# sum_exp_per_ct = sum_exp_per_ct.loc[num_genes_expressed_by_ct.sort_values(ascending=False).index]
# suplt.close()_ct.plot(kind='barh', color='purple', alpha=0.7, figsize=(12, 12))
# plt.title('Sum of expression by CellType (ordered by number of genes expressed)')
# plt.xlabel('Sum of expression')
# plt.ylabel('CellType')
# plt.tight_layout()
# plt.savefig('/scratch/mtn1n22/affinity-matrix/output/figures/sum_exp_per_cell_type_ordered.png')
# plt.clf()
plt.close()