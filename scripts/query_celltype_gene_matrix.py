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
    # Number of cell types that express each gene
    sum_exp_per_gene = ct_gene.sum(axis=1)
    
    # Number of genes expressed by each cell type
    sum_exp_per_ct = ct_gene.sum(axis=0)


    # Plot the number of cell types that express each gene
    sum_exp_per_gene.plot(kind='hist', bins=100, color='orange', alpha=0.9)
    plt.title('Sum of expression per gene')
    plt.xlabel('Sum of expression')
    plt.ylabel('Frequency')
    plt.tight_layout()
    plt.savefig('/scratch/mtn1n22/affinity-matrix/output/figures/query/celltype_gene_matrix/sum_exp_per_gene.png')
    plt.clf()


    # Plot the Sum of expression by CellType
    sum_exp_per_ct.sort_values(ascending=False).plot(kind='barh', color='purple', alpha=0.7, figsize=(12, 12))
    plt.title('Sum of expression by CellType')
    plt.xlabel('Sum of expression')
    plt.ylabel('CellType')
    plt.tight_layout()
    plt.savefig('/scratch/mtn1n22/affinity-matrix/output/figures/query/celltype_gene_matrix/sum_exp_per_cell_type.png')
    plt.clf()

# # load binarized celltype-gene matrix
# ct_gene_mat_bin = pd.read_csv('output/celltype_gene_matrix_binarized.csv', index_col=0)
# num_ct_expressing_gene = ct_gene_mat_bin.sum(axis=1)
# num_genes_expressed_by_ct.sort_values(ascending=False)

# # Plot the Sum of expression by CellType again,
# # this time ordered to match num_genes_expressed_by_ct.sort_values(ascending=False)
# sum_exp_per_ct = sum_exp_per_ct.loc[num_genes_expressed_by_ct.sort_values(ascending=False).index]
# sum_exp_per_ct.plot(kind='barh', color='purple', alpha=0.7, figsize=(12, 12))
# plt.title('Sum of expression by CellType (ordered by number of genes expressed)')
# plt.xlabel('Sum of expression')
# plt.ylabel('CellType')
# plt.tight_layout()
# plt.savefig('/scratch/mtn1n22/affinity-matrix/output/figures/sum_exp_per_cell_type_ordered.png')
# plt.clf()
