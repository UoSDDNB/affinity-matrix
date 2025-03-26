### LOAD PACKAGES
import os
import scanpy as sc
import json
import pandas as pd

# run the functions from the functions scripts
from scripts.build_cellType_function_mat import build_ct_fun_mat
from scripts.build_cellType_gene_matrix import build_ct_gene_mat
from scripts.plot_expression import plt_gene_exp_hist
from scripts.plot_PCA_UMAP import plot_pca, plot_umap

### SET DIRECTORY
os.chdir('/home/mtn1n22/affinity-matrix/')

### LOAD THE DATA
adata = sc.read_h5ad('data/human_cell_landscape.h5ad')

### FILTER THE DATA
#subset to only include the normal cells
#adata = adata[adata.obs['disease'] == 'normal']

### BUILD THE CELLTYPE GENE EXPRESSION MATRIX
# must be run on a login node.
ct_gene_mat = build_ct_gene_mat(cdata)
# Save the celltype-gene matrix to a CSV file
#ct_gene_mat.to_csv('output/celltype_gene_matrix.csv')

### PLOT THE GENE EXPRESSION
# run the plot_expression.py script
plt_gene_exp_hist(ct_gene_mat, 'output/figures/gene_expression_hist.png', bins=100)

plt_gene_exp_hist(ct_gene_mat,
                  'output/figures/gene_expression_hist.png',
                  xlim=3,
                  ylim=8000,
                  bins=1000)

# print the maximum gene expression value
print(ct_gene_mat.values.max())
# condsider not binarizing, instead normailizing the data?

### BINARIZE THE CELLTYPE GENE EXPRESSION MATRIX
ct_gene_mat_bin = (ct_gene_mat > 0).astype(int)
# save the binarized matrix to a CSV file
ct_gene_mat_bin.to_csv('output/celltype_gene_matrix_binarized.csv')

### BUILD THE CELLTYPE FUNCTION MATRIX
# load the go terms
with open('data/c5.go.v2024.1.Hs.json', 'r') as f:
    gene_functions = json.load(f)  # gene_functions is a dictionary

# build the celltype-function matrix
ct_fun_mat = build_ct_fun_mat(ct_gene_mat_bin, gene_functions)
# save the celltype-function matrix to a CSV file
ct_fun_mat.to_csv('output/celltype_function_matrix.csv')

# load the celltype-function matrix as a pd dataframe
ct_fun_mat = pd.read_csv('output/celltype_function_matrix.csv', index_col=0)

### PCA
# run pca on the celltype-function matrix
celltype_pca = plot_pca(ct_fun_mat)

### UMAP
plot_umap(celltype_pca, ct_fun_mat, num_pcs=12)

