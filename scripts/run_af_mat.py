### LOAD PACKAGES
import os
import scanpy as sc
import json
import pandas as pd
import matplotlib.pyplot as plt

### SET DIRECTORY
os.chdir('/home/mtn1n22/affinity-matrix/')

# run the functions from the functions scripts
from scripts.build_cellType_function_mat import build_ct_fun_mat
from scripts.build_cellType_gene_matrix import build_ct_gene_mat
from scripts.PCA import plot_pca
from scripts.query_celltype_gene_matrix import query_celltype_gene_matrix
from scripts.query_celltype_gene_matrix_binarized import query_celltype_gene_matrix_binarized
from scripts.query_cellType_function_matrix import query_celltype_function_matrix

### LOAD THE DATA
adata = sc.read_h5ad('data/human_cell_landscape.h5ad')

### FILTER THE DATA
#subset to only include the normal cells
#adata = adata[adata.obs['disease'] == 'normal']

# # Filter celltypes with >= 1000 cells and exclude those containing "Fetal"
adata = adata[adata.obs['author_cell_type'].value_counts()[adata.obs['author_cell_type']].ge(1000)]
adata = adata[~adata.obs['author_cell_type'].str.contains('Fetal')]

### BUILD THE CELLTYPE GENE EXPRESSION MATRIX
# must be run on a login node.
ct_gene_mat = build_ct_gene_mat(adata)
# Save the celltype-gene matrix to a CSV file
# ct_gene_mat.to_csv('output/celltype_gene_matrix.csv')
# load the celltype-gene matrix as a pd dataframe
# ct_gene_mat2 = pd.read_csv('output/celltype_gene_matrix.csv', index_col=0)

# query the celltype-gene matrix
query_celltype_gene_matrix(ct_gene_mat)

#### remove genes with low variance (e.g. variance < mean)
ct_gene_mat = ct_gene_mat[ct_gene_mat.var(axis=1) >= ct_gene_mat.var(axis=1).mean()]
#ct_gene_mat = ct_gene_mat[ct_gene_mat.var(axis=1) >= 0.2]


# consider not binarizing, instead normalizing the data?

### BINARIZE THE CELLTYPE GENE EXPRESSION MATRIX
# calculate the mean expression for the sum of each genes
bincut = ct_gene_mat.values.flatten().mean()
ct_gene_mat_bin = (ct_gene_mat > bincut).astype(int)
# # save the binarized matrix to a CSV file
# ct_gene_mat_bin.to_csv('output/celltype_gene_matrix_binarized.csv')
# # load the binarized celltype-gene matrix as a pd dataframe
# ct_gene_mat_bin = pd.read_csv('output/celltype_gene_matrix_binarized.csv', index_col=0)

# query the celltype-gene-bin matrix
query_celltype_gene_matrix_binarized(ct_gene_mat_bin)

## identify celltypes which express fewer than 100 genes, remove those columns from the matrix
ct_gene_mat_bin = ct_gene_mat_bin.loc[:, ct_gene_mat_bin.sum(axis=0) >= 100]

## identify the genes which have a variance below the mean and remove them
#ct_gene_mat_bin = ct_gene_mat_bin[ct_gene_mat_bin.var(axis=1) >= ct_gene_mat_bin.var(axis=1).mean()]

### BUILD THE CELLTYPE FUNCTION MATRIX
# load the go terms
with open('data/c5.go.v2024.1.Hs.json', 'r') as f:
    gene_functions = json.load(f)  # gene_functions is a dictionary

# build the celltype-function matrix
ct_fun_mat = build_ct_fun_mat(ct_gene_mat_bin, gene_functions)

# # save the celltype-function matrix to a CSV file
# ct_fun_mat.to_csv('output/celltype_function_matrix.csv')
# # load the celltype-function matrix as a pd dataframe
# ct_fun_mat = pd.read_csv('output/celltype_function_matrix.csv', index_col=0)

## query the celltype-function matrix
query_celltype_function_matrix(ct_fun_mat)

# ## remove functions with low variance
# # # remove functions with variance less than x* the mean
# ct_fun_mat = ct_fun_mat[ct_fun_mat.var(axis=1) >= 2 * ct_fun_mat.var(axis=1).mean()]
ct_fun_mat = ct_fun_mat[ct_fun_mat.var(axis=1) >= 0.002]

### PCA
# run pca on the celltype-function matrix
plot_pca(ct_fun_mat)


# print all celltypes
print(ct_fun_mat.columns.tolist())