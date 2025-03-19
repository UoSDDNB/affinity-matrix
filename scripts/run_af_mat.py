### LOAD PACKAGES
import os
import scanpy as sc
import pandas as pd
import json

### SET DIRECTORY
os.chdir('/home/mtn1n22/affinity-matrix/')

### LOAD THE DATA
adata = sc.read_h5ad('data/human_cell_landscape.h5ad')

### FILTER THE DATA
#subset to only include the normal cells
ndata = adata[adata.obs['disease'] == 'normal']


### BUILD THE CELLTYPE EXPRESSION MATRIX



### PLOT THE GENE EXPRESSION


### BINARIZE THE CELLTYPE GENE EXPRESSION MATRIX
celltype_gene_matrix[celltype_gene_matrix > 0] = 1


### BUILD THE CELLTYPE FUNCTION MATRIX


### PCA



### UMAP
