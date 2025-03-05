# Build a matrix with celltype by function
# load the necessary libraries
import scanpy
import pandas

# load the celltype_gene matrix
celltype_gene_matrix = pandas.read_csv("output/celltype_gene_matrix.csv", index_col = 0)

# load the gene_function json file
import json
with open("data/c5.go.v2024.Hs.json") as f:
    gene_function = json.load(f)

# Number of GO terms
print(len(gene_function))  

#there are 10454 terms.

## create a new matrix that is celltype x function
# the columns will be the celltypes
# the rows will be the functions
# the values will be the number of genes in that function that are expressed in that celltype.


