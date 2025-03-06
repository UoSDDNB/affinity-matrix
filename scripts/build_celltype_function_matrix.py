# Build a matrix with celltype by function
# load the necessary libraries
import pandas

# load the celltype_gene matrix
celltype_gene_matrix = pandas.read_csv("/scratch/mtn1n22/affinity-matrix/output/celltype_gene_matrix_binarized.csv", index_col = 0)

# head the gene names
print(celltype_gene_matrix.head())

# load the gene_function json file
import json
with open("/scratch/mtn1n22/affinity-matrix/data/c5.go.v2024.1.Hs.json") as f:
	gene_function = json.load(f)

# Number of GO terms
print(len(gene_function))  
#there are 10454 terms.


## create a new matrix that is celltype x function
# the columns will be the celltypes from the celltype_gene_matrix
# the rows will be the functions from the gene_function
# the values will be the number of genes in that function that are expressed in that celltype.

# Create a DataFrame to store the result
celltype_function_matrix = pandas.DataFrame(
	index=gene_function.keys(),
	columns=celltype_gene_matrix.keys(),
	dtype=int)

# Precompute expressed genes as sets for each cell type
celltype_expressed_sets = {
    celltype: set(expressed_genes.index[expressed_genes == 1.0])
    for celltype, expressed_genes in celltype_gene_matrix.items()
}

# Get the total number of functions
total_functions = len(gene_function)
# Loop through each function
for i, (function, function_data) in enumerate(gene_function.items(), start=1):
	if i % 500 == 0 or i == total_functions:
		print(f"Processing function {i}/{total_functions}...")
	genes = set(function_data['geneSymbols'])  # Convert genes to a set
    # Compute intersection with each cell type's expressed genes
	for celltype, expressed_gene_set in celltype_expressed_sets.items():
		celltype_function_matrix.loc[function, celltype] = len(genes & expressed_gene_set)

# Save the matrix
celltype_function_matrix.to_csv("/scratch/mtn1n22/affinity-matrix/output/celltype_function_matrix.csv", index=True)

