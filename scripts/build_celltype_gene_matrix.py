# load the necessary libraries
import scanpy
import pandas

# load the new adata object
adata = scanpy.read('/scratch/mtn1n22/affinity-matrix/data/human_cell_landscape.h5ad')

# extract the counts matrix 
# build a new matrix that is celltype x genes
# find the celltype information in the metadata

# for every celltype, find the average expression of each gene 
# from all the cells of that celltype
# this will give us a celltype x genes

# get the unique celltypes
celltypes = adata.obs['cell_type'].unique()
# create an empty matrix
celltype_gene_matrix = pandas.DataFrame(index=adata.var["feature_name"], columns=celltypes)
for celltype in celltypes:
    # subset the data to only include the cells of the current celltype
    cdata = adata[adata.obs['cell_type'] == celltype]
    # calculate the average of each gene
    avg = cdata.X.mean(axis = 0).A1
    # add the average to the matrix
    celltype_gene_matrix[celltype] = avg

# print the shape of the matrix
print(celltype_gene_matrix.shape)  

# remove the rows without a gene symbol feature name
celltype_gene_matrix = celltype_gene_matrix[~celltype_gene_matrix.index.astype(str).str.startswith("ENSG00000")]

# save the matrix as a csv file
celltype_gene_matrix.to_csv("/scratch/mtn1n22/affinity-matrix/output/celltype_gene_matrix.csv", index = True)

# # build another matrix that uses the raw counts data
# # get the unique celltypes
# celltypes = adata.obs['cell_type'].unique()
# # create an empty matrix
# RAW_celltype_gene_matrix = pandas.DataFrame(index=adata.var["feature_name"], columns=celltypes)
# for celltype in celltypes:
#     # subset the data to only include the cells of the current celltype
#     cdata = adata[adata.obs['cell_type'] == celltype]
#     # calculate the average of each gene
#     avg = cdata.raw.X.mean(axis = 0).A1
#     # add the average to the matrix
#     RAW_celltype_gene_matrix[celltype] = avg

# # print the shape of the matrix
# print(RAW_celltype_gene_matrix.shape)  

# # save the matrix as a csv file
# RAW_celltype_gene_matrix.to_csv("/scratch/mtn1n22/affinity-matrix/output/RAW_celltype_gene_matrix.csv", index = True)



######################################################################
# TO DO: perhaps in later scripts

# binarize the matrix in python !? what a pain... 
# maybe extract the matrix, convert to R and use the binarize function in R.
# done.


# convert the counts matrix from cells x genes to cell_type x genes
# done

# convert the cell type counts matrix to cell type by function matrix

# using the celltype by function matrix produce a cell type by cell type matrix
# where the values indicate a similarity between the cell types
# perhaps use the cosine similarity? or Euclidean, Manhattan, etc? 
# no need for k-nn? distance is enough?

