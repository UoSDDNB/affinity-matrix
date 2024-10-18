# load the necessary libraries
import scanpy
import pandas

# # Load the h5ad file
# adata = scanpy.read('d72106bd-d03b-45e6-a0fa-ca2a831ef092.h5ad')

# # view the headings in the metadata
# # print the col names
# print(adata.obs.columns)
# # head the "disease" column
# print(adata.obs['disease'].head())

# # subset to only include the normal cells
# bdata = adata[adata.obs['disease'] == 'normal']

# # check the shape of the new data 
# # compared to the original data
# print(adata.shape)
# print(bdata.shape)

# # # save the new adata object
# # adata.write('data/normal.h5ad')

# load the new adata object
adata = scanpy.read('data/normal.h5ad')


# extract the counts matrix 
# buld a new matrix that is celltype x genes
# find the celltype information in the metadata

# for every celltype, find the average expression of each gene 
# from all the cells of that celltype
# this will give us a celltype x genes

# get the unique celltypes
celltypes = adata.obs['cell_type'].unique()
# create an empty matrix
celltype_matrix = pandas.DataFrame()
for celltype in celltypes:
    # subset the data to only include the cells of the current celltype
    cdata = adata[adata.obs['cell_type'] == celltype]
    # calculate the average of each gene
    avg = cdata.X.mean(axis = 0).A1
    # add the average to the matrix
    celltype_matrix[celltype] = avg

#
print(celltype_matrix.shape)  


# convert the counts matrix from cells x genes to cell_type x genes

# convert the cell type counts matrix to cell type by function matrix

# using the celltype by function matrix produce a cell type by cell type matrix
# where the values indicate a similarity between the cell types
# perhaps use the cosine similarity? or Euclidean, Manhattan, etc? 
# no need for k-nn? distance is enough?