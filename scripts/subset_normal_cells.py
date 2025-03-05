# load the necessary libraries
import scanpy

# Load the h5ad file
adata = scanpy.read('d72106bd-d03b-45e6-a0fa-ca2a831ef092.h5ad')

# view the headings in the metadata
# print the col names
print(adata.obs.columns)
# head the "disease" column
print(adata.obs['disease'].head())

# subset to only include the normal cells
bdata = adata[adata.obs['disease'] == 'normal']

# check the shape of the new data 
# compared to the original data
print(adata.shape)
print(bdata.shape)

# # save the new adata object
adata.write('data/normal.h5ad')