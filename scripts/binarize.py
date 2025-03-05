# import a counts matrix, choose a cutoff, binarize it, save. 


# load the necessary libraries
import pandas
import matplotlib.pyplot as plt

# load 
celltype_gene_matrix = pandas.read_csv("/scratch/mtn1n22/affinity-matrix/output/celltype_gene_matrix.csv", index_col = 0)

# draw a histogram of gene expression from the celltype matrix
plt.hist(celltype_gene_matrix.values.flatten(), bins=100)

# set the x-axis limits
plt.xlim(0, 30)
# set the y-axis limits
plt.ylim(0, 2000)

# add a title and labels
plt.title("Histogram of gene expression")
plt.xlabel("Gene expression")
# save the plot to a file
plt.savefig("/scratch/mtn1n22/affinity-matrix/output/figures/HIST-cutoff.png")
# clear the plot
plt.clf()

## repeat for the RAW celltype matrix
# load
RAW_celltype_gene_matrix = pandas.read_csv("/scratch/mtn1n22/affinity-matrix/output/RAW_celltype_gene_matrix.csv", index_col = 0)

# draw a histogram of gene expression from the celltype matrix
plt.hist(RAW_celltype_gene_matrix.values.flatten(), bins=5000)
# set the x-axis limits
plt.xlim(0, 30)
# set the y-axis limits
plt.ylim(0, 2000)
# add a title and labels
plt.title("Histogram of RAW gene expression")
plt.xlabel("RAW Gene expression")
# save the plot to a file
plt.savefig("/scratch/mtn1n22/affinity-matrix/output/figures/HIST-RAW-cutoff.png")
# clear the plot
plt.clf()

# As both the Raw and the normalised data have the same distribution,
# and theres no sign of a bi-modal distribution, we can use a cutoff of >0 to binarize the data.

# binarize the data
celltype_gene_matrix[celltype_gene_matrix > 0] = 1

# check the summary of the matrix
print(celltype_gene_matrix.describe())

# save the matrix as a csv file
celltype_gene_matrix.to_csv("/scratch/mtn1n22/affinity-matrix/output/celltype_gene_matrix_binarized.csv", index = True)