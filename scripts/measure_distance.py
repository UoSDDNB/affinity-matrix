# Measure distance between celltypes in the celltype_funcition_matrix.
# perhaps use the cosine similarity to measure the distance between celltypes.
# or Euclidean distance.

# save the distances in a matrix, the affinity matrix.



# consider normalizing the values from celltype_function_matrix before measuring the distance.
# normalize by the sum of the values in each row.

# import necessary libraries
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

# load the celltype_function_matrix
celltype_function_matrix = pd.read_csv("/scratch/mtn1n22/affinity-matrix/output/celltype_function_matrix.csv", index_col = 0)

# two rows have NA valeus, im not sure why.
# for now, I will remove the rows with NA values.
# remove the rows with NA values
celltype_function_matrix = celltype_function_matrix.dropna()

# normalize the values in the matrix
celltype_function_matrix = celltype_function_matrix.div(celltype_function_matrix.sum(axis=1), axis=0)

# measure the distance between celltypes (columns) in the matrix.
# compute the cosine similarity between celltypes
celltype_similarity_matrix = cosine_similarity(celltype_function_matrix.T)
# convert to a distance matrix
#elltype_distance_matrix = 1 - celltype_similarity_matrix
# convert without generating nagative values
# generated due to floating point errors.
celltype_distance_matrix = np.maximum(0, 1 - celltype_similarity_matrix)


# save the distance matrix with column names.
np.savetxt("/scratch/mtn1n22/affinity-matrix/output/celltype_distance_matrix.csv",
    celltype_distance_matrix,
    delimiter=",",
    header=",".join(celltype_function_matrix.columns))


# plot the diatnaces using PCA
from sklearn.decomposition import PCA

# PCA
pca = PCA(n_components=2)
pca.fit(celltype_distance_matrix)
celltype_pca = pca.transform(celltype_distance_matrix)

# plot the PCA projection
import seaborn as sns
import matplotlib.pyplot as plt

# plot the PCA projection
plt.figure(figsize=(10, 10))
sns.scatterplot(x=celltype_pca[:, 0], y=celltype_pca[:, 1])
plt.title("PCA Projection of Cosine Distance Matrix")
plt.xlabel("PCA 1")
plt.ylabel("PCA 2")
# save the plot to a file
plt.savefig("/scratch/mtn1n22/affinity-matrix/output/figures/PCA.png")
# clear the plot
plt.close()


# plot the distance matrix
import umap

# UMAP embedding
reducer = umap.UMAP(metric="precomputed", random_state=42)
embedding = reducer.fit_transform(celltype_distance_matrix)

# plot the embedding
import seaborn as sns
import matplotlib.pyplot as plt

# plot the UMAP projection
plt.figure(figsize=(10, 10))
sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1])
plt.title("UMAP Projection of Cosine Distance Matrix")
plt.xlabel("UMAP 1")
plt.ylabel("UMAP 2")
# save the plot to a file
plt.savefig("/scratch/mtn1n22/affinity-matrix/output/figures/UMAP.png")
# clear the plot
plt.close()


