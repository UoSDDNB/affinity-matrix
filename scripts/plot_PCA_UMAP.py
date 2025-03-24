import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import umap
from sklearn.decomposition import PCA
from adjustText import adjust_text

def plot_pca(celltype_function_matrix, output_dir, label_points=True, num_pcs=20):
    """
    Perform PCA on celltype_function_matrix and plot PCA projection & elbow plot.

    Args:
        celltype_function_matrix (pd.DataFrame): Input matrix with genes as rows and cell types as columns.
        output_dir (str): Directory to save the plots.
        label_points (bool): Whether to label each point with the column name.
        num_pcs (int): Number of principal components to retain for UMAP.

    Returns:
        np.ndarray: PCA-transformed data (first num_pcs dimensions).
    """
    pca = PCA()
    pca.fit(celltype_function_matrix.T)
    celltype_pca = pca.transform(celltype_function_matrix.T)

    # PCA plot
    plt.figure(figsize=(10, 10))
    sns.scatterplot(x=celltype_pca[:, 0], y=celltype_pca[:, 1], alpha=0.7)

    if label_points:
        texts = [plt.text(celltype_pca[i, 0], celltype_pca[i, 1], label, fontsize=8, alpha=0.7)
                 for i, label in enumerate(celltype_function_matrix.columns)]
        adjust_text(texts)  # Avoid label overlap

    plt.title("PCA Projection of celltype_function_matrix")
    plt.xlabel("PCA 1")
    plt.ylabel("PCA 2")
    plt.savefig(f"{output_dir}/PCA_function.png")
    plt.clf()

    # Elbow plot
    plt.figure(figsize=(8, 5))
    plt.plot(np.cumsum(pca.explained_variance_ratio_), marker="o")
    plt.xlabel("Number of Principal Components")
    plt.ylabel("Cumulative Explained Variance")
    plt.title("Choosing Number of PCs for UMAP")
    plt.axhline(y=0.9, color='r', linestyle='--')  # 90% variance threshold
    plt.axhline(y=0.95, color='g', linestyle='--')  # 95% variance threshold
    plt.savefig(f"{output_dir}/ELBOW_function.png")
    plt.clf()

    return celltype_pca[:, :num_pcs]  # Return first num_pcs dimensions for UMAP


def plot_umap(pca_reduced, celltype_function_matrix, output_dir, label_points=True):
    """
    Perform UMAP on PCA-reduced data and plot the UMAP projection.

    Args:
        pca_reduced (np.ndarray): PCA-transformed data (first few principal components).
        celltype_function_matrix (pd.DataFrame): Input matrix (used for labeling).
        output_dir (str): Directory to save the plots.
        label_points (bool): Whether to label each point with the column name.
    """
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, metric='euclidean', random_state=42)
    celltype_umap = reducer.fit_transform(pca_reduced)

    # UMAP plot
    plt.figure(figsize=(10, 10))
    sns.scatterplot(x=celltype_umap[:, 0], y=celltype_umap[:, 1], alpha=0.7)

    if label_points:
        texts = [plt.text(celltype_umap[i, 0], celltype_umap[i, 1], label, fontsize=8, alpha=0.7)
                 for i, label in enumerate(celltype_function_matrix.columns)]
        adjust_text(texts)

    plt.title("UMAP Projection of celltype_function_matrix")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.savefig(f"{output_dir}/UMAP_function.png")
    plt.clf()


# Example usage:
celltype_function_matrix = pd.read_csv("/scratch/mtn1n22/affinity-matrix/output/celltype_function_matrix.csv", index_col=0)
output_dir = "/scratch/mtn1n22/affinity-matrix/output/figures"

# Run PCA
pca_reduced = plot_pca(celltype_function_matrix, output_dir, label_points=True, num_pcs=20)

# Run UMAP
plot_umap(pca_reduced, celltype_function_matrix, output_dir, label_points=True)
