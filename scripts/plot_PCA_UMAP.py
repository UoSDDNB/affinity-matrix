import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import umap
from sklearn.decomposition import PCA
from adjustText import adjust_text

def plot_pca(celltype_function_matrix, label_points=True):
    """
    Perform PCA on celltype_function_matrix and plot PCA projection & elbow plot.

    Args:
        celltype_function_matrix (pd.DataFrame): Input matrix with genes as rows and cell types as columns.
        label_points (bool): Whether to label each point with the column name.

    Returns:
        np.ndarray: PCA-transformed data (first num_pcs dimensions).
    """
    # Define groupings and colors
    groupings = {
        "Immune cells": ["macrophage", "monocyte", "T cell", "B cell", "dendritic cell",
                         "neutrophil", "mast cell", "inflammatory cell", "myeloid cell",
                         "professional antigen presenting cell"],
        "Barrier and structural cells": ["epithelial cell", "fibroblast", "mesothelial cell",
                                         "smooth muscle cell", "endothelial cell", "chondrocyte",
                                         "stromal cell", "basal cell", "mesenchymal stem cell"],
        "Secretory cells": ["acinar cell", "goblet cell", "endocrine cell", "plasma cell",
                            "peptic cell", "epithelial cell of exocrine pancreas",
                            "thyroid follicular cell", "Sertoli cell"],
        "Contractile cells": ["ventricular cardiac muscle cell", "cell of skeletal muscle",
                              "smooth muscle cell"],
        "Neural cells": ["neuron", "astrocyte", "oligodendrocyte"]
    }
    colors = {
        "Immune cells": "red",
        "Barrier and structural cells": "blue",
        "Secretory cells": "green",
        "Contractile cells": "purple",
        "Neural cells": "orange",
        "Other": "grey"  # Default color for ungrouped cell types
    }

    # Create a lookup dictionary
    celltype_to_group = {celltype: group for group, celltypes in groupings.items()
                         for celltype in celltypes}

    # Assign each column (cell type) to a group, defaulting to "Other"
    celltype_groups = [celltype_to_group.get(celltype, "Other")
                       for celltype in celltype_function_matrix.columns]
    celltype_colors = [colors[group] for group in celltype_groups]


    # Perform PCA
    pca = PCA()
    pca.fit(celltype_function_matrix.T)
    celltype_pca = pca.transform(celltype_function_matrix.T)

     # PCA plot
    plt.figure(figsize=(10, 10))
    sns.scatterplot(
        x=celltype_pca[:, 0], 
        y=celltype_pca[:, 1], 
        hue=celltype_groups,
        palette=colors,  # Now all categories exist in colors
        alpha=0.7
    )

    if label_points:
        texts = [plt.text(celltype_pca[i, 0], celltype_pca[i, 1], label, fontsize=8, alpha=0.7)
                 for i, label in enumerate(celltype_function_matrix.columns)]
        adjust_text(texts)  # Avoid label overlap

    plt.title("PCA Projection of celltype_function_matrix")
    plt.xlabel("PCA 1")
    plt.ylabel("PCA 2")
    plt.legend(title="Cell Groups", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig("output/figures/PCA_function_colored.png", bbox_inches="tight")
    plt.clf()

    # Elbow plot
    plt.figure(figsize=(8, 5))
    plt.plot(np.cumsum(pca.explained_variance_ratio_), marker="o")
    plt.xlabel("Number of Principal Components")
    plt.ylabel("Cumulative Explained Variance")
    plt.title("Choosing Number of PCs for UMAP")
    plt.axhline(y=0.9, color='r', linestyle='--')  # 90% variance threshold
    plt.axhline(y=0.95, color='g', linestyle='--')  # 95% variance threshold
    plt.savefig("output/figures/ELBOW_function.png")
    plt.clf()

    return celltype_pca  # Return PCA-transformed data for UMAP


def plot_umap(pca, celltype_function_matrix, num_pcs, label_points=True):
    """
    Perform UMAP on PCA-transformed data and plot the UMAP projection,
    colouring points based on predefined cell type groupings.

    Args:
        pca (np.ndarray): PCA-transformed data.
        celltype_function_matrix (pd.DataFrame): Input matrix (used for labeling).
        num_pcs (int): Number of principal components to use for UMAP.
        label_points (bool): Whether to label each point with the column name.
    """
    # Define groupings and colors
    groupings = {
        "Immune cells": ["macrophage", "monocyte", "T cell", "B cell", "dendritic cell",
                         "neutrophil", "mast cell", "inflammatory cell", "myeloid cell",
                         "professional antigen presenting cell"],
        "Barrier and structural cells": ["epithelial cell", "fibroblast", "mesothelial cell",
                                         "smooth muscle cell", "endothelial cell", "chondrocyte",
                                         "stromal cell", "basal cell", "mesenchymal stem cell"],
        "Secretory cells": ["acinar cell", "goblet cell", "endocrine cell", "plasma cell",
                            "peptic cell", "epithelial cell of exocrine pancreas",
                            "thyroid follicular cell", "Sertoli cell"],
        "Contractile cells": ["ventricular cardiac muscle cell", "cell of skeletal muscle",
                              "smooth muscle cell"],
        "Neural cells": ["neuron", "astrocyte", "oligodendrocyte"]
    }
    colors = {
        "Immune cells": "red",
        "Barrier and structural cells": "blue",
        "Secretory cells": "green",
        "Contractile cells": "purple",
        "Neural cells": "orange",
        "Other": "grey"  # Default color for ungrouped cell types
    }

    # Create lookup dictionary for cell type groups
    celltype_to_group = {celltype: group for group, celltypes in groupings.items()
                         for celltype in celltypes}
    
    # Assign each column (cell type) to a group, defaulting to "Other"
    celltype_groups = [celltype_to_group.get(celltype, "Other")
                       for celltype in celltype_function_matrix.columns]
    celltype_colors = [colors[group] for group in celltype_groups]

    # Reduce PCA dimensions
    pca_reduced = pca[:, :num_pcs]

    # UMAP
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, metric='euclidean', random_state=42)
    celltype_umap = reducer.fit_transform(pca_reduced)

    # UMAP plot
    plt.figure(figsize=(10, 10))
    sns.scatterplot(x=celltype_umap[:, 0], y=celltype_umap[:, 1], hue=celltype_groups,
                    palette=colors, alpha=0.7)

    if label_points:
        texts = [plt.text(celltype_umap[i, 0], celltype_umap[i, 1], label, fontsize=8, alpha=0.7)
                 for i, label in enumerate(celltype_function_matrix.columns)]
        adjust_text(texts)

    plt.title("UMAP Projection of celltype_function_matrix")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.legend(title="Cell Groups", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig("output/figures/UMAP_function.png", bbox_inches="tight")
    plt.clf()


