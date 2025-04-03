# imports
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
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
    # Define groupings
    groupings = {  
        "Immune Cells": [
            'Macrophage', 'Monocyte', 'Dendritic cell', 'B cell', 'T cell', 'Mast cell', 
            'Neutrophil', 'Myeloid cell', 'M2 Macrophage', 'B cell (Plasmocyte)', 
            'Antigen presenting cell (RPS high)', 'Adrenal gland inflammatory cell', 
            'Proliferating T cell', 'Immature Sertoli cell (Pre-Sertoli cell)'
        ],

        "Barrier and Structural Cells": [
            'Epithelial cell', 'Fibroblast', 'Stromal cell', 'Endothelial cell', 'Mesothelial cell', 
            'Basal cell', 'Sinusoidal endothelial cell', 'Endothelial cell (APC)', 
            'Endothelial cell (endothelial to mesenchymal transition)', 
            'Epithelial cell (intermediated)', 'Intermediated cell', 'Stratified epithelial cell', 
            'Goblet cell', 'Enterocyte', 'Pancreas exocrine cell', 'Hepatocyte/Endodermal cell', 
            'Loop of Henle', 'Chondrocyte', 'Fetal chondrocyte', 'Fetal fibroblast', 
            'Fetal acinar cell', 'Gastric chief cell'
        ],

        "Contractile Cells": [
            'Smooth muscle cell', 'Ventricle cardiomyocyte', 'Cell of skeletal muscle', 
            'Fetal skeletal muscle cell'
        ],

        "Neural Cells": [
            'Neuron', 'Astrocyte', 'Oligodendrocyte', 'Fetal neuron', 'Fetal astrocyte'
        ],

        "Blood and Hematopoietic Cells": [
            'Erythroid cell', 'Erythroid progenitor cell (RP high)', 'CB CD34+'
        ],

        "Endocrine Cells": [
            'Endocrine cell', 'Gastric endocrine cell', 'Fetal endocrine cell', 'Thyroid follicular cell', 
            'Pancreas exocrine cell'
        ],

        "Mesodermal and Connective Tissue Cells": [
            'Mesothelial cell', 'Fibroblast', 'Stromal cell', 'Chondrocyte', 'Fetal chondrocyte',
            'Fetal fibroblast'
        ],
    }

    # Define colors
    colors = {
        "Immune Cells": "red",
        "Barrier and Structural Cells": "blue",
        "Contractile Cells": "purple",
        "Neural Cells": "orange",
        "Blood and Hematopoietic Cells": "brown",
        "Endocrine Cells": "green",
        "Mesodermal and Connective Tissue Cells": "pink",
        "Other": "grey"  # Default color for ungrouped cell types
    }

    # Create a lookup dictionary
    celltype_to_group = {celltype: group for group, celltypes in groupings.items()
                         for celltype in celltypes}

    # Assign each column (cell type) to a group, defaulting to "Other"
    celltype_groups = [celltype_to_group.get(celltype, "Other")
                       for celltype in celltype_function_matrix.columns]
    # Ensure all groups have corresponding colors, defaulting to "grey" if missing
    celltype_colors = [colors.get(group, "grey") for group in celltype_groups]

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
    plt.close()

    # Elbow plot
    plt.figure(figsize=(8, 5))
    cumulative_variance = np.cumsum(pca.explained_variance_ratio_)
    plt.plot(range(1, len(cumulative_variance) + 1), cumulative_variance, marker="o")
    plt.xlabel("Number of Principal Components")
    plt.ylabel("Cumulative Explained Variance")
    plt.title("Choosing Number of PCs for UMAP")
    plt.axhline(y=0.9, color='r', linestyle='--')  # 90% variance threshold
    plt.axhline(y=0.95, color='g', linestyle='--')  # 95% variance threshold
    plt.ylim(0, 1)  # Ensure y-axis starts at 0 and ends at 1
    # Adjust x-axis limit to cutoff when cumulative variance is 97%, or use total components as fallback
    x_limit = np.argmax(cumulative_variance >= 0.97) + 1 if np.any(cumulative_variance >= 0.97) else len(cumulative_variance)
    plt.xlim(0, x_limit)
    plt.xticks(range(5, len(cumulative_variance) + 1, 5))
    plt.savefig("output/figures/ELBOW_function.png")
    plt.close()
