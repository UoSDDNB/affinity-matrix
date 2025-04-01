# load the necessary libraries
import matplotlib.pyplot as plt

def query_celltype_function_matrix(ct_fun):
    """
    Query the celltype-function matrix and return a histogram representing the variance of each function.

    Args:
        ct_fun (pd.DataFrame): Matrix with functions as rows and cell types as columns.

    Returns:
        None
    """

    # plot a histogram of the variance of each row
    ct_fun.var(axis=1).sort_values(ascending=True).plot(kind='hist', bins=100, color='purple', alpha=0.7)
    plt.title('Variance of each function')
    plt.xlabel('Variance')
    plt.ylabel('Frequency')
    plt.tight_layout()
    # draw a line at the mean
    mean = ct_fun.var(axis=1).mean()
    plt.axvline(mean, color='red', linestyle='dashed', linewidth=1)
    plt.savefig('/scratch/mtn1n22/affinity-matrix/output/figures/query/celltype_function_matrix/variance_of_each_function.png')
    plt.clf()