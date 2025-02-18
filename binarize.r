#binarize celltype_matrix 
# first I will:
# convert this to work with a matrix
# then I will:
# convert it to work in python!?


# test the function
library(GeneSwitches)
#load "data/celltype_matrix.csv"
celltype_matrix <- read.csv("data/celltype_matrix.csv", row.names = 1)

# plot a histogram of the gene expression
hist(as.numeric(as.matrix(celltype_matrix)), 
     breaks = 100, 
     main = "Histogram of gene expression", 
     xlab = "Gene expression")

# zoom in and set limits
hist(as.numeric(as.matrix(celltype_matrix)), 
     breaks = 1000000, 
     main = "Histogram of gene expression", 
     xlab = "Gene expression",
     xlim = c(0, 0.0002),
     ylim = c(0, 5000))
# add a vertical line at 0.00002
abline(v = 0.00002, col = "red")

# binarize the data
library(parallel)
library(mixtools)
binary_matrix <- binarize_exp(celltype_matrix,
                              fix_cutoff = TRUE,
                              binarize_cutoff = 0.00002,
                              ncores = 6)

# save the binary matrix as a csv
write.csv(binary_matrix, "data/binary_celltype_matrix.csv")

#' @title Binarize gene expression
#'
#' @description This function generates on/off binarized data for gene expression
#'
#' @param expdata Matrix containing gene expression data
#' @param fix_cutoff Logical. If TRUE, use a fixed global cutoff for binarization, default is FALSE
#' @param binarize_cutoff Numeric. Fixed global cutoff for binarization, default is 0.2
#' @param ncores Integer. Number of cores to use for parallel processing, default is 3
#' @return Matrix with binarized gene expression data
#'
#' @import parallel
#' @importFrom mixtools normalmixEM
#' @export
#'
binarize_exp <- function(expdata, fix_cutoff = FALSE, binarize_cutoff = 0.2, ncores = 4) {
  # Calculate the percentage of zero expression for each gene
  zerop_g <- c()
  for (i in 1:nrow(expdata)) {
    zp <- length(which(expdata[i, ] == 0)) / ncol(expdata)
    zerop_g <- c(zerop_g, zp)
  }

  if (fix_cutoff == TRUE) {
    # Use fixed global cutoff for binarization
    is.na(expdata) <- expdata == 0
    exp_reduced_binary <- as.matrix((expdata > binarize_cutoff) + 0)
    exp_reduced_binary[is.na(exp_reduced_binary)] = 0
    binary <- exp_reduced_binary
  } else {
    # Use mixture model for binarization
    # Add Gaussian noise to gene expression matrix with a standard deviation of 0.1
    LogCountsadd = expdata + matrix(rnorm(nrow(expdata) * ncol(expdata),
                                          mean = 0, sd = 0.1),
                                    nrow(expdata), ncol(expdata))
    # Fit mixture models for each gene in parallel
    oupBinary = do.call(
      rbind, mclapply(rownames(LogCountsadd), function(iGene){
        set.seed(42)   # Set seed for consistency
        tmpMix = normalmixEM(LogCountsadd[iGene, ], k = 2)
        if (tmpMix$mu[1] < tmpMix$mu[2]) {
          tmpOup = data.frame(geneID = iGene,
                              mu1 = tmpMix$mu[1],
                              mu2 = tmpMix$mu[2],
                              sigma1 = tmpMix$sigma[1],
                              sigma2 = tmpMix$sigma[2],
                              lambda1 = tmpMix$lambda[1],
                              lambda2 = tmpMix$lambda[2],
                              loglik = tmpMix$loglik)
        } else {
          tmpOup = data.frame(geneID = iGene,
                              mu1 = tmpMix$mu[2],
                              mu2 = tmpMix$mu[1],
                              sigma1 = tmpMix$sigma[2],
                              sigma2 = tmpMix$sigma[1],
                              lambda1 = tmpMix$lambda[2],
                              lambda2 = tmpMix$lambda[1],
                              loglik = tmpMix$loglik)
        }
        return(tmpOup)
      }, mc.cores = ncores))

    # Identify non-bimodal genes
    oupBinary$passBinary = TRUE
    oupBinary[oupBinary$lambda1 < 0.1, ]$passBinary = FALSE
    oupBinary[oupBinary$lambda2 < 0.1, ]$passBinary = FALSE
    oupBinary[(oupBinary$mu2 - oupBinary$mu1) < (oupBinary$sigma1 + oupBinary$sigma2), ]$passBinary = FALSE

    # Solve for intersection for remaining genes
    oupBinary$root = -1
    for(iGene in oupBinary[oupBinary$passBinary == TRUE, ]$geneID){
      tmpMix = oupBinary[oupBinary$geneID == iGene, ]
      tmpInt = uniroot(function(x, l1, l2, mu1, mu2, sd1, sd2) {
        dnorm(x, m = mu1, sd = sd1) * l1 -
          dnorm(x, m = mu2, sd = sd2) * l2},
        interval = c(tmpMix$mu1, tmpMix$mu2),
        l1 = tmpMix$lambda1, mu1 = tmpMix$mu1, sd1 = tmpMix$sigma1,
        l2 = tmpMix$lambda2, mu2 = tmpMix$mu2, sd2 = tmpMix$sigma2)
      oupBinary[oupBinary$geneID == iGene, ]$root = tmpInt$root
    }

    # Binarize expression data
    binLogCounts = expdata[oupBinary$geneID, ]
    binLogCounts = t(scale(t(binLogCounts), scale = FALSE, center = oupBinary$root))
    binLogCounts[binLogCounts >= 0] = 1
    binLogCounts[binLogCounts < 0] = 0
    binary <- binLogCounts
  }
  return(binary)
}
