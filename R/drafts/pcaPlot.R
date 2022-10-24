#' Plot Principle Components
#'
#' Makes scatterplot of samples using principal component scores.
#'
#' @param data_matrix The Data matrix output from readNULISAseq.R
#' or normalized data from normalization functions.
#
#'
#' @return Draws a plot.
#'
#' @examples
#' 
#'
#' @export
#' 
# NEEDS EDITING
### INPUTS
## countMatrix: data frame matchingReadCountsMatrix 
### from readNULIAseq.R, or normalized count matrix
## PCs: numeric vector of length 2 that specifies 
### PCs to plot on the x and y axes. Default is 
### first and second components, i.e., c(1, 2).
## 
## labels: (optional) character vector for labelling 
### points on plot. Length should match number of samples 
### and by in same order as countMatrix rows.
## col: (optional) vector specifying sample point colors.
### Length should match number of samples 
### and by in same order as countMatrix rows.
### OUTPUTS
## makes a scatter plot of sample points 
###############################################
PCAplot <- function(countMatrix,
                           PCs=c(1,2),
                           colors='black',
                           labels=FALSE,
                           label.colors='black'){
  log_counts <- log(t(countMatrix[,3:ncol(countMatrix)])+1, base=2)
  colnames(log_counts) <- countMatrix$BarcodeBName
  pca <- prcomp(log_counts)
  plot(pca$rotation[,PCs[1]], pca$rotation[,PCs[2]], col=colors, 
       xlab=paste0('PC', PCs[1]), ylab=paste0('PC', PCs[2]))
  if (labels==TRUE){
    par(xpd=TRUE)
    text(pca$rotation[,PCs[1]], pca$rotation[,PCs[2]], 
         labels=countMatrix$BarcodeBName,
         cex=0.5, pos=4,
         col=label.colors)
  }
}
