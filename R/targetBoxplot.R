#' Target Boxplot
#'
#' Makes boxplot of target counts.
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
targetBoxplot <- function(data_matrix){
  colors_ICs <- unlist(alamarColorPalette(n=1, nReps=2, interpolate=TRUE, interpolateIndex=6))
  boxplot(log(t(countMatrix1[,3:98])+1, base=2)[,order(apply(countMatrix1[,3:98], 1, median),
                                                       apply(countMatrix1[,3:98], 1, mean))],
          names=countMatrix1$BarcodeName[order(apply(countMatrix1[,3:98], 1, median),
                                               apply(countMatrix1[,3:98], 1, mean))],
          las=2, cex.axis=0.35,
          yaxt='n',
          ylab='log2(count+1)',
          main='NextSeq Run #1 -- target distributions')
  axis(side=2, las=1)
  
  boxplot_target <- function(countMatrix, title){
    boxplot(log(countMatrix[,3:ncol(countMatrix)] + 1, base=2),
            col='white',
            yaxt='n',
            ylab=expression('log'[2]*'(count+1)'),
            main=title,
            outline=FALSE,
            border='white',
            btn='n',
            las=3,
            cex.axis=0.5)
    axis(side=2, las=1)
    grid(nx=0, ny=NULL, lty=1)
    boxplot(log(countMatrix[,3:ncol(countMatrix)] + 1, base=2),
            col=colors_target,
            yaxt='n',
            ylab='',
            xaxt='n',
            main='',
            add=TRUE, 
            whisklty=1,
            staplewex=0)
  }
}
