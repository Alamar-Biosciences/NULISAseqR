#' Sample Boxplot
#'
#' Makes sample boxplots.
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
colors_ICs <- unlist(alamarColorPalette(n=1, nReps=2, interpolate=TRUE, interpolateIndex=6))
plot_ICs <- function(countMatrix,
                     ICindex,
                     PICindex,
                     well_order,
                     title,
                     names=NULL){
  orderedLogCountMatrix <- log(countMatrix[,3:ncol(countMatrix)] + 1, base=2)[,well_order]
  boxplot(orderedLogCountMatrix,
          col='lightgrey',
          yaxt='n',
          ylab=expression('log'[2]*'(count+1)'),
          main=title,
          outline=TRUE,
          border='grey',
          whisklty=1,
          staplewex=0,
          las=3,
          ylim=c(0, max(orderedLogCountMatrix)),
          cex.axis=0.5,
          names=names)
  axis(side=2, las=1)
  # plot NIC
  lines(1:ncol(orderedLogCountMatrix), orderedLogCountMatrix[ICindex,], lty=1, col=colors_ICs[1])
  # plot PIC
  lines(1:ncol(orderedLogCountMatrix), orderedLogCountMatrix[PICindex,], col=colors_ICs[2], lty=1)
  legend('topright', 
         legend=c('NIC', 'MCherry'),
         lty=c(1, 1), 
         col=colors_ICs, bty='n', cex=0.75)
}

boxplot_sample <- function(countMatrix, techReps, title, names=NULL, colors){
  boxplot(t(log(countMatrix[order(techReps),3:ncol(countMatrix)] + 1, base=2)),
          col='white',
          yaxt='n',
          ylab=expression('log'[2]*'(count+1)'),
          main=title,
          outline=FALSE,
          border='white',
          btn='n',
          las=3,
          ylim=c(0, max(log(countMatrix[order(techReps),3:ncol(countMatrix)] + 1, base=2))),
          cex.axis=0.5,
          names=names)
  axis(side=2, las=1)
  grid(nx=0, ny=NULL, lty=1)
  boxplot(t(log(countMatrix[order(techReps),3:ncol(countMatrix)] + 1, base=2)),
          col=colors,
          yaxt='n',
          ylab='',
          xaxt='n',
          main='',
          add=TRUE, 
          whisklty=1,
          staplewex=0)
}
