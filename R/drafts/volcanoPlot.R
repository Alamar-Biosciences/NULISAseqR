#' Volcano Plot
#'
#' Makes volcano plot using output of deTest function.
#'
#' @param 
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
volcanoPlot <- function(estimate, pval, title, xlimits=NULL, ylimits=NULL){
  par(mar=c(3, 3, 2, 0.3))
  truePos <- 41:50
  color <- (1:50 %in% truePos)*1 + 1
  if (is.null(xlimits)){
    xmin <- min(estimate)
    xmax <- max(estimate)
    xlimits <- c(xmin, xmax)
  } 
  if (is.null(ylimits)){
    ymin <- min(estimate)
    ymax <- max(estimate)
    ylimits <- c(ymin, ymax)
  }
  log10pval <- -log(pval, base=10)
  plot(estimate, log10pval,
       col=color,
       las=1,
       xlab='',
       ylab='',
       xlim=xlimits,
       ylim=ylimits,
       main=title,
       cex.main=1,
       cex.axis=0.9)
  abline(h=-log(0.05, base=10))
  abline(v=0)
  grid()
  title(xlab='estimate', line=2, cex.lab=0.9)
  title(ylab=expression('-log'[10]*'(Holm-corrected p-value)'), line=2, cex.lab=0.9)
}
