#' Coefficient of Variation Boxplot
#'
#' Makes boxplot of intra or inter-plate %CV by target.
#'
#' @param cv_table Output of intraCV or interCV function.
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
cv_boxplot <- function(cv_table, names, mar, title){
  par(mar=mar)
  boxplot(cv_table, col='white', border='white',
          las=3, names=names,
          ylab='', pars=list(font.lab=2),
          ylim=c(-5, max(cv_table, na.rm=TRUE)+2), xaxt='n')
  axis(side=1, at=1:8, labels=names, las=3, tick=FALSE)
  title(ylab=title, font.lab=2, line=2.5)
  par(xpd=FALSE)
  abline(h=(0:100)*2, col='lightgrey', lty='dotted')
  abline(h=0)
  boxplot(cv_table, col=alamarColorPalette(8),
          names=c('','','','','','','',''), 
          ann=FALSE, add=TRUE, xaxt='n')
  cv_medians <- sprintf(round(apply(cv_table, 2, median, na.rm=TRUE), 2), fmt = '%#.2f')
  text(c(0.5, 1:8),-4,labels=c('median', cv_medians), 
       srt=90, cex=0.75, font=2, col='darkred')
}

interCV_boxplot <- function(cv_table, names, mar, title){
  par(mar=mar)
  boxplot(cv_table, col='white', border='white',
          las=3, names=names,
          ylab='', pars=list(font.lab=2),
          ylim=c(-10, max(cv_table, na.rm=TRUE)+2), xaxt='n')
  axis(side=1, at=1:8, labels=names, las=3, tick=FALSE)
  title(ylab=title, font.lab=2, line=2.5)
  par(xpd=FALSE)
  abline(h=(0:100)*2, col='lightgrey', lty='dotted')
  abline(h=0)
  boxplot(cv_table, col=alamarColorPalette(8),
          names=c('','','','','','','',''), 
          ann=FALSE, add=TRUE, xaxt='n')
  cv_medians <- sprintf(round(apply(cv_table, 2, median, na.rm=TRUE), 2), fmt = '%#.2f')
  text(c(0.5, 1:8),-9,labels=c('median', cv_medians), 
       srt=90, cex=0.6, font=2, col='darkred')
}

colors_target_CV <- alamarColorPalette(n=62, interpolate=TRUE, interpolateIndex=1:8)

boxplot_target_CV <- function(avgCV, colors, title, ylimits=NULL){
  par(mar=c(6, 4, 2, 0.5))
  # avgCV <- avgCV[,c(1:51, 65:68, 61:64, 57:60, 53:56, 69, 52)]
  boxplot(avgCV[,2:ncol(avgCV)],
          names=colnames(avgCV)[2:ncol(avgCV)],
          col='white',
          yaxt='n',
          xaxt='n',
          main=title,
          ylab='Inter-Plate CV %', 
          xlab='',
          outline=FALSE,
          border='white',
          ylim=ylimits)
  abline(h=0:100*10, col='grey', lty='dotted')
  boxplot(avgCV[,2:ncol(avgCV)],
          names=colnames(avgCV)[2:ncol(avgCV)],
          col=colors,
          las=3,
          yaxt='n',
          cex.axis=0.4,
          add=TRUE,
          ylim=ylimits)
  axis(side=2, las=1)
}