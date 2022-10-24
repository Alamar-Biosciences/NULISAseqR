#' Fold Change Boxplot
#'
#' Makes boxplot of fold change or fold change error.
#'
#' @param fc_table Output of foldChange function.
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
fcBoxplot <- function(FCdata, ylabel=NULL, ylimits=NULL, label_yval){
  par(mar=c(15, 6, 2, 0.5), xpd=FALSE)
  FCdata[!is.finite(FCdata)] <- NA
  MAE <- apply(FCdata, 2, function(x) {
    x05 <- quantile(x, 0.05, na.rm=TRUE)
    x95 <- quantile(x, 0.95, na.rm=TRUE)
    x90 <- x[(x >= x05) & (x <= x95)]
    mean(abs(x90-2), na.rm=TRUE)})
  boxplot(FCdata[,order(MAE)],
          names=methodNames[order(MAE)],
          col='white',
          yaxt='n',
          xaxt='n',
          ylab=ylabel, 
          xlab='',
          outline=FALSE,
          border='white',
          ylim=ylimits)
  par(xpd=FALSE)
  abline(h=-30:50*0.5, col='grey', lty='dotted')
  boxplot(FCdata[,order(MAE)],
          names=methodNames[order(MAE)],
          col=colors_CV_medians[order(MAE)],
          las=3,
          yaxt='n',
          cex.axis=0.5,
          add=TRUE,
          ylim=ylimits)
  axis(side=2, las=1)
  abline(h=2)
  # MAE <- round(apply(FC2_IL28, 2, function(x) sqrt(mean((x-2)^2, na.rm=TRUE))), 2)
  MAE <- formatC(round(MAE, 3), format = 'f', flag='0', digits = 3)[order(MAE)]
  text(x=0:length(MAE), y=rep(label_yval, length(MAE)+1),
       labels=c("MAE", MAE), srt=90, cex=0.5, adj=c(0.25, 0.5),
       font=2, col='darkred')
}

FC_boxplot <- function(FCdata, ylabel=NULL, ylimits=NULL, label_yval){
  par(mar=c(15, 6, 2, 0.5), xpd=FALSE)
  FCdata[!is.finite(FCdata)] <- NA
  yminimum <- min(FCdata, na.rm=TRUE)
  ymaximum <- max(FCdata, na.rm=TRUE)
  MAE <- apply(FCdata, 2, function(x) {
    x05 <- quantile(x, 0.05, na.rm=TRUE)
    x95 <- quantile(x, 0.95, na.rm=TRUE)
    x90 <- x[(x >= x05) & (x <= x95)]
    mean(abs(x90-2), na.rm=TRUE)})
  boxplot(FCdata[,order(MAE)],
          names=methodNames[order(MAE)],
          col='white',
          yaxt='n',
          xaxt='n',
          ylab=ylabel, 
          xlab='',
          outline=FALSE,
          border='white',
          ylim=c(yminimum-2, ymaximum))
  par(xpd=FALSE)
  abline(h=-60:300*0.5, col='grey', lty='dotted')
  boxplot(FCdata[,order(MAE)],
          names=methodNames[order(MAE)],
          col=colors_CV_medians[order(MAE)],
          las=3,
          yaxt='n',
          cex.axis=0.5,
          add=TRUE,
          ylim=c(yminimum-2, ymaximum))
  axis(side=2, las=1)
  abline(h=2)
  # MAE <- round(apply(FC2_IL28, 2, function(x) sqrt(mean((x-2)^2, na.rm=TRUE))), 2)
  MAE <- formatC(round(MAE, 3), format = 'f', flag='0', digits = 3)[order(MAE)]
  text(x=0:length(MAE), y=rep(yminimum-2, length(MAE)+1),
       labels=c("MAE", MAE), srt=90, cex=0.5, adj=c(0.25, 0.5),
       font=2, col='darkred')
}
