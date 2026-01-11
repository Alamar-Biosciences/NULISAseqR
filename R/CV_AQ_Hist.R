#' Make a CV Histogram for AQ Data
#'
#' Generates a coefficient of variation histogram using output of 
#' \code{CV_AQ()} function.
#'
#' @param CV_data A vector of CV values output from \code{CV_AQ()}
#' function. 
#' @param plot_type Controls the tile of plot. Can be \code{'intra'} or 
#' {'inter'} for intra-plate or inter-plate CV, respectively. 
#' Default title will be "AQ Intra-plate CV (\%)" or "AQ Inter-plate CV (\%)".
#' @param title Optional title for the plot. Overrides default title.
#'
#'
#' @return Generates a histogram plot.
#'
#' @export
#' 
CV_AQ_Hist <- function(CV_data,
                       plot_type=NULL,
                       title=NULL){
  
  if(is.null(plot_type) & is.null(title)) stop('Must provide plot_type or title.')
  if(!is.null(plot_type)){
    if(plot_type=='intra') plot_title <- 'AQ Intra-plate CV (%)'
    if(plot_type=='inter') plot_title <- 'AQ Inter-plate CV (%)'
    if(!(plot_type %in% c('intra', 'inter')) & is.null(title)){
      stop('plot_type must be either \'intra\' or \'inter\' or NULL. If NULL then title must be given.')
    }
  } 
  if(!is.null(title)) plot_title <- title
  
  colors <- alamarColorPalette(2, nReps = 5, palette = 2)
  par(mar=c(4.5, 4.5, 2, 1))
  hist(CV_data,
       main=plot_title,
       xlab='CV (%)', las=1, breaks=50,
       xlim= c(0, max(CV_data, na.rm=TRUE)),
       col=colors[[2]][5],
       border=colors[[2]][2])
  abline(v=mean(CV_data, na.rm=TRUE), col=colors[[1]][3], lwd=2)
  abline(v=median(CV_data, na.rm=TRUE), col=colors[[1]][3], lty=2, lwd=2)
  legend('topright', c(paste0('mean CV% = ', format(round(mean(CV_data, na.rm=TRUE), 1), nsmall=1)), 
                       paste0('median CV% = ', format(round(median(CV_data, na.rm=TRUE), 1), nsmall=1))),
         lty=c(1,2), col=colors[[1]][3], bty='n', lwd=2, cex=1.25)
}

