#' Make a Quantifiability Histogram
#'
#' Generates a quantifiability histogram using output of 
#' \code{quantifiability()} function.
#'
#' @param quantifiability_data A list output from \code{quantifiability()}
#' function. The list has elements named \code{quant} which is a 
#' data frame with quantifiability values where targets are in rows and 
#' subgroups are in columns, and \code{n_samples} which is a named vector with 
#' sample sizes corresponding to the columns and column names of quant.
#' @param group_name The name of the subgroup for which to plot the data.
#' This name should match one of the columns of \code{quantifiability_data$quant}
#' and one of the entries of \code{quantifiability_data$n_samples}.
#' @param title Optional title. By default the title is generated using the 
#' group name and sample size. 
#'
#'
#' @return Generates a histogram plot.
#' 
#'
#'
#' @export
#' 
quantHist <- function(quantifiability_data,
                      group_name,
                      title=NULL){
  
  if(is.null(title)) title <- paste0('AQ Target quantifiability (%):\n', group_name, 
                                     ' (n = ', quantifiability_data$n_samples[group_name],')')
  
  colors <- alamarColorPalette(2, nReps = 5, palette = 2)
  par(mar=c(4.5, 4.5, 2, 1))
  hist(quantifiability_data$quant[,group_name],
       main=title,
       xlab='Quantifiability (%)', las=1, breaks=50, 
       col=colors[[2]][5],
       border=colors[[2]][2],
       xlim=c(0, 100))
  abline(v=50, col=colors[[1]][3], lwd=1)
  ydist <- par('usr')[4] - par('usr')[3]
  text(47, par('usr')[4] - 0.1 * ydist, 'quantifiability = 50%',
       col=colors[[1]][3], adj=c(1,1), cex=0.5, srt=90)
  legend('topleft', c(paste0('mean = ', format(round(mean(quantifiability_data$quant[,group_name], na.rm=TRUE), 1), nsmall=1)),
                      paste0('median = ', format(round(median(quantifiability_data$quant[,group_name], na.rm=TRUE), 1), nsmall=1))),
         bty='n', cex=1)
  
}

