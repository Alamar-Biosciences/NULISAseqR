#' Generate a Target QC Plot
#'
#' Generates a quantifiability histogram using output of 
#' \code{quantifiability()} function.
#'
#' @param runs A named list of run data output from \code{laodNULISAseq()} function 
#' or a list of these outputs for multiple runs.
#'
#' @return Generates a target QC plot.
#' 
#'
#'
#' @export
#' 
targetQCplot <- function(runs){
  
  # check if run data is not in a list and if so put into a list
  if('RunSummary' %in% names(runs)){
    runs <- list(runs)
    names(runs) <- 'Plate 01'
  } 
  
  targetQC <- lapply(runs, function(x) {
    target_QC_output <- QCFlagTarget(AQdata=x$AQ$Data_AQ_aM, 
                                     raw = x$Data,
                                     IPCnormed = x$normed$interNormData[[1]],
                                     detectability=x$detectability$all$detectability,
                                     aboveLOD=x$lod$aboveLOD,
                                     withinDR=x$AQ$withinDR, 
                                     absRun = !is.null(x$AQ$Data_AQ_aM),
                                     targets=x$targets, 
                                     samples=x$samples, 
                                     SCparams=x$AQ$targetAQ_param$SC_conc)
    target_QC_output <- target_QC_output[target_QC_output$flagName=='Target_Conc_Accuracy',]
    target_QC_output$percent_recovery <- as.numeric(target_QC_output$val) * 100
    target_QC_output <- target_QC_output[order(target_QC_output$target),]
    return(target_QC_output)
    
  })
  
  # define target labels
  labels <- targetQC[[1]]$target
  
  # define plate colors
  interpolate <- FALSE
  if(length(targetQC) > 13) interpolate <- TRUE
  plate_colors <- alamarColorPalette(n=length(targetQC), interpolate=interpolate)
  if(length(targetQC)==1) plate_colors <- alamarColorPalette(n=1, palette = 2)
  
  # get QC threshold
  QC_threshold <- as.numeric(unlist(strsplit(targetQC[[1]]$QCthreshold[1], ','))) * 100
  
  # plot first plate
  yval <- targetQC[[1]]$percent_recovery
  # replace values < -100 or > 100 
  yval[yval < -100 & !is.na(yval)] <- -100
  yval[yval > 100 & !is.na(yval)] <- 100
  # set point shape based on recovery error
  point_shape <- rep(1, length(yval))
  point_shape[yval < QC_threshold[1] & !is.na(yval)] <- 4
  point_shape[yval > QC_threshold[2] & !is.na(yval)] <- 4
  
  par(mfcol=c(1, 1), mar=c(4,4,5,0.5))
  
  plot(x=1:length(yval),
       y=as.numeric(yval), 
       ylab="Sample Control Recovery Error (%)", 
       xlab='',
       xaxt='n', 
       las=1, 
       cex=0.7, ylim=c(-100, 100), 
       col=plate_colors[1],
       pch=point_shape,
       main='')
  title(main='Target QC: Sample Control Recovery', line=4)
  axis(side=1, at=1:length(yval), labels=NA, tck=-0.01, 
       xpd=NA, las=2, adj=1, cex.axis=0.6)
  axis(side=1, at=1:length(yval), labels=labels, tck=NA, 
       xpd=NA, las=2, adj=1, cex.axis=0.4, line=-0.6, lwd=0)
  abline(v=1:length(yval), col='lightgray', lty='dotted')
  abline(h = QC_threshold[1], col='gray41')
  abline(h = QC_threshold[2], col='gray41')
  abline(h = 0, col='gray41', lty=2)
  
  # loop over remaining plates
  if(length(targetQC) > 1){
    for(i in 2:length(targetQC)){
      yval <- targetQC[[i]]$percent_recovery
      # replace values < -100 or > 100 
      yval[yval < -100 & !is.na(yval)] <- -100
      yval[yval > 100 & !is.na(yval)] <- 100
      # set point shape based on recovery error
      point_shape <- rep(1, length(yval))
      point_shape[yval < QC_threshold[1] & !is.na(yval)] <- 4
      point_shape[yval > QC_threshold[2] & !is.na(yval)] <- 4
      
      points(x=1:length(yval), y=as.numeric(yval), 
             col=plate_colors[i], cex=0.7, pch=point_shape)
    }
  }
  

  legend_names <- c(names(targetQC), 'Within +/- 30%', 'Outside +/- 30%')
  suppressWarnings(
    legend('bottom', 
           legend=legend_names, 
           col=c(plate_colors, 'black', 'black'), 
           pch = c(rep(15, length(targetQC)), 1, 4), cex=0.75, 
           bty='n', inset=c(0,1), xpd=TRUE, horiz=FALSE, ncol = 10)
  )
}

