#' Draw Plate Heatmap
#'
#' Makes heatmap for a specified target or total counts.
#'
#' @param target_data A numeric vector of length 96 that represents the 
#' counts that will be shown on the heatmap.
#' @param well_order Optional. A numeric vector used to order the wells 
#' from A01-A12, B01-B12, ..., H01-H12. If not provided, then 
#' it is assumed that the target_data is already in this order.
#' @param print_counts A TRUE / FALSE that specifies whether or not 
#' counts will be printed in each cell of the heatmap. Default is TRUE.
#' @param cex Character expansion factor that controls size of the 
#' printed counts. Default is 0.7.
#' @param cex.axis Character expansion factor that controls size of the 
#' axis labels. Default is 0.8.
#' @param digits number of digits after the decimal point to print. 
#' Default is 0.
#' @param title Optional. Plot title.
#'
#' @return Draws a heatmap.
#'
#' @examples
#' 
#'
#' @export
#' 
plateHeatmap <- function(target_data, 
                         well_order=NULL,
                         print_counts=TRUE,
                         cex=0.6,
                         cex.axis=0.8,
                         digits=0,
                         title=NULL,relative=F){
  if(!is.null(well_order)){
    target_data <- target_data[well_order]
  }
  if (relative){
    target_data <- (target_data - median(target_data))/median(target_data)*100
  }
  data_matrix <- t(matrix(target_data, nrow=8, ncol=12, byrow=TRUE)[8:1,])
  if(relative){
    paletteLength <- 50
    col <- colorRampPalette(c("royalblue1", "white", "red"))(paletteLength)
    myBreaks <- c(seq(min(data_matrix), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(data_matrix)/paletteLength, max(data_matrix), length.out=floor(paletteLength/2)))
    image.plot(data_matrix,
          xaxt='n', yaxt='n', main='', 
          col=col, breaks=myBreaks, axis.args=list(cex.axis=cex.axis, tck=-0.15, mgp=c(3,0.5,0)))
  }else{
    col <-  hcl.colors(12, "YlOrRd", rev = TRUE)
    image.plot(data_matrix,
          xaxt='n', yaxt='n', main='',
          col=col, axis.args=list(cex.axis=cex.axis, tck=-0.15, mgp=c(3,0.5,0)))
  }
  axis(2, (0:7)/7, c('H','G','F','E','D','C','B','A'), col=NA, col.ticks=NA, las=1, cex.axis=cex.axis, tck=-0.01, line=-0.5)
  axis(1, (0:11)/11, 1:12, col=NA, col.ticks=1, las=1, cex.axis=cex.axis, tck=-0.01, padj=-3)
  mtext(title, side=3, line=1, font=2, cex.axis=1)
  if (print_counts==TRUE){
    for (x in 1:nrow(data_matrix))
      for (y in 1:ncol(data_matrix))
        text(1/(nrow(data_matrix)-1)*(x-1), 
             1/(ncol(data_matrix)-1)*(y-1), 
             round(data_matrix[x,y], digits=digits), cex=cex)
  }
}
