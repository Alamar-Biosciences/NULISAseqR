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
                         cex=0.7,
                         title=NULL){
  data_matrix <- t(matrix(target_data, nrow=8, ncol=12, byrow=TRUE)[8:1,])
  image(data_matrix,
        xaxt='n', yaxt='n', main='', 
        col=hcl.colors(12, "YlOrRd", rev = TRUE))
  axis(2, (0:7)/7, c('H','G','F','E','D','C','B','A'), las=1, cex.axis=0.8)
  axis(1, (0:11)/11, 1:12, las=1, cex.axis=0.8)
  mtext(title, side=3, line=1, font=2, cex.axis=1)
  if (print_counts==TRUE){
    for (x in 1:nrow(data_matrix))
      for (y in 1:ncol(data_matrix))
        text(1/(nrow(data_matrix)-1)*(x-1), 
             1/(ncol(data_matrix)-1)*(y-1), 
             data_matrix[x,y], cex=cex)
  }
}
