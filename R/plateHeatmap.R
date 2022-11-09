#' Draw Plate Heatmap
#'
#' Makes heatmap for a specified target or total counts.
#'
#' @param data_matrix The Data matrix output from readNULISAseq.R
#' or normalized data from normalization functions.
#' @param target The row name or row index of the target.
#'
#' @return Draws a plot.
#'
#' @examples
#' 
#'
#' @export
#' 
# NEEDS EDITING
plateHeatmap <- function(data_matrix){
  m <- t(matrix(as.numeric(countMatrix1[NIC,3:98]), nrow=8, ncol=12, byrow=TRUE)[8:1,])
  image(m,
        xaxt='n', yaxt='n', main=paste0('NextSeq Run #1 ', countMatrix1$BarcodeName[NIC]), 
        col=hcl.colors(12, "YlOrRd", rev = TRUE))
  for (x in 1:nrow(m))
    for (y in 1:ncol(m))
      text(1/(nrow(m)-1)*(x-1), 1/(ncol(m)-1)*(y-1), m[x,y], cex=0.7)
}
