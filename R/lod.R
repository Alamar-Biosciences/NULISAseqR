#' Calculate Limits of Detection
#'
#' Calculates limit of detection (LoD) for each target based on the
#' negative controls. Designates data as either above or below LoD. 
#' Option to specify minimum count threshold for detectability.
#'
#' @param data_matrix The Data matrix output from readNULISAseq.R
#' or normalized data from normalization functions.
#
#'
#' @return 
#'
#' @examples
#' 
#'
#' @export
#' 
# NEEDS EDITING
lod <- function(plate, blanks){
  blank_mean <- rowMeans(plate[,blanks])
  blank_sd <- apply(plate[,blanks], 1, sd)
  blank_lod <- blank_mean + 3*blank_sd
  return(blank_lod)
}