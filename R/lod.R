# Convenience function to replace outliers using MAD
replace_outliers_mad <- function(col) {
  if (max(col, na.rm = TRUE) >= 100 && length(col) >= 4) { # only apply outlier detection if median > 100 reads
    outliers_mad(col, threshold = 2.5, maxRemoval = 1 , returnVal = "replaceNA")
  } else {
    col
  }
}

#' Calculate Limits of Detection
#'
#' Calculates limit of detection (LoD) for each target based on the
#' negative controls (blanks). LoD = mean(blanks) + 3*SD(blanks).
#' Designates data as either above or below LoD. 
#' Option to specify minimum count threshold for detectability.
#'
#' @param blanks Column indices or column names of the blanks in the
#' data_matrix.
#' @param data_matrix The Data matrix output from readNULISAseq.R
#' or normalized data from normalization functions.
#' @param min_count Optional count threshold to apply in addition
#' to the LoD. Default is 0.
#'
#'
#' @return A list.
#' @param LOD Vector of limits of detection.
#' @param aboveLOD Logical matrix indicating whether counts are 
#' above or below LoD for that target.
#'
#'
#' @export
#' 
lod <- function(data_matrix, blanks, min_count=0){
  if (!is.numeric(blanks)){
    blanks <- which(colnames(data_matrix) %in% blanks)
  }
  # replace NA with 0
  blank_data <- data_matrix[,blanks]
  blank_data[is.na(blank_data)] <- 0

  blank_mean <- colMeans(apply(blank_data, 1, replace_outliers_mad), na.rm=TRUE)
  blank_sd <- apply(apply(blank_data, 1, replace_outliers_mad), 2, sd, na.rm=TRUE)
  LOD <- blank_mean + 3*blank_sd
  if(min_count > 0){
    LOD[LOD < min_count] <- min_count
  }
  names(LOD) <- rownames(data_matrix)
  aboveLOD <- data_matrix
  aboveLOD <- apply(aboveLOD, 2, function(x){
    x > LOD
  })
  return(list(LOD=LOD, 
              aboveLOD=aboveLOD))
}
