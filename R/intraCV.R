#' Calculate Intra-plate Coefficients of Variation
#'
#' Calculates CV% for the specified replicates on a given plate. 
#'
#' @param data_matrix The Data matrix output from readNULISAseq.R
#' or normalized data from normalization functions.
#
#'
#' @return matrix with %CV values for each target for each 
#' multiple-replicate sample.
#'
#' @examples
#' 
#'
#' @export
#' 
# NEEDS EDITING
intraCV <- function(countMatrix, samples){
  countMatrix <- countMatrix[-c(NIC, MCherry),]
  unique_samples <- unique(samples)[order(unique(samples))][c(16:21,1,7:14,2:6,15,22,28:35,23:27,36:44)]
  cv_results <- data.frame(matrix(nrow=nrow(countMatrix), ncol=length(unique_samples)))
  rownames(cv_results) <- countMatrix$BarcodeName
  colnames(cv_results) <- unique_samples
  for (i in 1:length(unique_samples)){
    sample_data <- countMatrix[,3:98][,samples==unique_samples[i]]
    sample_means <- rowMeans(sample_data)
    sample_sds <- apply(sample_data,1,sd)
    sample_cv <- sample_sds/sample_means*100
    cv_results[,i] <- sample_cv
  }
  return(cv_results)
}