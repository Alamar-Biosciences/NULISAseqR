#' Calculate Intra-plate Coefficients of Variation
#'
#' Calculates CV% for multiple targets 
#' for the specified sample replicates on a given plate. 
#'
#' @param data_matrix The Data matrix output from readNULISAseq.R
#' or normalized data from normalization functions.
#' @param samples A vector of length equal to ncol(data_matrix) which encodes
#' the sets of sample replicates that will be used to calculate intra-plate CVs. 
#' Can be either a numeric or a character vector, where matching entries 
#' represent replicates for a given sample. 
#' Any NAs will be excluded (it may be desirable to 
#' exclude negative controls or inter-plate controls, for example). 
#' Entries of this vector will be used to identify unique sample replicate sets 
#' in the results output. 
#' @param aboveLOD Logical TRUE/FALSE matrix output by the NULISAseqR::lod()
#' function. Indicates which values to include / exclude in the CV calculation. 
#' Default is NULL, which includes all values in CV calculation.
#' @param exclude_targets A vector of row names or row indices representing the 
#' targets that should be excluded from intra-plate CV calculation. For example,
#' one might want to exclude internal controls. Default is NULL, which includes
#' all targets in the data_matrix.
#'
#'
#' @return A list.
#' @param cv_matrix matrix with %CV values for each target for each 
#' multiple-replicate sample.
#'
#' @examples
#' 
#'
#' @export
#' 
intraCV <- function(data_matrix, 
                    samples,
                    aboveLOD=NULL,
                    exclude_targets=NULL){
  # replace values below LOD with NA
  data_matrix[aboveLOD==FALSE] <- NA
  # remove excluded targets
  if (!is.null(exclude_targets)){
    if (!is.numeric(exclude_targets)){
      exclude_targets <- which(rownames(data_matrix) %in% exclude_targets)
    }
    data_matrix <- data_matrix[-exclude_targets,]
  }
  # loop through the unique replicate sets to get CVs
  unique_samples <- unique(na.exclude(samples))
  cv_matrix <- matrix(nrow=nrow(data_matrix), ncol=length(unique_samples))
  rownames(cv_matrix) <- rownames(data_matrix)
  colnames(cv_matrix) <- unique_samples
  for (i in 1:length(unique_samples)){
    sample_data <- data_matrix[, samples==unique_samples[i]]
    sample_means <- rowMeans(sample_data, na.rm=TRUE)
    sample_sds <- apply(sample_data, 1, sd, na.rm=TRUE)
    sample_cv <- sample_sds/sample_means*100
    cv_matrix[,i] <- sample_cv
  }
  return(cv_matrix)
}
