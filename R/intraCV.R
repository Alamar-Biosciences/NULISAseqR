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
#' @param method Use "count" (default) for unnormalized or 
#' normalized count data. Calculates CV using the formula 
#' \code{cv = sd / mean}. Use "log2" for log2 transformed data. Calculates CV using 
#' the formula \code{sigma = sd(x) log(2)} and 
#' \code{cv = 100 sqrt(exp(sigma^2)-1)}.
#'
#'
#' @return A list.
#' @param cv_matrix matrix with %CV values for each target for each 
#' multiple-replicate sample.
#'
#' 
#'
#' @export
#' 
intraCV <- function(data_matrix, 
                    samples,
                    aboveLOD=NULL,
                    exclude_targets=NULL,
                    method='count'){
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
  if(method=='count'){
    for (i in 1:length(unique_samples)){
      sample_data <- data_matrix[,!is.na(samples) & samples==unique_samples[i]]
      sample_means <- rowMeans(sample_data, na.rm=TRUE)
      sample_sds <- apply(sample_data, 1, sd, na.rm=TRUE)
      sample_cv <- sample_sds/sample_means*100
      cv_matrix[,i] <- sample_cv
    }
  } else if(method=='log2'){
    # cv formula
    cv_formula <- function(x){
      sigma <- apply(x, 1, function(x) sd(x, na.rm=TRUE)*log(2))
      100*sqrt(exp(sigma^2)-1)
    }
    for (i in 1:length(unique_samples)){
      sample_data <- data_matrix[,!is.na(samples) & samples==unique_samples[i]]
      cv_matrix[,i] <- cv_formula(sample_data)
    }
  }
  return(cv_matrix)
}
