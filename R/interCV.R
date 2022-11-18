#' Calculate Inter-plate Coefficients of Variation
#'
#' Calculates CV% for the specified replicates across multiple plates. 
#' This function only calculates inter-plate CV for 
#' targets that are present across all plates in the data_list. Any
#' targets that are present on only a subset of plates will be excluded.
#' In addition, individual targets, such as internal controls, can be 
#' explicitly excluded using the exclude_targets argument.
#'
#' @param data_list A list of data matrices that contain data for 
#' the sample replicates for calculating inter-plate CV.
#' @param samples A list of vectors of length equal to 
#' the number of columns in the corresponding data matrix which encodes
#' the sets of sample replicates that will be used to calculate inter-plate CVs. 
#' Can be a list of either numeric or character vectors, where matching entries 
#' represent replicates for a given sample. 
#' Any NAs in the sample vectors will be excluded (it may be desirable to 
#' exclude negative controls or inter-plate controls, for example). 
#' Unique entries of these vectors will be used to identify 
#' unique sample replicate sets in the results output. 
#' @param aboveLOD A list of logical TRUE/FALSE matrices output by the 
#' NULISAseqR::lod() function. Indicates which values to include / exclude 
#' in the CV calculation. 
#' Default is NULL, which includes all values in CV calculation.
#' @param exclude_targets A list of vectors of row names or row indices 
#' representing the 
#' targets that should be excluded from intra-plate CV calculation. For example,
#' one might want to exclude internal controls. Default is NULL, which includes
#' all common targets in the data_list matrices.
#' @param useMean A logical TRUE / FALSE that indicates whether or not 
#' to calculate inter-CV using the intra-plate means of replicates (TRUE, default), 
#' or to use the pooled replicate data (FALSE). 
#'
#' @return matrix with inter-plate %CV values for each target for each sample.
#'
#' @examples
#' 
#'
#' @export
#'
interCV <- function(data_list, 
                    samples,
                    aboveLOD=NULL,
                    exclude_targets=NULL,
                    useMean=TRUE){
  # replace values below LOD with NA
  for (i in 1:length(data_list)){
    data_list[[i]][aboveLOD[[i]]==FALSE] <- NA
  }
  # remove excluded targets
  if(!is.null(exclude_targets)){
    for (i in 1:length(data_list)){
      if (!is.numeric(exclude_targets[[i]])){
        exclude_targets[[i]] <- which(rownames(data_list[[i]]) %in% exclude_targets[[i]])
      }
      data_list[[i]] <- data_list[[i]][-exclude_targets[[i]],]
    }
  }
  # get the matching subset of targets across all plates
  target_list <- list()
  for (i in 1:length(data_list)){
    target_list[[i]] <- rownames(data_list[[i]])
  }
  matching_targets <- Reduce(intersect, target_list)
  all_targets <- unique(unlist(target_list))
  omitted_targets <- all_targets[!(all_targets %in% matching_targets)]
  if(length(omitted_targets) > 0){
    warning(paste0('The following targets are missing from some plates and will be omitted from inter-plate CV results:\n', 
                   paste(omitted_targets, collapse='\n')))
  }
  # remove the omitted targets
  for (i in 1:length(data_list)){
    data_list[[i]] <- data_list[[i]][rownames(data_list[[i]]) %in% matching_targets,]
    # make sure rows are sorted in the same order
    data_list[[i]] <- data_list[[i]][match(matching_targets, rownames(data_list[[i]])),]
  }
  # loop through the unique replicate sets to get CVs
  unique_samples <- unique(na.exclude(unlist(samples)))
  cv_matrix <- matrix(nrow=length(matching_targets), ncol=length(unique_samples))
  rownames(cv_matrix) <- matching_targets
  colnames(cv_matrix) <- unique_samples
  for (i in 1:length(unique_samples)){
    sample_data <- list()
    for (j in 1:length(data_list)){
      sample_data[[j]] <- data_list[[j]][,samples[[j]]==unique_samples[i]]
    }
    if (useMean==TRUE){
      sample_data <- lapply(sample_data, rowMeans, na.rm=TRUE)
    } 
    sample_data <- do.call(cbind, sample_data)
    sample_means <- rowMeans(sample_data, na.rm=TRUE)
    sample_sds <- apply(sample_data, 1, sd, na.rm=TRUE)
    sample_cv <- sample_sds/sample_means*100
    cv_matrix[,i] <- sample_cv
  }
  return(cv_matrix)
}
