#' Calculate Detectability for a Set of Targets and Samples
#'
#' Calculates detectability for each target based on the
#' aboveLOD matrix output of lod() function. Detectability is percent of samples
#' above LOD for that target.
#'
#' @param aboveLOD_matrix aboveLOD matrix output from lod() function.
#' Rows are targets and columns are samples. 
#' @param sample_subset A vector of column names or numeric column indices that
#' represent the sample subset that detectability will be calculated for. 
#' Default uses all columns. NCs and IPCs should probably be excluded.
#' @param sample_identity A string vector defining sample types/identities that
#' represent the sample type subsets that detectability should be calculated 
#' for. Default value NULL provides overall detectability assuming all samples
#' are of the same type (i.e., plasma).
#' @param exclude_targets A vector of row names or numeric row indices 
#' representing the 
#' targets that should be excluded from detectability calculation. For example,
#' one might want to exclude internal controls. Default is NULL, which includes
#' all targets in the aboveLOD_matrix.
#'
#'
#' @return A nested list containing detectability information.
#' The list has two high-level names: "all" and "sample".
#' The "all" name corresponds to detectability information for all samples
#' combined as a single unit. If the \code{sample_identity} parameter is not NULL,
#' the "sample" name provides detectability data for unique sample types.
#' The nested list has three names to access detectability and sample information.
#' @param sampleNumber The number of samples.
#' @param detectability A named vector of detectability values in percentages.
#' The names correspond to target row names.
#' @param detectable A named logical vector indicating whether detectability is
#' above 50% for each target. The names correspond to target row names.
#'
#'
#' @export
#' 
detectability <- function(aboveLOD_matrix, 
                          sample_subset=NULL,
                          sample_identity=NULL,
                          exclude_targets=NULL){
  # return named list
  detectability_data <- list()
  
  # remove excluded targets
  if (!is.null(exclude_targets)){
    if (!is.numeric(exclude_targets)){
      exclude_targets <- which(rownames(aboveLOD_matrix) %in% exclude_targets)
    }
    aboveLOD_matrix <- aboveLOD_matrix[-exclude_targets,]
  }
  
  # get detectability data for samples: Plasma, serum, etc.
  if (!is.null(sample_identity)){
    # make case insensitive
    sample_identity <- tolower(sample_identity)
    
    # set all samples if sample_subset == NULL
    if (is.null(sample_subset)){
      all_samples <- colnames(aboveLOD_matrix)
    } else{
      all_samples <- sample_subset
    }

    # collect sample type/identity specific detectability data
    for (i in unique(sample_identity)){
      sample_identity_inds <- which(sample_identity %in% i)
      sample_subset_inds <- which(colnames(aboveLOD_matrix) %in% all_samples[sample_identity_inds])
      sample_aboveLOD_matrix <- aboveLOD_matrix[,sample_subset_inds]
      
      detectability_data$sample$sampleNumber[[i]] <- length(sample_identity_inds)
      detectability_data$sample$detectability[[i]] <- apply(sample_aboveLOD_matrix, 1, function(x) sum(x)/length(x)*100)
      detectability_data$sample$detectable[[i]] <- detectability_data$sample$detectability[[i]] > 50
    }
  }
  
  # get overall detectability data: Apply sample_subset if applicable
  if (!is.null(sample_subset)){
    if (!is.numeric(sample_subset)){
      sample_subset <- which(colnames(aboveLOD_matrix) %in% sample_subset)
    }
    aboveLOD_matrix <- aboveLOD_matrix[,sample_subset]
  }
  detect <- apply(aboveLOD_matrix, 1, function(x) sum(x)/length(x)*100)
  detectable <- detect > 50
  names(detect) <- rownames(aboveLOD_matrix)
  names(detectable) <- rownames(aboveLOD_matrix)
  
  detectability_data$all$sampleNumber <- length(colnames(aboveLOD_matrix))
  detectability_data$all$detectability <- detect
  detectability_data$all$detectable <- detectable
  
  return(detectability_data)
}
