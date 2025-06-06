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
#' @param sample_groups A string vector defining sample types / identities that
#' represent the sample type subsets that detectability should be calculated 
#' for. Default value NULL provides overall detectability assuming all samples
#' are of the same type (e.g., plasma). NOTE that sample groups must be the 
#' same length and in the same order as \code{sample_subset}. 
#' @param exclude_targets A vector of row names or numeric row indices 
#' representing the 
#' targets that should be excluded from detectability calculation. For example,
#' one might want to exclude internal controls. Default is NULL, which includes
#' all targets in the aboveLOD_matrix.
#'
#'
#' @return A nested list containing detectability information.
#' The list has two high-level names: "all" and "sample_group".
#' The "all" name corresponds to detectability information for all samples
#' combined as a single unit. If the \code{sample_groups} parameter is not NULL,
#' the "sample_group" name provides detectability data for unique sample groups.
#' The nested list has three names to access detectability and sample information.
#' \item{sampleNumber}{The number of samples for each group.}
#' \item{detectability}{A named vector of detectability values in percentages.
#' The names correspond to target row names.}
#' \item{detectable}{A named logical vector indicating whether detectability is
#' above 50\% for each target. The names correspond to target row names.}
#'
#'
#' @export
#' 
detectability <- function(aboveLOD_matrix, 
                          sample_subset=NULL,
                          sample_groups=NULL,
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
  
  # get sample subset
  if (!is.null(sample_subset)){
    if (!is.numeric(sample_subset)){
      sample_subset <- which(colnames(aboveLOD_matrix) %in% sample_subset)
    }
    aboveLOD_matrix <- aboveLOD_matrix[,sample_subset]
  }
  
  # get detectability data for samples: Plasma, serum, etc.
  if (!is.null(sample_groups)){
    # make case insensitive
    sample_groups <- tolower(sample_groups)
    
    # collect sample group specific detectability data
    for (i in unique(sample_groups)){
      if(is.null(dim(aboveLOD_matrix))){
        aboveLOD_matrix <- matrix(aboveLOD_matrix, ncol=1)
      }
      sample_aboveLOD_matrix <- as.matrix(aboveLOD_matrix[,sample_groups %in% i], nrow=nrow(aboveLOD_matrix))

      detectability_data$sample_group$sampleNumber[[i]] <- ncol(sample_aboveLOD_matrix)
      detectability_data$sample_group$detectability[[i]] <- apply(as.matrix(sample_aboveLOD_matrix), 1, function(x) sum(x)/length(x)*100)
      detectability_data$sample_group$detectable[[i]] <- detectability_data$sample_group$detectability[[i]] > 50
    }
  }
  
  # get overall detectability data
  detect <- apply(as.matrix(aboveLOD_matrix), 1, function(x) sum(x)/length(x)*100)
  detectable <- detect > 50
  if(is.null(dim(aboveLOD_matrix))){
    names(detect) <- rownames(as.matrix(aboveLOD_matrix))
    names(detectable) <- rownames(as.matrix(aboveLOD_matrix))
  } else{
    names(detect) <- rownames(aboveLOD_matrix)
    names(detectable) <- rownames(aboveLOD_matrix)
  }
  
  detectability_data$all$sampleNumber <- ncol(aboveLOD_matrix)
  detectability_data$all$detectability <- detect
  detectability_data$all$detectable <- detectable
  
  return(detectability_data)
}
