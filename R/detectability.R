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
#' @param exclude_targets A vector of row names or numeric row indices 
#' representing the 
#' targets that should be excluded from detectability calculation. For example,
#' one might want to exclude internal controls. Default is NULL, which includes
#' all targets in the aboveLOD_matrix.
#'
#'
#' @return A list.
#' @param detectability A named vector of detectability values as percentages. 
#' Names correspond to target row names.
#' @param detectable A named TRUE/FALSE vector that indicates whether 
#' detectability is above 50% for that target. 
#' Names correspond to target row names.
#'
#' @examples
#' 
#'
#' @export
#' 
detectability <- function(aboveLOD_matrix, 
                          sample_subset=NULL,
                          exclude_targets=NULL){
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
  detect <- apply(aboveLOD_matrix, 1, function(x) sum(x)/length(x)*100)
  detectable <- detect > 50
  names(detect) <- rownames(aboveLOD_matrix)
  names(detectable) <- rownames(aboveLOD_matrix)
  return(list(detectability=detect,
              detectable=detectable))
}
