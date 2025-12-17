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
    aboveLOD_matrix <- aboveLOD_matrix[,sample_subset, drop=FALSE]
  }
  
  # get detectability data for samples: Plasma, serum, etc.
  if (!is.null(sample_groups)){
    # make case insensitive
    sample_groups <- toupper(sample_groups)
    
    # collect sample group specific detectability data
    for (i in unique(sample_groups)){
      if(is.null(dim(aboveLOD_matrix))){
        aboveLOD_matrix <- matrix(aboveLOD_matrix, ncol=1)
      }
      sample_aboveLOD_matrix <- as.matrix(aboveLOD_matrix[,sample_groups %in% i], nrow=nrow(aboveLOD_matrix))
      
      detectability_data$sample_group$sampleNumber[[i]] <- ncol(sample_aboveLOD_matrix)
      detectability_data$sample_group$detectability[[i]] <- apply(as.matrix(sample_aboveLOD_matrix), 1, function(x) sum(x, na.rm=TRUE) / sum(!is.na(x))*100)
      detectability_data$sample_group$detectable[[i]] <- detectability_data$sample_group$detectability[[i]] > 50
    }
  }
  
  # get overall detectability data
  detect <- apply(as.matrix(aboveLOD_matrix), 1, function(x) sum(x, na.rm=TRUE) / sum(!is.na(x))*100)
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




#' Summarize Detectability Across Multiple Runs and Sample Groups
#'
#' Summarizes detectability across multiple runs and if applicable 
#' broken down by sample groups based on the sample_group_covar input
#' to loadNULISAseq.
#'
#' @param runs A named list of run data output from \code{laodNULISAseq()} function 
#' or a list of these outputs for multiple runs. To make output more interpretable,
#' it is recommended to name each run according to the plate ID for that run.
#'
#'
#' @return A list with items "sample_group", if applicable, "all", and "run_detectability".
#' Items sample_group and all each have a data frame with a row for each target 
#' and detectability in columns. "all" has only one column summarizing detectability for all 
#' samples. "sample_group" has a column summarizing detectability for each 
#' sample group. Each of there items also has a sampleNumber vector which gives the 
#' sample size for the groups. The run_detectability output is a list of the 
#' original individual run detectability objects from loadNULISAseq.
#' 
#'
#'
#' @export
#' 
detectability_summary <- function(runs){
  
  # check if run data is not in a list and if so put into a list
  if('RunSummary' %in% names(runs)){
    runs <- list(runs)
    names(runs) <- 'Plate 01'
  } 
  
  # get detectability data for each run
  detect <- lapply(runs, function(x) x$detectability)
  
  # summarize detectability for all samples 
  all <- list()
  all$sampleNumber <- sum(sapply(detect, function(x) x$all$sampleNumber))
  all_detectability <- lapply(detect, function(x) {
    all_detect <- data.frame(Target=names(x$all$detectability), 
                             detectability=x$all$detectability / 100 * x$all$sampleNumber)
    return(all_detect)
  })
  all_detectability <- suppressWarnings(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Target", all = TRUE),
                                               all_detectability))
  all_detectability$total_detectability <- rowSums(as.matrix(all_detectability[,2:ncol(all_detectability)])) / all$sampleNumber * 100
  all$detectability <- data.frame(Target=all_detectability$Target, 
                                  detectability=all_detectability$total_detectability)
  rownames(all$detectability) <- all$detectability$Target
  colnames(all$detectability)[2] <- paste0('detectability (n = ', all$sampleNumber, ')')
  
  output <- list(all=all,
                 run_detectability=detect)
  
  # get sample group names for each run
  group_names <- unique(unlist(lapply(detect, function(x) names(x$sample_group$sampleNumber))))
  names(group_names) <- group_names
  if(length(group_names) > 0){
    sample_group_detect <- lapply(group_names, function(x) {
      sampleNumber <- sum(unlist(lapply(detect, function(y) {
        sampleNumber <- y$sample_group$sampleNumber[[x]]
      })))
      sample_group_detectability <- lapply(detect, function(y) {
        if(x %in% names(y$sample_group$detectability)){
          data.frame(Target=names(y$sample_group$detectability[[x]]), 
                     detectability=y$sample_group$detectability[[x]] / 100 * y$sample_group$sampleNumber[[x]])
        } else {
          NULL
        }
      })
      
      sample_group_detectability <- sample_group_detectability[sapply(sample_group_detectability, function(x) !is.null(x))]
      sample_group_detectability <- suppressWarnings(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Target", all = TRUE),
                                                            sample_group_detectability))
      sample_group_detectability$total_detectability <- rowSums(as.matrix(sample_group_detectability[,2:ncol(sample_group_detectability)])) / sampleNumber * 100
      sample_group_detectability <- data.frame(Target=sample_group_detectability$Target, 
                                               detectability=sample_group_detectability$total_detectability)
      rownames(sample_group_detectability) <- sample_group_detectability$Target
      colnames(sample_group_detectability)[2] <- paste0(x, ' (n = ', sampleNumber, ')')
      
      return(list(sampleNumber=sampleNumber,
                  sample_group_detectability=sample_group_detectability))
    })
    
    # combine into a single sampleNumber vector and detectability data frame
    sampleNumber <- sapply(sample_group_detect, function(x) x$sampleNumber)
    sample_group_detectability <- lapply(sample_group_detect, function(x) x$sample_group_detectability)
    sample_group_detectability <- suppressWarnings(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Target", all = TRUE),
                                                          sample_group_detectability))
    rownames(sample_group_detectability) <- sample_group_detectability$Target
    
    sample_group <- list()
    sample_group$sampleNumber <- sampleNumber
    sample_group$detectability <- sample_group_detectability
    
    output <- list(sample_group=sample_group,
                   all=all,
                   run_detectability=detect)
  }
  
  return(output)
}
