#' NULISAseq Inter-Plate Normalization
#'
#' Does inter-plate normalization using one or both of the following methods:
#' inter-plate control (IPC; divides each target 
#' by the median IPC count for that target), intensity normalization (IN; 
#' sets median of each target on each plate to equal the global median). Default
#' is to do IPC only. If both options are chosen, IPC is done first followed
#' by IN.
#'   
#' Input is a list of data matrices, one for each plate. 
#' These would typically be normData matrix from the intraPlateNorm function 
#' or the Data matrix of readNULISAseq.R (if no intra-plate normalization is 
#' done). Other required input depends on methods used. 
#' 
#' Output is a single matrix.
#'
#' @param data_list A list of data matrices, one for each plate, 
#' to be normalized together. 
#' @param IPC Logical TRUE/FALSE indicating whether a IPC 
#' normalization should be done. Default is TRUE.
#' @param IN Logical TRUE/FALSE indicating whether intensity normalization
#' should be done. Recommended only when samples are randomized 
#' across the plates. If targets do not match across all plates, only the 
#' matching targets will be returned in the output. Default is FALSE.
#' @param IPC_wells List of vectors of either the column indices (numeric) 
#' or the column names (character string) of the inter-plate control wells 
#' for each plate. IPC wells are omitted from the 
#' intensity normalization median calculation. 
#' Required for IPC normalization and recommended for IN.  
#' @param NC_wells Recommended for intensity normalization. 
#' List of vectors of either the row indices (numeric) or the 
#' row names (character string) 
#' of the negative control wells for each plate. NC wells are omitted from 
#' intensity normalization. 
#' @param IPC_method 'median' is the default. Other options include 
#' 'mean' (arithmetic mean) and 'geom_mean' (geometric mean). 
#' Determines how the counts are summarized 
#' across the IPC wells on a given plate.
#' @param IN_samples Optional argument. A list of column names or
#' indices specifying which subset of samples to use for 
#' intensity normalization step for each plate in data_list. 
#' Will over-ride the IPC_wells and 
#' NC_wells arguments for IN. 
#' @param scaleFactor. Optional numeric value used to rescale all data 
#' after normalizing. Default is 1. This may be desirable to avoid
#' normalized quantities between 0 and 1 (which will be negative
#' in the log scale).
#'
#' @return 
#' @param interNormData A matrix of normalized count data (not log-transformed).
#' @param plate A vector of integers 1-(number of plates) 
#' corresponding to the interNormData matrix columns
#' that indicates which plate samples are from.
#'
#' @examples
#' 
#'
#' @export
#' 
interPlateNorm - function(data_list,
                          IPC=TRUE,
                          IN=FALSE,
                          IPC_wells=NULL,
                          NC_wells=NULL,
                          IPC_method='median',
                          IN_samples=NULL,
                          scaleFactor=1){
  # inter-plate control normalization
  if (IPC==TRUE){
    if(is.null(IPC_wells)){
      stop('Must provide IPC_wells for IPC normalization.')
    }
    IPC_factors <- list()
    data_list_IPC <- list()
    for (i in 1:length(data_list)){
      IPC_cols <- data_list[[i]][,IPC_wells[[i]]]
      if (IPC_method=='median'){
        IPC_factors_i <- apply(IPC_cols, 1, median, na.rm=TRUE)
      } else if (IPC_method=='mean'){
        IPC_factors_i <- rowMeans(IPC_cols, median, na.rm=TRUE)
      } else if (IPC_method=='geom_mean'){
        IPC_factors_i <- apply(IPC_cols, 1, function(x){
          x <- x + 1
          prod(x, na.rm=TRUE)^(1/length(x))
        })
      }
      # replace 0s or NAs with 1s
      # ADD WARNING HERE SAYING WHICH TARGETS THIS WAS DONE FOR
      IPC_factors_i[IPC_factors_i==0] <- 1
      IPC_factors_i[is.na(IPC_factors_i)] <- 1
      IPC_factors[[i]] <- IPC_factors_i
      data_list_IPC[[i]] <- data_list[[i]] / IPC_factors_i
    }
    data_list <- data_list_IPC
  } # end IPC normalization
  
  # intensity normalization
  if (IN==TRUE){
    # get the matching subset of targets
    target_list <- list()
    for (i in 1:length(data_list)){
      target_list[[i]] <- rownames(data_list[[i]])
    }
    matching_targets <- Reduce(intersect, target_list)
    all_targets <- unique(unlist(target_list))
    omitted_targets <- all_targets[!(all_targets %in% matching_targets)]
    if(length(omitted_targets) > 0){
      warning(paste0('The following targets are missing from some plates and will be omitted from intensity normalization:\n', 
                    paste(omitted_targets, collapse='\n')))
    }
    # remove the omitted targets
    for (i in 1:length(data_list)){
      data_list[[i]] <- data_list[[i]][rownames(data_list[[i]]) %in% matching_targets,]
      # make sure rows are sorted in the same order
      data_list[[i]] <- data_list[[i]][match(matching_targets, rownames(data_list[[i]])),]
    }
    # if IN_samples is undefined
    # omit the IPC and / or NC wells if specified
    if (is.null(IN_samples)){
      IN_samples <- list()
      for (i in 1:length(data_list)){
        if (!is.null(IPC_wells)){
          if (!is.numeric(IPC_wells[[i]][1])){
            IPC_wells[[i]] <- which(colnames(data_list[[i]]) %in% IPC_wells[[i]])
          }
        }
        if (!is.null(NC_wells)){
          if (!is.numeric(NC_wells[[i]][1])){
            NC_wells[[i]] <- which(colnames(data_list[[i]]) %in% NC_wells[[i]])
          }
        }
        IN_omit <- c(IPC_wells[[i]], NC_wells[[i]]) 
        # use all sample columns except IPCs and NCs
        if (!is.null(IN_omit)){
          col_index <- 1:ncol(data_list[[i]])
          IN_samples[[i]] <- col_index[!(col_index %in% IN_omit)]
        } else {
          IN_samples[[i]] <- col_index
        }
      }
    } # end defining IN_samples
    # get the subset of data used for IN
    # and get target medians
    data_list_IN_samples <- list()
    IN_medians <- matrix(nrow=length(matching_targets), ncol=length(data_list))
    rownames(IN_medians) <- matching_targets
    for (i in 1:length(data_list)){
      data_list_IN_samples[[i]] <- data_list[[i]][,IN_samples[[i]]]
      IN_medians[,i] <- apply(data_list_IN_samples[[i]], 1, median, na.rm=TRUE)
    }
    global_medians <- apply(IN_medians, 1, median, na.rm=TRUE)
    IN_factors <- (1/IN_medians)*global_medians
    # IN rescale data
    for (i in 1:length(data_list)){
      data_list[[i]] <- data_list[[i]]*IN_factors[,i]
    }
  } # end intensity normalization
  
  # apply the scale factor to the data
  data_list <- lapply(data_list, function(x) x*scaleFactor)
  # generate the plate assignment vector
  plateNs <- unlist(lapply(data_list, ncol))
  plate <- NULL
  for (i in 1:length(plateNs)){
    plate <- c(plate, rep(i, plateNs[i]))
  }
  
  # return output
  return(interNormData=interNormData,
         plate=plate)
}
