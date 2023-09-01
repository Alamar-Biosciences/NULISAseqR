#' NULISAseq Inter-Plate Normalization
#'
#' Does inter-plate normalization using one of the following methods:
#' inter-plate control (IPC; divides each target 
#' by the median IPC count for that target), intensity normalization (IN; 
#' sets median of each target on each plate to equal the global median). Default
#' is to do IPC only. Any intensity normalization or bridge normalization steps
#' should typically be done after IPC normalization.
#' 
#' If the IPC median for a target is zero, it is replaced with 1. Similarly, if 
#' intensity normalization plate or global target median is zero, 
#' it is replaced with 1.
#' 
#' The intensity normalization option can also be used to do bridge sample-based
#' normalization by specifying \code{IN_samples}.
#'   
#' Input is a list of data matrices, one for each plate. 
#' These would typically be normData matrix from the intraPlateNorm function. 
#' Other required input depends on methods used. 
#' 
#' @description Use option \code{dataScale = log} when data is already on 
#' a log scale.
#' 
#' @description Output is a list.
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
#' @param dataScale 'count' is the default and interplate 
#' normalization is multiplicative. Use option 'log' for log-transformed
#' data; normalization is additive on the log scale. 
#' @param scaleFactor Optional numeric value used to rescale all data 
#' after IPC normalization. Default is 10^4. This shifts the data distribution 
#' to larger positive values, and helps prevent negative values after 
#' adding 1 and log-transforming. Only useful for count scale data.
#' @param transformReverse A vector of target names that use reverse curve
#' quantification. The reverse curve transformation will be applied to these targets.
#' @param transformReverse_scaleFactor The scaling factor used in the 
#' reverse curve transformation. Default is 1e4. Reverse curve transformation is 
#' \code{transformReverse_scaleFactor / (IPC normalized count + 1)}. Then the log2
#' tranformation is applied to this value, after adding 1, to obtain NPQ.
#' @return A list.
#' \item{interNormData}{A list of matrices of normalized count data (not 
#' log-transformed, unless input data was log-transformed, which should use `dataScale='log'`).}
#' \item{plate}{A vector of integers 1-(number of plates) 
#' corresponding to the interNormData matrix columns
#' that indicates which plate samples are from.} 
#' \item{normFactors}{A list of vectors with the normalization factor for each 
#' target in the corresponding plate}
#' \item{scaleFactor}{The scale factor that was used in the final scaling step.}
#'
#' 
#'
#' @export
#' 
interPlateNorm <- function(data_list,
                           IPC=TRUE,
                           IN=FALSE,
                           IPC_wells=NULL,
                           NC_wells=NULL,
                           IPC_method='median',
                           IN_samples=NULL,
                           dataScale='count',
                           scaleFactor=10^4,
                           transformReverse=NULL,
                           transformReverse_scaleFactor=10^4){

  # inter-plate control normalization
  if (IPC==TRUE){
    if(is.null(IPC_wells)){
      stop('Must provide IPC_wells for IPC normalization.')
    }
    IPC_factors <- list()
    data_list_IPC <- list()
    normFactors <- list()
    for (i in 1:length(data_list)){
      IPC_cols <- data_list[[i]][,IPC_wells[[i]]]
      if (IPC_method=='median'){
        IPC_factors_i <- apply(IPC_cols, 1, median, na.rm=TRUE)
      } else if (IPC_method=='mean'){
        IPC_factors_i <- rowMeans(IPC_cols, median, na.rm=TRUE)
      } else if (IPC_method=='geom_mean'){
        IPC_factors_i <- apply(IPC_cols, 1, function(x){
          # replace zeros or NA with 1s
          x[x==0] <- 1
          x[is.na(x)] <- 1
          prod(x, na.rm=TRUE)^(1/length(x))
        })
      }
      if (dataScale=='count'){
        # replace 0s or NAs with 1s
        # ADD WARNING HERE SAYING WHICH TARGETS THIS WAS DONE FOR
        IPC_factors_i[IPC_factors_i==0] <- 1
        IPC_factors_i[is.na(IPC_factors_i)] <- 1
        IPC_factors[[i]] <- IPC_factors_i
        data_list_IPC[[i]] <- data_list[[i]] / IPC_factors_i
        normFactors[[i]] <- 1 / IPC_factors_i
      } else if (dataScale=='log'){
        # replace NAs with 0s
        IPC_factors_i[is.na(IPC_factors_i)] <- 0
        IPC_factors[[i]] <- IPC_factors_i
        data_list_IPC[[i]] <- data_list[[i]] - IPC_factors_i
        normFactors[[i]] <- 1 / IPC_factors_i
      }
    } # end loop over data list
    data_list <- data_list_IPC
    # apply the scale factor to the data
    data_list <- lapply(data_list, function(x) x*scaleFactor)
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
          col_index <- 1:ncol(data_list[[i]])
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
    # replace zeros with 1s
    IN_medians[IN_medians==0] <- 1
    # combine data to calculate global medians
    all_IN_sample_data <- do.call(cbind, data_list_IN_samples)
    global_medians <- apply(all_IN_sample_data, 1, median, na.rm=TRUE)
    # replace zeros with 1s
    global_medians[global_medians==0] <- 1
    # save normFactors
    normFactors <- list()
    if (dataScale=='count'){
      # intensity normalization -- count scale
      IN_factors <- (1/IN_medians)*global_medians
      # IN rescale data
      for (i in 1:length(data_list)){
        data_list[[i]] <- data_list[[i]]*IN_factors[,i]
        normFactors[[i]] <- IN_factors[,i]
      } 
    } else if (dataScale=='log'){
      # intensity normalization -- log scale
      IN_factors <- global_medians - IN_medians
      # IN shift data
      for (i in 1:length(data_list)){
        data_list[[i]] <- data_list[[i]] + IN_factors[,i]
        normFactors[[i]] <- IN_factors[,i]
      }
    }
  } # end intensity normalization
  
  # generate the plate assignment vector
  plateNs <- unlist(lapply(data_list, ncol))
  plate <- NULL
  for (i in 1:length(plateNs)){
    plate <- c(plate, rep(i, plateNs[i]))
  }
  # apply the reverse curve transformation to specified targets
  if(!is.null(transformReverse)){
    reverses <- which(names(data_list[[i]][,1]) %in% transformReverse) 
    for (i in 1:length(data_list)){
      data_list[[i]][reverses, ] <- (transformReverse_scaleFactor) / (data_list[[i]][reverses,]/scaleFactor + 1)
    } 
  } 
  # log2 transform the output
  log2_interNormData <- lapply(data_list, function(x) log2(x + 1))
  
  # return output
  return(list(interNormData=data_list,
              log2_interNormData=log2_interNormData,
              plate=plate,
              normFactors=normFactors,
              scaleFactor=scaleFactor,
              transformReverse=transformReverse,
              transformReverse_scaleFactor=transformReverse_scaleFactor))
}
