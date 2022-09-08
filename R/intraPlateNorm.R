#' NULISAseq Intra-Plate Normalization
#'
#' Does intra-plate normalization either using single internal control, 
#' multiple internal controls (using geometric mean), or total data_matrix. 
#' Input is the output of readNULISAseq.R Data matrix.
#'
#' @param data_matrix The Data matrix output from readNULISAseq.R.
#' @param method \code{'IC'} (internal control; default) divides 
#' each target count by the specified IC count for a given sample well; 
#' if multiple ICs are specified uses geometric mean of several ICs.
#' \code{'TC'} (total count)  scales each sample well so data_matrix sum to 
#' one million (except NC wells and omitted target data_matrix).
#' @param IC Required for 'single' and 'geom_mean' methods. 
#' Vector of either the row(s) (numeric) or the 
#' row name(s) (character string) of the internal control(s) to be used 
#' for normalization. Raw data_matrix will be divided by IC data_matrix or 
#' by the geometric mean of IC data_matrix.
#' @param NC Recommended for 'total_count' method when negative control 
#' wells are present. Vector of either 
#' the column(s) (numeric) or the column name(s) (character string) 
#' of the negative control(s). For total count method, NCs will be scaled such
#' that the (total non-NC count : total NC count) ratio is the same 
#' in both the un-normalized and the normalized data.
#' @param TC_omit Option for 'total_count' method only.
#' Vector of either the row indices (numeric) or the 
#' row names (character string) of any targets, if any, 
#' that should be omitted
#' from the total count normalization. For example, it may be desirable to
#' omit IC data_matrix if IC data_matrix are high. 
#' Omitted targets will be rescaled but will not be used in calculating
#' the scaling factors.
#' @param scaleFactor. Optional numeric value used to rescale data 
#' after normalizing. Default is 1. This may be desirable if there are 
#' many targets with read counts less than the IC counts, to avoid
#' normalized quantities between 0 and 1 (which will be negative
#' in the log scale).
#'
#' @return 
#' @param normData A matrix of normalized count data (not log-transformed).
#' @param normFactor A vector of the normalization factors used for 
#' each sample (column) in the data matrix.
#'
#' @examples
#' 
#'
#' @export
#' 
intraPlateNorm <- function(data_matrix, 
                           method='IC',
                           IC=NULL,
                           NC=NULL,
                           TC_omit=NULL,
                           scaleFactor=1){
  
  if (method=='IC'){
    if (is.null(IC)){
      stop('Must specify ICs.')
    }
    IC_missing <- sum(is.na(data_matrix[IC,])) + sum(data_matrix[IC,]==0, na.rm=TRUE)
    if (IC_missing > 0) {
      cat(paste0('Warning: ' , IC_missing, ' missing or zero values in internal control data. 
                 Normalized data for these samples will all be missing values.'))
    }
    if (length(IC) == 1){
      # normalize using single internal control
      cat('Using single IC to normalize data.')
      normFactor <- scaleFactor/data_matrix[IC,]
      normCounts <- t(t(data_matrix)*normFactor)
    } else if (length(IC) > 1){
      # normalize using geometric mean of internal controls
      cat('Using geometric mean of multiple ICs to normalize data.')
      normFactor <- scaleFactor/apply(data_matrix[IC,], 2, function(x) {
        prod(x)^(1/length(x))
      })
      normCounts <- t(t(data_matrix)*normFactor)
    }
  } else if (method=='TC'){
    # rescale each row to sum to 10^6 (excluding "TC_omit" columns)
    if (!is.null(TC_omit)){
      if (!is.numeric(TC_omit[1])){
        TC_omit <- which(rownames(data_matrix) %in% TC_omit)
      }
      data_matrix_TC <- data_matrix[-TC_omit,]
    } else {
      data_matrix_TC <- data_matrix
    }
    totalCounts <- colSums(data_matrix_TC, na.rm=TRUE)
    if (!is.null(NC)){
      # calculate non-NC to NC count ratio
      if (!is.numeric(NC[1])){
        NC <- which(colnames(data_matrix) %in% NC)
      }
      NC_counts <- sum(data_matrix_TC[,NC], na.rm=TRUE)
      non_NC_counts <- sum(data_matrix_TC[,-NC], na.rm=TRUE)
      NC_count_ratio <- NC_counts/non_NC_counts
      n_NC <- length(NC)
      n_non_NC <- ncol(data_matrix) - n_NC
      NC_scaleFactor <- rep(1, length(totalCounts))
      NC_scaleFactor[NC] <- NC_count_ratio*(n_non_NC/n_NC)
    } else {
      NC_scaleFactor <- rep(1, length(totalCounts))
    }
    normFactor <- scaleFactor*(10^6/totalCounts)*NC_scaleFactor
    normCounts <- t(t(data_matrix)*normFactor)
  }
  # return normalized count matrix
  return(list(
    normData=normCounts,
    normFactor=normFactor
    ))
}
