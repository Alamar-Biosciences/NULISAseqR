#' NULISAseq Intra-Plate Normalization
#'
#' Does intra-plate normalization either using single internal control, 
#' multiple internal controls (using geometric mean), or total data_matrix. 
#' Input is the output of readNULISAseq.R Data matrix.
#'
#' @param data_matrix The Data matrix output from readNULISAseq.R.
#' @param method \code{'IC'} (default) divides each target count by 
#' the specified IC count for a given sample well; 
#' if multiple ICs are specified uses geometric mean of several ICs.
#' \code{'total_count'} scales each sample well so data_matrix sum to 
#' one million (except NC wells and omitted target data_matrix).
#' @param IC Required for 'single' and 'geom_mean' methods. 
#' Vector of either the row(s) (numeric) or the 
#' row name(s) (character string) of the internal control(s) to be used 
#' for normalization. Raw data_matrix will be divided by IC data_matrix or 
#' by the geometric mean of IC data_matrix.
#' @param NC Required for 'total_count' method. Vector of either 
#' the column(s) (numeric) or the column name(s) (character string) 
#' of the internal control(s).
#' @param TC_omit Required for 'total_count' method only.
#' Vector of either the row indices (numeric) or the 
#' row names (character string) of any targets that should be omitted
#' from the total count normalization. For example, it may be desirable to
#' omit IC data_matrix if IC data_matrix are high. 
#' Omitted targets will be rescaled but will not be used in calculating
#' the scaling factors.
#' @param scaleFactor. Optional numeric value used to rescale data 
#' after normalizing. Default is 1.
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
                           method='single',
                           IC=NULL,
                           TC_omit=NULL,
                           scaleFactor=1){
  if (method=='single' | method=='geom_mean'){
    if (is.null(IC)){
      stop('Must specify ICs.')
    }
    IC_missing <- sum(is.na(data_matrix[IC,])) + sum(data_matrix[IC,]==0, na.rm=TRUE)
    if (IC_missing > 0) {
      cat(paste0('Warning: ' , IC_missing, ' missing or zero values in internal control data. 
                 Normalized data for these samples will all be missing values.'))
    }
  }
  if (method=='single'){
    # normalize using single internal control
    if (length(IC) > 1){
      stop('Must specify only one IC for method single.')
    }
    normFactor <- scaleFactor/data_matrix[IC,]
    normCounts <- t(t(data_matrix)*normFactor)
  } else if (method=='geom_mean'){
    # normalize using geometric mean of internal controls
    normFactor <- scaleFactor/apply(data_matrix[IC,], 2, function(x) {
      prod(x)^(1/length(x))
    })
    normCounts <- t(t(data_matrix)*normFactor)
  } else if (method=='total_count'){
    # rescale each row to sum to 10^6 (excluding "TC_omit" columns)
    if (!is.null(TC_omit)){
      data_matrix_TC <- data_matrix[-TC_omit,]
    } else {
      data_matrix_TC <- data_matrix
    }
    totalCounts <- colSums(data_matrix_TC, na.rm=TRUE)
    normFactor <- scaleFactor*(10^6/totalCounts)
    normCounts <- t(t(data_matrix)*normFactor)
  }
  # return normalized count matrix
  return(normCounts)
}
