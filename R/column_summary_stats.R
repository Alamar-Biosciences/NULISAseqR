#' Generate Summary Statistics for Columns of a Matrix or Data Frame
#'
#' Calculates mean, standard deviation, median, min, max, 
#' and number of missing values
#' for the columns of a matrix or for a vector. 
#'
#' @param x A matrix, data frame, or vector.
#' @param n_digits An optional vector of length 5 that indicates the number
#' of digits after the decimal place that will be used for rounding and 
#' formatting 
#' the mean, sd, median, min, and max, respectively. Default is
#' \code{c(2, 2, 2, 3, 1)]}. If NULL, no rounding or formatting will be done.
#' @param rowNames Names for the rows of output. Default will use 
#' column names of \code{x} as rows.
#'
#' @return A data frame with columns mean, sd, median, min, max, and 
#' number of missing values. Each row corresponds to the columns of the input. 
#'
#' 
#'
#' @export

column_summary_stats <- function(x, 
                                 n_digits=c(2,2,2,3,1),
                                 rowNames=colnames(x)){
  
  if(!is.array(x)) x <- data.frame(x)
  
  output <- cbind(colMeans(x, na.rm=TRUE),
                  apply(x, 2, sd, na.rm=TRUE),
                  apply(x, 2, median, na.rm=TRUE),
                  apply(x, 2, min, na.rm=TRUE),
                  apply(x, 2, max, na.rm=TRUE),
                  apply(x, 2, function(y){
                    sum(is.na(y))
                  }))
  colnames(output) <- c('mean', 'sd', 'median', 'min', 'max', 'missing')
  rownames(output) <- rowNames
  
  # Convert to data frame to allow mixed types
  output <- as.data.frame(output)

  if(!is.null(n_digits)){
    output[,1] <- format(round(as.numeric(output[,1]), n_digits[1]), nsmall=n_digits[1])
    output[,2] <- format(round(as.numeric(output[,2]), n_digits[2]), nsmall=n_digits[2])
    output[,3] <- format(round(as.numeric(output[,3]), n_digits[3]), nsmall=n_digits[3])
    output[,4] <- format(round(as.numeric(output[,4]), n_digits[4]), nsmall=n_digits[4])
    output[,5] <- format(round(as.numeric(output[,5]), n_digits[5]), nsmall=n_digits[5])
  }

  # Always format missing column as integers (use character to prevent decimal formatting in reactable)
  output[,6] <- as.character(as.integer(output[,6]))

  return(output)
}
