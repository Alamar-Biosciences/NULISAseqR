#' Detect Outliers using MAD (Median Absolute Deviation)
#'
#' This function detects outliers in a numeric dataset using the Median Absolute Deviation (MAD) method.
#'
#' @param data A numeric vector or data frame containing the data from which outliers should be detected.
#' @param threshold A numeric value representing the threshold for identifying outliers. Default is 2.5.
#' @param maxRemoval An integer specifying the maximum number of the worst outliers to remove. Default is 0 (no removal).
#' @param returnVal A character string specifying the type of result to return. Options are "indices" (default), "outliers", or "keeps".
#'
#' @details The function calculates the MAD, identifies outliers based on the threshold, and can optionally remove a specified number of the worst outliers.
#'
#' @return Depending on the 'returnVal' parameter, the function returns:
#'   - "indices": A numeric vector of the indices of the detected outliers.
#'   - "outliers": A numeric vector of the values of the detected outliers.
#'   - "keeps": A numeric vector of the values of the data points that are not identified as outliers.
#'   - "replaceNA": A numeric vector of the values of the data points with outliers replaced with NA.
#'
#' @examples
#' # Example usage
#' your_data <- c(10, 12, 13, 15, 9, 11, 100, 14, 16, 8, 7)
#' outlier_indices <- outliers_mad(your_data, threshold = 2, maxRemoval = 2, returnVal = "indices")
#' cat("Outlier indices: ", outlier_indices, "\n")
#' removed_outliers <- outliers_mad(your_data, threshold = 2, maxRemoval = 2, returnVal = "keeps")
#' cat("Data after removing outliers: ", removed_outliers, "\n")
#'
#' @seealso
#' 'mad' function from the 'stats' package for calculating the MAD.
#'
#' @export
outliers_mad <- function(data, threshold = 2.5, maxRemoval = 0, returnVal = "indices") {
  if(maxRemoval + 1 > length(data)){
    stop("Error: maxRemoval must be less than the number of data points provided!")
  }
  # Calculate the median of the data
  median_value <- median(data)
  
  # Calculate the MAD (Median Absolute Deviation)
  mad_value <- mad(data)
  
  # Calculate the modified Z-score
  z_scores <- 0.6745 * (data - median_value) / mad_value

  # Identify outliers based on the threshold
  outliers <- abs(z_scores) > threshold
 
  # Limit the number of outliers to remove based on maxRemoval
  if (maxRemoval > 0) {
    sorted_indices <- order(abs(z_scores), decreasing = TRUE)
    removal_indices <- sorted_indices[(maxRemoval+1):length(data)]
    outliers[removal_indices] <- FALSE
  }

  # Return from the function
  if(returnVal == "indices"){ # Return the indices of the outliers
    return(which(outliers))
  } else if(returnVal == "outliers"){ # Return the values of outliers
    return(data[which(outliers)])
  } else if(returnVal == "keeps"){  # Return the values of data to keep
    return(data[-which(outliers)])
  } else if(returnVal == "replaceNA"){ # Return the original dataset but with outliers replaced by NA
    data[which(outliers)] = NA
    return(data)
  } else{ # Return an error
    stop(paste0("Error: Invalid returnVal function argument \"", returnVal, "\"!"))
  }
}

# Example usage
# Replace 'your_data' with your actual dataset
your_data <- c(10, 12, 13, 15, 9, 11, 100, 14, 16, 8, 7)
#your_data <- c(1000,126,282,198,-1000)

# Detect outliers in the dataset
outlier_indices <- outliers_mad(your_data)
cat("Outlier indices: ", outlier_indices, "\n")
outlier_indices <- outliers_mad(your_data, returnVal="outliers")
cat("Outlier values: ", outlier_indices, "\n")
outlier_indices <- outliers_mad(your_data, returnVal="keeps")
cat("Keep values: ", outlier_indices, "\n")
#outlier_indices <- detect_outliers_mad(your_data, threshold = 2, returnVal="kees")
#cat("Outlier values: ", outlier_indices, "\n")
