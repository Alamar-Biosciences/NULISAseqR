# Convenience function to replace outliers using MAD
replace_outliers_mad <- function(col, min_blanks = 4, threshold = 2.5) {
  if (length(col) >= min_blanks) {
    outliers_mad(col, threshold = threshold, maxRemoval = 1, returnVal = "replaceNA")
  } else {
    col
  }
}

# Convenience function to return outlier index using MAD
outliers_index_mad <- function(col, min_blanks = 4, threshold = 2.5) {
  if (length(col) >= min_blanks) {
    ind <- outliers_mad(col, threshold = threshold, maxRemoval = 1, returnVal = "indices")
    if(length(ind) > 0) return(ind)
    if(length(ind) == 0) return(NA)
  } else {
    return(NA)
  }
}

#' Calculate Limits of Detection
#'
#' Calculates limit of detection (LoD) for each target based on the
#' negative controls (blanks). LoD = mean(blanks) + 3*SD(blanks).
#' Designates data as either above or below LoD. 
#' Option to specify minimum count threshold for detectability.
#'
#' @param blanks Column indices or column names of the blanks in the
#' data_matrix.
#' @param data_matrix The Data matrix output from readNULISAseq.R
#' or normalized data from normalization functions.
#' @param min_count Optional count threshold to apply in addition
#' to the LoD. Default is 0.
#' @param min_blank_no Optional numeric parameter defining the minimum number of
#' blanks required to enable MAD outlier detection and removal. Default is 4.
#' @param mad_threshold Optional numeric parameter defining the threshold used for
#' MAD outlier identification. Default is 2.5.
#' @param ignore_target_blank List of targets and corresponding blanks/NCs to exclude.
#' Names of the list represent targets and
#' values which are arrays represent blank/NC names to ignore during LoD calculation.
#' @param targetNoOutlierDetection Option to provide targets which should NOT have 
#' outlier detection applied
#' Defaults to NULL.
#' @param match_matrix Matrix of indices provided by calcSampleTargetNAs. 
#' Lists samples/targets that should not be reported 
#'
#'
#' @return A list.
#' @param LOD Vector of limits of detection.
#' @param aboveLOD Logical matrix indicating whether counts are 
#' above or below LoD for that target.
#'
#' @export
#' 
lod <- function(data_matrix, blanks, min_count = 0, min_blank_no = 4, mad_threshold = 2.5, ignore_target_blank = NULL, targetNoOutlierDetection = NULL, match_matrix = NULL) {
  # Determine blank names if blank indices are provided
  if (is.numeric(blanks)) {
    blank_names <- colnames(data_matrix)[blanks]
  }

  # Determine blank indices if blank names are provided
  # Store blank names in a temp var for later use
  if (!is.numeric(blanks)) {
    blank_names <- blanks
    blanks <- which(colnames(data_matrix) %in% blanks)
  }

  blank_data <- data_matrix[, blanks]

  # Exclude negative controls/blanks by targets
  # This helps handling NC spike outs in a target specific manner
  # Set the count value to NA based on user provided target/blank combinations
  # Only existing blanks will be updated
  if (!is.null(ignore_target_blank) && !is.null(names(ignore_target_blank))) {
    for (target in names(ignore_target_blank)) {
      ignore_blank_names <- ignore_target_blank[[target]]

      # Fail safe step: Avoid trying to find an blank that's not within original blanks
      ignore_blank_names <- intersect(blank_names, ignore_blank_names)

      # Determine the remaining number of blanks
      # Remaining number of blanks should be greater than or equal to 2
      # Otherwise an LoD cannot be calculated
      # Therefore we should not replace values with NA
      remaining_blank_no <- length(blank_names) - length(ignore_blank_names)

      # Replacing blank count values with NA in a target specific manner
      if (length(ignore_blank_names) > 0 && remaining_blank_no >= 2) {
        ignore_blank_names_str <- paste(ignore_blank_names, collapse = ", ")
        message(paste("Replacing", target, "NC", ignore_blank_names_str, "value with NA!"))
        blank_data[target, ignore_blank_names] <- NA
      }
    }
  }

  # Compute blank mean
  # Apply MAD outlier removal: Gets applied only when blanks > min_blank_no
  blank_mean <- colMeans(
    apply(blank_data, 1, function(row) replace_outliers_mad(row, min_blanks = min_blank_no, threshold = mad_threshold)),
    na.rm = TRUE
  )
  
  # create table with NC name of outliers by target
  blank_outliers <- apply(blank_data, 1, function(row) outliers_index_mad(row, min_blanks = min_blank_no, threshold = mad_threshold))
  blank_outliers <- blank_outliers[!is.na(blank_outliers)]
  if(length(blank_outliers) > 0){
    blank_outlier_table <- data.frame(Target=names(blank_outliers),
                                      NC_outlier=blank_names[blank_outliers],
                                      NC_outlier_value=round(blank_data[cbind(names(blank_outliers),blank_names[blank_outliers])], 1),
                                      NC_mean=round(blank_mean[names(blank_outliers)], 1))
  } else {
    blank_outlier_table <- NULL
  }

  # Compute blank SD
  blank_sd <- apply(
    apply(blank_data, 1, function(row) replace_outliers_mad(row, min_blanks = min_blank_no, threshold = mad_threshold)),
    2, sd, na.rm = TRUE
  )

  # Calculate the blank_mean and blank_sd assuming no outlier detection (blank_meanCount, blank_sdCount)
  # replace blank_mean and blank_sd with blank_meanCount / blank_sdCount for targets where no outlier 
  # deteection is desired
  if(!is.null(targetNoOutlierDetection)){
    blank_meanCount <- rowMeans(blank_data, na.rm = TRUE )
    blank_sdCount <- apply(blank_data, 1, sd, na.rm = TRUE)
    blank_mean[targetNoOutlierDetection] <- blank_meanCount[targetNoOutlierDetection]
    blank_sd[targetNoOutlierDetection] <- blank_sdCount[targetNoOutlierDetection]
    # remove these targets from blank outlier table
    if(!is.null(blank_outlier_table)){
      blank_outlier_table <- blank_outlier_table[!(blank_outlier_table$Target %in% targetNoOutlierDetection), ]
    }
  }

  # Calculate LoD
  LOD <- blank_mean + 3 * blank_sd
  if(min_count > 0){
    LOD[LOD < min_count] <- min_count
  }
  names(LOD) <- rownames(data_matrix)
  aboveLOD <- data_matrix
  aboveLOD <- apply(aboveLOD, 2, function(x){
    result <- x > LOD
    result[is.na(result)] <- FALSE
    return(result)
  })

  # Do not report values for a specific target / sample_matrix combination if specified in barcodeA
  if(!is.null(match_matrix) && nrow(match_matrix) > 0 ) {
    aboveLOD[cbind(match_matrix[,"row"], match_matrix[,"col"])] <- NA
  }
  return(list(LOD=LOD, 
              aboveLOD=aboveLOD,
              blank_outlier_table=blank_outlier_table))
}
