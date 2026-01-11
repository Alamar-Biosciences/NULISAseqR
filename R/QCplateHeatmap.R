#' Order samples by well position
#'
#' Orders samples either by well position alone or by sample type and well position.
#'
#' @param samples Data frame containing sample information with AUTO_WELLCOL and AUTO_WELLROW columns
#' @param sampleTypeFactor Optional factor to sort by sample type first
#' @param index Logical indicating whether to return indices or sample names
#' @return Ordered vector of sample names or indices
#' @keywords internal
wellorder <- function(samples, sampleTypeFactor = NULL, index = FALSE) {
  cols <- formatC(as.numeric(samples$AUTO_WELLCOL), width = 2, flag = 0)
  rows <- samples$AUTO_WELLROW
  inds <- if(is.null(sampleTypeFactor)) {
    sort(paste0(rows, cols), index.return = TRUE)
  } else {
    sort(paste0(as.numeric(sampleTypeFactor), rows, cols), index.return = TRUE)
  }
  if(index) {
    return(inds$ix)
  } else {
    return(samples$sampleName[inds$ix])
  }
}

#' Assign colors to wells based on their type
#'
#' @param data Data matrix with well names as column names
#' @param well_types List of well type patterns to match
#' @param boxplot_colors Vector of colors to assign
#' @return Vector of colors for each well
#' @keywords internal
wellcolors <- function(data, well_types, boxplot_colors) {
  colors <- rep("black", ncol(data))
  for (j in seq_along(well_types)) {
    indices <- which(grepl(paste0(well_types[[j]], collapse = "|"), colnames(data), fixed = TRUE))
    if(length(indices) > 0) {
      colors[indices] <- boxplot_colors[j + 1]
    }
  }
  return(colors)
}

#' Convert sample data to plate matrix format
#'
#' @param runSamples Sample data from a run
#' @param plate Logical indicating whether to fill empty wells with ""
#' @return Matrix representing the plate layout
#' @keywords internal
matrixify <- function(runSamples, plate = FALSE) {
  val <- matrix(rep(NA, 96), nrow = 8)
  colnames(val) <- 1:12
  rownames(val) <- LETTERS[1:8]
  for(j in LETTERS[1:8]) {
    for(k in 1:12) {
      sample_jk <- runSamples$sampleName[runSamples$AUTO_WELLROW == j & 
                                           as.numeric(runSamples$AUTO_WELLCOL) == k]
      if(length(sample_jk) != 0) {
        val[rownames(val) == j, colnames(val) == k] <- sample_jk
      }
      if(plate && length(sample_jk) == 0) {
        val[rownames(val) == j, colnames(val) == k] <- ""
      }
    }
  }
  return(val)
}

#' Generate plate heatmaps showing values relative to plate median
#'
#' Creates heatmaps of NULISAseq plate layouts showing:
#' 1) log2(total counts) as percent relative to plate median
#' 2) Internal control values as percent relative to plate median
#'
#' @param runs A named list of run data output from \code{loadNULISAseq()} function 
#' or a list of these outputs for multiple runs.
#' @param heatMapRel Logical indicating whether to show values as percent relative to 
#'                  plate median (TRUE) or absolute values (FALSE) (default: TRUE).
#' @param digitsIC Number of digits after the decimal points to display for IC values 
#' (default: 1 for relative, 0 for absolute). 
#' @param cex Numeric character expansion factor for text in heatmaps (default: 0.5).
#' @param cex.axis Numeric expansion factor for axis text (default: 0.5).
#' @param ic_name Name to use for internal control in plot titles (default: "IC").
#' @param save_plots Logical whether to save plots as PDF (default: TRUE).
#' @param output_dir Directory to save PDF files (NULL for no saving).
#'
#' @return Generates plate heatmaps.
#'
#' @examples
#' \dontrun{
#' # Show values as percent relative to plate median
#' plot_qc_plate_heatmaps_relative(runs)
#' 
#' # Show absolute values instead
#' plot_qc_plate_heatmaps_relative(runs, heatMapRel = FALSE)
#' }
#'
#' @export
#' @importFrom graphics par
QCplateHeatmap <- function(runs, heatMapRel = TRUE, 
                           digitsIC = if(heatMapRel) 1 else 0,
                           cex = 0.5, cex.axis = 0.5,
                           ic_name = "IC", 
                           save_plots = FALSE, 
                           output_dir = NULL) {
  
  # Input validation
  if(!is.list(runs)) stop("runs must be a list of run objects")
  if(!is.logical(heatMapRel)) stop("heatMapRel must be logical")
  
  # Create output directory if needed
  if (save_plots && !is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  for(i in seq_along(runs)) {
    # Check run structure
    if(!all(c("targets", "Data", "plateID", "samples") %in% names(runs[[i]]))) {
      stop(paste("Run", i, "doesn't have the required components"))
    }
    
    ICs <- which(tolower(runs[[i]]$targets$targetType) == "control")
    val <- matrixify(runs[[i]]$samples)
    vals <- unlist(as.list(t(val)))
    
    # Set up plotting parameters
    par(mfrow = c(ceiling((length(ICs) + 1)/2), 2), 
        oma = c(1, 1, 1, 1), 
        mar = c(2, 2, 2, 1))
    
    # Calculate and plot log2(total counts) relative to median
    logVals <- log2(colSums(runs[[i]]$Data, na.rm = TRUE))
    well_order <- wellorder(runs[[i]]$samples)
    val <- vals
    val[which(!is.na(vals))] <- logVals[well_order]
    
    # Save the plot if requested
    if (save_plots && !is.null(output_dir)) {
      pdf_file <- file.path(output_dir, paste0(runs[[i]]$plateID, "_log2_total_counts.pdf"))
      grDevices::pdf(pdf_file, width = 6, height = 4)
      on.exit(dev.off(), add = TRUE)  # Close the device after plotting
    }
    
    # Plot the heatmap for log2 values
    plateHeatmap(as.numeric(val), 
                 title = paste0(runs[[i]]$plateID, ": log2(total counts)\n(% relative to plate median)"), 
                 cex = cex, digits = 1,
                 relative = heatMapRel, cex.axis = cex.axis)
    
    # Plot ICs relative to median
    for (j in seq_along(ICs)) {
      val <- vals
      val[which(!is.na(vals))] <- runs[[i]]$Data[ICs[j],][well_order]
      current_ic_name <- if(length(ICs) == 1) ic_name else runs[[i]]$IC[j]
      
      # Save the plot for each IC if requested
      if (save_plots && !is.null(output_dir)) {
        pdf_file_ic <- file.path(output_dir, paste0(runs[[i]]$plateID, "_", current_ic_name, "_relative_to_median.pdf"))
        grDevices::pdf(pdf_file_ic, width = 6, height = 4)
        on.exit(dev.off(), add = TRUE)  # Close the device after plotting
      }
      
      # Plot the heatmap for the ICs
      plateHeatmap(as.numeric(val), 
                   title = paste0(runs[[i]]$plateID, ": ", current_ic_name, 
                                  "\n(% relative to plate median)"), 
                   cex = cex, digits = digitsIC, 
                   relative = heatMapRel, cex.axis = cex.axis)
    }
  }
}