#' Generate Sample Distribution Boxplots
#'
#' Creates horizontal or vertical boxplots showing sample distributions (normalized and/or unnormalized) with internal 
#' control reference lines. Options to save as PDF.
#'
#' @param runs A named list of run data output from \code{loadNULISAseq()} function 
#' or a list of these outputs for multiple runs.
#' @param metrics Character specifying which data to plot: "unnormed", "normed", or "all" (default: "normed").
#' @param plateIDs Optional vector of plate names (default: extract from runs).
#' @param output_dir Directory to save PDF files (default: NULL).
#' @param plot_name Optional base name for output files, uses plateIDs automatically (default: NULL).
#' @param plot_title Optional custom plot title in additional to plateIDs (default: NULL).
#' @param horizontal Logical indicating horizontal (TRUE) or vertical (FALSE) orientation (default: TRUE).
#' @param internal_control_label Label for internal control line (default: "Internal Control").
#' @param plot_width Optional plot width in inches (default: NULL uses automatic sizing).
#' @param plot_height Optional plot height in inches (default: NULL uses automatic sizing).
#' @param title_cex Numeric value controlling title size (default: 0.65).
#' @param axis_label_cex Numeric value controlling axis label size (default: 0.5).
#' @param data_label_cex Numeric value controlling data label size (default: 0.3).
#' @param IC_line_width Numeric value controlling internal control line width (default: 1).
#' @param sample_subset_list Optional list of sample subsets for each plate. Can be either:
#'   \itemize{
#'     \item{A list of character vectors - each vector contains sample names to include for the corresponding plate 
#'           (e.g., `list(c('A_01_C001 P3', 'A_02_C044 P1'), c('A_01_C040 P1'))` for two plates)}
#'     \item{A single character vector - only subsets samples from the first plate}
#'   }
#'   Note: Control samples (IPC, NC, SC, Bridge) are always included regardless of subsetting.
#'   (default: NULL)
#'   
#' @param target_subset_list Optional list of target subsets for each plate. Can be either:
#'   \itemize{
#'     \item{A list of character vectors - each vector contains target names to include for the corresponding plate}
#'     \item{A single character vector - only subsets targets from the first plate}
#'   }
#'   (default: NULL)
#'   
#' @param sample_exclude_list Optional list of samples to exclude for each plate. Can be either:
#'   \itemize{
#'     \item{A list of character vectors - each vector contains sample names to exclude from the corresponding plate}
#'     \item{A single character vector - only excludes samples from the first plate}
#'   }
#'   Note: Control samples (IPC, NC, SC, Bridge) are never excluded.
#'   (default: NULL)
#'   
#' @param target_exclude_list Optional list of targets to exclude for each plate. Can be either:
#'   \itemize{
#'     \item{A list of character vectors - each vector contains target names to exclude from the corresponding plate}
#'     \item{A single character vector - only excludes targets from the first plate}
#'   }
#'   (default: NULL)
#'
#' @return A list of plot data and saved PDF files in the specified output directory.
#' @export
#' @importFrom graphics boxplot axis lines mtext legend abline par
#' @importFrom grDevices hcl.colors pdf dev.off
sampleBoxplot <- function(runs,
                          metrics = c("unnormed", "normed", "all"),
                          plateIDs = NULL,
                          output_dir = NULL,
                          plot_name = NULL,
                          plot_title = NULL,
                          horizontal = TRUE,
                          internal_control_label = "Internal Control",
                          plot_width = NULL,
                          plot_height = NULL,
                          title_cex=0.65,       
                          axis_label_cex=0.5,   
                          data_label_cex=0.3,   
                          IC_line_width=1,
                          sample_subset_list = NULL,
                          target_subset_list = NULL,
                          sample_exclude_list = NULL,
                          target_exclude_list = NULL) {
  
  # Input validation checks
  if (!is.list(runs)) {
    stop("runs must be a list of run objects")
  }
  
  if (length(runs) == 0) {
    stop("runs list cannot be empty")
  }
  
  required_elements <- c("samples", "Data", "normed", "plateID", "IC", "ExecutionDetails")
  missing_elements <- setdiff(required_elements, unique(unlist(lapply(runs, names))))
  if (length(missing_elements) > 0) {
    stop("Each run object must contain: ", paste(missing_elements, collapse = ", "))
  }
  
  
  if(length(metrics) > 1 || !metrics %in% c("unnormed", "normed", "all")) {
    stop("Invalid metrics argument. Use 'unnormed', 'normed', or 'all'")
  }
  
  # If not saving to files, plotting will go to the current device.
  # (Removed graphics.off() to avoid interfering with user's devices)
  
  # Create output directory if needed
  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get plate names if not provided
  if (is.null(plateIDs)) {
    plateIDs <- vapply(runs, function(x) x$plateID, character(1))
  }
  # Get IC name (assuming same across runs)
  IC_name <- Reduce(intersect, lapply(runs, function(x) x$IC))
  
  # Set default colors 
  boxplot_colors <- grDevices::hcl.colors(5, palette = "Set3")[c(4, 3, 1, 5, 2)]
  label_colors <- grDevices::hcl.colors(5, palette = "Set3")[c(4, 3, 5, 2, 1)]
  label_colors[1] <- 'black'
  
  # Set default plot dimensions if not provided
  if (is.null(plot_width)) {
    plot_width <- if (metrics == "all") {
      if (horizontal == TRUE) 9 else 8
    } else {
      if (horizontal == TRUE) 9 else 4
    }
  }
  if (is.null(plot_height)) {
    plot_height <- if (metrics == "all") {
      if (horizontal == TRUE) 8 else 9
    } else {
      if (horizontal == TRUE) 4 else 9
    }
  }
  
  # Process each plate
  for (i in seq_along(runs)) {
    
    # Get filters for current plate
    current_sample_subset <- if (is.list(sample_subset_list)) {
      if (i <= length(sample_subset_list)) sample_subset_list[[i]] else NULL
    } else if (i == 1) {
      sample_subset_list
    } else {
      NULL
    }
    
    current_sample_exclude <- if (is.list(sample_exclude_list)) {
      if (i <= length(sample_exclude_list)) sample_exclude_list[[i]] else NULL
    } else if (i == 1) {
      sample_exclude_list
    } else {
      NULL
    }
    
    current_target_subset <- if (is.list(target_subset_list)) {
      if (i <= length(target_subset_list)) target_subset_list[[i]] else NULL
    } else if (i == 1) {
      target_subset_list
    } else {
      NULL
    }
    
    current_target_exclude <- if (is.list(target_exclude_list)) {
      if (i <= length(target_exclude_list)) target_exclude_list[[i]] else NULL
    } else if (i == 1) {
      target_exclude_list
    } else {
      NULL
    }
    
    # Apply filters
    filtered_run <- filter_run_data(
      runs[[i]],
      sample_subset = current_sample_subset,
      sample_exclude = current_sample_exclude,
      target_subset = current_target_subset,
      target_exclude = current_target_exclude
    )
    
    if (nrow(filtered_run$samples) == 0 || nrow(filtered_run$targets) == 0) {
      warning("No data remaining after filtering for Plate ", filtered_run$plateID, ". Skipping the plate.")
      next
    }
    
    # Determine if absolute quantification data
    AbsData <- "Abs" %in% names(filtered_run$ExecutionDetails)
    CAL_IPC <- ifelse(AbsData, "CAL", "IPC")
    
    # Get sample types and ordering
    sample_type <- filtered_run$samples$sampleType
    sample_type <- factor(sample_type, levels = c('Sample', 'IPC', 'SC', 'Bridge', 'NC'))
    
    # Apply colors based on sample type
    boxplot_colors2 <- boxplot_colors[c(1, 2, 4, 5, 3)][as.numeric(sample_type)]
    label_colors2 <- label_colors[as.numeric(sample_type)]
    
    # Get well ordering
    ordering2 <- wellorder(filtered_run$samples, index = TRUE)
    ordering <- NULL
    types <- sort(unique(sample_type))
    for (j in seq_along(types)) {
      inds <- which(sample_type[ordering2] == types[j])
      ordering <- c(ordering, ordering2[inds])
    }
    ordering_v <- rev(ordering)  # Reverse so special wells are at bottom
    
    if(horizontal==TRUE) {
      final_ordering <- rev(ordering_v)
    } else {
      final_ordering <- ordering_v
    }
    
    # Handle absolute quantification naming
    if (AbsData) {
      colnames(filtered_run$Data) <- gsub("^IPC_([1234])$", "CAL_\\1", colnames(filtered_run$Data))
      colnames(filtered_run$Data) <- gsub("^SC_([1234])$", "AQSC_\\1", colnames(filtered_run$Data))
      colnames(filtered_run$normed$interNormData[[1]]) <- gsub("^IPC_([1234])$", "CAL_\\1", colnames(filtered_run$normed$interNormData[[1]]))
      colnames(filtered_run$normed$interNormData[[1]]) <- gsub("^SC_([1234])$", "AQSC_\\1", colnames(filtered_run$normed$interNormData[[1]]))
    }
    
    # Set up PDF output if requested
    if (!is.null(output_dir)) {
      if(!is.null(plot_name)) {
        pdf_file <- file.path(output_dir, paste0("sample_boxplot_", metrics, "_", plateIDs[i], "_", plot_name, ".pdf"))
      } else {
        pdf_file <- file.path(output_dir, paste0("sample_boxplot_", metrics, "_", plateIDs[i], ".pdf"))
      }
      grDevices::pdf(pdf_file, width = plot_width, height = plot_height)
      on.exit(grDevices::dev.off())
    }
    
    # Set up plot parameters based on metrics
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar), add = TRUE)
    
    title_pre <- ifelse(!is.null(plot_name), 
                        paste0(plot_name, ", ", plateIDs[i]), 
                        plateIDs[i])
    
    if (metrics == "all") {
      # Force vertical for side-by-side plots
      if (horizontal==TRUE) {
        par(mfrow = c(2, 1), mar = c(6, 4, 2, 1))
      } else {
        par(mfrow = c(1, 2), mar = c(2, 2, 2, 1))
      }
      
      sampleboxplot(filtered_run$Data, data_order=final_ordering, 
                    title=paste0(title_pre, ": unnormalized"), 
                    IC_name=IC_name,
                    boxplot_colors=boxplot_colors2,
                    label_colors=label_colors2,
                    data_axis_label = "log2(count + 1)",
                    horiz = !horizontal,
                    title_cex = title_cex,        
                    axis_label_cex = axis_label_cex,   
                    data_label_cex = data_label_cex,   
                    IC_line_width = IC_line_width)     
      
      sampleboxplot(filtered_run$normed$interNormData[[1]], data_order=final_ordering, 
                    title=paste0(title_pre, ": IC + ", CAL_IPC, " normalized"), 
                    IC_name=IC_name,
                    boxplot_colors=boxplot_colors2,
                    label_colors=label_colors2,
                    data_axis_label = "NPQ",
                    horiz = !horizontal,
                    title_cex = title_cex,        
                    axis_label_cex = axis_label_cex,   
                    data_label_cex = data_label_cex,   
                    IC_line_width = IC_line_width)  
    } 
    else if (metrics == "unnormed"){
      plot_suffix <- ": Sample Distributions Before Normalization"
      if (horizontal==TRUE) {
        par(mfrow = c(1, 1), mar = c(6, 4, 2, 1))
      } else {
        par(mfrow = c(1, 1), mar = c(6, 2, 2, 2))
      }
      suppressWarnings({
        sampleboxplot(filtered_run$Data, data_order=final_ordering, 
                      title=paste0(title_pre, plot_suffix), 
                      IC_name=IC_name,
                      boxplot_colors=boxplot_colors2,
                      label_colors=label_colors2,
                      horiz=!horizontal, 
                      data_axis_label = "log2(count + 1)",
                      title_cex = title_cex,        
                      axis_label_cex = axis_label_cex,   
                      data_label_cex = data_label_cex,   
                      IC_line_width = IC_line_width)  
      })
    }
    else if (metrics == "normed"){
      plot_suffix <- ": Sample Distributions After Normalization"
      if (horizontal==TRUE) {
        par(mfrow = c(1, 1), mar = c(6, 4, 2, 1))
      } else {
        par(mfrow = c(1, 1), mar = c(6, 2, 2, 2))
      }
      suppressWarnings({
        sampleboxplot(filtered_run$normed$interNormData[[1]], data_order=final_ordering, 
                      title=paste0(title_pre, plot_suffix), 
                      IC_name=IC_name,
                      boxplot_colors=boxplot_colors2,
                      label_colors=label_colors2,
                      data_axis_label = "NPQ",
                      horiz=!horizontal,
                      title_cex = title_cex,        
                      axis_label_cex = axis_label_cex,   
                      data_label_cex = data_label_cex,   
                      IC_line_width = IC_line_width)  
      })
    }
    # Save plot as PDF for each plate
    if (!is.null(output_dir)) {
      dev.off()  # Close the current PDF device after saving the plot
    }
  }
  
  invisible(NULL)
}


#' Sample Boxplot Visualization base function
#'
#' Internal function to creates boxplots of sample data with optional internal control reference lines. 
#' Used internally by \code{\link{sampleBoxplot}}.
#' 
#' @param data Matrix or data frame containing the data to plot (samples in columns).
#' @param data_order Numeric vector specifying the order of samples to plot.
#' @param title Character string for the plot title.
#' @param IC_name Character string specifying which column represents the internal control.
#' @param boxplot_colors Vector of colors for the boxplots (one per sample).
#' @param label_colors Vector of colors for the sample labels (one per sample).
#' @param horiz Logical indicating whether to plot horizontally (default: TRUE).
#' @param data_axis_limits Optional numeric vector of length 2 specifying axis limits.
#' @param internal_control_label Label for the internal control legend (default: "Internal Control").
#' @param data_axis_label Label for the data axis (default: "NPQ").
#' @param title_cex Numeric value controlling title size (default: 0.65).
#' @param axis_label_cex Numeric value controlling axis label size (default: 0.5).
#' @param data_label_cex Numeric value controlling data label size (default: 0.3).
#' @param IC_line_width Numeric value controlling internal control line width (default: 1).
#'
#' @return Invisibly returns NULL. Produces a boxplot.
#' @keywords internal
#' @noRd
sampleboxplot <- function(data, 
                          data_order, 
                          title, 
                          IC_name,
                          boxplot_colors, 
                          label_colors, 
                          horiz=TRUE,
                          data_axis_limits=NULL,
                          internal_control_label='Internal Control',
                          data_axis_label='NPQ',
                          title_cex=0.65,       
                          axis_label_cex=0.5,   
                          data_label_cex=0.3,   
                          IC_line_width=1) {   
  suppressWarnings({
    if(horiz==TRUE){par(mar=c(2,6,2,0.5))} else{par(mar=c(6,2,2,0.5))}
    if(horiz==TRUE){
      sample_axis_side <- 2
      data_axis_side <- 1
      xlimits <- NULL
      ylimits <- data_axis_limits
    } else if(horiz==FALSE){
      sample_axis_side <- 1
      data_axis_side <- 2
      xlimits <- NULL
      ylimits <- data_axis_limits
    }
    
    boxplot(log2(data[,data_order] + 1),
            xlim=xlimits,
            ylim=ylimits,
            las=1, 
            yaxt='n',
            xaxt='n',
            ylab='',
            main=title,
            outcex=0.5,
            col=boxplot_colors[data_order], 
            horizontal=horiz, 
            cex=0.5,
            cex.main=title_cex,        
            cex.lab=0.75)
    # draw axes
    for (j in 1:ncol(data)){
      axis(sample_axis_side, at=j, labels=FALSE, cex.axis=data_label_cex, las=2, 
           col.axis=label_colors[data_order[j]], tck=-0.01)
      axis(sample_axis_side, at=j, labels=colnames(data)[data_order[j]], 
           cex.axis=data_label_cex, las=2, col.axis=label_colors[data_order[j]], 
           tick=FALSE, line=-0.75)
    }
    axis(side=data_axis_side, labels=FALSE, tck=-0.01)
    axis(side=data_axis_side, labels=TRUE, cex.axis=axis_label_cex, line=-1, tick=FALSE)
    # add internal control line
    if (horiz==TRUE) {
      lines(cbind(log2(data[IC_name, data_order] + 1), 1:ncol(data)), 
            col='red', las=1, lwd=IC_line_width)  
    } else {
      lines(cbind(1:ncol(data), log2(data[IC_name, data_order] + 1)), 
            col = "red", las=1, lwd=IC_line_width) 
    }
    # add data axis label
    mtext(data_axis_label, side = if (horiz==TRUE) 1 else 2, line=0.75, cex=axis_label_cex)
    # legend
    legend('topleft', internal_control_label, col='red', lty=1, cex=0.4, bty='n')
    # add grid lines
    if(horiz==TRUE){
      abline(v = seq(-50, 50, by = 5), col = 'grey', lty = 3)
    } else {
      abline(h = seq(-50, 50, by = 5), col = 'grey', lty = 3)
    }
  })
}