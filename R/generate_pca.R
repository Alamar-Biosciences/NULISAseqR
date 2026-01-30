#' Generate PCA Biplot for NULISAseq Data with PCAtools
#'
#' Performs Principal Component Analysis (PCA) and generates a biplot for visualizing
#' sample relationships based on gene expression data. Uses PCAtools for analysis
#' and ggplot2 for visualization. Colors are automatically generated from RColorBrewer
#' palettes.
#'
#' @param data A matrix with targets in rows, samples in columns. 
#' Row names should be the target names, and column names are the sample names.
#' It is assumed that data has already been transformed using \code{log2(x + 1)} 
#' for each NULISAseq normalized count value \code{x}, i.e. NPQ.
#' @param sampleInfo A data frame with sample metadata. Rows are samples, 
#' columns are sample metadata variables.
#' @param sampleName_var Character string specifying the name of the column in \code{sampleInfo} 
#' that matches the column names of \code{data}.
#' @param annotate_sample_by Character string specifying the column name from \code{sampleInfo}
#' to use for coloring points. Only one variable is allowed; defaults to \code{NULL}.
#' @param label_points Logical indicating whether to add sample labels to the plot;
#' defaults to \code{FALSE}.
#' @param sample_subset Vector of sample names for selected samples to include in PCA, 
#' should match the existing column names of \code{data}; defaults to \code{NULL} (all samples).
#' @param target_subset Vector of target names for selected targets to include in PCA, 
#' should match the existing row names of \code{data}; defaults to \code{NULL} (all targets).
#' @param shape_by Character string specifying the column name from \code{sampleInfo}
#' to use for point shapes; defaults to \code{NULL}.
#' @param ellipse Logical indicating whether to draw ellipses around groups;
#' defaults to \code{TRUE}.
#' @param ellipseType Character string specifying the type of ellipse. Options include
#' \code{"t"} for t-distribution and \code{"norm"} for normal distribution;
#' defaults to \code{"t"}.
#' @param ellipseAlpha Numeric value between 0 and 1 for ellipse transparency;
#' defaults to \code{0.15}.
#' @param ellipseFill Logical indicating whether to fill the ellipses;
#' defaults to \code{TRUE}.
#' @param ellipseLineSize Numeric value for the ellipse border line width;
#' defaults to \code{0} (no border).
#' @param sample_colors Named vector of custom colors for sample groups.
#' Names should match the levels in \code{annotate_sample_by}; defaults to \code{NULL}.
#' @param components Integer vector of length 2 specifying which principal components to plot.
#' For example, \code{c(1, 2)} plots PC1 vs PC2, \code{c(2, 3)} plots PC2 vs PC3;
#' defaults to \code{c(1, 2)}.
#' @param output_dir Character string specifying the directory path to save the plot. 
#' If \code{NULL}, the plot is not saved; defaults to \code{NULL}. If provided without
#' \code{plot_name}, a default filename with timestamp will be generated.
#' @param plot_name Character string specifying the filename for the saved plot, 
#' including file extension (.pdf, .png, .jpg, or .svg). If \code{NULL} and \code{output_dir}
#' is provided, a default filename with timestamp will be used; defaults to \code{NULL}.
#' @param plot_title Character string for the title of the PCA plot; defaults to \code{NULL}.
#' @param plot_width Numeric value for the width of the saved plot in inches; defaults to 10.
#' @param plot_height Numeric value for the height of the saved plot in inches; defaults to 8.
#' @param ... Additional arguments passed to \code{PCAtools::biplot} function.
#'
#' @return A list containing:
#' \describe{
#'   \item{targets_used}{Character vector of target names used in the PCA after filtering.}
#'   \item{pca_results}{The PCAtools PCA object containing all PCA results.}
#'   \item{rotated}{Data frame containing the PC scores (PC1, PC2, etc.) for each sample.}
#'   \item{plot}{The ggplot2 object of the PCA biplot.}
#'   \item{output_path}{Character string of the full path to the saved file, or \code{NULL} if not saved.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Filters data to specified samples and targets
#'   \item Removes targets with all zero values
#'   \item Scales data by row (Z-score transformation)
#'   \item Removes rows with NA, NaN, or Inf values after scaling
#'   \item Performs PCA using PCAtools
#'   \item Generates biplot with specified aesthetics
#'   \item Optionally saves to file
#' }
#'
#' @section Custom Colors:
#' To specify custom colors for sample groups, use the \code{sample_colors} parameter:
#' \preformatted{
#' my_colors <- c("Control" = "#FF0000", "Treatment" = "#0000FF")
#' }
#'
#' @examples
#' \dontrun{
#' # Basic PCA plot
#' result <- generate_pca(
#'   data = Data_NPQ,
#'   sampleInfo = sample_metadata,
#'   sampleName_var = "SampleName",
#'   annotate_sample_by = "Group"
#' )
#'
#' # PCA with sample labels and custom shapes
#' result <- generate_pca(
#'   data = Data_NPQ,
#'   sampleInfo = sample_metadata,
#'   sampleName_var = "SampleName",
#'   annotate_sample_by = "Group",
#'   shape_by = "Batch",
#'   label_points = TRUE
#' )
#'
#' # PCA with custom colors
#' custom_colors <- c("Control" = "blue", "Treatment" = "red")
#' result <- generate_pca(
#'   data = Data_NPQ,
#'   sampleInfo = sample_metadata,
#'   sampleName_var = "SampleName",
#'   annotate_sample_by = "Group",
#'   sample_colors = custom_colors
#' )
#'
#' # Save PCA plot to file
#' result <- generate_pca(
#'   data = Data_NPQ,
#'   sampleInfo = sample_metadata,
#'   sampleName_var = "SampleName",
#'   annotate_sample_by = "Group",
#'   output_dir = "output/figures",
#'   plot_name = "pca_analysis.pdf",
#'   plot_title = "PCA Analysis of Expression Data",
#'   plot_width = 12,
#'   plot_height = 10
#' )
#' }
#'
#' @importFrom PCAtools pca biplot
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr filter distinct select everything
#' @importFrom tibble column_to_rownames
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 aes ggtitle theme element_text element_rect scale_fill_manual scale_color_manual scale_shape_manual
#' @importFrom grDevices pdf png jpeg dev.off
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom svglite svglite
#' @importFrom rlang sym
#'
#' @export
generate_pca <- function(data,
                         sampleInfo,
                         sampleName_var,
                         annotate_sample_by = NULL,
                         label_points = FALSE,
                         sample_subset = NULL,
                         target_subset = NULL,
                         shape_by = NULL,
                         ellipse = TRUE,
                         ellipseType = "t",
                         ellipseAlpha = 0.15,
                         ellipseFill = TRUE,
                         ellipseLineSize = 0,
                         sample_colors = NULL,
                         components = c(1, 2),
                         output_dir = NULL,
                         plot_name = NULL,
                         plot_title = NULL,
                         plot_width = 10,
                         plot_height = 8,
                         ... ) {
  
  # Ensure valid subsetting
  sample_subset <- intersect(sample_subset, colnames(data))
  target_subset <- intersect(target_subset, rownames(data))
  
  if (!is.null(sample_subset) && length(sample_subset) == 0) {
    stop("Error: No valid samples in 'sample_subset' found in data.")
  }
  
  if (!is.null(target_subset) && length(target_subset) == 0) {
    stop("Error: No valid targets in 'target_subset' found in data.")
  }
  
  # Subset data correctly
  filtered_data <- data
  if (!is.null(sample_subset)) filtered_data <- filtered_data[, sample_subset, drop = FALSE]
  if (!is.null(target_subset)) filtered_data <- filtered_data[target_subset, , drop = FALSE]
  
  # Check if filtering removed all data
  if (nrow(filtered_data) == 0 || ncol(filtered_data) == 0) {
    stop("Error: After filtering, data is empty. Check 'sample_subset' and 'target_subset'.")
  }
  
  # Remove rows with all zero values
  target_remove <- rownames(filtered_data[rowSums(filtered_data != 0, na.rm = TRUE) == 0, ])
  if (length(target_remove) != 0) {
    filtered_data <- filtered_data[!rownames(filtered_data) %in% target_remove, , drop = FALSE]
    warning("Warning: Some targets contain only zeros and were removed from the PCA:\n", 
            paste(target_remove, collapse = ", "))
  }
  
  # Scaling data
  scaled_data <- t(scale(t(filtered_data), center = TRUE, scale = TRUE))
  
  # Identify rows with any NA, NaN, or Inf
  bad_rows <- rownames(scaled_data)[apply(scaled_data, 1, function(x) any(is.na(x) | is.nan(x) | is.infinite(x)))]
  
  # Check for scaling issues (NA, NaN or Inf values)
  if (any(is.na(scaled_data)) || any(is.nan(scaled_data)) || any(is.infinite(scaled_data))) {
    warning("Warning: Scaling resulted in NA or NaN or Inf values. Remove rows/columns with missing values.\n",
            "Check your data for extreme values or missing entries.\n",
            "Targets excluded: ", paste(bad_rows, collapse = ", "))
    scaled_data <- scaled_data[!rownames(scaled_data) %in% bad_rows, , drop = FALSE]
  }
  
  target_used <- rownames(scaled_data)
  
  # Check if annotation variable is provided
  if (is.null(annotate_sample_by)) {
    stop("Error: 'annotate_sample_by' must be provided for PCA visualization.")
  }
  
  # Only 1 annotation variable allowed
  if (length(annotate_sample_by) > 1) {
    stop("Error: Only 1 annotation variable is allowed for coloring points in PCA.")
  }
  
  # Prepare sample annotations
  sample_annotation <- sampleInfo %>% 
    dplyr::filter(get(sampleName_var) %in% colnames(scaled_data)) %>%  
    dplyr::distinct(get(sampleName_var), .keep_all = TRUE) 
  
  # Ensure annotation variable exists
  missing_annotations <- setdiff(annotate_sample_by, colnames(sample_annotation))
  if (length(missing_annotations) > 0) {
    stop("Error: The annotation variable '", annotate_sample_by, "' does not exist in sampleInfo.")
  }
  
  sample_annotation <- sample_annotation %>%
    # Reorder columns before setting rownames
    dplyr::select(
      !!rlang::sym(annotate_sample_by),  # Put target column first
      dplyr::everything()
    ) %>%
    tibble::column_to_rownames(sampleName_var)
  
  # Generate or use provided colors
  if (!is.null(sample_colors)) {
    # Use user-provided colors
    unique_vals <- sort(unique(sample_annotation[[annotate_sample_by]]))
    
    # Check if all values have colors assigned
    missing_vals <- setdiff(unique_vals, names(sample_colors))
    if (length(missing_vals) > 0) {
      stop("Error: Missing color assignments for '", annotate_sample_by, "' values: ", 
           paste(missing_vals, collapse = ", "))
    }
    
    pca_annotation_colors <- sample_colors[unique_vals]
  } else {
    # Generate colors automatically using the helper function
    different_sets <- "Set1"
    
    suppressWarnings(
      pca_annotation_colors <- generate_cov_colors(
        cov = annotate_sample_by, 
        data = sample_annotation,
        set = different_sets
      )
    )
  }
  
  # Perform PCA
  pca_res <- PCAtools::pca(scaled_data, 
                           metadata = sample_annotation[colnames(scaled_data), ],
                           center = TRUE,
                           scale = TRUE)
  
  # Generate biplot with specified components
  pc_x <- paste0("PC", components[1])
  pc_y <- paste0("PC", components[2])

  plot <- PCAtools::biplot(pca_res,
                           x = pc_x,
                           y = pc_y,
                           colby = annotate_sample_by,
                           lab = NULL,
                           shape = shape_by,
                           hline = 0,
                           vline = 0,
                           pointSize = 2,
                           legendPosition = 'right',
                           ellipse = ellipse,
                           ellipseType = ellipseType,
                           ellipseAlpha = ellipseAlpha,
                           ellipseFill = ellipseFill,
                           ellipseLineSize = ellipseLineSize,
                           ...) 
  
  # Add title if specified
  if (!is.null(plot_title)) {
    plot <- plot + ggplot2::ggtitle(plot_title)
  }
  
  # Add labels if specified
  if (label_points) {
    pc_coords <- as.data.frame(pca_res$rotated[, c(pc_x, pc_y)])
    pc_coords$Label <- rownames(pca_res$rotated)

    plot <- plot +
      ggrepel::geom_text_repel(data = pc_coords,
                               ggplot2::aes(x = .data[[pc_x]], y = .data[[pc_y]], label = Label),
                               size = 3,
                               box.padding = 0.35,
                               point.padding = 0.3,
                               max.overlaps = 40)
  }
  
  # Apply theme and colors
  plot <- plot + 
    ggplot2::theme(
      text = ggplot2::element_text(size = 25),
      panel.background = ggplot2::element_rect(fill = 'white'),
      plot.background = ggplot2::element_rect(fill = 'transparent', color = NA),
      legend.background = ggplot2::element_rect(fill = 'transparent'),
      legend.key = ggplot2::element_rect(fill = "transparent", color = NA)
    ) + 
    ggplot2::scale_fill_manual(values = pca_annotation_colors) +   
    ggplot2::scale_color_manual(values = pca_annotation_colors) 
  
  # Add shape scale if specified
  if (!is.null(shape_by)) {
    # set number of shapes
    n_shape <- length(unique(sample_annotation[colnames(scaled_data), shape_by])) 
    plot <- plot + ggplot2::scale_shape_manual(values = c(1:n_shape))
  }
  
  # Save plot if output directory is provided
  output_path <- NULL
  if (!is.null(output_dir)) {
    # Check if plot_name is provided, if not use default
    if (is.null(plot_name)) {
      plot_name <- paste0("pca_", format(Sys.time(), "%Y%m%d"), ".pdf")
      message("No plot_name provided. Using default filename: ", plot_name)
    }
    
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Full path for the output file
    output_path <- file.path(output_dir, plot_name)
    
    # Determine file type from extension
    file_ext <- tools::file_ext(plot_name)
    
    # Save based on file extension
    if (file_ext == "pdf") {
      grDevices::pdf(output_path, width = plot_width, height = plot_height)
      print(plot)
      grDevices::dev.off()
    } else if (file_ext %in% c("png", "PNG")) {
      grDevices::png(output_path, width = plot_width * 100, height = plot_height * 100, res = 100)
      print(plot)
      grDevices::dev.off()
    } else if (file_ext %in% c("jpg", "jpeg", "JPG", "JPEG")) {
      grDevices::jpeg(output_path, width = plot_width * 100, height = plot_height * 100, res = 100)
      print(plot)
      grDevices::dev.off()
    } else if (file_ext %in% c("svg", "SVG")) {
      svglite::svglite(output_path, width = plot_width, height = plot_height)
      print(plot)
      grDevices::dev.off()
    } else {
      warning("Unsupported file extension: ", file_ext, ". Saving as PDF instead.")
      output_path <- paste0(tools::file_path_sans_ext(output_path), ".pdf")
      grDevices::pdf(output_path, width = plot_width, height = plot_height)
      print(plot)
      grDevices::dev.off()
    }
    
    message("PCA plot saved to: ", output_path)
  } else {
    # Just print the plot
    print(plot)
  }
  
  # Return results
  return(list(
    targets_used = target_used,
    pca_results = pca_res,
    rotated = pca_res$rotated,
    plot = plot,
    output_path = output_path
  ))
}