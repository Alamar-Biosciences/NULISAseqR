#' Generate Heatmap for NULISAseq Data with ComplexHeatmap
#'
#' Draws a heatmap for a set of samples and targets based on sample and target metadata. 
#' Uses ComplexHeatmap for generating plots. Colors of covariates are automatically generated 
#' from palettes in RColorBrewer, or can be user-specified. Supports both standard and transposed 
#' orientations, row and column annotations, and various clustering options.
#'
#' @param data A matrix with targets in rows, samples in columns. 
#' Row names should be the target names, and column names are the sample names.
#' It is assumed that data has already been transformed using \code{log2(x + 1)} 
#' for each NULISAseq normalized count value \code{x}, i.e. NPQ.
#' @param sampleInfo A data frame with sample metadata. Rows are samples, 
#' columns are sample metadata variables. Heatmap will only show the samples in 
#' \code{sample_subset} from \code{sampleInfo}.
#' @param targetInfo A data frame with target metadata. Rows are targets,
#' columns are target metadata variables. Required if \code{annotate_target_by} 
#' or \code{row_split_by} are specified; defaults to \code{NULL}.
#' @param sampleName_var Character string specifying the name of the column in \code{sampleInfo} 
#' that matches the column names of \code{data}. This variable will be used for subsetting 
#' samples from \code{sampleInfo}.
#' @param targetName_var Character string specifying the name of the column in \code{targetInfo}
#' that matches the row names of \code{data}; defaults to \code{"Target"}.
#' @param sample_subset Vector of sample names for selected samples to show in the heatmap, 
#' should match the existing column names of \code{data}; defaults to \code{NULL} (all samples).
#' @param target_subset Vector of target names for selected targets to show in the heatmap, 
#' should match the existing row names of \code{data}; defaults to \code{NULL} (all targets).
#' @param annotate_sample_by Character vector of column names from \code{sampleInfo} that will
#' be used for sample annotations (shown as colored bars); defaults to \code{NULL}.
#' @param annotate_target_by Character vector of column names from \code{targetInfo} that will
#' be used for target annotations (shown as colored bars on the left or top depending on orientation); 
#' defaults to \code{NULL}.
#' @param column_split_by Character string specifying the name of the column from \code{sampleInfo} 
#' that will be used for supervised clustering of columns (samples). This creates separate column 
#' slices in the heatmap; defaults to \code{NULL}.
#' @param row_split_by Character string specifying the name of the column from \code{targetInfo} 
#' that will be used for supervised clustering of rows (targets). This creates separate row 
#' slices in the heatmap; defaults to \code{NULL}.
#' @param row_fontsize Numeric value for the text size of the row labels in heatmap; 
#' defaults to 4.
#' @param col_fontsize Numeric value for the text size of the column labels in heatmap; 
#' defaults to 4.
#' @param column_split_num Integer specifying the number of slices that the columns are split into
#' via unsupervised clustering; defaults to \code{NULL}. Ignored if \code{column_split_by} is specified.
#' @param name Character string used as the title of the heatmap legend; 
#' defaults to \code{"Z-Score"}.
#' @param row_title_rot Numeric value for rotation of row titles in degrees. 
#' Only 0, 90, 270 are allowed; defaults to 0.
#' @param clustering_method_columns Character string specifying the method to perform 
#' hierarchical clustering, passed to \code{hclust}; defaults to \code{"ward.D2"}.
#' @param row_split Integer specifying the number of slices that the rows are split into
#' via unsupervised clustering; defaults to \code{NULL}. Ignored if \code{row_split_by} is specified.
#' When \code{NULL} and \code{targetInfo} is not provided, no row splitting is performed.
#' @param cluster_rows Logical indicating whether to cluster rows (targets). If \code{NULL} (default),
#' clustering is enabled when \code{targetInfo} is provided, and disabled when \code{targetInfo} is \code{NULL}.
#' @param cluster_column_slices Logical indicating whether to perform clustering on column slices 
#' if columns are split; defaults to \code{FALSE}.
#' @param cluster_row_slices Logical indicating whether to perform clustering on row slices
#' if rows are split; defaults to \code{TRUE}.
#' @param transpose Logical indicating whether to transpose the heatmap (samples in rows, 
#' targets in columns); defaults to \code{FALSE}.
#' @param numeric_color_palette Character vector of colors to use for numeric annotations.
#' Each color will be used with a white-to-color gradient for continuous variables;
#' defaults to \code{c("blue", "red", "green", "purple", "orange", "brown")}.
#' @param sample_colors Named list of custom colors for categorical sample annotations.
#' List names should match column names in \code{annotate_sample_by}. Each element should be
#' a named vector where names are category levels and values are color codes; defaults to \code{NULL}.
#' @param output_dir Character string specifying the directory path to save the plot. 
#' If \code{NULL}, the plot is not saved; defaults to \code{NULL}.
#' @param plot_name Character string specifying the filename for the saved plot, 
#' including file extension (.pdf, .png, .jpg, or .svg). Required if \code{output_dir} is specified; 
#' defaults to \code{NULL}.
#' @param plot_title Character string for the title of the heatmap; defaults to \code{NULL}.
#' @param plot_width Numeric value for the width of the saved plot in inches; defaults to 14.
#' @param plot_height Numeric value for the height of the saved plot in inches; defaults to 7.
#' @param ... Additional arguments passed to \code{ComplexHeatmap::Heatmap} function.
#'
#' @return A list containing:
#' \describe{
#'   \item{targets_used}{Character vector of target names used in the heatmap after filtering.}
#'   \item{heatmap}{The ComplexHeatmap object.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Filters data to specified samples and targets
#'   \item Removes targets with all zero values
#'   \item Scales data by row (Z-score transformation)
#'   \item Removes rows with NA, NaN, or Inf values after scaling
#'   \item Optionally transposes the matrix
#'   \item Generates or uses custom colors for annotations
#'   \item Creates the heatmap with specified annotations and clustering
#'   \item Optionally saves to file
#' }
#'
#' @section Custom Colors:
#' To specify custom colors for categorical sample annotations, use the \code{sample_colors} parameter:
#' \preformatted{
#' my_colors <- list(
#'   Group = c("Control" = "#FF0000", "Treatment" = "#0000FF"),
#'   Batch = c("Batch1" = "#FFA500", "Batch2" = "#800080")
#' )
#' }
#'
#' @examples
#' \dontrun{
#' # Basic heatmap with sample annotations (no row clustering by default)
#' result <- generate_heatmap(
#'   data = Data_NPQ,
#'   sampleInfo = sample_metadata,
#'   sampleName_var = "SampleName",
#'   annotate_sample_by = c("Group", "Batch")
#' )
#' 
#' # Heatmap with row clustering enabled
#' result <- generate_heatmap(
#'   data = Data_NPQ,
#'   sampleInfo = sample_metadata,
#'   sampleName_var = "SampleName",
#'   annotate_sample_by = c("Group", "Batch"),
#'   cluster_rows = TRUE,
#'   row_split = 3
#' )
#'
#' # Heatmap with target annotations and custom colors
#' custom_colors <- list(
#'   Group = c("Control" = "blue", "Treatment" = "red")
#' )
#' 
#' result <- generate_heatmap(
#'   data = Data_NPQ,
#'   sampleInfo = sample_metadata,
#'   targetInfo = target_metadata,
#'   sampleName_var = "SampleName",
#'   targetName_var = "Target",
#'   annotate_sample_by = c("Group", "Batch"),
#'   annotate_target_by = c("Pathway", "Domain"),
#'   sample_colors = custom_colors,
#'   row_split_by = "Pathway",
#'   column_split_by = "Group"
#' )
#'
#' # Save heatmap to file (supports PDF, PNG, JPG, SVG)
#' result <- generate_heatmap(
#'   data = expr_matrix,
#'   sampleInfo = sample_metadata,
#'   sampleName_var = "Sample_ID",
#'   annotate_sample_by = c("Group"),
#'   output_dir = "output/figures",
#'   plot_name = "expression_heatmap.svg",
#'   plot_title = "Gene Expression Heatmap",
#'   plot_width = 12,
#'   plot_height = 8
#' )
#' }
#' @import showtext
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation draw
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom grid gpar
#' @importFrom dplyr filter distinct
#' @importFrom tibble column_to_rownames
#' @importFrom showtext showtext_auto
#' @importFrom grDevices colorRampPalette pdf png jpeg dev.off
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom svglite svglite
#'
#' @export
generate_heatmap <- function(data, 
                             sampleInfo,
                             targetInfo = NULL,
                             sampleName_var,
                             targetName_var = "Target",
                             sample_subset = NULL, 
                             target_subset = NULL, 
                             annotate_sample_by = NULL,
                             annotate_target_by = NULL,
                             column_split_by = NULL,
                             row_split_by = NULL,
                             row_fontsize = 4,
                             col_fontsize = 4,
                             column_split_num = NULL,
                             name = "Z-Score",
                             row_title_rot = 0,
                             clustering_method_columns = "ward.D2",
                             row_split = NULL,
                             cluster_rows = TRUE,
                             cluster_column_slices = FALSE,
                             cluster_row_slices = TRUE,
                             transpose = FALSE,
                             numeric_color_palette = c("blue", "red", "green", "purple", "orange", "brown"),
                             sample_colors = NULL,
                             output_dir = NULL,
                             plot_name = NULL,
                             plot_title = NULL,
                             plot_width = 14,
                             plot_height = 7,
                             ...){
  
  showtext::showtext_auto() ## able to output beta symbol and others in pdf
  
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
    warning("Warning: Some targets contain only zeros and were removed from the heatmap:\n", 
            paste(target_remove, collapse = ", "))
  }
  
  # Scaling data (always scale by targets/rows)
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
  
  # Transpose if requested
  if (transpose) {
    scaled_data <- t(scaled_data)
  }
  
  # Set row_split based on targetInfo if not specified
  if (is.null(row_split) && !is.null(targetInfo)) {
    row_split <- 2  # Default to 2 splits when targetInfo is provided
  }
  
  # Remove annotation columns that contain only NA values
  valid_annotate_sample_by <- annotate_sample_by[sapply(sampleInfo[annotate_sample_by], function(x) any(!is.na(x)))]
  
  if (length(valid_annotate_sample_by) == 0 && !is.null(annotate_sample_by)) {
    stop("Error: No valid annotation columns found after removing NA-only columns.")
  }
  
  # Prepare sample annotations
  sample_annotation <- sampleInfo %>% 
    dplyr::filter(get(sampleName_var) %in% (if (transpose) rownames(scaled_data) else colnames(scaled_data))) %>%  
    dplyr::distinct(get(sampleName_var), .keep_all = TRUE) %>% 
    tibble::column_to_rownames(sampleName_var)
  
  # Ensure annotation variables exist
  if (!is.null(annotate_sample_by)) {
    missing_annotations <- setdiff(annotate_sample_by, colnames(sample_annotation))
    if (length(missing_annotations) > 0) {
      stop("Error: The following annotation variables are missing from sampleInfo: ", 
           paste(missing_annotations, collapse = ", "))
    }
    
    # Subset annotation columns based on valid annotations
    sample_annotation <- sample_annotation[, valid_annotate_sample_by, drop = FALSE]
  }
  
  sample_annotation_df <- sample_annotation[(if (transpose) rownames(scaled_data) else colnames(scaled_data)), annotate_sample_by, drop = FALSE]
  
  # Prepare column split dataframe (keep as dataframe, not extract as vector)
  sample_annotation_col <- NULL
  column_split_factor <- NULL
  if (!is.null(column_split_by)) {
    sample_annotation_col <- sample_annotation[(if (transpose) rownames(scaled_data) else colnames(scaled_data)), column_split_by, drop = FALSE]
    # Convert to factor for ComplexHeatmap column_split parameter to display proper labels
    column_split_factor <- factor(sample_annotation_col[[column_split_by]])
  }
  
  # Identify numeric and categorical columns
  is_numeric_col <- sapply(sample_annotation_df, function(x) is.numeric(x) && length(unique(x[!is.na(x)])) > 5)
  
  # Generate annotation colors for sample annotations
  different_sets <- rep_len(c("Set1", "Set2", "Set3", "Dark2"), length.out = ncol(sample_annotation_df))
  
  # Track which numeric color to use
  numeric_color_index <- 1
  
  heatmap_annotation_colors <- list()
  
  for (i in seq_along(colnames(sample_annotation_df))) {
    col_name <- colnames(sample_annotation_df)[i]
    
    if (is_numeric_col[col_name]) {
      # For numeric columns, use a continuous color scale from white to a color
      current_color <- numeric_color_palette[numeric_color_index]
      heatmap_annotation_colors[[col_name]] <- circlize::colorRamp2(
        breaks = seq(min(sample_annotation_df[[col_name]], na.rm = TRUE),
                     max(sample_annotation_df[[col_name]], na.rm = TRUE),
                     length.out = 100),
        colors = grDevices::colorRampPalette(c("white", current_color))(100)
      )
      # Cycle through available colors
      numeric_color_index <- (numeric_color_index %% length(numeric_color_palette)) + 1
    } else {
      # Check if user provided custom colors for this annotation
      if (!is.null(sample_colors) && col_name %in% names(sample_colors)) {
        # Use user-provided colors
        user_colors <- sample_colors[[col_name]]
        
        # Get unique values in the data
        unique_vals <- sort(unique(sample_annotation_df[[col_name]]))
        
        # Check if all values have colors assigned
        missing_vals <- setdiff(unique_vals, names(user_colors))
        if (length(missing_vals) > 0) {
          stop("Error: Missing color assignments for '", col_name, "' values: ", 
               paste(missing_vals, collapse = ", "))
        }
        
        # Use only the colors for values present in the data
        heatmap_annotation_colors[[col_name]] <- user_colors[unique_vals]
      } else {
        # For categorical columns, use discrete colors
        suppressWarnings(
          heatmap_annotation_colors[[col_name]] <- generate_cov_colors(
            cov = col_name,
            set = different_sets[i],
            data = sample_annotation_df
          )
        )
      }
    }
  }
  
  # Initialize row annotation and split variables
  row_ha <- NULL
  row_split_vector <- NULL
  col_ha <- NULL
  col_split_vector <- NULL
  
  # Handle row/column splitting and annotations based on transpose
  if (!is.null(targetInfo) && (!is.null(row_split_by) || !is.null(annotate_target_by))) {
    # Check if targetName_var exists in targetInfo
    if (!targetName_var %in% colnames(targetInfo)) {
      stop("Error: '", targetName_var, "' column not found in targetInfo.")
    }
    
    # Determine which dimension targets are in
    target_dimension <- if (transpose) colnames(scaled_data) else rownames(scaled_data)
    
    # Prepare target data
    target_data <- targetInfo %>% 
      dplyr::filter(get(targetName_var) %in% target_dimension) %>%  
      dplyr::distinct(get(targetName_var), .keep_all = TRUE) %>% 
      tibble::column_to_rownames(targetName_var)
    
    # Match to scaled_data
    target_data <- target_data[match(target_dimension, rownames(target_data)), , drop = FALSE]
    
    # Handle row splitting
    if (!is.null(row_split_by)) {
      if (!row_split_by %in% colnames(targetInfo)) {
        stop("Error: '", row_split_by, "' not found in targetInfo columns.")
      }
      
      split_vector <- target_data[, row_split_by]
      split_vector <- ifelse(is.na(split_vector), "NA", as.character(split_vector))
      
      if (transpose) {
        col_split_vector <- split_vector
      } else {
        row_split_vector <- split_vector
      }
    }
    
    # Handle target annotations
    if (!is.null(annotate_target_by)) {
      # Check if annotation columns exist in targetInfo
      missing_target_annotations <- setdiff(annotate_target_by, colnames(targetInfo))
      if (length(missing_target_annotations) > 0) {
        stop("Error: The following annotation variables are missing from targetInfo: ", 
             paste(missing_target_annotations, collapse = ", "))
      }
      
      # Remove annotation columns that contain only NA values
      valid_annotate_target_by <- annotate_target_by[sapply(target_data[annotate_target_by], function(x) any(!is.na(x)))]
      
      if (length(valid_annotate_target_by) > 0) {
        # Subset to valid annotation columns
        target_annotation_df <- target_data[, valid_annotate_target_by, drop = FALSE]
        
        # Replace NA values with "NA" string for visualization
        target_annotation_df[] <- lapply(target_annotation_df, function(x) {
          ifelse(is.na(x), "NA", as.character(x))
        })
        
        # Generate colors for target annotations
        different_sets_target <- rep_len(c("Set3", "Dark2", "Set1", "Set2"), length.out = ncol(target_annotation_df))
        
        suppressWarnings(
          target_annotation_colors <- mapply(generate_cov_colors, cov = colnames(target_annotation_df), 
                                             set = different_sets_target, 
                                             MoreArgs = list(data = target_annotation_df), SIMPLIFY = FALSE)
        )
        
        # Create annotation based on orientation
        if (transpose) {
          col_ha <- ComplexHeatmap::HeatmapAnnotation(
            df = target_annotation_df,
            col = target_annotation_colors,
            show_annotation_name = TRUE
          )
        } else {
          row_ha <- ComplexHeatmap::rowAnnotation(
            df = target_annotation_df,
            col = target_annotation_colors,
            show_annotation_name = TRUE
          )
        }
      } else {
        warning("Warning: All target annotation columns contain only NA values. No target annotation will be added.")
      }
    }
  }
  
  # Generate heatmap
  if (transpose) {
    # Transposed view: samples in rows, targets in columns
    suppressMessages(
      h1 <- ComplexHeatmap::Heatmap(
        name = name,
        matrix = scaled_data,
        left_annotation = if (!is.null(annotate_sample_by)) 
          ComplexHeatmap::rowAnnotation(df = sample_annotation_df, col = heatmap_annotation_colors),
        top_annotation = col_ha,
        row_names_gp = grid::gpar(fontsize = row_fontsize),
        column_names_gp = grid::gpar(fontsize = col_fontsize),
        cluster_columns = cluster_rows,  # Swapped for transpose
        cluster_rows = TRUE,  # Always cluster samples (columns in original)
        cluster_column_slices = cluster_row_slices,  # Swapped for transpose
        cluster_row_slices = cluster_column_slices,  # Swapped for transpose
        column_split = if (!is.null(col_split_vector)) col_split_vector else if (!is.null(row_split_by)) NULL else column_split_num,
        row_split = if (!is.null(column_split_by)) sample_annotation_col else NULL,
        clustering_method_columns = clustering_method_columns,  # Apply to targets now
        clustering_method_rows = clustering_method_columns,  # Apply to samples
        row_title_rot = row_title_rot,
        row_title_gp = grid::gpar(fontsize = 12),  # Set font size for row split labels
        column_title_gp = grid::gpar(fontsize = 12),  # Set font size for column split labels
        ...
      )
    )
  } else {
    # Normal view: targets in rows, samples in columns
    suppressMessages(
      h1 <- ComplexHeatmap::Heatmap(
        name = name,
        matrix = scaled_data,
        left_annotation = row_ha,
        top_annotation = if (!is.null(annotate_sample_by)) 
          ComplexHeatmap::HeatmapAnnotation(df = sample_annotation_df, col = heatmap_annotation_colors),
        row_names_gp = grid::gpar(fontsize = row_fontsize),
        column_names_gp = grid::gpar(fontsize = col_fontsize),
        cluster_rows = cluster_rows,
        cluster_columns = TRUE,  # Always cluster samples
        cluster_column_slices = cluster_column_slices,
        cluster_row_slices = cluster_row_slices,
        column_split = if (!is.null(column_split_by)) column_split_factor else column_split_num,
        row_split = if (!is.null(row_split_vector)) row_split_vector else row_split,
        clustering_method_columns = clustering_method_columns,
        row_title_rot = row_title_rot,
        column_title_side = "top",
        column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
        show_column_dend = TRUE,
        show_column_names = TRUE,
        ...
      )
    )
  }
  
  # Save plot if output directory is provided
  if (!is.null(output_dir)) {
    # Check if plot_name is provided, if not use default
    if (is.null(plot_name)) {
      plot_name <- paste0("heatmap_", format(Sys.time(), "%Y%m%d"), ".pdf")
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
      ComplexHeatmap::draw(h1, background = "transparent", padding = grid::unit(c(2, 10, 2, 2), "mm"))
      grDevices::dev.off()
    } else if (file_ext %in% c("png", "PNG")) {
      grDevices::png(output_path, width = plot_width * 100, height = plot_height * 100, res = 100)
      ComplexHeatmap::draw(h1, background = "transparent", padding = grid::unit(c(2, 10, 2, 2), "mm"))
      grDevices::dev.off()
    } else if (file_ext %in% c("jpg", "jpeg", "JPG", "JPEG")) {
      grDevices::jpeg(output_path, width = plot_width * 100, height = plot_height * 100, res = 100)
      ComplexHeatmap::draw(h1, background = "transparent", padding = grid::unit(c(2, 10, 2, 2), "mm"))
      grDevices::dev.off()
    } else if (file_ext %in% c("svg", "SVG")) {
      svglite::svglite(output_path, width = plot_width, height = plot_height)
      ComplexHeatmap::draw(h1, background = "transparent", padding = grid::unit(c(2, 10, 2, 2), "mm"))
      grDevices::dev.off()
    } else {
      warning("Unsupported file extension: ", file_ext, ". Saving as PDF instead.")
      output_path <- paste0(tools::file_path_sans_ext(output_path), ".pdf")
      grDevices::pdf(output_path, width = plot_width, height = plot_height)
      ComplexHeatmap::draw(h1, background = "transparent", padding = grid::unit(c(2, 10, 2, 2), "mm"))
      grDevices::dev.off()
    }
    
    message("Heatmap saved to: ", output_path)
  } else {
    # Just draw to current device
    ComplexHeatmap::draw(h1, background = "transparent", padding = grid::unit(c(2, 10, 2, 2), "mm"))
  }
  
  return(list(
    targets_used = target_used,
    heatmap = h1,
    output_path = if (!is.null(output_dir) && !is.null(plot_name)) output_path else NULL
  ))
}
#' Generate Colors for Categorical Annotations
#'
#' Helper function to generate color palettes for categorical variables in heatmap annotations.
#' Uses RColorBrewer palettes and falls back to additional palettes if needed.
#'
#' @param cov Character string specifying the column name of the annotation variable.
#' @param data Data frame containing the annotation data.
#' @param set Character string specifying the primary RColorBrewer palette to use.
#' @param fallback_palettes Vector of character strings specifying fallback RColorBrewer 
#' palettes if the primary set doesn't have enough colors; defaults to a predefined list.
#'
#' @return A named vector of colors where names correspond to unique levels of the categorical variable.
#'
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#'
#' @keywords internal
generate_cov_colors <- function(cov, data, set, 
                                fallback_palettes = c("Set1", "Set2", "Set3", "Dark2", "Pastel1", 
                                                      "Pastel2", "Purples", "Oranges",
                                                      "Blues", "Paired", "Set1")) {
  
  # Extract unique groups
  unique_groups <- sort(unique(data[,cov]))
  unique_groups_length <- length(unique_groups)
  
  # Define empty color list
  all_colors <- c()
  
  # Loop through provided color sets (including fallback)
  for (color_set in c(set, fallback_palettes)) {
    # Get colors from each set using brewer.pal
    colors <- RColorBrewer::brewer.pal(min(unique_groups_length, RColorBrewer::brewer.pal.info[color_set, "maxcolors"]), 
                                       name = color_set)
    all_colors <- c(all_colors, colors)
  }
  
  # Limit colors to unique groups length and assign names
  color_subset <- all_colors[1:unique_groups_length]
  names(color_subset) <- unique_groups
  
  # Return the modified color list
  color_subset
}
