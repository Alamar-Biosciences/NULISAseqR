#' Plot Plate Layout with Well Type Coloring 
#' 
#' This function generates plots and HTML tables for plate layouts, highlighting various well types (IPC, NC, SC, Bridge).
#' 
#' @param runs runs A named list of run data output from \code{loadNULISAseq()} function 
#' or a list of these outputs for multiple runs.
#' @param plate_names A vector of plate names, default is NULL and will be inferred from `runs`.
#' @param colors A vector of colors for different well types (IPC, NC, etc.), default is NULL to use a predefined color palette.
#' @param col_width The column width for the HTML table.
#' @param font_size The font size for the HTML table.
#' 
#' @return A list of HTML elements with tables and plots for each plate layout.
#' @export
plot_plateLayout <- function(runs,
                             output_dir = NULL,
                             plate_names = NULL,
                             colors = NULL,
                             col_width = "1.75cm",
                             font_size = 8) {
  # Input validation checks
  if (!is.list(runs)) {
    stop("runs must be a list of run objects")
  }
  
  if (length(runs) == 0) {
    stop("runs list cannot be empty")
  }
  
  if (!all(c("samples", "plateID", "IPC", "NC") %in% unique(unlist(lapply(runs, names))))) {
    stop("Each run object must contain samples, plateID, IPC, and NC elements")
  }
  
  # Set default colors if not provided
  if (is.null(colors)) {
    colors <- grDevices::hcl.colors(5, palette = "Set3")[c(4, 3, 1, 5, 2)]
  }
  
  # Create output directory if needed
  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get plate names if not provided
  if (is.null(plate_names)) {
    plate_names <- vapply(runs, function(x) x$plateID, character(1))
  }
  
  # Extract control wells
  IPC_wells <- lapply(runs, function(x) x$IPC)
  NC_wells <- lapply(runs, function(x) x$NC)
  SC_wells <- lapply(runs, function(x) x$SC)
  Bridge_wells <- lapply(runs, function(x) x$Bridge)
  
  # Set SC indicator
  indicatorSC <- all(sapply(runs, function(x) length(x$SC) > 0))
  # Set Bridge indicator
  indicatorBridge <- if(sum(sapply(runs, function(x) length(x$Bridge))) == 0) FALSE else TRUE
  
  # Output list to store HTML elements
  html_output <- list()
  
  # Process each plate
  for (i in seq_along(runs)) {
    # Get plate matrix using internal matrixify function
    val <- matrixify(runs[[i]]$samples, plate = TRUE)
    well_types <- list(IPC_wells[[i]], NC_wells[[i]])
    if (indicatorSC) { well_types <- append(well_types, list(SC_wells[[i]])) }
    if (indicatorBridge) { well_types <- append(well_types, list(Bridge_wells[[i]])) }
    
    # Apply colors based on well types
    for (j in 1:length(well_types)) {
      inds <- which(grepl(paste0(well_types[[j]], collapse = "|"), val, fixed = TRUE))
      if (length(inds) > 0) {
        val[inds] <- kableExtra::cell_spec(val[inds], "html", background = colors[j + 1])
      }
    }
    
    # Replace underscores with spaces for text wrapping
    val <- gsub('_', ' ', val)
    
    # Plot the plate layout as a matrix
    val_matrix <- matrix(0, nrow = nrow(val), ncol = ncol(val))
    for (j in 1:length(well_types)) {
      inds <- which(grepl(paste0(well_types[[j]], collapse = "|"), val))
      if (length(inds) > 0) {
        val_matrix[inds] <- j  # Assign numeric identifiers for well types
      }
    }
    
    # Generate plot
    plot_title <- paste0(plate_names[i], ": Plate Layout")
    plot_file <- file.path(output_dir, paste0("plate_", i, "_layout.png"))
    grDevices::png(plot_file)
    image(val_matrix, main = plot_title, col = colors, axes = FALSE)  # Plot plate layout
    dev.off()
    
    # Create table
    table.attr <- paste0("id=\"plate-summary", i, "\"")
    table_html <- kableExtra::kbl(val, caption = paste0(runs[[i]]$plateID, ": Plate Layout"),
                                  align = "c", escape = FALSE, table.attr = table.attr) %>%
      kableExtra::kable_styling(font_size = font_size, bootstrap_options = c("striped", "hover", "condensed"), 
                                full_width = TRUE) %>%
      kableExtra::column_spec(column = 1:13, width = col_width, border_left = TRUE, border_right = TRUE) %>%
      kableExtra::row_spec(0:8, extra_css = "border-bottom: 1px solid")
    
    # Store the generated table HTML and plot image file path in the output list
    html_output[[i]] <- list(
      table_html = table_html,
      plot_file = plot_file
    )
  }
  
  # Return HTML output (tables and plot image file paths)
  return(html_output)
}