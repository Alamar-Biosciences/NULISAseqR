#' Volcano plot for NULISAseq differential expression test
#'
#' Draws a volcano plot for a set of coefficients and p-values. Uses ggplot
#' for plotting and ggrepel package for labels. By default, plot labeling
#' assumes FDR-corrected p-values. Can plot both unadjusted and FDR-adjusted
#' p-values with different colors.
#'
#' @param coefs Vector of coefficients.
#' @param p_vals Vector of p-values OR a named list with 'unadj' and 'fdr' p-values.
#' If a list, the y-axis will plot unadjusted p-values with light colors for unadjusted
#' significant targets and darker colors for FDR significant targets.
#' @param target_labels Vector of target names.
#' @param title Plot title.
#' @param xlabel X axis label.
#' @param ylabel Y axis label; defaults to \code{expression('-log'[10]*'(FDR-adjusted p-value)')}.
#' Automatically overridden to show unadjusted p-value when in dual mode.
#' @param xlimits X axis limits
#' @param ylimits Y axis limits
#' @param sig_threshold P-value significance threshold OR a named list with 'unadj' and 
#' 'fdr' thresholds. Must match the structure of p_vals.
#' @param sig_label Label to mark significance threshold on plot; defaults to \code{"FDR = 5\%"}.
#' Automatically overridden to \code{"Unadj = X\%"} when in dual mode (where X is the unadj threshold percentage).
#' Always a single string.
#' @param label_all_targets Logical TRUE or FALSE. 
#' Should all targets be labelled? Default is FALSE 
#' and will only label significant targets.
#' @param upper_log2FC_threshold Coefficients/fold change upper threshold cutoff for 
#' labeling insignificant targets. Insignificant targets above the upper threshold will be labelled.
#' Default is \code{NULL}
#' @param lower_log2FC_threshold Coefficients/fold change lower threshold cutoff for 
#' labeling insignificant targets. Insignificant targets below the lower threshold will be labelled.
#' Default is \code{NULL}.
#' @param target_labels_off Logical TRUE or FALSE. Should no targets be labelled? 
#' Default is \code{FALSE}.
#' @param target_label_colors  A vector of the same length as number of targets
#' and same order as targets with colors for the target labels. Default is to 
#' color upregulated target labels red and downregulated target labels blue.
#' @param target_point_colors A vector of the same length as number of targets
#' and same order as targets with colors for the points. Default is to 
#' make significant upregulated target points red, significant downregulated 
#' target points blue, and nonsignificant target points grey.
#' @param target_label_size Font size of target labels.
#' @param target_label_segment_color Color of line segments for target labels.
#' Can be a single value or a vector same length as number of targets.
#' @param max.overlaps Passed to \code{ggrepel::geom_text_repel}. Integer which determines how many targets 
#' will have labels. "Inf" (default) labels all targets. 
#' @param force Passed to \code{ggrepel::geom_text_repel}. Force of repulsion between overlapping text labels. Defaults to 1.
#' @param force_pull Passed to \code{ggrepel::geom_text_repel}. Force of attraction between a text label and its corresponding data point. Defaults to 1.
#' @param plot_title_font_size Font size of plot title.
#' @param axis_label_font_size Font size of axis title.
#' @param tick_label_font_size Font size of axis tick marks.
#' @param plot_aspect_ratio Aspect ratio for plot. Default is 1 (square plot).
#' @param log_y Logical \code{TRUE} (default) or \code{FALSE}. Should y-axis be 
#' -log10 transformed? (\code{TRUE} recommended for plotting p-values.)
#' @param output_dir Directory path where plot should be saved; defaults to \code{NULL} 
#' (plot not saved to file).
#' @param plot_name File name for saved plot including extension (.pdf, .png, .jpg, .svg); 
#' defaults to \code{NULL}. If \code{output_dir} is specified but \code{plot_name} is \code{NULL}, 
#' the plot will be saved as "volcano_plot_YYYYMMDD.pdf" with current date.
#' @param plot_width Numeric value for the width of the saved plot in inches; defaults to 7.
#' @param plot_height Numeric value for the height of the saved plot in inches; defaults to 7.
#'
#' @return Outputs a ggplot object. If \code{output_dir} is specified, also saves the plot 
#' in the specified format (PDF, PNG, JPG, or SVG based on file extension in \code{plot_name}).
#'
#' 
#'
#' @export
#'
#'
volcanoPlot <- function(coefs,
                        p_vals,
                        target_labels,
                        title = NULL,
                        xlabel = expression("log"[2] * "(fold change)"),
                        ylabel = expression("-log"[10] * "(FDR-adjusted p-value)"),
                        xlimits = NULL,
                        ylimits = NULL,
                        sig_threshold = 0.05,
                        sig_label = "FDR = 5%",
                        label_all_targets = FALSE,
                        upper_log2FC_threshold = NULL,
                        lower_log2FC_threshold = NULL,
                        target_labels_off = FALSE,
                        target_label_colors = NULL,
                        target_point_colors = NULL,
                        target_label_size = 2,
                        target_label_segment_color = 'grey',
                        max.overlaps = Inf,
                        force = 1,
                        force_pull = 1,
                        plot_title_font_size = 14,
                        axis_label_font_size = 12,
                        tick_label_font_size = 12,
                        plot_aspect_ratio = 1,
                        log_y = TRUE,
                        output_dir = NULL,
                        plot_name = NULL,
                        plot_width = 7,
                        plot_height = 7){
  
  # Check if p_vals is a list (dual mode) or vector (single mode)
  is_dual_mode <- is.list(p_vals) && !is.null(names(p_vals))
  
  if(is_dual_mode){
    # Validate dual mode inputs
    if(!all(c("unadj", "fdr") %in% names(p_vals))){
      stop("Error: When p_vals is a list, it must contain names 'unadj' and 'fdr'")
    }
    
    # Check sig_threshold
    if(is.list(sig_threshold)){
      if(!all(c("unadj", "fdr") %in% names(sig_threshold))){
        stop("Error: When sig_threshold is a list, it must contain names 'unadj' and 'fdr'")
      }
      sig_thresh_unadj <- sig_threshold$unadj
      sig_thresh_fdr <- sig_threshold$fdr
    } else {
      # Use same threshold for both
      sig_thresh_unadj <- sig_threshold
      sig_thresh_fdr <- sig_threshold
    }
    
    # Auto-override ylabel for dual mode if it's still the single mode default
    if(identical(ylabel, expression("-log"[10] * "(FDR-adjusted p-value)"))){
      ylabel <- expression("-log"[10] * "(unadjusted p-value)")
    }
    
    # Validate sig_label is not a list
    if(is.list(sig_label)){
      stop("Error: sig_label must be a single string, not a list. In dual mode, only the unadjusted threshold line is displayed.")
    }
    
    # Auto-override sig_label for dual mode if it's still the single mode default
    if(identical(sig_label, "FDR = 5%")){
      sig_label <- paste0("Unadj = ", sig_thresh_unadj * 100, "%")
    }
    
    # Extract p-values
    p_vals_unadj <- p_vals$unadj
    p_vals_fdr <- p_vals$fdr
    p_vals_plot <- p_vals_unadj  # Plot unadjusted on y-axis
    
  } else {
    # Single mode (original behavior)
    if(is.list(sig_threshold)){
      stop("Error: sig_threshold cannot be a list when p_vals is a vector")
    }
    
    # Validate sig_label is not a list
    if(is.list(sig_label)){
      stop("Error: sig_label must be a single string, not a list")
    }
    
    # ylabel and sig_label already have appropriate defaults in function signature
    p_vals_plot <- p_vals
    sig_thresh_single <- sig_threshold
  }
  
  # Handle labeling logic
  if(is_dual_mode){
    # Dual mode: label based on either unadj or fdr significance
    if(label_all_targets==FALSE){
      # Label only significant targets (either unadj or fdr significant)
      target_labels[p_vals_unadj > sig_thresh_unadj & p_vals_fdr > sig_thresh_fdr] <- ''
    }
  } else {
    # Single mode: original behavior
    if(label_all_targets==FALSE){
      target_labels[p_vals_plot > sig_thresh_single] <- ''
    }
  }
  
  if(target_labels_off==TRUE){
    target_labels <- rep('', length(target_labels))
  }
  
  # Handle fold change threshold labeling (works for both modes)
  if(!is.null(upper_log2FC_threshold) && !is.null(lower_log2FC_threshold)){
    if(is_dual_mode){
      cutoff <- (coefs > upper_log2FC_threshold) | (coefs < lower_log2FC_threshold) | 
        (p_vals_unadj < sig_thresh_unadj) | (p_vals_fdr < sig_thresh_fdr)
    } else {
      cutoff <- (coefs > upper_log2FC_threshold) | (coefs < lower_log2FC_threshold) | 
        (p_vals_plot < sig_thresh_single)
    }
    target_labels[cutoff == 0] <- ''
  }
  if(!is.null(upper_log2FC_threshold) && is.null(lower_log2FC_threshold)){
    if(is_dual_mode){
      cutoff <- (coefs > upper_log2FC_threshold) | (p_vals_unadj < sig_thresh_unadj) | 
        (p_vals_fdr < sig_thresh_fdr)
    } else {
      cutoff <- (coefs > upper_log2FC_threshold) | (p_vals_plot < sig_thresh_single)
    }
    target_labels[cutoff == 0] <- ''
  }
  if(!is.null(lower_log2FC_threshold) && is.null(upper_log2FC_threshold)){
    if(is_dual_mode){
      cutoff <- (coefs < lower_log2FC_threshold) | (p_vals_unadj < sig_thresh_unadj) | 
        (p_vals_fdr < sig_thresh_fdr)
    } else {
      cutoff <- (coefs < lower_log2FC_threshold) | (p_vals_plot < sig_thresh_single)
    }
    target_labels[cutoff == 0] <- ''
  }
  # create point colors
  red_blue <- alamarColorPalette(n=2, palette=2)
  
  if(is.null(target_point_colors)){
    target_point_colors <- rep('grey', length(target_labels))
    
    if(is_dual_mode){
      # Light colors for unadjusted significance only
      # Darker colors for FDR significance
      light_red <- scales::alpha(red_blue[2], 0.4)
      light_blue <- scales::alpha(red_blue[1], 0.4)
      dark_red <- red_blue[2]
      dark_blue <- red_blue[1]
      
      # Unadjusted significant only (light colors)
      target_point_colors[p_vals_unadj < sig_thresh_unadj & p_vals_fdr >= sig_thresh_fdr & coefs > 0] <- light_red
      target_point_colors[p_vals_unadj < sig_thresh_unadj & p_vals_fdr >= sig_thresh_fdr & coefs < 0] <- light_blue
      
      # FDR significant (dark colors, overrides unadjusted)
      target_point_colors[p_vals_fdr < sig_thresh_fdr & coefs > 0] <- dark_red
      target_point_colors[p_vals_fdr < sig_thresh_fdr & coefs < 0] <- dark_blue
      
      # Create significance category for legend
      sig_category <- rep('Not significant', length(target_labels))
      sig_category[p_vals_unadj < sig_thresh_unadj & p_vals_fdr >= sig_thresh_fdr & coefs > 0] <- 'Unadj sig (upregulated)'
      sig_category[p_vals_unadj < sig_thresh_unadj & p_vals_fdr >= sig_thresh_fdr & coefs < 0] <- 'Unadj sig (downregulated)'
      sig_category[p_vals_fdr < sig_thresh_fdr & coefs > 0] <- 'FDR sig (upregulated)'
      sig_category[p_vals_fdr < sig_thresh_fdr & coefs < 0] <- 'FDR sig (downregulated)'
      
    } else {
      # Single mode: original behavior
      target_point_colors[p_vals_plot < sig_thresh_single & coefs > 0] <- red_blue[2]
      target_point_colors[p_vals_plot < sig_thresh_single & coefs < 0] <- red_blue[1]
      sig_category <- NULL
    }
  } else {
    sig_category <- NULL
  }
  
  point_colors <- unique(target_point_colors)
  names(point_colors) <- point_colors
  # create label colors
  if(is.null(target_label_colors)){
    target_label_colors <- rep('grey', length(target_labels))
    
    if(is_dual_mode){
      # Use same logic as points for consistency
      light_red <- scales::alpha(red_blue[2], 0.6)
      light_blue <- scales::alpha(red_blue[1], 0.6)
      dark_red <- red_blue[2]
      dark_blue <- red_blue[1]
      
      target_label_colors[p_vals_unadj < sig_thresh_unadj & p_vals_fdr >= sig_thresh_fdr & coefs > 0] <- light_red
      target_label_colors[p_vals_unadj < sig_thresh_unadj & p_vals_fdr >= sig_thresh_fdr & coefs < 0] <- light_blue
      target_label_colors[p_vals_fdr < sig_thresh_fdr & coefs > 0] <- dark_red
      target_label_colors[p_vals_fdr < sig_thresh_fdr & coefs < 0] <- dark_blue
    } else {
      target_label_colors[p_vals_plot < sig_thresh_single & coefs > 0] <- red_blue[2]
      target_label_colors[p_vals_plot < sig_thresh_single & coefs < 0] <- red_blue[1]
    }
  }
  label_colors <- unique(target_label_colors)
  names(label_colors) <- label_colors
  if(is.null(xlimits)){
    xlimits <- c(min(coefs), max(coefs))
    xmin <- min(coefs)
  } else if(!is.null(xlimits)){
    xmin <- xlimits[1]
  }
  if(is.null(ylimits)){
    if (log_y == TRUE) {
      max_val <- max(-log10(p_vals_plot))
      if (is.finite(max_val)) {
        ymax <- max_val
      } else {
        ymax <- max(-log10(p_vals_plot)[is.finite(-log10(p_vals_plot))])
        warning("Some targets have p-values that are numerically equal to zero.")
      }
    }
    if (log_y==FALSE) ymax <- max(p_vals_plot)
    ylimits <- c(0, ymax)
  }
  # organize plot data
  if (log_y==TRUE) y_vals <- -log10(p_vals_plot)
  if (log_y==FALSE) y_vals <- (p_vals_plot)
  
  plot_data <- data.frame(coefs=coefs,
                          minus_log10_p_vals=y_vals,
                          target_labels=target_labels,
                          target_point_colors=target_point_colors)
  
  # Add significance category for legend in dual mode
  if(is_dual_mode && !is.null(sig_category)){
    plot_data$sig_category <- factor(sig_category, 
                                     levels = c('Not significant', 
                                                'Unadj sig (downregulated)',
                                                'Unadj sig (upregulated)',
                                                'FDR sig (downregulated)', 
                                                'FDR sig (upregulated)'))
  }
  
  # Create base plot
  if(is_dual_mode && !is.null(sig_category)){
    # Dual mode with legend
    light_red <- scales::alpha(red_blue[2], 0.4)
    light_blue <- scales::alpha(red_blue[1], 0.4)
    dark_red <- red_blue[2]
    dark_blue <- red_blue[1]
    
    color_mapping <- c('Not significant' = 'grey',
                       'Unadj sig (downregulated)' = light_blue,
                       'Unadj sig (upregulated)' = light_red,
                       'FDR sig (downregulated)' = dark_blue,
                       'FDR sig (upregulated)' = dark_red)
    
    volcano_plot <- ggplot2::ggplot(plot_data, 
                                    ggplot2::aes(x=coefs, 
                                                 y=minus_log10_p_vals,
                                                 label=target_labels)) + 
      ggplot2::xlim(xlimits) +
      ggplot2::ylim(ylimits) +
      ggplot2::geom_vline(xintercept = 0, color='grey') +
      ggplot2::geom_hline(yintercept = 0, color='grey') +
      ggplot2::geom_point(shape=20, ggplot2::aes(color=sig_category)) +
      ggplot2::scale_color_manual(values=color_mapping, name="Significance") +
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel) +
      ggplot2::labs(title=title) +
      ggplot2::theme_light() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face='bold', 
                                                        size=plot_title_font_size),
                     axis.title = ggplot2::element_text(size=axis_label_font_size),
                     axis.text = ggplot2::element_text(size=tick_label_font_size),
                     aspect.ratio = plot_aspect_ratio,
                     plot.margin = ggplot2::margin(4,4,4,4),
                     legend.position='right')
    
    # Add significance threshold line for unadjusted p-value (only one line in dual mode)
    volcano_plot <- volcano_plot +
      ggplot2::geom_hline(yintercept=-log10(sig_thresh_unadj), color='red') +
      ggplot2::annotate('text', x=xmin, y=-log10(sig_thresh_unadj), label=sig_label, 
                        color='red', size=2, hjust=0.25, vjust=-0.5,
                        fontface='italic')
    
  } else {
    # Single mode (original)
    volcano_plot <- ggplot2::ggplot(plot_data, 
                                    ggplot2::aes(x=coefs, 
                                                 y=minus_log10_p_vals,
                                                 label=target_labels)) + 
      ggplot2::xlim(xlimits) +
      ggplot2::ylim(ylimits) +
      ggplot2::geom_vline(xintercept = 0, color='grey') +
      ggplot2::geom_hline(yintercept = 0, color='grey') +
      ggplot2::geom_point(shape=20, ggplot2::aes(color=target_point_colors)) +
      ggplot2::scale_color_manual(values=point_colors) +
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel) +
      ggplot2::labs(title=title) +
      ggplot2::theme_light() +
      ggplot2::geom_hline(yintercept=-log10(sig_thresh_single), color='red') +
      ggplot2::annotate('text', x=xmin, y=-log10(sig_thresh_single), label=sig_label, 
                        color='red', size=2, hjust=0.25, vjust=-0.5,
                        fontface='italic') +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face='bold', 
                                                        size=plot_title_font_size),
                     axis.title = ggplot2::element_text(size=axis_label_font_size),
                     axis.text = ggplot2::element_text(size=tick_label_font_size),
                     aspect.ratio = plot_aspect_ratio,
                     plot.margin = ggplot2::margin(4,4,4,4),
                     legend.position='none')
  }
  
  # Add text labels (common to both modes)
  volcano_plot <- volcano_plot +
    ggrepel::geom_text_repel(size=target_label_size, 
                             max.overlaps = max.overlaps, 
                             force=force,
                             force_pull=force_pull,
                             segment.color=target_label_segment_color,
                             color=target_label_colors) +
    ggplot2::coord_cartesian(clip='off')
  
  # Add fold change threshold lines if specified
  if(!is.null(upper_log2FC_threshold) && !is.null(lower_log2FC_threshold)){
    volcano_plot <- volcano_plot + 
      ggplot2::geom_vline(xintercept=upper_log2FC_threshold, color='grey', linetype = "longdash") +
      ggplot2::annotate('text', x=upper_log2FC_threshold, y=0, 
                        label=paste0('log2FC = ', round(upper_log2FC_threshold, 2)), 
                        color='grey', size=2, hjust=0.25, vjust=-0.5,
                        fontface='italic') + 
      ggplot2::geom_vline(xintercept=lower_log2FC_threshold, color='grey', linetype = "longdash") +
      ggplot2::annotate('text', x=lower_log2FC_threshold, y=0, 
                        label=paste0('log2FC = ', round(lower_log2FC_threshold, 2)), 
                        color='grey', size=2, hjust=0.25, vjust=-0.5,
                        fontface='italic')
  }
  if(!is.null(upper_log2FC_threshold) && is.null(lower_log2FC_threshold)){
    volcano_plot <- volcano_plot + 
      ggplot2::geom_vline(xintercept=upper_log2FC_threshold, color='grey', linetype = "longdash") +
      ggplot2::annotate('text', x=upper_log2FC_threshold, y=0, 
                        label=paste0('log2FC = ', round(upper_log2FC_threshold, 2)), 
                        color='grey', size=2, hjust=0.25, vjust=-0.5,
                        fontface='italic') 
  }
  if(!is.null(lower_log2FC_threshold) && is.null(upper_log2FC_threshold)){
    volcano_plot <- volcano_plot + 
      ggplot2::geom_vline(xintercept=lower_log2FC_threshold, color='grey', linetype = "longdash") +
      ggplot2::annotate('text', x=lower_log2FC_threshold, y=0, 
                        label=paste0('log2FC = ', round(lower_log2FC_threshold, 2)), 
                        color='grey', size=2, hjust=0.25, vjust=-0.5,
                        fontface='italic')
  }
  
  # Save plot if output_dir is specified
  if(!is.null(output_dir)){
    # Check if plot_name is provided, if not use default
    if(is.null(plot_name)){
      plot_name <- paste0("volcano_plot_", format(Sys.time(), "%Y%m%d"), ".pdf")
      message("No plot_name provided. Using default filename: ", plot_name)
    }
    
    # Create output directory if it doesn't exist
    if(!dir.exists(output_dir)){
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Full path for the output file
    file_path <- file.path(output_dir, plot_name)
    
    # Determine file type from extension
    file_ext <- tools::file_ext(plot_name)
    
    # Save based on file extension
    if(file_ext == "pdf"){
      grDevices::pdf(file_path, width = plot_width, height = plot_height)
      print(volcano_plot)
      grDevices::dev.off()
    } else if(file_ext %in% c("png", "PNG")){
      grDevices::png(file_path, width = plot_width * 100, height = plot_height * 100, res = 100)
      print(volcano_plot)
      grDevices::dev.off()
    } else if(file_ext %in% c("jpg", "jpeg", "JPG", "JPEG")){
      grDevices::jpeg(file_path, width = plot_width * 100, height = plot_height * 100, res = 100)
      print(volcano_plot)
      grDevices::dev.off()
    } else if(file_ext %in% c("svg", "SVG")){
      svglite::svglite(file_path, width = plot_width, height = plot_height)
      print(volcano_plot)
      grDevices::dev.off()
    } else {
      warning("Unsupported file extension: ", file_ext, ". Saving as PDF instead.")
      file_path <- paste0(tools::file_path_sans_ext(file_path), ".pdf")
      grDevices::pdf(file_path, width = plot_width, height = plot_height)
      print(volcano_plot)
      grDevices::dev.off()
    }
    
    message("Volcano plot saved to: ", file_path)
  }
  
  return(volcano_plot)
}