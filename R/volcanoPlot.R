#' Volcano plot for NULISAseq differential expression test
#'
#' Draws a volcano plot for a set of coefficients and p-values. Uses ggplot
#' for plotting and ggrepel package for labels. By default, plot labeling
#' assumes FDR-corrected p-values. 
#'
#' @param coefs Vector of coefficients.
#' @param p_vals Vector of p-values.
#' @param target_labels Vector of target names.
#' @param title Plot title.
#' @param xlabel X axis label.
#' @param ylabel Y axis label.
#' @param xlimits X axis limits
#' @param ylimits Y axis limits
#' @param sig_threshold P-value significance threshold.
#' @param sig_label Label to mark significance threshold on plot.
#' @param label_all_targets Logical TRUE or FALSE. 
#' Should all targets be labelled? Default is FALSE 
#' and will only label significant targets.
#' @param target_labels_off Logical TRUE or FALSE. Should no targets be labelled? 
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
#'
#' @return Outputs a ggplot object. 
#'
#' 
#'
#' @export
#'
#'
volcanoPlot <- function(coefs,
                        p_vals,
                        target_labels,
                        title=NULL,
                        xlabel=expression('log'[2]*'(fold change)'),
                        ylabel=expression('-log'[10]*'(FDR-adjusted p-value)'),
                        xlimits=NULL,
                        ylimits=NULL,
                        sig_threshold=0.05,
                        sig_label='FDR = 5%',
                        label_all_targets=FALSE,
                        target_labels_off=FALSE,
                        target_label_colors=NULL,
                        target_point_colors=NULL,
                        target_label_size=2,
                        target_label_segment_color='grey',
                        max.overlaps=Inf,
                        force=1,
                        force_pull=1,
                        plot_title_font_size=14,
                        axis_label_font_size=12,
                        tick_label_font_size=12,
                        plot_aspect_ratio=1,
                        log_y=TRUE){
  
  if(label_all_targets==FALSE){
    # do not label insignificant targets
    target_labels[p_vals > sig_threshold] <- ''
  }
  if(target_labels_off==TRUE){
    target_labels <- rep('', length(target_labels))
  }
  # create point colors
  red_blue <- alamarColorPalette(n=2, palette=2)
  if(is.null(target_point_colors)){
    target_point_colors <- rep('grey', length(target_labels))
    target_point_colors[p_vals < sig_threshold & coefs > 0] <- red_blue[2]
    target_point_colors[p_vals < sig_threshold & coefs < 0] <- red_blue[1]
  }
  point_colors <- unique(target_point_colors)
  names(point_colors) <- point_colors
  # create label colors
  if(is.null(target_label_colors)){
    target_label_colors <- rep('grey', length(target_labels))
    target_label_colors[p_vals < sig_threshold & coefs > 0] <- red_blue[2]
    target_label_colors[p_vals < sig_threshold & coefs < 0] <- red_blue[1]
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
      max_val <- max(-log10(p_vals))
      if (is.finite(max_val)) {
        ymax <- max_val
      } else {
        ymax <- max(-log10(p_vals)[is.finite(-log10(p_vals))])
        warning("Some targets have p-values that are numerically equal to zero.")
      }
    }
    if (log_y==FALSE) ymax <- max(p_vals)
    ylimits <- c(0, ymax)
  }
  # organize plot data
  if (log_y==TRUE) y_vals <- -log10(p_vals)
  if (log_y==FALSE) y_vals <- (p_vals)
  
  plot_data <- data.frame(coefs=coefs,
                          minus_log10_p_vals=y_vals,
                          target_labels=target_labels,
                          target_point_colors=target_point_colors)
  
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
    ggplot2::geom_hline(yintercept=-log10(sig_threshold), color='red') +
    ggplot2::annotate('text', x=xmin, y=-log10(sig_threshold), label=sig_label, 
                      color='red', size=2, hjust=0.25, vjust=-0.5,
                      fontface='italic') +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face='bold', 
                                                      size=plot_title_font_size),
                   axis.title = ggplot2::element_text(size=axis_label_font_size),
                   axis.text = ggplot2::element_text(size=tick_label_font_size),
                   aspect.ratio = plot_aspect_ratio,
                   plot.margin = ggplot2::margin(4,4,4,4),
                   legend.position='none') +
    ggrepel::geom_text_repel(size=target_label_size, 
                             max.overlaps = max.overlaps, 
                             force=force,
                             force_pull=force_pull,
                             segment.color=target_label_segment_color,
                             color=target_label_colors) +
    ggplot2::coord_cartesian(clip='off')
  return(volcano_plot)
}






