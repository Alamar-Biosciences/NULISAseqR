#' Spaghetti plot for visualizing longitudinal NULISAseq data
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
#' @param max.overlaps Passed to ggrepel. Integer which determines how many targets 
#' will have labels. "Inf" (default) labels all targets. 
#' @param plot_title_font_size Font size of plot title.
#' @param axis_label_font_size Font size of axis title.
#' @param tick_label_font_size Font size of axis tick marks.
#' @param target_label_size Font size of target labels.
#' @param target_label_segment_color Color of line segments for target labels.
#' @param plot_aspect_ratio Aspect ratio for plot. Default is 1 (square plot).
#'
#' @return Outputs a ggplot object. 
#'
#' 
#'
#' @export
#'
#'
spaghetti_plot <- function(target, 
                           data, 
                           sampleInfo,
                           time_var_name,
                           id_var_name,
                           group_var_name=NULL,
                           axis_label_font_size=0.55, 
                           legend=TRUE){
  target_data <- data.frame(sampleName=colnames(data),
                            target_data=data[target,])
  plot_data <- merge(sampleInfo, target_data,
                      all.x=TRUE, all.y=FALSE,
                      by.x=sampleName_var, by.y='sampleName')

  # define colors
  colors <- alamarColorPalette(n=length(unique(plot_data$group_var_name)), nReps=5)
  

  mean_color <- 'black'
  control_color <- colors[[1]][2]
  

  # plot first patient
  plot(as.numeric(plot_data$visit[plot_data$id==unique_patients[1]]),
       plot_data$log2count[plot_data$id==unique_patients[1]],
       ylim=c(ymin, ymax2),
       xlim=c(0, 6),
       type='o',
       las=1,
       xlab='',
       ylab='',
       main=paste0(panel, ' - ', target),
       col=patient_colors[plot_data$traj_group[plot_data$id==unique_patients[1]][1]],
       xaxt='n')
  axis(side=1, at=0:6, labels=c('control', levels(plot_data$visit)), cex.axis=cex_x_axis)
  mtext('log2(normalized count)', side=2, line=2, cex=0.7)
  # loop over patients and draw lines
  for(i in 2:length(unique_patients)){
    lines(as.numeric(plot_data$visit[plot_data$id==unique_patients[i]]),
          plot_data$log2count[plot_data$id==unique_patients[i]],
          col=patient_colors[plot_data$traj_group[plot_data$id==unique_patients[i]][1]],
          type='o')
  }
  # plot controls
  points(rep(0, length(control_data)), control_data, col=colors[[1]][4])
  points(0, mean(control_data), col=colors[[1]][2], pch=16)
  # plot COVID mean fit line
  fixed_estimates <- fixef(target_data$full_model_visit)
  COVID_mean <- c(fixed_estimates[1],
                  fixed_estimates[1] + fixed_estimates[2],
                  fixed_estimates[1] + fixed_estimates[3],
                  fixed_estimates[1] + fixed_estimates[4],
                  fixed_estimates[1] + fixed_estimates[5],
                  fixed_estimates[1] + fixed_estimates[6])
  
  # plot the mean line
  lines(1:6, COVID_mean, col=mean_color, lwd=2, lty=1, pch=16, type='o')
  
  
  # add legend
  if(legend==TRUE){
    legend('topleft',
           legend=c('COVID mean',
                    'COVID trajectory group 1-3',
                    'COVID trajectory group 4-5',
                    'Control mean',
                    'Control sample'),
           col=c('black',
                 patient_colors[1],
                 patient_colors[2],
                 colors[[1]][2],
                 colors[[1]][4]),
           lwd=c(2, 1, 1, NA, NA), 
           pch=c(16, 1, 1, 16, 1),
           bty='n', cex=0.5, ncol=2) 
  }
  # add statistics
  # visit
  if (KR_result[target,]$visit_pval_FDR < 0.001){
    visit_pval_FDR <- 'p_FDR < 0.001'
  } else {
    visit_pval_FDR <- paste0('p_FDR = ', 
                             format(round(KR_result[target,]$visit_pval_FDR, 3), nsmall=3))
  }
  Fstat <- paste0('Visit: F = ', 
                  format(round(KR_result[target,]$visit_Fstat, 2), nsmall=2))
  if(stats==TRUE){
    legend('topright',
           legend=c(Fstat, visit_pval_FDR),
           bty='n', cex=0.5)
  }
}

