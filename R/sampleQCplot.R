#' Creates Plots Showing Sample-Specific Quality Control Metrics
#'
#' This function generates QC plots for the following metrics:
#' \describe{
#'   \item{Detectability}{Percentage of targets with reads above the Limit of Detection (LOD). 
#'     Minimum thresholds vary by sample type. 
#'     \itemize{
#'       \item{Plasma: 90\%}
#'       \item{Serum: 90\%}
#'       \item{CSF: 70\%}
#'       \item{Urine: 65\%}
#'       \item{Cell culture: 30\%}
#'       \item{NHP plasma: 55\%}
#'       \item{NHP CSF: 35\%}
#'       \item{Dried blood spot: 75\%}
#'       \item{Control: 90\%}
#'       \item{Other: 0\%}
#'     }
#'   }
#'   \item{IC Reads}{Number of Internal Control (IC) reads within a sample. Minimum threshold: 1,000 reads}
#'   \item{Reads}{Total number of reads within a sample. Minimum threshold: 500,000 reads}
#'   \item{IC Median}{Sample IC reads relative to the plate median. Acceptable range: within +/-40\% of plate median}
#' }
#'   
#' 
#' For "all" option, plots all four QC metrics in a consolidated view.
#' @param runs A named list of run data output from \code{loadNULISAseq()} function 
#' or a list of these outputs for multiple runs.
#' @param qc_metric The QC metric to plot (default: "Detectability"; select other QC metrics: "IC Reads", "Reads", "IC Median", "all" (showing all 4 metrics)).
#' @param plateIDs Optional vector of plate names. (default: extract from runs).
#' @param output_dir Directory to save PDF files (NULL for no saving).
#' @param plot_name Optional name for the plot file (default: NULL).
#' @param plot_title Optional title for the plot (default: NULL).
#' @param plot_width PDF width in inches (default: 14 for single QC metric plot, 8 for all four QC metric plots).
#' @param plot_height PDF height in inches (default: 7 for single QC metric plot, 6 for all four QC metric plots).
#' @param label_cex Text size for labels (default: 0.45).
#' @param point_cex Point size (default: 0.75).
#' @param legend_cex Legend text size (default: 0.75).
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
#' @importFrom graphics plot axis text abline points legend par
#' @importFrom grDevices pdf dev.off
sampleQCplot <- function(runs,
                         qc_metric = "Detectability",
                         plateIDs = NULL,
                         output_dir = NULL,
                         plot_name = NULL,
                         plot_title = NULL,
                         plot_width = ifelse(qc_metric == "all", 8, 14),
                         plot_height = ifelse(qc_metric == "all", 6, 7),
                         label_cex = 0.45,
                         point_cex = 0.75,
                         legend_cex = 0.75,
                         sample_subset_list = NULL,
                         target_subset_list = NULL,
                         sample_exclude_list = NULL,
                         target_exclude_list = NULL) {
  #############
  ## Checks
  #############
  # Input validation checks
  if (!is.list(runs)) {
    stop("runs must be a list of run objects")
  }
  
  if (length(runs) == 0) {
    stop("runs list cannot be empty")
  }
  
  if(!qc_metric %in% c("Detectability", "IC Reads", "Reads", "IC Median", "all")) {
    stop("qc_metric must be one of 'Detectability', 'IC Reads', 'Reads', 'IC Median', or 'all'")
  }
  
  # Set control indicators
  indicatorSC <- all(vapply(runs, function(x) length(x$SC) > 0, logical(1)))
  indicatorBridge <- any(vapply(runs, function(x) length(x$Bridge) > 0, logical(1)))
  
  # Set default colors
  colors <- if (qc_metric == "all") {
    grDevices::hcl.colors(5, palette = "Set3")[c(4,3,1,5,2)]
  } else {
    unlist(lapply(alamarColorPalette(n = 5, nReps = 5), function(x) x[2]))
  }
  
  # Get plate IDs if not provided
  if (is.null(plateIDs)) {
    plateIDs <- vapply(runs, function(x) x$plateID, character(1))
  }
  
  # Create output directory if needed
  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  # (Removed graphics.off() to avoid interfering with user's graphics devices)
  
  if (qc_metric == "all") QCplot_metrics <- generate_QCplot_metrics(runs) 
  
  # Process each plate
  for (i in seq_along(runs)) {
    current_run <- runs[[i]]
    
    ###############
    ## Filter data
    ###############
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
    
    filtered_run <- filter_run_data(
      current_run,
      sample_subset = current_sample_subset,
      sample_exclude = current_sample_exclude,
      target_subset = current_target_subset,
      target_exclude = current_target_exclude
    )
    
    if (nrow(filtered_run$samples) == 0 || nrow(filtered_run$targets) == 0) {
      warning("No data remaining after filtering for Plate ", filtered_run$plateID, ". Skipping the plate.")
      next
    }
    
    # Sort samples by type
    sampleType_factor <- factor(filtered_run$samples$sampleType,
                                levels = c('Sample', 'IPC', 'SC', 'Bridge', 'NC'))
    well_order <- wellorder(filtered_run$samples,
                            sampleTypeFactor = sampleType_factor)
    
    # Perform QC flagging
    qcSample <- QCFlagSample(filtered_run$Data, 
                             filtered_run$lod$aboveLOD, 
                             filtered_run$samples,
                             filtered_run$targets, 
                             well_order,
                             TAP = TRUE, 
                             QCS = filtered_run$QCS, 
                             SN = filtered_run$SN)
    
    # Handle absolute quantification naming
    AbsData <- "Abs" %in% names(filtered_run$ExecutionDetails)
    if (AbsData) {
      qcSample$sampleName <- gsub("^IPC_([1234])$", "CAL_\\1", qcSample$sampleName)
      qcSample$sampleName <- gsub("^SC_([1234])$", "AQSC_\\1", qcSample$sampleName)
    }
    
    # Define sample types
    if (indicatorSC) {
      inds <- which(grepl(paste(filtered_run$SC, collapse = "|"), qcSample$sampleName))
      qcSample[inds, ]$sampleType <- "SC"
    }
    if (indicatorBridge) {
      inds <- which(grepl(paste(filtered_run$Bridge, collapse = "|"), qcSample$sampleName))
      qcSample[inds, ]$sampleType <- "Bridge"
    }
    
    # Set up PDF output
    if (!is.null(output_dir)) {
      if(is.null(plot_name)) {
        plot_name_full <- if (qc_metric == "all") {
          paste0("sample_QC_all_metrics_", plateIDs[i])
        } else {
          paste0("sample_QC_", gsub(" ", "_", qc_metric), "_", plateIDs[i])
        }
      } else {
        plot_name_full <- if (qc_metric == "all") {
          paste0("sample_QC_all_metrics_", plateIDs[i], "_", plot_name)
        } else {
          paste0("sample_QC_", gsub(" ", "_", qc_metric), "_", plateIDs[i], "_", plot_name)
        }
      }
      
      pdf_file <- file.path(output_dir, paste0(plot_name_full, ".pdf"))
      grDevices::pdf(pdf_file, width = plot_width, height = plot_height)
      on.exit(grDevices::dev.off())
    }
    
    # Branch based on qc_metric
    if (qc_metric == "all") {
      # All metrics plot 
      criteria <- QCplot_metrics$criteria
      
      # Set up plot layout
      par(mfcol = c(1, length(criteria$thresholdNames) + 1), 
          mar = c(5, 0.05, 2, 0.08))
      
      # Get all sample names for the first column (sample names)
      all_sample_names <- qcSample$sampleName
      if (AbsData) {
        all_sample_names <- gsub("^IPC_([1234])$", "CAL_\\1", all_sample_names)
        all_sample_names <- gsub("^SC_([1234])$", "AQSC_\\1", all_sample_names)
      }
      
      # Create empty plot for sample names (using all samples, not just inds)
      dotchart(x = rep(0, nrow(qcSample)),
               xlab = '',
               xlim = c(50, 100),
               labels = all_sample_names,
               color = "white",
               lcolor = 'white',
               gcolor = 'black',
               xaxt = 'n', 
               frame.plot = FALSE,
               offset = -1,
               cex = 0.45)
      title(main = filtered_run$plateID)
      
      # Create plots for each QC metric
      for(j in 1:length(names(criteria$thresholdNames))) {
        inds <- which(qcSample$flagName == names(criteria$thresholdNames)[j])
        
        # Set up plot labels and values based on metric type
        xlab <- ""
        xval <- NULL
        if(names(criteria$thresholdNames)[j] == "Detectability") {
          xlab <- "Detectability %"
          xval <- as.numeric(qcSample[inds, ]$val) * 100
        } else if(names(criteria$thresholdNames)[j] == "IC_Median") {
          xlab <- "% of IC Median"
          xval <- as.numeric(qcSample[inds, ]$val) * 100
        } else if(names(criteria$thresholdNames)[j] == "QCS") {
          xlab <- "QCS"
          xval <- as.numeric(qcSample[inds, ]$val)
        } else if(names(criteria$thresholdNames)[j] == "SN") {
          xlab <- "SN / QCS"
          xval <- as.numeric(qcSample[inds, ]$val)
        } else {
          xlab <- "Number of Reads"
          xval <- as.numeric(qcSample[inds, ]$val)
          zs <- which(xval == 0)
          if(length(zs) > 0) {
            xval[zs] <- xval[zs] + 1
          }
        }
        
        # Determine if log scale should be used
        log <- if(names(criteria$thresholdNames)[j] != "Detectability" && 
                  names(criteria$thresholdNames)[j] != "IC_Median" && 
                  names(criteria$thresholdNames)[j] != "QCS" && 
                  names(criteria$thresholdNames)[j] != "SN") "x" else ""
        
        # Set labels (only for first plot)
        if(j == 1) {
          labels <- qcSample[inds, ]$sampleName
          if(AbsData) {
            labels <- gsub("^IPC_([1234])$", "CAL_\\1", labels)
            labels <- gsub("^SC_([1234])$", "AQSC_\\1", labels)
          }
        } else {
          labels <- rep('', length(xval))
        }
        
        # Set colors based on sample type
        color <- rep("black", length(inds))
        color[which(qcSample[inds,]$sampleType == "IPC")] <- colors[2]
        color[which(qcSample[inds,]$sampleType == "NC")] <- colors[3]
        if(length(which(qcSample[inds,]$sampleType == "SC")) > 0) {
          color[which(qcSample[inds,]$sampleType == "SC")] <- colors[4]
        }
        if(length(which(qcSample[inds,]$sampleType == "Bridge")) > 0) {
          color[which(qcSample[inds,]$sampleType == "Bridge")] <- colors[5]
        }
        
        # Get axis limits from pre-calculated metrics
        valMin <- if(is.finite(QCplot_metrics$QCplot_axis_limits[[i]]$rmin[j])) 
          QCplot_metrics$QCplot_axis_limits[[i]]$rmin[j] else 0
        valMax <- if(is.finite(QCplot_metrics$QCplot_axis_limits[[i]]$rmax[j])) 
          QCplot_metrics$QCplot_axis_limits[[i]]$rmax[j] else 1
        
        # Create the dotchart plot
        dotchart(x = as.numeric(xval),
                 xlab = xlab,
                 labels = rep('', length(xval)), 
                 las = 1,
                 cex = 0.45, 
                 xlim = c(valMin, valMax),
                 color = color, 
                 log = log)
        
        # Add sample names to first plot
        if(j == 1) {
          for(k in 1:length(xval)) {
            axis(side = 2, 
                 at = k,
                 labels = labels[k], 
                 col.axis = color[k], 
                 las = 1,
                 cex.axis = 0.5)
          }
        }
        
        # Set warning/fail colors and text
        warnFailTxt <- "Warning"
        warnFailColor <- '#ebc634'
        if(names(criteria$thresholdNames)[j] == "QCS" || names(criteria$thresholdNames)[j] == "SN") {
          warnFailTxt <- "Fail"
          warnFailColor <- 'red'
        }
        
        # Plot points with appropriate colors
        xvalFail <- xval
        xvalFail[which(qcSample[inds,]$status != TRUE)] <- NA
        xvalPass <- xval
        xvalPass[which(qcSample[inds,]$status != FALSE)] <- NA
        points(as.numeric(xvalFail), 1:length(xvalFail), pch = 19, col = warnFailColor)
        points(as.numeric(xvalPass), 1:length(xvalPass), pch = 19, col = "green")
        
        # Add threshold lines
        if(names(criteria$thresholdNames)[j] == "IC_Median") {
          minmaxVal <- unlist(strsplit(criteria$thresholds["IC_Median"], ","))
          abline(v = as.numeric(minmaxVal[[1]]) * 100, col = 'brown')
          abline(v = as.numeric(minmaxVal[[2]]) * 100, col = 'brown')
          abline(v = 0, col = 'brown', lty = 2)
        } else if(names(criteria$thresholdNames)[j] == "Detectability") {
          value <- if(criteria$format[j] == "percentage") as.numeric(QCplot_metrics$val[[j]]) * 100 else as.numeric(QCplot_metrics$val[[j]])
          detectVals <- as.numeric(qcSample$QCthreshold[which(qcSample$flagName == "Detectability")]) * 100
          for(k in 1:length(detectVals)) {
            if(!is.na(detectVals[k])) {
              segments(detectVals[k], k - 0.5, detectVals[k], k + 0.5, col = "brown")
            }
          }
        } else {
          thresh <- criteria$thresholds[which(names(criteria$format)[j] == names(criteria$thresholds))]
          value <- if(criteria$format[j] == "percentage") as.numeric(thresh) * 100 else as.numeric(thresh)
          abline(v = value, col = 'brown')
        }
        
        # Add legend
        legend('bottomright', 
               legend = c('Pass', warnFailTxt), 
               col = c('green', warnFailColor), 
               pch = 19, 
               cex = 0.4, 
               bty = 'n', 
               inset = c(0, 1), 
               xpd = TRUE, 
               horiz = TRUE)
        
        title(main = criteria$properNames[j], cex.main = 1)
      }
      
    } else {
      # Single metric plot 
      # Set metric-specific parameters
      if(qc_metric == "Detectability") {
        ylabel <- "Detectability %"
        qc <- qc_metric_text <- "Detectability"
        y_limits <- c(0, 100)
      } else if(qc_metric == "IC Reads") {
        ylabel <- "Number of Reads"
        qc_metric_text <- "Internal Control Reads"
        qc <- "ICReads"
      } else if(qc_metric == "Reads") {
        ylabel <- "Number of Reads"
        qc_metric_text <- "Total Reads"
        qc <- "NumReads"
      } else if(qc_metric == "IC Median") {
        ylabel <- "% of IC Median"
        qc_metric_text <- "Internal Control Median"
        qc <- "IC_Median"
        y_limits <- c(-100, 100)
      }
      
      # Get QC metrics values
      inds <- which(qcSample$flagName == qc)
      if(qc_metric %in% c("IC Reads", "Reads")){
        xval <- as.numeric(qcSample[inds, ]$val)
        zs <- which(xval == 0)
        if(length(zs) > 0){
          xval[zs] <- xval[zs]+1
        }
        xval <- log10(xval)
        rmax <- max(xval, na.rm = TRUE) 
        rmin <- min(xval, na.rm = TRUE) 
        rlim_delta <- log10(rmax) - log10(rmin)
        rmax <- 10^(log10(rmax) + 0.025 * rlim_delta)
        y_limits <- c(1, rmax)
      } else {
        xval <- as.numeric(qcSample[inds, ]$val) * 100
        xval[xval < -100] <- -100
        xval[xval > 100] <- 100
      }
      
      # Define labels and colors
      labels <- qcSample[inds, ]$sampleName
      color <- rep("black", length(inds))
      color[which(qcSample[inds, ]$sampleType == "IPC")] <- colors[2]
      color[which(qcSample[inds, ]$sampleType == "NC")] <- colors[3]
      if (indicatorSC) {
        color[which(qcSample[inds, ]$sampleType == "SC")] <- colors[4]
      }
      if (indicatorBridge) {
        color[which(qcSample[inds, ]$sampleType == "Bridge")] <- colors[5]
      }
      
      yaxis_label <- ifelse(ylabel == "Number of Reads", "n", "s")
      
      title <- ifelse(is.null(plot_title), 
                      paste0(plateIDs[i], ' Sample QC: ', qc_metric_text), 
                      paste0(plot_title, ' , ', plateIDs[i], ' Sample QC: ', qc_metric_text))
      
      # Set up plot parameters
      opar <- par(no.readonly = TRUE)
      on.exit(par(opar), add = TRUE)
      par(mfcol = c(1, 1), mar = c(5, 4, 2, 0.5))
      
      # Create base plot
      plot(x = seq_along(xval),
           y = rev(xval),
           ylab = ylabel,
           xlab = '',
           xaxt = 'n',
           yaxt = yaxis_label,
           las = 1,
           cex = point_cex,
           ylim = y_limits,
           col = 'white',
           main = title)
      
      if (qc_metric %in% c("Reads", "IC Reads")) {
        axis(2, at = seq_along(xval), labels = format(10^rev(xval), scientific = TRUE, digits = 1), las = 1, cex.axis = 0.5)
      } 
      
      # Add custom x-axis labels
      axis(side = 1, at = seq_along(xval), labels = FALSE, tck = -0.01)
      for (k in seq_along(xval)) {
        text(x = k,
             y = par('usr')[3] - diff(par('usr')[3:4]) * 0.05,
             labels = rev(labels)[k],
             xpd = NA,
             col = rev(color)[k],
             srt = 45,
             adj = 0.9,
             cex = label_cex)
      }
      
      # Add grid lines
      abline(v = seq_along(xval), col = 'lightgray', lty = 'dotted')
      
      # Add points colored by QC status
      xvalFail <- xval
      xvalFail[which(qcSample[inds, ]$status != TRUE)] <- NA
      xvalPass <- xval
      xvalPass[which(qcSample[inds, ]$status != FALSE)] <- NA
      
      points(seq_along(xvalPass), rev(xvalPass), pch = 19, col = "green", cex = point_cex)
      points(seq_along(xvalFail), rev(xvalFail), pch = 19, col = "red", cex = point_cex)
      
      # Add threshold lines
      if(qc_metric == "IC Median") {
        threshold <- 40
        abline(h = c(-threshold, threshold), col = 'gray41')
        abline(h = 0, col = 'gray41', lty = 2)
      } else if(qc_metric == "Detectability") {
        threshold <- as.numeric(qcSample[which(qcSample$flagName == qc),]$QCthreshold[[1]]) * 100
        abline(h = threshold, col = 'gray41')
      } else {
        threshold <- log10(as.numeric(qcSample[which(qcSample$flagName == qc),]$QCthreshold[[1]]))
        abline(h = threshold, col = 'gray41')
      }
      
      # Add legend
      legend('bottomright',
             legend = c('Pass', "Warning"),
             col = c('green', 'red'),
             pch = 19,
             cex = legend_cex,
             bty = 'n',
             inset = c(0, 1),
             xpd = TRUE,
             horiz = TRUE)
    }
    
    # Close the PDF device if output_dir was specified
    if (!is.null(output_dir)) {
      dev.off()
    }
  }
  
  invisible(NULL)
}

#' Filter Run Data for Each Plate
#'
#' Internal function to filter sample and target data from each plate in the \code{runs} list.
#' Filters applied to \code{samples}, \code{targets}, \code{Data}, \code{lod$aboveLOD}, and \code{normed$interNormData} based on sample and target names.
#'
#' @param run A single plate or item from the named list of run data output from \code{loadNULISAseq()} function 
#' or a list of these outputs for multiple runs.
#' @param sample_subset Optional vector of sample names to include
#' @param sample_exclude Optional vector of sample names to exclude
#' @param target_subset Optional vector of target names to include
#' @param target_exclude Optional vector of target names to exclude
#'
#' @return Filtered run object
#' @keywords internal
filter_run_data <- function(run, sample_subset = NULL, sample_exclude = NULL, 
                            target_subset = NULL, target_exclude = NULL) {
  
  # Validate filter inputs
  validate_filter <- function(filter, name) {
    if (!is.null(filter)) {
      if (!is.character(filter)) {
        stop(paste(name, "must be NULL or a character vector"))
      }
    }
  }
  
  validate_filter(sample_subset, "sample_subset")
  validate_filter(sample_exclude, "sample_exclude")
  validate_filter(target_subset, "target_subset")
  validate_filter(target_exclude, "target_exclude")
  
  # Always include controls
  control_samples <- c(run$IPC, run$NC, run$SC)
  if (!is.null(run$Bridge)) control_samples <- c(control_samples, run$Bridge)
  
  # Apply sample filters if provided
  if (!is.null(sample_subset)) {
    samples_to_keep <- unique(c(sample_subset, control_samples))
  } else {
    samples_to_keep <- run$samples$sampleName
  }
  
  if (!is.null(sample_exclude)) {
    samples_to_keep <- setdiff(samples_to_keep, sample_exclude)
  }
  
  # Apply target filters if provided
  targets_to_keep <- if (!is.null(target_subset)) {
    target_subset
  } else {
    run$targets$targetName
  }
  
  if (!is.null(target_exclude)) {
    targets_to_keep <- setdiff(targets_to_keep, target_exclude)
  }
  
  # Apply the filters
  run$samples <- run$samples[run$samples$sampleName %in% samples_to_keep, ]
  run$targets <- run$targets[run$targets$targetName %in% targets_to_keep, ]
  
  # Only subset data if we actually filtered
  if (!is.null(sample_subset) || !is.null(sample_exclude)) {
    run$Data <- run$Data[, colnames(run$Data) %in% samples_to_keep, drop = FALSE]
    run$lod$aboveLOD <- run$lod$aboveLOD[, colnames(run$lod$aboveLOD) %in% samples_to_keep, drop = FALSE]
    run$normed$interNormData[[1]] <- run$normed$interNormData[[1]][, colnames(run$normed$interNormData[[1]]) %in% samples_to_keep, drop = FALSE]
  }
  
  if (!is.null(target_subset) || !is.null(target_exclude)) {
    run$Data <- run$Data[rownames(run$Data) %in% targets_to_keep, , drop = FALSE]
    run$lod$aboveLOD <- run$lod$aboveLOD[rownames(run$lod$aboveLOD) %in% targets_to_keep, , drop = FALSE]
    run$normed$interNormData[[1]] <- run$normed$interNormData[[1]][rownames(run$normed$interNormData[[1]]) %in% targets_to_keep, , drop = FALSE]
  }
  
  return(run)
}


#' Generate QC Plot Metrics
#'
#' Internal function to calculate QC plot axis limits for sample QC and target QC plots and generate QC summary.
#'
#' @param runs A named list of run data output from \code{loadNULISAseq()} function 
#' or a list of these outputs for multiple runs.
#'
#' @return List containing QC metrics and axis limits
#' @keywords internal
generate_QCplot_metrics <- function(runs) {
  # make empty sample QC summary table
  criteria <- QCSampleCriteria()
  rowNameQC <- c()
  QCSummary <- matrix(nrow=length(runs) + 1, ncol=length(criteria$thresholds) - length(which(startsWith(names(criteria$thresholds), "Detectability"))) + 1 + 1)
  
  for (i in 1:length(runs)){
    QCSummary[i,1] <- length(runs[[i]]$SampleNames)
    rowNameQC <- c(rowNameQC,  runs[[i]]$plateID) 
  }
  
  val <- valT <- list()
  QCplot_axis_limits <- vector(mode='list', length=length(runs))
  QCplot_target_axis_limits <- vector(mode='list', length=length(runs))
  criteria <- QCSampleCriteria()
  
  for (i in 1:length(runs)){
    AbsData <- "Abs" %in% names(runs[[i]]$ExecutionDetails)
    if(AbsData){
      Tcriteria <- QCTargetCriteria(AQ=AbsData, advancedQC=FALSE)
    } else Tcriteria = NULL
    
    well_order <- wellorder(runs[[i]]$samples)
    qcSample <- runs[[i]]$qcSample
    
    rmin <- rmax <- rminT <- rmaxT <- c()
    for (j in 1:length(criteria$thresholdNames)){
      inds <- which(startsWith(names(criteria$thresholdNames)[j], qcSample$flagName))
      QC_criterion_j_data <- qcSample[inds, ]
      if (startsWith(names(criteria$thresholdNames)[j], "Detectability")){
        val[j] <- criteria$thresholds[1]
      }else{
        val[j] <- criteria$thresholds[which(names(criteria$thresholds) == names(criteria$thresholdNames)[j])]
      }
      
      val[j] <- criteria$thresholds[j]
      if(names(criteria$thresholdNames)[j] == "Detectability"){
        xval <- as.numeric(QC_criterion_j_data$val)*100
        rmax[j] <- 100.0
        rmin[j] <- 0.0 #min(as.numeric(val[[j]])*100, as.numeric(xval), na.rm=TRUE)
      }else if (names(criteria$thresholdNames)[j] == "IC_Median"){
        xval <- as.numeric(QC_criterion_j_data$val)*100
        minmaxVal <- unname(unlist(strsplit(val[[j]], ",")))
        rmin[j] <- min(as.numeric(minmaxVal[1])*100, as.numeric(xval), na.rm=TRUE)
        rmax[j] <- max(as.numeric(minmaxVal[2])*100, as.numeric(xval), na.rm=TRUE)
      }else{
        xval <- as.numeric(QC_criterion_j_data$val)
        rmin[j] <- min(as.numeric(val[[j]]), as.numeric(xval), na.rm=TRUE)
        rmax[j] <- max(as.numeric(val[[j]]), as.numeric(xval), na.rm=TRUE)
      }
      
      # add a small amount of padding to x axis limits
      if(names(criteria$thresholdNames)[j] %in% c("ICReads", "NumReads")){
        rlim_delta <- log10(rmax[j]) - log10(rmin[j])
        rmin[j] <- 10^(log10(rmin[j]) - 0.025 * rlim_delta)
        rmax[j] <- 10^(log10(rmax[j]) + 0.025 * rlim_delta)
      } else {
        rlim_delta <- rmax[j] - rmin[j]
        rmin[j] <- rmin[j] - 0.025 * rlim_delta
        rmax[j] <- rmax[j] + 0.025 * rlim_delta
      }
      # count how many samples fail QC (excludes controls)
      QCSummary[i, j+1] <- length(which(QC_criterion_j_data[QC_criterion_j_data$sampleName %in% runs[[i]]$SampleNames,]$status=="TRUE"))
    }
    QCplot_axis_limits[[i]] <- list(rmin=rmin, rmax=rmax)
    
    if(length(Tcriteria) > 0){
      for (j in 1:length(Tcriteria$thresholdNames)){
        # Calculate limits for Targets
        valT[j] <- Tcriteria$thresholds[which(names(Tcriteria$thresholds) == names(Tcriteria$thresholdNames)[j])]
        indsT <- which(startsWith(names(Tcriteria$thresholdNames)[j], runs[[i]]$qcTarget$flagName))  
        QC_criterion_Target_j_data <- runs[[i]]$qcTarget[indsT, ]
        
        xval <- as.numeric(QC_criterion_Target_j_data$val) * 100
        if(names(Tcriteria$thresholdNames)[j] == "Target_Conc_Accuracy"){
          minmaxValT <- unname(unlist(strsplit(valT[[j]], ",")))
          rminT[j] <- min(as.numeric(minmaxValT[1])*100, as.numeric(xval), na.rm=TRUE)
          rmaxT[j] <- max(as.numeric(minmaxValT[2])*100, as.numeric(xval), na.rm=TRUE)
        } else if(names(Tcriteria$thresholdNames)[j] %in%  c("Target_Conc_CV", "Target_Detectability", "Target_Min_Reads", "Target_Conc_CV_RQ")){
          rminT[j] <- min(as.numeric(valT[[j]])*100, as.numeric(xval), na.rm=TRUE)
          rmaxT[j] <- max(as.numeric(valT[[j]])*100, as.numeric(xval), na.rm=TRUE)
        } else if(names(Tcriteria$thresholdNames)[j] %in% c("Target_IPC_Min_Reads")){
          rminT[j] <- 0
          rmaxT[j] <- log2(max(as.numeric(valT[[j]]), as.numeric(xval), na.rm=TRUE))
        }
      }
      QCplot_target_axis_limits[[i]] <- list(rmin=rminT, rmax=rmaxT)
    }
  }
  
  if(length(Tcriteria) > 0){
    output <- list(
      criteria = criteria,
      Tcriteria = Tcriteria,
      QCplot_axis_limits = QCplot_axis_limits,
      QCplot_target_axis_limits = QCplot_target_axis_limits,
      QCSummary = QCSummary,
      val = val,
      Tval = valT
    )
  } else {
    output <- list(
      criteria = criteria,
      QCplot_axis_limits = QCplot_axis_limits,
      QCplot_target_axis_limits = QCplot_target_axis_limits,
      QCSummary = QCSummary,
      val = val
    )
  }
  # Return the three objects as a list
  return(output)
}
