#' Get Reverse Curve Target Names
#'
#' Returns the target names for reverse curve targets, identified by
#' \code{Curve_Quant} starting with "R". These targets have no LOD and
#' should be excluded from detectability computation.
#'
#' @param targets A data frame with at least \code{targetName} and
#'   \code{Curve_Quant} columns (e.g., \code{data$targets}).
#'
#' @return A character vector of reverse curve target names. Returns
#'   \code{character(0)} if none are found or if \code{Curve_Quant}
#'   is not present.
#'
#' @keywords internal
get_reverse_curve_targets <- function(targets) {
  if (!"Curve_Quant" %in% names(targets)) return(character(0))
  ind <- which(substr(targets$Curve_Quant, 1, 1) == "R")
  if (length(ind) > 0) targets$targetName[ind] else character(0)
}


#' Get noDetectability Target Names
#'
#' Returns the target names that have the \code{noDetectability} modifier set
#' in the XML configuration. These include both reverse curve targets and
#' rare-case targets (e.g., variant-specific or PTM targets).
#'
#' @param targets A data frame with at least \code{targetName} and
#'   \code{noDetectability} columns (e.g., \code{data$targets}).
#'
#' @return A character vector of noDetectability target names. Returns
#'   \code{character(0)} if the column is not present or none are flagged.
#'
#' @keywords internal
get_noDetectability_targets <- function(targets) {
  if (!"noDetectability" %in% names(targets)) return(character(0))
  ind <- which(targets$noDetectability)
  if (length(ind) > 0) targets$targetName[ind] else character(0)
}


#' Compute Summary Statistics Row for Detectability Values
#'
#' Creates a single-row character matrix of formatted summary statistics
#' for a set of detectability values. Used internally by
#' \code{format_detectability_report()}.
#'
#' @param detectability_values Numeric vector of detectability percentages,
#'   already filtered for any noDetectability targets. Names are optional.
#' @param type_label Character string label for this row (e.g., \code{"all"}
#'   or a sample group name).
#' @param sample_number Integer number of samples.
#'
#' @return A 1-row character matrix with columns: type, \# samples, mean, sd,
#'   median, min, max, \# of detectable targets (\%). When all values are NA or
#'   the vector is empty, stat columns contain empty strings.
#'
#' @keywords internal
format_detectability_stats_row <- function(detectability_values, type_label, sample_number) {
  col_names <- c("type", "# samples", "mean", "sd", "median",
                 "min", "max", "# of detectable targets (%)")

  n_valid <- sum(!is.na(detectability_values))

  if (n_valid > 0) {
    n_detectable <- sum(detectability_values > 50, na.rm = TRUE)
    detectable_pct <- format(round(n_detectable / n_valid * 100, 1), nsmall = 1)
    mean_val <- format(round(mean(detectability_values, na.rm = TRUE), 1), nsmall = 1)
    sd_val <- format(round(sd(detectability_values, na.rm = TRUE), 1), nsmall = 1)
    # round(..., 1) for median/min/max normalizes per-plate formatting to match
    # the Overall section, which already rounded these stats
    median_val <- formatC(round(median(detectability_values, na.rm = TRUE), 1), width = 2)
    min_val <- formatC(round(min(detectability_values, na.rm = TRUE), 1), width = 2)
    max_val <- formatC(round(max(detectability_values, na.rm = TRUE), 1), width = 2)
  } else {
    n_detectable <- 0
    detectable_pct <- ""
    mean_val <- sd_val <- median_val <- min_val <- max_val <- ""
  }

  row <- cbind(
    type_label,
    sample_number,
    mean_val,
    sd_val,
    median_val,
    min_val,
    max_val,
    paste0(n_detectable, " / ", n_valid, " (", detectable_pct, "%)")
  )
  colnames(row) <- col_names
  row
}


#' Format Detectability Report for a Single Plate
#'
#' Takes the output of \code{detectability()} for a single plate and returns
#' formatted summary and target tables ready for rendering.
#'
#' @param detect_result Return value from \code{detectability()}.
#' @param noDetectability_targets Character vector of target names to exclude
#'   from summary statistics but keep in the target table. Default
#'   \code{character(0)}.
#' @param reverse_curve_targets Character vector of reverse curve target names
#'   to append as "High Abundance" rows for PLASMA or SERUM matrix types. Other
#'   matrix types will show as NA. Default \code{character(0)}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{summary}{Character matrix with columns: type, \# samples, mean, sd,
#'     median, min, max, \# of detectable targets (\%). One row per type
#'     ("all" + sample groups).}
#'   \item{targets}{Character matrix with rows = targets, columns = "all" +
#'     groups, values formatted with \code{format(round(x, 1), nsmall = 1)},
#'     plus "High Abundance" rows for reverse curve targets for PLASMA / SERUM.}
#' }
#'
#' Label Detectability Table for Display
#'
#' Converts a numeric detectability data frame into a display version where
#' reverse curve targets show "High Abundance" in plasma/serum columns and
#' NA elsewhere. Non-RC rows are formatted as character with one decimal place.
#'
#' @param detect_table A data frame with a \code{Target} column and numeric
#'   detectability columns (as returned in \code{merged$detectability}).
#' @param reverse_curve_targets Character vector of reverse curve target names.
#'
#' @return A data frame with character value columns suitable for display.
#' @keywords internal
label_detectability_for_display <- function(detect_table, reverse_curve_targets = character(0)) {
  display <- detect_table
  for (col in colnames(display)[-1]) {
    vals <- display[[col]]
    display[[col]] <- ifelse(is.na(vals), NA_character_,
                             format(round(vals, 1), nsmall = 1, trim = TRUE))
  }
  if (length(reverse_curve_targets) > 0) {
    is_rc <- display$Target %in% reverse_curve_targets
    for (col in colnames(display)[-1]) {
      is_plasma_serum <- grepl("^(plasma|serum)\\s*\\(", col, ignore.case = TRUE)
      if (is_plasma_serum) {
        display[[col]][is_rc] <- "High Abundance"
      }
    }
  }
  display
}

#' @keywords internal
format_detectability_report <- function(detect_result,
                                        noDetectability_targets = character(0),
                                        reverse_curve_targets = character(0)) {
  summary_matrix <- NULL
  targets_matrix <- NULL
  col_names_targets <- NULL

  # "all" row
  raw <- detect_result$all$detectability
  summary_vals <- raw[!names(raw) %in% noDetectability_targets]
  summary_matrix <- rbind(summary_matrix, format_detectability_stats_row(
    summary_vals, "all", detect_result$all$sampleNumber
  ))
  col_names_targets <- c(col_names_targets, "all")
  targets_matrix <- cbind(targets_matrix, ifelse(is.na(raw), "", format(round(raw, 1), nsmall = 1, trim=TRUE)))

  # sample group rows
  if (!is.null(detect_result$sample_group)) {
    all_names <- names(detect_result$all$detectability)
    for (j in names(detect_result$sample_group$detectability)) {
      raw <- detect_result$sample_group$detectability[[j]]
      # align to "all" target order (fills NA for targets missing from this group)
      raw <- raw[all_names]
      names(raw) <- all_names
      summary_vals <- raw[!names(raw) %in% noDetectability_targets]
      summary_matrix <- rbind(summary_matrix, format_detectability_stats_row(
        summary_vals, j, detect_result$sample_group$sampleNumber[[j]]
      ))
      col_names_targets <- c(col_names_targets, j)
      targets_matrix <- cbind(targets_matrix, ifelse(is.na(raw), "", format(round(raw, 1), nsmall = 1, trim=TRUE)))
    }
  }

  colnames(targets_matrix) <- col_names_targets

  # Add "High Abundance" for reverse curve targets for PLASMA or SERUM columns only
  # Otherwise put NA
  if (length(reverse_curve_targets) > 0) {
    is_plasma_serum <- grepl("^(PLASMA|SERUM)$", colnames(targets_matrix), ignore.case = TRUE)
    rc_matrix <- matrix(ifelse(is_plasma_serum, "High Abundance", NA_character_),
                        nrow = length(reverse_curve_targets),
                        ncol = ncol(targets_matrix),
                        byrow = TRUE)
    rownames(rc_matrix) <- reverse_curve_targets
    colnames(rc_matrix) <- colnames(targets_matrix)
    targets_matrix <- rbind(targets_matrix, rc_matrix)
  }

  # Sort target rows alphabetically by name
  targets_matrix <- targets_matrix[order(tolower(rownames(targets_matrix))), , drop = FALSE]

  list(summary = summary_matrix, targets = targets_matrix)
}


#' Calculate Detectability for a Set of Targets and Samples
#'
#' Calculates detectability for each target based on the
#' aboveLOD matrix output of lod() function. Detectability is percent of samples
#' above LOD for that target.
#'
#' @param aboveLOD_matrix aboveLOD matrix output from lod() function.
#' Rows are targets and columns are samples. 
#' @param sample_subset A vector of column names or numeric column indices that
#' represent the sample subset that detectability will be calculated for. 
#' Default uses all columns. NCs and IPCs should probably be excluded.
#' @param sample_groups A string vector defining sample types / identities that
#' represent the sample type subsets that detectability should be calculated 
#' for. Default value NULL provides overall detectability assuming all samples
#' are of the same type (e.g., plasma). NOTE that sample groups must be the 
#' same length and in the same order as \code{sample_subset}. 
#' @param exclude_targets A vector of row names or numeric row indices 
#' representing the 
#' targets that should be excluded from detectability calculation. For example,
#' one might want to exclude internal controls. Default is NULL, which includes
#' all targets in the aboveLOD_matrix.
#'
#'
#' @return A nested list containing detectability information.
#' The list has two high-level names: "all" and "sample_group".
#' The "all" name corresponds to detectability information for all samples
#' combined as a single unit. If the \code{sample_groups} parameter is not NULL,
#' the "sample_group" name provides detectability data for unique sample groups.
#' The nested list has three names to access detectability and sample information.
#' \item{sampleNumber}{The number of samples for each group.}
#' \item{detectability}{A named vector of detectability values in percentages.
#' The names correspond to target row names.}
#' \item{detectable}{A named logical vector indicating whether detectability is
#' above 50\% for each target. The names correspond to target row names.}
#'
#'
#' @export
#' 
detectability <- function(aboveLOD_matrix, 
                          sample_subset=NULL,
                          sample_groups=NULL,
                          exclude_targets=NULL){
  # return named list
  detectability_data <- list()
  
  # remove excluded targets
  if (!is.null(exclude_targets)){
    if (!is.numeric(exclude_targets)){
      exclude_targets <- which(rownames(aboveLOD_matrix) %in% exclude_targets)
    }
    aboveLOD_matrix <- aboveLOD_matrix[-exclude_targets,]
  }
  
  # get sample subset
  if (!is.null(sample_subset)){
    if (!is.numeric(sample_subset)){
      sample_subset <- which(colnames(aboveLOD_matrix) %in% sample_subset)
    }
    aboveLOD_matrix <- aboveLOD_matrix[,sample_subset, drop=FALSE]
  }
  
  # get detectability data for samples: Plasma, serum, etc.
  if (!is.null(sample_groups)){
    # make case insensitive
    sample_groups <- toupper(sample_groups)
    
    # collect sample group specific detectability data
    for (i in unique(sample_groups)){
      if(is.null(dim(aboveLOD_matrix))){
        aboveLOD_matrix <- matrix(aboveLOD_matrix, ncol=1)
      }
      sample_aboveLOD_matrix <- as.matrix(aboveLOD_matrix[,sample_groups %in% i], nrow=nrow(aboveLOD_matrix))
      
      detectability_data$sample_group$sampleNumber[[i]] <- ncol(sample_aboveLOD_matrix)
      detectability_data$sample_group$detectability[[i]] <- apply(as.matrix(sample_aboveLOD_matrix), 1, function(x) sum(x, na.rm=TRUE) / sum(!is.na(x))*100)
      detectability_data$sample_group$detectable[[i]] <- detectability_data$sample_group$detectability[[i]] > 50
    }
  }
  
  # get overall detectability data
  detect <- apply(as.matrix(aboveLOD_matrix), 1, function(x) sum(x, na.rm=TRUE) / sum(!is.na(x))*100)
  detectable <- detect > 50
  if(is.null(dim(aboveLOD_matrix))){
    names(detect) <- rownames(as.matrix(aboveLOD_matrix))
    names(detectable) <- rownames(as.matrix(aboveLOD_matrix))
  } else{
    names(detect) <- rownames(aboveLOD_matrix)
    names(detectable) <- rownames(aboveLOD_matrix)
  }
  
  detectability_data$all$sampleNumber <- ncol(aboveLOD_matrix)
  detectability_data$all$detectability <- detect
  detectability_data$all$detectable <- detectable
  
  return(detectability_data)
}




#' Summarize Detectability Across Multiple Runs and Sample Groups
#'
#' Summarizes detectability across multiple runs and if applicable
#' broken down by sample groups based on the sample_group_covar input
#' to loadNULISAseq. When \code{format = TRUE}, returns formatted summary
#' and target tables ready for rendering, similar to
#' \code{quantifiability()$summary_tables}.
#'
#' @param runs A named list of run data output from \code{loadNULISAseq()} function
#' or a list of these outputs for multiple runs. To make output more interpretable,
#' it is recommended to name each run according to the plate ID for that run.
#' @param exclude_targets Optional targets to exclude from the aggregation. Either
#'   a character vector (applied to all plates) or a list of character vectors
#'   (one per plate, matching the order of \code{runs}).
#' @param exclude_noDetect_targets Logical. When \code{TRUE} (default),
#'   auto-detects reverse curve targets (excluded from computation entirely)
#'   and non-RC noDetectability targets (excluded from summary statistics only,
#'   kept in target tables). Emits \code{message()} listing affected targets.
#' @param format Logical. When \code{TRUE} (default), returns additional
#'   \code{run_summary} and \code{run_targets} fields containing formatted
#'   tables ready for rendering. When \code{FALSE}, returns only the original
#'   raw fields for backward compatibility.
#'
#'
#' @return A list with items "all", optionally "sample_group",
#' "run_detectability", and when \code{format = TRUE}, additional formatted
#' output fields. The structure matches \code{detectability()} output:
#' \item{all}{A list with \code{sampleNumber} (integer), \code{detectability}
#'   (named numeric vector of percentages), and \code{detectable} (named
#'   logical vector, TRUE when detectability > 50\%).}
#' \item{sample_group}{If applicable, a list with \code{sampleNumber} (named
#'   list of integers), \code{detectability} (named list of named numeric
#'   vectors), and \code{detectable} (named list of named logical vectors).}
#' \item{run_detectability}{A list of the original individual run detectability
#'   objects from \code{loadNULISAseq()}.}
#' \item{reverse_curve_targets}{Character vector of reverse curve target names
#'   detected across all runs (when \code{exclude_noDetect_targets = TRUE}).}
#' \item{noDetect_nonRC_targets}{Character vector of non-RC noDetectability
#'   target names detected across all runs (when
#'   \code{exclude_noDetect_targets = TRUE}).}
#' \item{run_summary}{Named list of formatted summary matrices per plate (and
#'   "Overall" for multi-plate), with columns: type, \# samples, mean, sd,
#'   median, min, max, \# of detectable targets (\%). Only when
#'   \code{format = TRUE}.}
#' \item{run_targets}{Named list of formatted target matrices per plate (and
#'   "Overall" for multi-plate), with rows = targets and columns = "all" +
#'   sample groups. Only when \code{format = TRUE}.}
#'
#'
#'
#' @export
#'
detectability_summary <- function(runs, exclude_targets = NULL,
                                  exclude_noDetect_targets = TRUE,
                                  format = TRUE){
  
  # check if run data is not in a list and if so put into a list
  if('RunSummary' %in% names(runs)){
    runs <- list(runs)
    names(runs) <- 'Plate 01'
  }

  # auto-detect reverse curve and noDetectability targets
  if (exclude_noDetect_targets) {
    all_reverseCurve <- unique(unlist(lapply(runs, function(x) get_reverse_curve_targets(x$targets))))
    all_noDetectability <- unique(unlist(lapply(runs, function(x) get_noDetectability_targets(x$targets))))
    all_noDetect_nonRC <- setdiff(all_noDetectability, all_reverseCurve)
    if (length(all_reverseCurve) > 0)
      message("Reverse curve targets excluded from detectability: ", paste(all_reverseCurve, collapse = ", "))
    if (length(all_noDetect_nonRC) > 0)
      message("Rare case targets excluded from summary statistics: ", paste(all_noDetect_nonRC, collapse = ", "))
  } else {
    all_reverseCurve <- character(0)
    all_noDetect_nonRC <- character(0)
  }

  # get detectability data for each run
  detect <- lapply(runs, function(x) x$detectability)

  # filter out excluded targets from per-run detectability before aggregation
  if (!is.null(exclude_targets)) {
    for (i in seq_along(detect)) {
      excl <- if (is.list(exclude_targets)) exclude_targets[[i]] else exclude_targets
      if (!is.null(excl) && length(excl) > 0) {
        detect[[i]]$all$detectability <- detect[[i]]$all$detectability[!names(detect[[i]]$all$detectability) %in% excl]
        detect[[i]]$all$detectable <- detect[[i]]$all$detectable[!names(detect[[i]]$all$detectable) %in% excl]
        if (!is.null(detect[[i]]$sample_group)) {
          for (g in names(detect[[i]]$sample_group$detectability)) {
            detect[[i]]$sample_group$detectability[[g]] <- detect[[i]]$sample_group$detectability[[g]][!names(detect[[i]]$sample_group$detectability[[g]]) %in% excl]
            detect[[i]]$sample_group$detectable[[g]] <- detect[[i]]$sample_group$detectable[[g]][!names(detect[[i]]$sample_group$detectable[[g]]) %in% excl]
          }
        }
      }
    }
  }

  # aggregate per-target detectability across plates for all samples
  all <- list()
  all$sampleNumber <- sum(sapply(detect, function(x) x$all$sampleNumber))
  all_detectability <- lapply(detect, function(x) {
    all_detect <- data.frame(Target=names(x$all$detectability), 
                             detectability=x$all$detectability / 100 * x$all$sampleNumber)
    return(all_detect)
  })
  all_detectability <- suppressWarnings(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Target", all = TRUE),
                                               all_detectability))
  if (ncol(all_detectability) < 2) {
    all$detectability <- numeric(0)
    all$detectable <- logical(0)
  } else {
    all_detectability$total_detectability <- rowSums(as.matrix(all_detectability[,2:ncol(all_detectability)]), na.rm = TRUE) / all$sampleNumber * 100
    detect_vec <- all_detectability$total_detectability
    names(detect_vec) <- all_detectability$Target
    all$detectability <- detect_vec
    all$detectable <- detect_vec > 50
  }
  
  output <- list(all=all,
                 run_detectability=detect)
  
  # get sample group names for each run
  group_names <- unique(unlist(lapply(detect, function(x) names(x$sample_group$sampleNumber))))
  names(group_names) <- group_names
  if(length(group_names) > 0){
    sample_group_detect <- lapply(group_names, function(x) {
      sampleNumber <- sum(unlist(lapply(detect, function(y) {
        sampleNumber <- y$sample_group$sampleNumber[[x]]
      })))
      sample_group_detectability <- lapply(detect, function(y) {
        if(x %in% names(y$sample_group$detectability)){
          data.frame(Target=names(y$sample_group$detectability[[x]]), 
                     detectability=y$sample_group$detectability[[x]] / 100 * y$sample_group$sampleNumber[[x]])
        } else {
          NULL
        }
      })
      
      sample_group_detectability <- sample_group_detectability[sapply(sample_group_detectability, function(x) !is.null(x))]
      sample_group_detectability <- suppressWarnings(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Target", all = TRUE),
                                                            sample_group_detectability))
      if (ncol(sample_group_detectability) < 2) {
        return(list(sampleNumber=sampleNumber,
                    detectability=numeric(0),
                    detectable=logical(0)))
      }
      sample_group_detectability$total_detectability <- rowSums(as.matrix(sample_group_detectability[,2:ncol(sample_group_detectability)]), na.rm = TRUE) / sampleNumber * 100
      detect_vec <- sample_group_detectability$total_detectability
      names(detect_vec) <- sample_group_detectability$Target

      return(list(sampleNumber=sampleNumber,
                  detectability=detect_vec,
                  detectable=detect_vec > 50))
    })

    # restructure into detectability()-compatible format
    sample_group <- list()
    sample_group$sampleNumber <- lapply(sample_group_detect, function(x) x$sampleNumber)
    sample_group$detectability <- lapply(sample_group_detect, function(x) x$detectability)
    sample_group$detectable <- lapply(sample_group_detect, function(x) x$detectable)

    output <- list(sample_group=sample_group,
                   all=all,
                   run_detectability=detect)
  }

  output$reverse_curve_targets <- all_reverseCurve
  output$noDetect_nonRC_targets <- all_noDetect_nonRC

  # build formatted tables when format = TRUE
  if (format) {
    run_summary <- list()
    run_targets <- list()

    # per-plate tables
    plate_names <- names(detect)
    if (is.null(plate_names)) plate_names <- paste("Plate", sprintf("%02d", seq_along(detect)))
    for (i in seq_along(detect)) {
      # per-plate RC and noDetect targets
      plate_rc <- if (exclude_noDetect_targets) {
        get_reverse_curve_targets(runs[[i]]$targets)
      } else {
        character(0)
      }
      plate_noDetect_nonRC <- if (exclude_noDetect_targets) {
        setdiff(get_noDetectability_targets(runs[[i]]$targets), plate_rc)
      } else {
        character(0)
      }

      plate_report <- format_detectability_report(
        detect_result = detect[[i]],
        noDetectability_targets = plate_noDetect_nonRC,
        reverse_curve_targets = plate_rc
      )
      run_summary[[plate_names[i]]] <- plate_report$summary
      run_targets[[plate_names[i]]] <- plate_report$targets
    }

    # overall tables (multi-plate)
    if (length(detect) > 1) {
      overall_report <- format_detectability_report(
        detect_result = output,
        noDetectability_targets = all_noDetect_nonRC,
        reverse_curve_targets = all_reverseCurve
      )
      run_summary[["Overall"]] <- overall_report$summary
      run_targets[["Overall"]] <- overall_report$targets
    }

    output$run_summary <- run_summary
    output$run_targets <- run_targets
  }

  return(output)
}
