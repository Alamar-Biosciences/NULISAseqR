#' Calculate Quantifiability for Each Target for a Set of AQ Runs
#'
#' Calculates quantifiability for each target for a single AQ run or 
#' a set of AQ runs. Quantifiability is the percent of samples
#' within the dynamic range for that target.
#' Function input is the output of loadNULISAseq or a list
#' of loadNULISAseq outputs. Includes option to calculate quantifiability for 
#' sample subsets using a column in the samples data frame. 
#'
#' @param runs A named list of run data output from \code{laodNULISAseq()} function 
#' or a list of these outputs for multiple runs. To make output more interpretable,
#' it is recommended to name each run according to the plate ID for that run.
#' @param sampleGroupCovar Character string. Optional name of column in the 
#' samples data frame 
#' of \code{runs} which represents subgroups for which to calculate 
#' quantifiability. Default is \code{'SAMPLE_MATRIX'}.
#' @param sample_subset A list of vectors in the same order as runs 
#' specifying column names or numeric column indices that
#' represent the sample subset that quantifiability will be calculated for each run. 
#' If a vector is provided instead of a list, it will be applied to all plates in runs list. 
#' The list length must equal to the number of runs. Default uses all sample columns. 
#' SCs, NCs and IPCs are excluded.
#' @param exclude_targets A list of vectors of row names or numeric row indices 
#' representing the targets that should be excluded from quantifiability calculation for each run.
#' For example, one might want to exclude internal controls. If a vector is provided instead 
#' of a list, it will be applied to all plates in the run. The list length must equal to the number of runs. 
#' Default is \code{NULL}, which includes all targets in the AQ data matrix.
#'
#' @return A nested list containing quantifiability information. List elements 
#' include:
#' \item{run_quantifiability}{A list with the same names as \code{runs}. Each run 
#'   has a sublist that includes a \code{quant} data frame with targets 
#'   in rows and overall sample quantifibility as 
#'   well as any subgroups as columns, and \code{n_samples} vector which gives
#'   number of samples corresponding to each column in \code{quant}.}
#'   
#' \item{combined_quantifiability}{Combines quantifiability across runs. 
#'   This sublist also includes a \code{quant} data frame with targets 
#'   in rows and overall sample quantifibility as 
#'   well as any subgroups as columns, and \code{n_samples} vector which gives
#'   number of samples corresponding to each column in \code{quant}. If only a 
#'   single run is input, \code{combined_quantifiability} will be identical 
#'   to \code{run_quantifiability}.}
#'   
#' \item{summary_tables}{A table or tables with quantifiability summary statistics
#'   for all samples and by sample group, if \code{sampleGroupCovar} is defined.}
#'
#'
#' @export
#' 
quantifiability <- function(runs,
                            sampleGroupCovar='SAMPLE_MATRIX',
                            sample_subset=NULL,
                            exclude_targets=NULL){
  
  # validate input parameters
  if (missing(runs) || is.null(runs)) {
    stop("The 'runs' parameter is required and cannot be NULL")
  }
  if (length(runs) == 0) {
    stop("The runs list is empty - no data to process")
  }
  if (!is.list(runs) || !("RunSummary" %in% names(runs) || any(sapply(runs, function(x) is.list(x) && "RunSummary" %in% names(x))))) {
    stop("The 'runs' parameter must be a list or contain a RunSummary element")
  }
  # check if run data is not in a list and if so put into a list
  if('RunSummary' %in% names(runs)){
    runs <- list(runs)
    names(runs) <- 'Plate 01'
  } 
  # validate exclude_targets
  if (!is.null(exclude_targets)) {
    if (is.list(exclude_targets)) {
      if (length(exclude_targets) != length(runs)) {
        stop("exclude_targets list must be the same length as the number of runs")
      }
    } else {
      # Vector input - apply the same targets to all plates
      exclude_targets <- rep(list(exclude_targets), length(runs))
    }
  } 
  
  # validate sample_subset
  if (!is.null(sample_subset)) {
    if (is.list(sample_subset)) {
      if (length(sample_subset) != length(runs)) {
        stop("sample_subset list must be the same length as the number of runs")
      }
    } else {
      # Vector input - apply the same samples to all plates
      sample_subset <- rep(list(sample_subset), length(runs))
    }
  } 
  
  run_quant <- lapply(seq_along(runs), function(k) {
    x <- runs[[k]]
    
    tryCatch(
      {
        if (is.null(x$AQ) || is.null(x$AQ$targetAQ_param)) {
          stop("AQ data or targetAQ_param is missing")
        }
        # check if LLOQ and ULOQ are present
        if(!("LLOQ" %in% colnames(x$AQ$targetAQ_param) & "ULOQ" %in% colnames(x$AQ$targetAQ_param))){
          stop('LLOQ and ULOQ are not present in run data for one or more runs. Quantifiability cannot be calculated.')
        }
        if(!is.null(exclude_targets)){
          # exclude targets from various data structures 
          exclude_list <- exclude_targets[[k]]
          if (!is.null(exclude_list) && length(exclude_list) > 0) {
            existing_targets <- exclude_list[exclude_list %in% x$AQ$targetAQ_param$targetName]
            if (length(existing_targets) < length(exclude_list)) {
              warning("Some exclude_targets not found in run ", k, ": ", 
                      paste(setdiff(exclude_list, existing_targets), collapse=", "))
            }
            
            x$AQ$targetAQ_param <- x$AQ$targetAQ_param[!x$AQ$targetAQ_param$targetName %in% exclude_list, ]
            x$AQ$Data_AQ_aM <- x$AQ$Data_AQ_aM[!rownames(x$AQ$Data_AQ_aM) %in% exclude_list, ]
            x$AQ$Data_AQ <- x$AQ$Data_AQ[!rownames(x$AQ$Data_AQ) %in% exclude_list, ]
            x$AQ$withinDR <- x$AQ$withinDR[!names(x$AQ$withinDR) %in% exclude_list]
          }
        }
        if(!is.null(sample_subset)){
          # sample subset from various data structures 
          sample_list <- sample_subset[[k]]
          if (!is.null(sample_list) && length(sample_list) > 0) {
            existing_samples <- sample_list[sample_list %in% x$SampleNames]
            if (length(existing_samples) < length(sample_list)) {
              warning("Some sample_subset not found in run ", current_run_index, ": ", 
                      paste(setdiff(sample_list, existing_samples), collapse=", "))
            }
            
            x$SampleNames <- x$SampleNames[x$SampleNames %in% sample_list]
          }
        }
        
        # get the LLOQ and ULOQ
        DR <- x$AQ$targetAQ_param[,c("targetName","LLOQ","ULOQ")]
        # merge LOD
        DR <- merge(DR, data.frame(targetName=names(x$lod$LOD_aM), LOD_aM=x$lod$LOD_aM),
                    all.x=TRUE, all.y=FALSE)
        # replace LLOQ with LOD if the LOD > LLOQ
        DR$LLOQ[DR$LOD_aM > DR$LLOQ & !is.na(DR$LOD_aM) & !is.na(DR$LLOQ)] <- DR$LOD_aM[DR$LOD_aM > DR$LLOQ & !is.na(DR$LOD_aM) & !is.na(DR$LLOQ)]
        
        # calculate overall quantifiability
        intersect_samples <- intersect(colnames(x$AQ$Data_AQ_aM), x$SampleNames)
        AQ_quant <- x$AQ$Data_AQ_aM[,intersect_samples, drop=FALSE] 
        AQ_quant <- merge(AQ_quant, DR, by.x='row.names', by.y='targetName')
        rownames(AQ_quant) <- AQ_quant[,1]
        AQ_quant <- AQ_quant[,2:ncol(AQ_quant)]
        AQ_quant$overall <- apply(AQ_quant, 1, function(y){
          target_data <- as.numeric(y[1:(length(y) - 3)])
          LLOQ <- as.numeric(y[(length(y) - 2)])
          ULOQ <- as.numeric(y[(length(y) - 1)])
          target_quant <- sum(target_data > LLOQ & target_data < ULOQ & !is.na(target_data)) / length(target_data) * 100
          return(target_quant)
        })
        
        AQ_quant_output_columns <- c('overall')
        n_samples <- c(overall=length(intersect_samples)) 
        
        # calculate subgroup quantifiability
        if(!is.null(sampleGroupCovar)){
          if(sampleGroupCovar %in% colnames(x$samples)){
            subgroup_names <- unique(x$samples[x$samples$sampleType=='Sample', sampleGroupCovar])
            
            for(i in 1:length(subgroup_names)){
              # Filter samples to those in intersect_samples AND belonging to this subgroup
              samples_in_subgroup <- x$samples$sampleName[x$samples[,sampleGroupCovar] == subgroup_names[i]]
              subgroup_samples <- intersect(samples_in_subgroup, intersect_samples)
              subgroup_sample_data <- AQ_quant[,subgroup_samples, drop=FALSE]
              AQ_quant[,as.character(subgroup_names[i])] <- NA
              AQ_quant_output_columns <- c(AQ_quant_output_columns, as.character(subgroup_names[i]))
              n_samples <- c(n_samples, length(subgroup_samples))
              names(n_samples) <- AQ_quant_output_columns
              
              for(j in 1:nrow(AQ_quant)){
                AQ_quant[j, as.character(subgroup_names[i])] <- sum(subgroup_sample_data[j,] > AQ_quant$LLOQ[j] & subgroup_sample_data[j,] < AQ_quant$ULOQ[j] & !is.na(subgroup_sample_data[j,])) / length(subgroup_sample_data[j,]) * 100
              }
            }
          }
        }
        quant_output <- AQ_quant[,AQ_quant_output_columns]
        if(!is.data.frame(quant_output)){
          quant_output <- data.frame(overall=quant_output)
          rownames(quant_output) <- rownames(AQ_quant)
        }
        
        if (nrow(quant_output) == 0) {
          warning("Run ", k, " produced no quantifiability results. Returning NULL.")
          return(NULL)
        }
        
        target_quant_output <- list(quant=quant_output,
                                    n_samples=n_samples)
        return(target_quant_output)
      }, error = function(e) {
        warning("Error processing run: ", conditionMessage(e))
        return(NULL)
      }
    )
  })

  # assign run names from input runs list
  names(run_quant) <- names(runs)

  # filter out failed runs
  run_quant <- run_quant[!sapply(run_quant, is.null)]
  if (length(run_quant) == 0) {
    stop("All runs failed during quantifiability calculation")
  }
  
  # combine the quantifiability across runs
  if(length(run_quant) > 1){
    unique_quant_groups <- unique(unlist(lapply(run_quant, function(x) colnames(x$quant))))
    combined_quant <- vector(mode='list', length=length(unique_quant_groups))
    names(combined_quant) <- unique_quant_groups
    for(i in 1:length(combined_quant)){
      runs_with_group <- which(sapply(run_quant, function(x) unique_quant_groups[i] %in% colnames(x$quant)))
      
      quant_group_i <- lapply(run_quant[runs_with_group], function(x){
        quant_group_data <- data.frame(target=rownames(x$quant), 
                                       quant=x$quant[,unique_quant_groups[i]])
        quant_group_n <- x$n_samples[unique_quant_groups[i]]
        quant_group_data$n_quantifiable <- quant_group_data$quant * quant_group_n / 100
        return(list(quant_group_data=quant_group_data,
                    quant_group_n=quant_group_n))
      })
      
      n_quantifiable <- suppressWarnings(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "target", all = TRUE),
                                                lapply(quant_group_i, function(x) x$quant_group_data[,c('target', 'n_quantifiable')])))
      n_quantifiable$total_quantifiable <- rowSums(n_quantifiable[,2:ncol(n_quantifiable), drop=FALSE], na.rm=TRUE)
      total_n <- sum(sapply(quant_group_i, function(x) x$quant_group_n))
      n_quantifiable[,unique_quant_groups[i]] <- n_quantifiable$total_quantifiable / total_n * 100
      
      combined_quant[[unique_quant_groups[i]]]$quantifiability <- n_quantifiable[,c("target",unique_quant_groups[i])]
      combined_quant[[unique_quant_groups[i]]]$n_samples <- total_n
    }
    
    # combine the subgroup quantifiability into a single data frame
    combined_quant_matrix <- suppressWarnings(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "target", all = TRUE),
                                                     lapply(combined_quant, function(x) x$quantifiability)))
    target_names <- combined_quant_matrix$target
    combined_quant_matrix <- as.data.frame(combined_quant_matrix[,2:ncol(combined_quant_matrix)])
    colnames(combined_quant_matrix)[1] <- 'overall'
    rownames(combined_quant_matrix) <- target_names
    combined_quant_n_samples <- sapply(combined_quant, function(x) x$n_samples)
    
    if (nrow(combined_quant_matrix) == 0) {
      warning("Combined quantifiability is empty across all runs. Returning NULL.")
      return(NULL)
    }
    
    combined_quant_output <- list(quant=combined_quant_matrix,
                                  n_samples=combined_quant_n_samples)
  }
  
  if(length(runs) == 1){
    function_output <- list(run_quantifiability=run_quant, 
                            combined_quantifiability=run_quant[[1]])
  } else if(length(runs) > 1){
    function_output <- list(run_quantifiability=run_quant, 
                            combined_quantifiability=combined_quant_output)
  }
  
  # make summary tables
  sample_groups <- colnames(function_output$combined_quantifiability$quant)
  names(sample_groups) <- sample_groups
  
  summary_tables <- lapply(sample_groups, function(x) {
    plate_stats <- column_summary_stats(do.call(cbind, lapply(function_output$run_quantifiability, function(y) {
      if(x %in% colnames(y$quant)){
        y$quant[,x]
      } else {
        NULL
      }
    })))[,1:5]
    
    n_samples <- unlist(lapply(function_output$run_quantifiability, function(y) {
      if(x %in% names(y$n_samples)) {
        y$n_samples[x]
      } else {
        NULL
      }
    }))
    names(n_samples) <- rownames(plate_stats)
    
    n_pct_quantifiabile <- unlist(lapply(function_output$run_quantifiability, function(y){
      if(x %in% colnames(y$quant)){
        n_quant <- sum(y$quant[,x] > 50, na.rm=TRUE)
        tot <- length(y$quant[,x])
        pct_quant <- format(round(n_quant / tot * 100, 1), nsmall=1)
        output <- paste0(n_quant, ' / ', tot, ' (', pct_quant, '%)')
        return(output)
      }
    }))
    
    # if only one plate
    if(is.null(nrow(plate_stats))){
      plate_stats <- matrix(c(n_samples, plate_stats, n_pct_quantifiabile), nrow=1)
      rownames(plate_stats) <- names(function_output$run_quantifiability)
    } else {
      # if multiple plates
      plate_stats <- data.frame(n_samples, plate_stats, n_pct_quantifiabile)
    }
    
    if(length(runs) > 1){
      combined_quant_stats <- column_summary_stats(function_output$combined_quantifiability$quant[,x])[,1:5]
      combined_quant_n_samples <- function_output$combined_quantifiability$n_samples[x]
      combined_n_quant <- sum(function_output$combined_quantifiability$quant[,x] > 50, na.rm=TRUE)
      combined_tot <- length(function_output$combined_quantifiability$quant[,x])
      combined_pct_quant <- format(round(combined_n_quant / combined_tot * 100, 1), nsmall=1)
      combined_n_pct_quantifiabile <- paste0(combined_n_quant, ' / ', combined_tot, ' (', combined_pct_quant, '%)')
      combined_stats <- matrix(c(combined_quant_n_samples, combined_quant_stats, combined_n_pct_quantifiabile), nrow=1)
      colnames(combined_stats) <- colnames(plate_stats)
      plate_stats <- rbind(plate_stats, combined_stats)
      rownames(plate_stats)[nrow(plate_stats)] <- 'Overall'
    }
    
    colnames(plate_stats) <- c('# samples', 'mean', 'sd', 'median', 'min', 'max', '# of quantifiable targets (%)')
    return(plate_stats)
  })
  
  function_output <- c(function_output, list(summary_tables=summary_tables))
  
  return(function_output)
}
