# helper function to summarize special well types:
# IPCs
# NCs
# SCs
# Bridge samples
typeSummary <- function(well_type, plate_data, total_plate_reads){    
  well_type_data <- plate_data$Data[,well_type]
  # total counts and percent
  well_type_total <- sum(well_type_data, na.rm=TRUE)
  well_type_total_perc <- format(round(well_type_total/total_plate_reads*100, 1), nsmall=1)
  well_type_totals <- colSums(well_type_data, na.rm=TRUE)
  well_type_totals_perc <- format(round(well_type_totals/total_plate_reads*100, 1), nsmall=1)
  well_type_medians <- apply(well_type_data, 1, median, na.rm=TRUE)
  well_type_medians_total <- sum(well_type_medians, na.rm=TRUE)
  well_type_medians_total_perc <- format(round(well_type_medians_total/total_plate_reads*100, 1), nsmall=1)
  total_well_type_counts <- c(well_type_total, well_type_totals, well_type_medians_total)
  total_well_type_perc <- c(well_type_total_perc, well_type_totals_perc, well_type_medians_total_perc)
  total_well_type_count_perc <- paste0(format(total_well_type_counts, big.mark=","), ' (',
                                       total_well_type_perc, '%)')
  # mean count per target
  well_type_total_mean <- mean(well_type_data, na.rm=TRUE)
  well_type_target_mean <- colMeans(well_type_data, na.rm=TRUE)
  well_type_medians_mean <- mean(well_type_medians, na.rm=TRUE)
  well_type_means_per_target <- c(well_type_total_mean, well_type_target_mean, well_type_medians_mean)
  well_type_means_per_target <- format(round(well_type_means_per_target, 1), big.mark=",", nsmall=1)
  # missing n %
  well_type_total_missing <- sum(is.na(well_type_data)) + sum(well_type_data==0, na.rm=TRUE)
  well_type_total_missing_perc <- well_type_total_missing/(nrow(well_type_data)*ncol(well_type_data))*100
  well_type_missing <- apply(well_type_data, 2, function(x){
    sum(is.na(x)) + sum(x==0, na.rm=TRUE)
  })
  well_type_missing_perc <- well_type_missing/nrow(well_type_data)*100
  well_type_medians_missing <- sum(is.na(well_type_medians)) + sum(well_type_medians==0, na.rm=TRUE)
  well_type_medians_missing_perc <- well_type_medians_missing/length(well_type_medians)*100
  well_type_missing_totals <- c(well_type_total_missing, well_type_missing, well_type_medians_missing)
  well_type_missing_percents <- c(well_type_total_missing_perc,
                                  well_type_missing_perc,
                                  well_type_medians_missing_perc)
  well_type_missing_percents <- format(round(well_type_missing_percents, 1), nsmall=1, big.mark=",")
  well_type_missing <- paste0(well_type_missing_totals, ' (',
                              well_type_missing_percents, '%)')
  # calculate CV
  well_type_means <- rowMeans(well_type_data, na.rm=TRUE)
  well_type_sds <- apply(well_type_data, 1, sd, na.rm=TRUE)
  well_type_cvs <- well_type_sds/well_type_means*100
  well_type_cv_mean <- c(paste0(format(round(mean(well_type_cvs, na.rm=TRUE), 1), nsmall=1), '%'),
                         rep('', length(well_type)),
                         '')
  well_type_cv_median <- c(paste0(format(round(median(well_type_cvs, na.rm=TRUE), 1), nsmall=1), '%'),
                           rep('', length(well_type)),
                           '')
  # create table
  well_type_table <- cbind(total_well_type_count_perc,
                           well_type_means_per_target,
                           well_type_missing,
                           well_type_cv_mean,
                           well_type_cv_median)
  colnames(well_type_table) <- c('Total count (%)',
                                 'Mean count per target',
                                 'Zeros n (%)',
                                 'CV% mean',
                                 'CV% median')
  rownames(well_type_table) <- c('Total',
                                 colnames(well_type_data),
                                 'Median')
  return(well_type_table)
}

#' NULISAseq Plate Summary
#'
#' Summarizes plate reads. Input is the output of readNULISAseq.R.
#'
#' @param plate_data The list output from readNULISAseq.R.
#' @param ICs Optional. Vector of either the row indices (numeric) or the 
#' row names (character string) of the internal controls.
#' @param IPCs Optional. Vector of either the column indices (numeric) or the 
#' column names (character string) of the inter-plate controls.
#' @param NCs Optional. Vector of either the column indices (numeric) or the 
#' column names (character string) of the negative controls.
#'
#' @return A table (character matrix) or list of tables (character matrices).
#' @param readsTable A table that summarizes read counts and percent of 
#' total read counts by count type. Also summarizes missing data.
#' @param IC_table If ICs are defined, outputs another table with number (%)
#' of samples that have missing or zero IC counts, total and percent reads (out of all
#' parseable match reads), mean, SD, and %CV reads across all samples.
#' @param IPC_table If IPCs are defined, outputs another table with total (%)
#' IPC reads, mean IPC count per target, number (%) of targets with missing or zero
#' IPC reads, and the mean and median CV% across targets. Data is shown for 
#' total IPC wells / targets, each IPC well individually, and the targetwise
#' medians across the IPC wells (potentially used for interplate normalization).
#' @param NC_table If NCs are defined, outputs another table with total (%)
#' NC reads, mean NC count per target, maximum NC target count, 
#' number (%) of targets with missing or zero
#' NC reads. Data is shown for 
#' total NC wells / targets, each NC well individually, and the targetwise
#' means across the NC wells (potentially used for LOD calculation).
#'
#' @examples
#' plate1summary <- plateSummary(plate1)
#' plate1summary
#'
#' @export
#' 
plateSummary <- function(plate_data, ICs=NULL, IPCs=NULL, NCs=NULL, SCs=NULL, Bridges=NULL){
  ##############################
  # main output
  ##############################
  # get reads
  parseable <- as.numeric(plate_data$RunSummary$Parseable)
  parseable_match <- as.numeric(plate_data$RunSummary$ParseableMatch)
  parseable_nonmatch <- parseable - parseable_match
  unparseable <- as.numeric(plate_data$RunSummary$Unparseable)
  balancers <- sum(as.numeric(plate_data$RunSummary$Balancers))
  totalReads <- as.numeric(plate_data$RunSummary$TotalReads) + balancers
  # calculate total samples*targets and missing
  total_samples <- ncol(plate_data$Data)
  total_targets <- nrow(plate_data$Data)
  samples_targets <- nrow(plate_data$Data)*ncol(plate_data$Data)
  missing_samples_targets <- sum(is.na(plate_data$Data)) + sum(plate_data$Data==0, na.rm=TRUE)
  # calculate percent
  parseable_perc <- format(round(parseable/totalReads*100, 1), nsmall=1)
  parseable_match_perc <- format(round(parseable_match/totalReads*100, 1), nsmall=1)
  parseable_nonmatch_perc <- format(round(parseable_nonmatch/totalReads*100, 1), nsmall=1)
  unparseable_perc <- format(round(unparseable/totalReads*100, 1), nsmall=1)
  balancers_perc <- format(round(balancers/totalReads*100, 1), nsmall=1)
  missing_perc <- paste0(' (', 
                         format(round(missing_samples_targets/samples_targets*100, 1), nsmall=1),
                         '%)')
  # paste reads and percents together
  reads <- c(totalReads, parseable, 
             parseable_match, parseable_nonmatch,
             unparseable, 
             total_samples, total_targets,
             samples_targets, missing_samples_targets)
  reads <- formatC(reads, format='d', big.mark=',')
  percents <- paste0(' (', 
                     c(parseable_perc, parseable_match_perc, 
                       parseable_nonmatch_perc,
                       unparseable_perc), '%)')
  reads_percents <- paste0(reads, c('', percents, '', '', '', missing_perc))
  reads_percents <- matrix(reads_percents, nrow=9)
  row.names(reads_percents) <- c('Total reads', 'Parseable',
                                 'Parseable Match', 'Parseable Non-match',
                                 'Unparseable',
                                 'Total samples',
                                 'Total targets',
                                 'Total samples * targets',
                                 'Zero values')
  output <- reads_percents
  # check if total reads in Data equals parseable match reads in Run Summary
  total_plate_reads <- sum(plate_data$Data, na.rm=TRUE)
  if (total_plate_reads != parseable_match){
    warning('Total plate reads in Data does not equal total parseable match reads in Run Summary.')
  }
  ##############################
  # summarize controls
  ##############################
  # ICs
  ############
  ICinds <- if(!is.null(ICs)) ICs else which(plate_data$targets$targetType == "Control")
  if (length(ICinds) > 0){
    IC_data <- plate_data$Data[ICinds,, drop=F]
    IC_totals <- rowSums(IC_data, na.rm=TRUE)
    IC_percents <- format(round(IC_totals/total_plate_reads*100, 1), nsmall=1)
    IC_total_percent <- paste0(format(IC_totals, big.mark=","), ' (', IC_percents, '%)')
    IC_missing_n <- apply(IC_data, 1, function(x) sum(is.na(x)) + sum(x==0, na.rm=TRUE))
    IC_missing_perc <- format(round(IC_missing_n/ncol(IC_data)*100, 1), nsmall=1)
    IC_missing_n_perc <- paste0(format(IC_missing_n, big.mark=","), ' (', IC_missing_perc, '%)')
    IC_means <- format(round(rowMeans(IC_data, na.rm=TRUE), 1), nsmall=1, big.mark=",")
    IC_sd <- format(round(apply(IC_data, 1, sd, na.rm=TRUE), 1), nsmall=1, big.mark=",")
    IC_cv <- paste0(format(round(apply(IC_data, 1, sd, na.rm=TRUE)/rowMeans(IC_data, na.rm=TRUE)*100, 1), nsmall=1),
                    '%')
    IC_table <- cbind(IC_missing_n_perc, IC_total_percent, IC_means, IC_sd, IC_cv)
    colnames(IC_table) <- c('Zeros n (%)', 'Total reads (%)', 'Mean', 'SD', 'CV%')
    rownames(IC_table) <- rownames(IC_data)
    output <- list(readsTable=output,
                   IC_table=IC_table)
  }
  ############
  # IPCs, SCs, Bridges
  ############
  IPCinds <- if(!is.null(IPCs)) IPCs else which(plate_data$samples$sampleType == "IPC")
  if(length(IPCinds) > 0 ){
    IPC_table <- typeSummary(IPCinds, plate_data, total_plate_reads)
    output <- if (is.list(output)==TRUE) c(output, list(IPC_table=IPC_table)) else list(readsTable=output, IPC_table=IPC_table)
  }
  SCinds <- if(!is.null(SCs)) SCs else which(plate_data$samples$sampleType == "SC")
  if(length(SCinds) > 0){
    SC_table <- typeSummary(SCinds, plate_data, total_plate_reads)
    output <- if (is.list(output)==TRUE) c(output, list(SC_table=SC_table)) else list(readsTable=output, SC_table=SC_table)
  }
  Bridgeinds <- if(!is.null(Bridges)) Bridges else which(plate_data$samples$sampleType == "Bridge")
  if(length(Bridgeinds)>0){
    Bridge_table <- typeSummary(Bridges, plate_data, total_plate_reads)
    output <- if (is.list(output)==TRUE) c(output, list(Bridge_table=Bridge_table)) else list(readsTable=output, Bridge_table=Bridge_table)
  }
  
  ############
  # NCs
  ############
  NCinds <- if(!is.null(NCs)) NCs else which(plate_data$samples$sampleType == "NC")
  if (length(NCinds)> 0){
    NC_data <- plate_data$Data[,NCinds]
    # total counts and percent
    NC_total <- sum(NC_data, na.rm=TRUE)
    NC_total_perc <- format(round(NC_total/total_plate_reads*100, 1), nsmall=1)
    NC_totals <- colSums(NC_data, na.rm=TRUE)
    NC_totals_perc <- format(round(NC_totals/total_plate_reads*100, 1), nsmall=1)
    NC_means <- apply(NC_data, 1, mean, na.rm=TRUE)
    NC_means_total <- sum(NC_means, na.rm=TRUE)
    NC_means_total_perc <- format(round(NC_means_total/total_plate_reads*100, 1), nsmall=1)
    total_NC_counts <- c(NC_total, NC_totals, NC_means_total)
    total_NC_perc <- c(NC_total_perc, NC_totals_perc, NC_means_total_perc)
    total_NC_count_perc <- paste0(format(total_NC_counts, big.mark=",", nsmall=1), ' (',
                                  total_NC_perc, '%)')
    # mean count per target
    NC_total_mean <- mean(NC_data, na.rm=TRUE)
    NC_target_mean <- colMeans(NC_data, na.rm=TRUE)
    NC_means_mean <- mean(NC_means, na.rm=TRUE)
    NC_means_per_target <- c(NC_total_mean, NC_target_mean, NC_means_mean)
    NC_means_per_target <- format(round(NC_means_per_target, 1), nsmall=1)
    # max target counts
    NC_max <- apply(NC_data, 2, max, na.rm=TRUE)
    NC_mean_max <- format(round(max(NC_means, na.rm=TRUE), 1), big.mark=",", nsmall=1)
    NC_maxs <- c('', format(NC_max, big.mark=","), NC_mean_max)
    # missing n %
    NC_total_missing <- sum(is.na(NC_data)) + sum(NC_data==0, na.rm=TRUE)
    NC_total_missing_perc <- NC_total_missing/(nrow(NC_data)*ncol(NC_data))*100
    NC_missing <- apply(NC_data, 2, function(x){
      sum(is.na(x)) + sum(x==0, na.rm=TRUE)
    })
    NC_missing_perc <- NC_missing/nrow(NC_data)*100
    NC_means_missing <- sum(is.na(NC_means)) + sum(NC_means==0, na.rm=TRUE)
    NC_means_missing_perc <- NC_means_missing/length(NC_means)*100
    NC_missing_totals <- c(NC_total_missing, NC_missing, NC_means_missing)
    NC_missing_percents <- c(NC_total_missing_perc,
                             NC_missing_perc,
                             NC_means_missing_perc)
    NC_missing_percents <- format(round(NC_missing_percents, 1), nsmall=1)
    NC_missing <- paste0(format(NC_missing_totals, big.mark=","), ' (',
                         NC_missing_percents, '%)')
    # create table
    NC_table <- cbind(total_NC_count_perc,
                      NC_means_per_target,
                      NC_maxs,
                      NC_missing)
    colnames(NC_table) <- c('Total count (%)',
                            'Mean count per target',
                            'Max target count',
                            'Missing n (%)')
    rownames(NC_table) <- c('NC total',
                            colnames(NC_data),
                            'NC means')
    if (is.list(output)==TRUE){
      output <- c(output, list(NC_table=NC_table))
    } else {
      output <- list(readsTable=output,
                     NC_table=NC_table)
    }
  }
  
  ##############################
  # return output
  ##############################
  return(output)
}
