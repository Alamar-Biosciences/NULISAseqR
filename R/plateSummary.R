typeSummary <- function(IPCs, plate_data, total_plate_reads){    
  IPC_data <- plate_data$Data[,IPCs]
  # total counts and percent
  IPC_total <- sum(IPC_data, na.rm=TRUE)
  IPC_total_perc <- format(round(IPC_total/total_plate_reads*100, 1), nsmall=1)
  IPC_totals <- colSums(IPC_data, na.rm=TRUE)
  IPC_totals_perc <- format(round(IPC_totals/total_plate_reads*100, 1), nsmall=1)
  IPC_medians <- apply(IPC_data, 1, median, na.rm=TRUE)
  IPC_medians_total <- sum(IPC_medians, na.rm=TRUE)
  IPC_medians_total_perc <- format(round(IPC_medians_total/total_plate_reads*100, 1), nsmall=1)
  total_IPC_counts <- c(IPC_total, IPC_totals, IPC_medians_total)
  total_IPC_perc <- c(IPC_total_perc, IPC_totals_perc, IPC_medians_total_perc)
  total_IPC_count_perc <- paste0(format(total_IPC_counts, big.mark=","), ' (',
                                 total_IPC_perc, '%)')
  # mean count per target
  IPC_total_mean <- mean(IPC_data, na.rm=TRUE)
  IPC_target_mean <- colMeans(IPC_data, na.rm=TRUE)
  IPC_medians_mean <- mean(IPC_medians, na.rm=TRUE)
  IPC_means_per_target <- c(IPC_total_mean, IPC_target_mean, IPC_medians_mean)
  IPC_means_per_target <- format(round(IPC_means_per_target, 1), big.mark=",", nsmall=1)
  # missing n %
  IPC_total_missing <- sum(is.na(IPC_data)) + sum(IPC_data==0, na.rm=TRUE)
  IPC_total_missing_perc <- IPC_total_missing/(nrow(IPC_data)*ncol(IPC_data))*100
  IPC_missing <- apply(IPC_data, 2, function(x){
    sum(is.na(x)) + sum(x==0, na.rm=TRUE)
  })
  IPC_missing_perc <- IPC_missing/nrow(IPC_data)*100
  IPC_medians_missing <- sum(is.na(IPC_medians)) + sum(IPC_medians==0, na.rm=TRUE)
  IPC_medians_missing_perc <- IPC_medians_missing/length(IPC_medians)*100
  IPC_missing_totals <- c(IPC_total_missing, IPC_missing, IPC_medians_missing)
  IPC_missing_percents <- c(IPC_total_missing_perc,
                            IPC_missing_perc,
                            IPC_medians_missing_perc)
  IPC_missing_percents <- format(round(IPC_missing_percents, 1), nsmall=1, big.mark=",")
  IPC_missing <- paste0(IPC_missing_totals, ' (',
                        IPC_missing_percents, '%)')
  # calculate CV
  IPC_means <- rowMeans(IPC_data, na.rm=TRUE)
  IPC_sds <- apply(IPC_data, 1, sd, na.rm=TRUE)
  IPC_cvs <- IPC_sds/IPC_means*100
  IPC_cv_mean <- c(paste0(format(round(mean(IPC_cvs, na.rm=TRUE), 1), nsmall=1), '%'),
                   rep('', length(IPCs)),
                   '')
  IPC_cv_median <- c(paste0(format(round(median(IPC_cvs, na.rm=TRUE), 1), nsmall=1), '%'),
                     rep('', length(IPCs)),
                     '')
  # create table
  IPC_table <- cbind(total_IPC_count_perc,
                     IPC_means_per_target,
                     IPC_missing,
                     IPC_cv_mean,
                     IPC_cv_median)
  colnames(IPC_table) <- c('Total count (%)',
                           'Mean count per target',
                           'Missing n (%)',
                           'CV% mean',
                           'CV% median')
  rownames(IPC_table) <- c('Total',
                           colnames(IPC_data),
                           'Median')
  return(IPC_table)
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
                                 'Missing/zero values')
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
  if (!is.null(ICs)){
    IC_data <- plate_data$Data[ICs,, drop=F]
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
    colnames(IC_table) <- c('Missing n (%)', 'Total reads (%)', 'Mean', 'SD', 'CV%')
    rownames(IC_table) <- rownames(IC_data)
    output <- list(readsTable=output,
                   IC_table=IC_table)
  }
  ############
  # IPCs, SCs, Bridges
  ############
  if(!is.null(IPCs)){
    IPC_table <- typeSummary(IPCs, plate_data, total_plate_reads)
    output <- if (is.list(output)==TRUE) c(output, list(IPC_table=IPC_table)) else list(readsTable=output, IPC_table=IPC_table)
  }
  if(!is.null(SCs)){
    SC_table <- typeSummary(SCs, plate_data, total_plate_reads)
    output <- if (is.list(output)==TRUE) c(output, list(SC_table=SC_table)) else list(readsTable=output, SC_table=SC_table)
  }
  if(!is.null(Bridges)){
    Bridge_table <- typeSummary(Bridges, plate_data, total_plate_reads)
    output <- if (is.list(output)==TRUE) c(output, list(Bridge_table=Bridge_table)) else list(readsTable=output, Bridge_table=Bridge_table)
  }

  ############
  # NCs
  ############
  if (!is.null(NCs)){
    NC_data <- plate_data$Data[,NCs]
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
