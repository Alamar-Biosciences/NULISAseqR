#' Calculate Coefficient of Variation for Each Target for a Set of AQ Runs
#'
#' Calculates intra-plate coefficient of variation (CV) 
#' for each target for a single AQ run or intra- and inter-plate CV for
#' a set of AQ runs. By default the function will output CVs for sample types 
#' SC (AQSC) and IPC (CAL). In the future, we might add the ability to also output 
#' CVs for samples that are grouped using a provided sample name pattern. 
#' Function input is the output of loadNULISAseq or a list
#' of loadNULISAseq outputs. CVs are calculated directly on the AQ data in 
#' aM units. Values outside of the dynamic range are excluded from CV calculation.
#'
#' @param runs A named list of a single run data output from \code{laodNULISAseq()} 
#' function or a named list of these outputs for multiple runs.
#' @param exclude_outside_DR Logical. Default is TRUE. Should values outside the 
#' dynamic range (i.e., values below LLOQ or above ULOQ) be excluded 
#' before calculating CV?
#'
#'
#' @return A nested list containing CV information. List elements 
#' include:
#' \item{run_intraCV}{A list with the same names as \code{runs}. Each run 
#'   has a sublist that includes a \code{intraCV} data frame with targets 
#'   in rows and samples as columns, and \code{n_samples} vector which gives
#'   number of sample replicates corresponding to each column in \code{intraCV}.}
#'   
#' \item{average_intraCV}{Summarizes intraCV across runs by taking the average 
#'   of the run-specific CVs. 
#'   This sublist also includes a \code{intraCV} data frame with targets 
#'   in rows and samples as columns, and \code{n_samples} vector which gives
#'   average number of sample replicates per run corresponding to each column 
#'   in \code{intraCV}. If only a 
#'   single run is input, \code{combined_intraCV} will be identical 
#'   to \code{run_intraCV}.}
#' 
#' \item{interCV}{If multiple runs are input, the \code{interCV} list is output.
#'   This gives the inter-plate CV calculated on the pooled replicates. Note this 
#'   is not a component inter-plate CV but a total CV across all runs. 
#'   This list includes a \code{interCV} data frame with targets 
#'   in rows and samples as columns, and \code{n_samples} vector which gives
#'   number of sample replicates corresponding to each column in \code{interCV}.}
#'
#'
#' @export
#' 
CV_AQ <- function(runs,
                  exclude_outside_DR=TRUE){
  
  # check if run data is not in a list and if so put into a list
  if('RunSummary' %in% names(runs)){
    runs <- list(runs)
    names(runs) <- 'Plate 01'
  } 
  
  # obtain the SC/AQSC and IPC/CAL AQ data
  SC_IPC_data <- lapply(runs, function(x) {
    
    # check if LLOQ and ULOQ are present
    if(!("LLOQ" %in% colnames(x$AQ$targetAQ_param) & "ULOQ" %in% colnames(x$AQ$targetAQ_param))){
      warning('LLOQ and ULOQ are not present in run data for one or more runs. Values outside dynamic range were not excluded for these run(s); all values were used to calculate AQ CVs.')
      x$AQ$targetAQ_param$LLOQ <- rep(0, nrow(x$AQ$targetAQ_param))
      x$AQ$targetAQ_param$ULOQ <- rep(Inf, nrow(x$AQ$targetAQ_param))
    }
    
    if(exclude_outside_DR==FALSE){
      warning('Values outside dynamic range were not excluded; all values were used to calculate AQ CVs.')
      x$AQ$targetAQ_param$LLOQ <- rep(0, nrow(x$AQ$targetAQ_param))
      x$AQ$targetAQ_param$ULOQ <- rep(Inf, nrow(x$AQ$targetAQ_param))
    }
    
    # get DR info to set values to missing 
    # get the LLOQ and ULOQ
    DR <- x$AQ$targetAQ_param[,c("targetName","LLOQ","ULOQ")]
    # merge LOD
    DR <- merge(DR, data.frame(targetName=names(x$lod$LOD_aM), LOD_aM=x$lod$LOD_aM),
                all.x=TRUE, all.y=FALSE)
    # replace LLOQ with LOD if the LOD > LLOQ
    DR$LLOQ[DR$LOD_aM > DR$LLOQ & !is.na(DR$LOD_aM) & !is.na(DR$LLOQ)] <- DR$LOD_aM[DR$LOD_aM > DR$LLOQ & !is.na(DR$LOD_aM) & !is.na(DR$LLOQ)]
    
    # extract SC and IPC data
    SC <- x$AQ$Data_AQ_aM[,x$samples$sampleType=='SC']
    IPC <- x$AQ$Data_AQ_aM[,x$samples$sampleType=='IPC']
    
    # set values outside DR to NA
    SC[SC < DR$LLOQ | SC > DR$ULOQ] <- NA
    IPC[IPC < DR$LLOQ | IPC > DR$ULOQ] <- NA
    
    return(list(SC=SC,
                IPC=IPC))
  })
  
  # calculate intra CVs
  intraCVs <- lapply(SC_IPC_data, function(x){
    
    run_intraCV <- lapply(seq_along(x), function(i){
      intraCV <- apply(x[[i]], 1, sd, na.rm=TRUE) / rowMeans(x[[i]], na.rm=TRUE) * 100
      output <- data.frame(target=names(intraCV), intraCV=intraCV)
      colnames(output)[2] <- names(x)[i]
      return(output)
    })
    
    # make one data frame for each run
    run_intraCV <- suppressWarnings(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "target", all = TRUE), x=run_intraCV))
    rownames(run_intraCV) <- run_intraCV$target
    run_intraCV <- run_intraCV[,2:ncol(run_intraCV)]
    
    # get number of replicates for each run and sample
    n_samples <- unlist(lapply(x, ncol))
    
    return(list(intraCV=run_intraCV,
                n_samples=n_samples))
  })
  
  # calculate average intraCV
  if(length(runs) == 1) average_intraCV <- intraCVs[[1]]
  if(length(runs) > 1) {
    
    SC_intraCVs <- lapply(intraCVs, function(x) {
      data.frame(target=rownames(x$intraCV), run_SC_CV=x$intraCV[,'SC'])
    })
    SC_intraCVs <- suppressWarnings(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "target", all = TRUE), x=SC_intraCVs))
    SC_intraCVs$SC <- rowMeans(SC_intraCVs[,2:ncol(SC_intraCVs)], na.rm=TRUE)
    
    IPC_intraCVs <- lapply(intraCVs, function(x) {
      data.frame(target=rownames(x$intraCV), run_IPC_CV=x$intraCV[,'IPC'])
    })
    IPC_intraCVs <- suppressWarnings(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "target", all = TRUE), x=IPC_intraCVs))
    IPC_intraCVs$IPC <- rowMeans(IPC_intraCVs[,2:ncol(IPC_intraCVs)], na.rm=TRUE)
    
    avg_intraCV <- merge(SC_intraCVs[,c("target","SC")], 
                         IPC_intraCVs[,c("target","IPC")], all=TRUE)
    rownames(avg_intraCV) <- avg_intraCV$target
    avg_intraCV <- avg_intraCV[,2:ncol(avg_intraCV)]
    
    avg_n_samples <- rowMeans(sapply(intraCVs, function(x) x$n_samples), na.rm=TRUE)
    
    average_intraCV <- list(intraCV=avg_intraCV,
                            n_samples=avg_n_samples)
  }
  
  # calculate inter CV
  if(length(runs) > 1){
    SC_data <- lapply(SC_IPC_data, function(x) {
      dat <- as.data.frame(x$SC)
      dat$target <- rownames(dat)
      return(dat)
    })
    SC_data <- suppressWarnings(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "target", all = FALSE), x=SC_data))
    SC_interCV <- data.frame(target=SC_data[,1],
                             SC=apply(SC_data[,2:ncol(SC_data)], 1, function(x) {
                               sd(x, na.rm=TRUE) / mean(x, na.rm=TRUE) * 100
                             }))
    
    IPC_data <- lapply(SC_IPC_data, function(x) {
      dat <- as.data.frame(x$IPC)
      dat$target <- rownames(dat)
      return(dat)
    })
    IPC_data <- suppressWarnings(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "target", all = FALSE), x=IPC_data))
    IPC_interCV <- data.frame(target=IPC_data[,1],
                              IPC=apply(IPC_data[,2:ncol(IPC_data)], 1, function(x) {
                                sd(x, na.rm=TRUE) / mean(x, na.rm=TRUE) * 100
                              }))
    
    interCV <- merge(SC_interCV, IPC_interCV, by='target')
    rownames(interCV) <- interCV[,1]
    interCV <- interCV[,2:ncol(interCV)]
    
    interCV_n_samples <- c(SC=ncol(SC_data) - 1,
                           IPC=ncol(IPC_data) - 1)
    
    interCV <- list(interCV=interCV,
                    n_samples=interCV_n_samples)
    
  }
  
  
  CV_AQ_output <- list(run_intraCV=intraCVs,
                       average_intraCV=average_intraCV)
  if(length(runs) > 1){
    CV_AQ_output <- list(run_intraCV=intraCVs,
                         average_intraCV=average_intraCV,
                         interCV=interCV)
  }
  
  return(CV_AQ_output)
}


