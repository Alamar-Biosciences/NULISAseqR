##' Process loadNULISAseq output for downstream analysis
#'
#' This function takes the output from `loadNULISAseq` and performs additional processing and cleaning to prepare the data for downstream analysis (e.g., NAS). It standardizes sample and plate identifiers, cleans covariate names, and generates additional normalized and transformed data matrices.
#'
#' @param data A list, the output from `loadNULISAseq` (from the NULISAseqR package), representing a single NULISAseq plate/run.
#'
#' @return A processed list with the following elements:
#'   \item{samples}{A cleaned data frame of sample metadata, with standardized column names and identifiers.}
#'   \item{covariateNames}{A named character vector mapping cleaned covariate names to their original names.}
#'   \item{Data_raw}{Matrix of raw counts, with cleaned sample names.}
#'   \item{Data_rawlog2}{Matrix of log2-transformed raw counts, with cleaned sample names.}
#'   \item{Data_IC}{Matrix of IC-normalized data, with cleaned sample names.}
#'   \item{Data_IClog2}{Matrix of log2-transformed IC-normalized data, with cleaned sample names.}
#'   \item{Data_Reverse}{Matrix of IC-normalized data for reverse curve targets, with cleaned sample names.}
#'   \item{Data_Reverselog2}{Matrix of log2-transformed IC-normalized data for reverse curve targets, with cleaned sample names.}
#'   \item{plateID}{Character, standardized plate ID.}
#'   \item{xmlFile}{Character, xml file name.}
#'   \item{ExecutionDetails}{List of run/plate metadata.}
#'   \item{numericCovariates}{Character vector of numeric covariate names.}
#'   \item{NPQ}{Matrix of NPQ values, with cleaned sample names.}
#'   \item{normed}{List of normalized data.}
#'   \item{normed_untransformedReverse}{List of normalized data for reverse curve targets.}
#'   \item{targets}{Data frame of target information.}
#'   \item{qcSample}{Data frame of sample-level QC metrics.}
#'   \item{qcPlate}{Data frame of plate-level QC metrics.}
#'   \item{aqParams}{Data frame of AQ target parameters (if present).}
#'   \item{detectability}{List of target detectability for all samples and by sample group/matrix.}
#'   \item{quantifiability}{List of target quantifiability for all samples and by sample group/matrix (if AQ data present).}
#'   \item{...}{Other elements from the original input, as needed.}
#'
#' @details
#' - Removes unnecessary columns and standardizes identifiers for compatibility with downstream tools.
#' - Cleans covariate (column) names to lowerCamelCase.
#' - Adds or updates normalized and log2-transformed data matrices.
#' - Used internally as a pre-processing step before merging or further analysis.
#'
#' @examples
#' # processed <- process_loadNULISAseq(loadNULISAseq_output)
#' @keywords internal
process_loadNULISAseq <- function(data) {
  processed_data <- data
  
  processed_data$samples <- processed_data$samples %>%
    select(!any_of(c("sampleName.1", "plateGUID")))
  
  logger::log_info("Processing samples dataframe")
  processed_data$samples <- process_samples_df(
    df = processed_data$samples,
    numericCovariates = processed_data$numericCovariates
  )
  
  logger::log_info("Clean columns for samples dataframe")
  defaultNames <- colnames(processed_data$samples)
  cleanNames <- clean_covariate_names(defaultNames, case = "lower_camel")
  names(defaultNames) <- cleanNames
  processed_data$covariateNames <- defaultNames

  # Preserve original sampleName from XML before any transformations
  processed_data$samples <- processed_data$samples %>%
    mutate(sampleName_original = sampleName)

  names_df <- processed_data$samples %>%
    filter(sampleType != "Sample") %>%
    select(sampleName, sampleID)

  processed_data$samples <- processed_data$samples %>%
    mutate(sampleName = dplyr::if_else(sampleType != "Sample",sampleID, sampleName))
  
  # save raw/normalized data as matrices 
  # raw counts 
  processed_data$NPQ <- processed_data$NPQ %>%
    rename_cols(., names_df = names_df)
  
  processed_data$Data_raw <- processed_data$Data %>%
    rename_cols(., names_df = names_df)
  # raw logged
  processed_data$Data_rawlog2 <- log2(processed_data$Data + 1) %>%
    rename_cols(., names_df = names_df)
  # IC norm
  processed_data$Data_IC <- processed_data$normed$interNormData[[1]] %>%
    rename_cols(., names_df = names_df)
  
  processed_data$Data_IClog2 <- processed_data$normed$log2_interNormData[[1]] %>%
    rename_cols(., names_df = names_df)
  # Create the Data_Reverse and Data_Reverselog2 with Reverse Curve Targets only
  processed_data$Data_Reverse <- processed_data$normed_untransformedReverse$interNormData[[1]] %>%
    rename_cols(., names_df = names_df)
  
  processed_data$Data_Reverselog2 <- processed_data$normed_untransformedReverse$log2_interNormData[[1]] %>%
    rename_cols(., names_df = names_df)
  
  if("AQ" %in% names(processed_data)){
    # AQ values can be available in aM and pg/mL units
    # In this event we will use the aM values
    processed_data$Data_AQ <- processed_data$AQ$Data_AQ_aM %>%
      rename_cols(., names_df = names_df)
    
    processed_data$Data_AQlog2 <- processed_data$Data_AQ %>%
      apply(., 2, function(x) log2(x))
    processed_data$AQ_unit <- "aM" 
    
    if("withinDR" %in% names(processed_data$AQ)){
      processed_data$withinDR <- processed_data$AQ$withinDR %>%
        rename_cols(., names_df = names_df)
    }
    if(!is.null(processed_data$AQ[["targetAQ_param"]])){
      if("MW_kDa" %in% names(processed_data$AQ$targetAQ_param)){
        processed_data$aqParams <- processed_data$AQ$targetAQ_param %>%
          tibble::as_tibble() %>%
          select(!any_of(c("Encrypted", "Decrypted")))
        
        logger::log_info("aqParams found, creating pg/mL matrix")
        processed_data$Data_AQ_pgmL <- processed_data[["Data_AQ"]] %>%
          tibble::as_tibble(rownames = "targetName") %>%
          left_join(
            processed_data$aqParams %>%
              select("targetName", "MW_kDa"),
            by = "targetName") %>%
          tidyr::pivot_longer(
            cols = !c(targetName, MW_kDa),
            values_to = "raw",
            names_to = "sampleName") %>%
          mutate(
            raw = dplyr::if_else(is.na(raw),
                          raw, 
                          NULISAseqAQ::unit_convert_am_conc(raw, MW_kDa))) %>%
          select(-MW_kDa) %>%
          tidyr::pivot_wider(
            id_cols = "targetName",
            names_from = "sampleName",
            values_from = "raw") %>%
          column_to_rownames("targetName") %>%
          as.matrix()
        
        logger::log_info("Data_AQlog2_pg/mL matrix created")
        
        processed_data$Data_AQlog2_pgmL <- log2(processed_data$Data_AQ_pgmL)
        processed_data$AQ_unit <- "pg/mL" 
      } else{
        logger::log_info("Molecular weights in kDa not available, only aM units available")
      }
    }
  }
  
  processed_data$NC <- processed_data$samples %>%
    filter(sampleType == "NC") %>%
    pull(sampleID)
  
  processed_data$IPC <-processed_data$samples %>%
    filter(sampleType == "IPC") %>%
    pull(sampleID)
  
  processed_data$Calibrator <-processed_data$samples %>%
    filter(sampleType == "Calibrator") %>%
    pull(sampleID)
  
  processed_data$SC <-processed_data$samples %>%
    filter(sampleType == "SC") %>%
    pull(sampleID)
  
  plate_samples <- processed_data$SampleNames
  
  processed_data$aboveLOD <- processed_data$lod$aboveLOD %>%
    rename_cols(., names_df = names_df)
  
  plate_detect <- processed_data$detectability$all$detectability %>%
    tibble::enframe(
      name = "targetName",
      value = "targetDetectability") %>%
    left_join(
      processed_data$lod$LOD %>%
        tibble::enframe(
          name = "targetName",
          value = "targetLOD"),
      by = "targetName") %>%
    left_join(
      processed_data$lod$LODNPQ %>%
        tibble::enframe(
          name = "targetName",
          value = "logged_LOD"),
      by = "targetName") 
  
  # TODO reproducible with `aq_example_P1.xml` & `aq_example_P2.xml`
  if ("AQ" %in% names(processed_data) &&
      !is.null(processed_data$quantifiability$all) &&
      length(processed_data$quantifiability$all$quantifiability) > 0) {
    plate_quant <- processed_data$quantifiability$all$quantifiability %>%
      tibble::enframe(
        name = "targetName",
        value = "targetQuantifiability") %>%
      left_join(
        processed_data$AQ$targetAQ_param %>%
          select(targetName, ULOQ, LLOQ, ULOQ_pg_ml, LLOQ_pg_ml), 
        by = "targetName") %>%
      rename(
        "targetULOQ_aM" = "ULOQ",
        "targetLLOQ_aM" = "LLOQ",
        "targetULOQ_pg_ml" = "ULOQ_pg_ml",
        "targetLLOQ_pg_ml" = "LLOQ_pg_ml"
      )
  }
  
  if("untransformedReverse_LODNPQ" %in% names(processed_data$lod)){
    plate_detect <- plate_detect %>%
      left_join(
        processed_data$lod$untransformedReverse_LODNPQ %>%
          tibble::enframe(
            name = "targetName",
            value = "rev_logged_LOD"),
        by = "targetName"
      )
  }
  
  if("AQ" %in% names(processed_data)){
    if(!is.null(processed_data$lod$LOD_pgmL)){
      pgmL <- tibble::enframe(processed_data$lod$LOD_pgmL, name = "targetName", value = "LOD_pgmL")
      plate_detect <- plate_detect %>%
        left_join(pgmL, by = "targetName") %>%
        mutate(logged_LOD_pgmL = log2(LOD_pgmL))
    } else{
      logger::log_info("pg/mL LOD not found")
    }
    
    if(!is.null(processed_data$lod$LOD_aM)){
      aM <- tibble::enframe(processed_data$lod$LOD_aM, name = "targetName", value = "LOD_aM")
      plate_detect <- plate_detect %>%
        left_join(aM, by = "targetName") %>%
        mutate(logged_LOD_aM = log2(LOD_aM))
    } else{
      logger::log_info("aM LOD not found")
    }
  }
  
  processed_data$targets <- left_join(processed_data$targets, plate_detect,
                                      by = "targetName")
  
  if("quantifiability" %in% names(processed_data) && 
     !is.null(processed_data$quantifiability$all) &&
     length(processed_data$quantifiability$all$quantifiability) > 0){
    processed_data$targets <- left_join(processed_data$targets, plate_quant,
                                        by = "targetName")
  }
  
  if("hide" %in% colnames(processed_data$targets)){
    rmTarget <- processed_data$targets %>%
      filter(hide & targetType == "target") %>%
      pull(targetName)
    processed_data$targets <- subset(processed_data$targets, !targetName %in% rmTarget)
  }
  
  # processed_data$samples$sampleID <- processed_data$samples$sampleName
  
  # convert wellCol to numric 
  processed_data$samples$wellCol <- as.numeric(processed_data$samples$wellCol)
  
  if(is.null(processed_data$plateID)) {
    logger::log_info("null plate id")
    processed_data$plateID <- unique(processed_data$samples$plateID)
  }
  
  # drop AUTO_Plate column 
  #processed_data$samples <- processed_data$samples[,!colnames(processed_data$samples) %in% c("AUTO_PLATE")]
  
  # restructure samples qc dataframe
  names(processed_data$qcSample)[names(processed_data$qcSample) == "type"] <- "sampleType"
  processed_data$qcSample <- processed_data$qcSample[,c("sampleName", "flagName", "val", "status","sampleBarcode", "sampleType")] 
  
  # convert val to be numeric
  processed_data$qcSample$val <- as.numeric(processed_data$qcSample$val)
  processed_data$qcPlate$val <- as.numeric(processed_data$qcPlate$val)
  
  # get updated sampleNames for control samples and add sampleName_original
  processed_data$qcSample <- processed_data$qcSample %>%
    left_join(names_df, by = "sampleName") %>%
    mutate(sampleName = dplyr::if_else(sampleType =="Sample", sampleName, sampleID)) %>%
    left_join(
      processed_data$samples %>% select(sampleName, sampleName_original),
      by = "sampleName"
    )
  
  # reorder samples df columns
  column_order <- c("plateID","plateGUID", "sampleType", "sampleID","sampleBarcode",
                    "sampleName", "sampleName_original", "wellRow", "wellCol", "matching", "non-matching")
  
  processed_data$samples <- processed_data$samples[,c(column_order, 
                                                      colnames(processed_data$samples)[!colnames(processed_data$samples) %in% column_order ])]
  
  # convert zero range values to factor
  to_convert <- processed_data$samples %>%
    select(where(is.numeric)) %>%
    purrr::map(.,zero_range) %>%
    do.call("c",.)
  to_convert <- names(which(to_convert == T))
  
  if(length(to_convert) > 0){
    processed_data$samples <- processed_data$samples %>%
      mutate(across(all_of(to_convert), function(x) forcats::fct(as.character(x))))
  }
  
  # remove unnecessary elements from the processed data 
  processed_data["normed"] <- NULL # now called Data_IC
  processed_data["Data"] <- NULL # now called Data_raw
  processed_data["IC_normed"] <- NULL # we use inter-plate normalized values now (npq)
  processed_data["normed_untransformedReverse"] <- NULL # now called Data_Reverse
  processed_data["AQ"] <- NULL
  
  matrices <-  c("Data_IC","Data_IClog2","Data_raw","Data_rawlog2","aboveLOD", "Data_AQ", "Data_AQlog2", "withinDR", "Data_Reverse", "Data_Reverselog2","Data_AQ_pgmL","Data_AQlog2_pgmL")
  matrices <- matrices[matrices %in% names(processed_data)]
  
  for (i in matrices) {
    logger::log_info("Removing samples with all NaN or NA values -- ",i)
    missing <- apply(processed_data[[i]], 2, function(x) all(is.nan(x))  || all(is.na(x)))
    
    if(i %in% c("Data_AQ", "Data_AQlog2", "Data_AQ_pgmL","Data_AQlog2_pgmL", "withinDR")){
      targets <- processed_data[["aqParams"]] %>%
        pull(targetName)
    } else{
      targets <- processed_data[["targets"]] %>%
        pull(targetName) 
    }
    
    processed_data[[i]] <- processed_data[[i]][targets, !missing]
  }
  
  processed_data[["samples"]] <- processed_data[["samples"]] %>%
    filter(sampleName %in% colnames(processed_data[["Data_IC"]])) %>%
    mutate(sampleType = forcats::fct_drop(sampleType))
  
  processed_data[["qcSample"]] <- processed_data[["qcSample"]] %>%
    filter(sampleName %in% colnames(processed_data[["Data_IC"]])) 
  
  return(processed_data)
}


##' Merge a list of NULISAseq plate data objects
#'
#' This function merges the outputs from multiple NULISAseq plate runs (as produced by process_loadNULISAseq), combining their metadata, sample, target, QC, and data matrices into a single unified object for downstream analysis.
#'
#' @param dataList A list of data objects, each produced by process_loadNULISAseq, representing individual NULISAseq plates.
#' @param fileNameList A list of filenames corresponding to each plate in dataList.
#' @param sample_group_covar Character, the covariate name used for sample grouping in recalculate multiple runs quantifiability. Default is 'SAMPLE_MATRIX', Function will check first to
#' be sure that the variable is present in the column names of the samples matrix.
#' Can be set to NULL to not use this feature.
#' 
#' @return A named list containing:
#'   \item{plateID}{Character vector of plate IDs}
#'   \item{fileNames}{Character vector of file names}
#'   \item{covariateNames}{Character vector of covariate (column) names}
#'   \item{ExecutionDetails}{List of execution details for each plate}
#'   \item{RunSummary}{Data frame summarizing each run}
#'   \item{IC}{Character, the internal control used}
#'   \item{targets}{Data frame of target information}
#'   \item{samples}{Data frame of sample information}
#'   \item{qcSample}{Data frame of sample-level QC metrics}
#'   \item{qcTarget}{Data frame of target-level QC metrics (if present)}
#'   \item{qcPlate}{Data frame of plate-level QC metrics}
#'   \item{aqParams}{Data frame of AQ target parameters (if present)}
#'   \item{excluded}{Vector of excluded targets (if any)}
#'   \item{Data_IC, Data_IClog2, Data_raw, Data_rawlog2, aboveLOD, Data_AQ, Data_AQlog2, Data_Reverse, Data_Reverselog2, Data_AQ_pgmL, Data_AQlog2_pgmL}{Matrices of merged data (if present)}
#'   \item{unit}{Character, AQ unit (if present)}
#'   \item{detectability}{Data frame of target detectability by sample group/matrix of combined runs.}
#'   \item{quantifiability}{Data frame of target quantifiability by sample group/matrix of combined runs (if AQ data present).}
#'
#' @details
#' - Checks for consistent internal control (IC) across plates.
#' - Handles duplicate sample names by appending plate IDs.
#' - Cleans covariate names and merges all available data matrices.
#' - Excludes targets not present in all plates.
#' - Data_IC: IC-normalized data
#' - Data_raw: un-normalized data
#' - Data_AQ: AQ data in aM
#' - Data_AQ_pgmL: AQ data in pgmL
#' - Data_Reverse: Data matrix for reverse curve targets
#' - aboveLOD: Above LOD matrix
#' @examples
#' # merged <- mergeNULISAseq(list_of_plate_data, list_of_filenames)
#'
#' @keywords internal
mergeNULISAseq <- function(dataList, fileNameList, sample_group_covar = "SAMPLE_MATRIX") {
  # get plate IDs
  plateID <- lapply(dataList, function(x){
    x$plateID
  })
  plateID <- unlist(plateID)

  # get Filenames
  fileNames <- do.call("c", fileNameList)
  # Use plateID for naming dataList to ensure unique names
  # (fileNames can be duplicated when same XML is loaded multiple times)
  names(dataList) <- plateID
  
  # merge Execution Details
  ExecutionDetails <- lapply(dataList, function(x){
    x$ExecutionDetails
  })
  names(ExecutionDetails) <- plateID
  
  check_assay_type(ExecutionDetails = ExecutionDetails, fileNames = fileNames)
  
  # merge Run Summary 
  process_RunSummary <- function(x, plateID) {
    x$RunSummary <- lapply(x$RunSummary, replace_empty_numeric)
    x$RunSummary <- data.frame(as.list(x$RunSummary))
    x$RunSummary$plateID <- plateID
    x$RunSummary
  }
  
  RunSummary <- mapply(process_RunSummary, dataList, plateID, SIMPLIFY = F)
  RunSummary <- do.call(rbind, RunSummary)
  
  # merge targets, samples, samples_combined
  # check if targets match across plates
  # (return all targets), returns a list of targets dataframe and vector of
  # inconsistent targets
  targets <- combine_targets(dataList, plateID)
  if("list" == class(targets)){
    excluded <- targets$excluded
    targets <- targets$targets
  }
  
  # merge ICs 
  IC <- unique(do.call("c", lapply(dataList, "[[", "IC")))
  if(length(IC) > 1){
    logger::log_error("Non-Unique IC normalization- ", paste(IC, collapse = ", "))
    stop("Non-Unique IC normalization")
  }
  
  # merge samples
  samples <- lapply(dataList, function(x){
    x$samples
  })
  
  # extract shared sample columns between plates
  # TODO https://github.com/Alamar-Biosciences/NULISA-Analysis-Software/issues/2192
  shared_sample_columns <- Reduce(union, lapply(samples, colnames))
  
  # get classes of all columns in the all the xmls
  check_class <- lapply(samples, function(x){
    do.call('c', lapply(x, class))
  }) %>%
    dplyr::bind_rows() %>%
    tidyr::pivot_longer(cols = everything()) %>%
    dplyr::rename(columns = name) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::select(columns, value) %>%
    unique()
  
  # check if any of the columns has more than one class
  change_class_of <- check_class %>%
    dplyr::summarize(n = dplyr::n(), .by = columns) %>%
    dplyr::filter(n > 1) %>%
    dplyr::pull(columns)
  
  if(length(change_class_of) != 0){
    samples <- lapply(samples, function(x){
      # Only convert columns that actually exist in this dataframe
      cols_to_convert_here <- intersect(change_class_of, colnames(x))
      if(length(cols_to_convert_here) > 0) {
        x[cols_to_convert_here] <- lapply(x[cols_to_convert_here], function(x) forcats::fct(as.character(x)))
      }
      x
    }) %>%
      purrr::map(., ~ {.x %>%
          dplyr::select(dplyr::any_of(shared_sample_columns))
      }) %>%
      dplyr::bind_rows() %>%
      tibble::as_tibble()
  } else{
    samples <- samples %>%
      purrr::map(., ~ {.x %>%
          dplyr::select(dplyr::any_of(shared_sample_columns))
      }) %>%
      dplyr::bind_rows() %>%
      tibble::as_tibble()
  }
  
  if(class(samples$plateID) != "factor") {
    samples$plateID <- forcats::fct(samples$plateID)
  }
  
  # clean column names of punctuation
  default_names <- colnames(samples)
  clean_names <- clean_covariate_names(default_names, case = "lower_camel")
  
  if(all(default_names == clean_names)){
    logger::log_info("No sanitation performed")
  }else{
    logger::log_info("Sanitation perfromed")
  }
  names(default_names) <- clean_names
  
  # Get a list of duplicated sampleNames
  dup_sampleName <- samples %>%
    group_by(sampleName) %>%
    filter(dplyr::n()>1) %>%
    pull(sampleName)
  
  if(length(dup_sampleName)>0){
    # create a map of the old and new names for the duplicated samples
    rename_map <- samples %>%
      filter(sampleName %in% dup_sampleName) %>%
      mutate(sampleID = dplyr::if_else(sampleName %in% dup_sampleName, 
                                paste(sampleName, plateID, sep = "_"), sampleName)) %>%
      select(plateID, sampleName, sampleID)
    
    samples <- samples %>% 
      select(-sampleID) %>%
      left_join(rename_map, by = c("plateID", "sampleName")) %>%
      mutate(sampleName = dplyr::if_else(is.na(sampleID), sampleName, sampleID),
             sampleID = sampleName)
  }else{
    rename_map <- data.frame("plateID"=NA, "sampleName" = NA, "sampleID" = NA)
  }
  
  matrices <- c("Data_IC","Data_IClog2","Data_raw","Data_rawlog2","aboveLOD","Data_AQ", "Data_AQlog2","Data_Reverse","Data_Reverselog2",
                "Data_AQ_pgmL","Data_AQlog2_pgmL")
  
  for (i in seq_along(dataList)) {
    pID <- dataList[[i]]$plateID
    if(pID %in% rename_map$plateID){

      plt_map <- rename_map %>%
        filter(plateID == pID) %>%
        select(-plateID)

      # Update samples dataframe with renamed sampleNames
      dataList[[i]]$samples <- dataList[[i]]$samples %>%
        left_join(plt_map, by = "sampleName") %>%
        mutate(sampleName = dplyr::if_else(!is.na(sampleID.y), sampleID.y, sampleName),
               sampleID.x = sampleName) %>%
        select(-sampleID.y) %>%
        rename(sampleID = sampleID.x)

      dataList[[i]]$qcSample <- dataList[[i]]$qcSample %>%
        select(-sampleID)%>%
        left_join(plt_map, by = "sampleName") %>%
        mutate(sampleName = dplyr::if_else(is.na(sampleID), sampleName, sampleID))

      # Update control sample name vectors (IPC, SC, NC, Bridge, Calibrator, SampleNames)
      if(!is.null(dataList[[i]]$IPC)) {
        dataList[[i]]$IPC <- sapply(dataList[[i]]$IPC, function(name) {
          ifelse(name %in% plt_map$sampleName, plt_map$sampleID[plt_map$sampleName == name], name)
        }, USE.NAMES = FALSE)
      }
      if(!is.null(dataList[[i]]$SC)) {
        dataList[[i]]$SC <- sapply(dataList[[i]]$SC, function(name) {
          ifelse(name %in% plt_map$sampleName, plt_map$sampleID[plt_map$sampleName == name], name)
        }, USE.NAMES = FALSE)
      }
      if(!is.null(dataList[[i]]$NC)) {
        dataList[[i]]$NC <- sapply(dataList[[i]]$NC, function(name) {
          ifelse(name %in% plt_map$sampleName, plt_map$sampleID[plt_map$sampleName == name], name)
        }, USE.NAMES = FALSE)
      }
      if(!is.null(dataList[[i]]$Bridge)) {
        dataList[[i]]$Bridge <- sapply(dataList[[i]]$Bridge, function(name) {
          ifelse(name %in% plt_map$sampleName, plt_map$sampleID[plt_map$sampleName == name], name)
        }, USE.NAMES = FALSE)
      }
      if(!is.null(dataList[[i]]$Calibrator)) {
        dataList[[i]]$Calibrator <- sapply(dataList[[i]]$Calibrator, function(name) {
          ifelse(name %in% plt_map$sampleName, plt_map$sampleID[plt_map$sampleName == name], name)
        }, USE.NAMES = FALSE)
      }
      if(!is.null(dataList[[i]]$SampleNames)) {
        dataList[[i]]$SampleNames <- sapply(dataList[[i]]$SampleNames, function(name) {
          ifelse(name %in% plt_map$sampleName, plt_map$sampleID[plt_map$sampleName == name], name)
        }, USE.NAMES = FALSE)
      }

      mats <- matrices [matrices %in% names(dataList[[i]])]
      ## TODO the two detectability files fail here
      for (j in mats) {
        dataList[[i]][[j]] <- dataList[[i]][[j]] %>%
          rename_cols(., names_df = plt_map)
      }
    }
  }

  # merge qcSample
  qcSample <- mapply(add_plateID, dataList, "qcSample", plateID, SIMPLIFY = F)
  qcSample <- do.call(rbind, qcSample)
  
  # merge qcPlate
  qcPlate <- mapply(add_plateID, dataList, "qcPlate", plateID, SIMPLIFY = F)
  qcPlate <- do.call(rbind, qcPlate)
  
  # merge qcTarget if it exists
  qcTargets_exists <- "qcTarget" %in% unique(unlist(lapply(dataList, names)))
  if(qcTargets_exists){
    logger::log_info("Target QC dataframe found")
    qcTarget <- mapply(add_plateID, dataList, "qcTarget", plateID, SIMPLIFY = F)
    qcTarget <- do.call(rbind, qcTarget)
  } else{
    logger::log_info("Target QC dataframe not found")
    qcTarget <- NULL
  }
  
  # merge AQ_params if it exists
  target_params_exists <- "aqParams" %in% unique(unlist(lapply(dataList, names)))
  if(target_params_exists){
    logger::log_info("AQ Target params dataframe found")
    aqParams <- mapply(add_plateID, dataList, "aqParams", plateID, SIMPLIFY = F)
    aqParams <- bind_rows(aqParams)
  } else{
    logger::log_info("AQ Target params dataframe not found")
    aqParams <- NULL 
  }
  
  # Merge all data matrices
  dataMatrix <- list()
  # get all matrices available in the dataList object
  mats <- do.call("c",lapply(dataList, lapply, class)) %>%
    stack() %>%
    filter(values %in% c("matrix","array")) %>%
    mutate(mats = sub(".*\\.", "", ind)) %>%
    filter(grepl("_|LOD", mats)) %>%
    pull(mats) %>%
    unique()
  
  unit <- NULL
  if("Data_AQ" %in% mats){
    unit <- unique(do.call("c", lapply(dataList,"[[","AQ_unit")))
    if(length(unit) > 1){
      stop("Non unique AQ units")
    }
  }
  
  for (i in mats) {
    logger::log_info("Merging NULISAseq Data Matrix -- ",i)
    
    dataObj <- lapply(dataList, function(x){
      x[[i]] %>%
        tibble::as_tibble(rownames = "targetName")
    })
    
    dataMatrix[[i]] <- dataObj %>%
      purrr::reduce(dplyr::full_join, by = "targetName") %>%
      filter(targetName %in% targets$targetName) %>%
      tibble::column_to_rownames("targetName") %>%
      as.matrix()
  }
  
  # Detetability by covariate
  detect_summary <- detectability_summary(dataList)
  
  if ("sample_group" %in% names(detect_summary)) {
    detect_table <- detect_summary$sample_group$detectability
  } else {
    detect_table <- detect_summary$all$detectability
  }
  detect_table[, 2:ncol(detect_table)] <- round(detect_table[,2:ncol(detect_table)], 1)
  colnames(detect_table) <- ifelse(colnames(detect_table) == "Target", "Target", tolower(colnames(detect_table)))
  
  # Quantifiabiltiy by covariate
  if (any(sapply(dataList, function(x) "quantifiability" %in% names(x)))) {
    
    aq_re <- lapply(dataList, function(x){
      x$AQ <- list(
        targetAQ_param = x$aqParams,
        Data_AQ_aM = x$Data_AQ,
        Data_AQ = x$Data_AQ_pgmL,
        withinDR = x$withinDR
      )
      return(x)
    })
    
    quant <- quantifiability(runs = aq_re, 
                             sampleGroupCovar = sample_group_covar)
    
    quant_table <- round(quant$combined_quantifiability$quant, 1)
    quant_table$target <- rownames(quant_table)
    colnames(quant_table)[colnames(quant_table) == "overall"] <- "quantifiability"
    colnames(quant_table) <- tolower(colnames(quant_table))
    colnames(quant_table)[1:(ncol(quant_table) - 1)] <- paste0(colnames(quant_table)[1:(ncol(quant_table) - 
                                                                                          1)], " (n = ", quant$combined_quantifiability$n_samples, 
                                                               ")")
    colnames(quant_table)[colnames(quant_table) == "target"] <- "Target"
    if (ncol(quant_table) < 3) {
      quant_table <- quant_table[, c(ncol(quant_table), 
                                     1)]
    } else {
      quant_table <- quant_table[, c(ncol(quant_table), 
                                     2:(ncol(quant_table) - 1))]
    }
  }
  
  # Extract control sample names and SampleNames from merged samples dataframe
  IPC_samples <- samples$sampleName[samples$sampleType == 'IPC']
  SC_samples <- samples$sampleName[samples$sampleType == 'SC']
  NC_samples <- samples$sampleName[samples$sampleType == 'NC']
  Bridge_samples <- samples$sampleName[samples$sampleType == 'Bridge']
  Calibrator_samples <- samples$sampleName[samples$sampleType == 'Calibrator']
  SampleNames_samples <- samples$sampleName[samples$sampleType == 'Sample']

  # Create base return list (always include IPC, SC, NC, SampleNames)
  return_list <- list(
    plateID = plateID,
    fileNames = fileNames,
    covariateNames = default_names,
    ExecutionDetails = ExecutionDetails,
    RunSummary = RunSummary,
    IC = IC,
    targets = targets,
    samples = samples,
    qcSample = qcSample,
    qcTarget = qcTarget,
    qcPlate = qcPlate,
    aqParams = aqParams,
    inconsistent_targets = excluded,
    detectability = detect_table,
    IPC = IPC_samples,
    SC = SC_samples,
    NC = NC_samples,
    SampleNames = SampleNames_samples
  )

  # Only include Bridge and Calibrator if they have samples (optional sample types)
  if (length(Bridge_samples) > 0) {
    return_list$Bridge <- Bridge_samples
  }
  if (length(Calibrator_samples) > 0) {
    return_list$Calibrator <- Calibrator_samples
  }
  if (exists("quant_table")) {
    return_list$quantifiability <- quant_table
  }
  # Add dataMatrix and unit
  return_list <- c(return_list, dataMatrix, unit = unit)
  
  # return the output
  return(return_list)
}

#' Import and Process Multiple NULISAseq Runs
#'
#' This function loads multiple NULISAseq data files, processes them, and optionally
#' merges them into a single dataset. It provides extensive control over quality
#' control parameters and data processing options for NULISAseq data.
#'
#' @param files Character vector of path and name of the NULISAseq data xml files.
#' @param plateName Optional character vector of names for each file/run/plate. If \code{NULL},
#' will use AUTO_PLATE variable within xml files or generate default names (Plate_01, Plate_02, etc)
#' if there are duplicate or missing AUTO_PLATE names. 
#' Length must match `files`.
#' @param return_type Character specifying output format: "all" (list with both individual runs 
#' and merged data), "run" (individual runs only), or "merged" 
#' (merged data only) (Default: "all").
#' @param IC Default is `NULL`. Optional character string giving internal 
#' control target name that  
#' will override pre-defined target type for all plates. 
#' @param IPC Default is `NULL`. Override default Inter-plate Control sample specification.
#' Optional vector (applied to all plates) 
#' or named list with `plateName` 
#' names and vectors of character string(s) that match the IPC sample names or 
#' substrings of sample names. These will override pre-defined sample type 
#' and use these samples for inter-plate control normalization. 
#' (e.g., `list("Plate_01" = c("Plate01_IPC_sampleName_1", "Plate01_IPC_sampleName_2", "Plate01_IPC_sampleName_3"), "Plate_02" = c("Plate_02_IPC_sampleName_1", "Plate_02_IPC_sampleName_2", "Plate_02_IPC_sampleName_3"))`). 
#' @param SC Default is `NULL`. Override default Sample Control sample specification. Same format as `IPC`.
#' @param NC Default is `NULL`. Override default Negative Control sample specification. Same format as `IPC`.
#' @param Bridge Default is `NULL`. Override default Bridge control sample specification. Same format as `IPC`.
#' @param Calibrator Default is `NULL`. Override default calibrator sample specification. Same format as `IPC`.
#' @param sample_group_covar Sample group covariate. 
#' Optional column name in the samples data matrix (from Barcode B file) for each
#' plate that represents subgroups for which detectability and quantifiability (if AQ available) will be calculated separately 
#' in addition to overall detectability and overall quantifiability. Default is `SAMPLE_MATRIX`. 
#' Function will first check to be sure that the variable is present in the 
#' column names of the samples matrix.
#' Can be set to NULL to not use this feature.
#' @param excludeSamples Samples to exclude from analysis. Can be: NULL (no exclusions),
#'   a vector of sample names (applied to all plates), or a named list 
#'   (e.g., `list("Plate_01" = c("Sample1", "Sample2"), "Plate_02" = c("Sample3", "Sample4"))`).
#' @param excludeTargets Targets to exclude from analysis. Can be: NULL (no exclusions),
#'   a vector of target names (applied to all plates), or a named list
#'   (e.g., `list("Plate_01" = c("Target1"), "Plate_02" = c("Target2", "Target3"))`).
#' @param include_qc Logical indicating whether to include QC (Quality Control) columns
#'   in the long format output (`Data_NPQ_long` and `Data_AQ_long`). When `TRUE` (default),
#'   Sample QC and Target QC metrics are added to the long format data. Target QC columns are
#'   automatically filtered based on the data mode: RQ data (`Data_NPQ_long`) excludes AQ-specific
#'   metrics like concentration accuracy and CV calculated on AQ values, while AQ data
#'   (`Data_AQ_long`) includes all Target QC metrics. Set to `FALSE` to exclude all QC
#'   information and reduce output size (Default: TRUE).
#' @param verbose Logical indicating whether to display progress messages (Default: TRUE).
#' @param ... Additional arguments passed to \code{\link{loadNULISAseq}}.
#'
#' @return Depending on `return_type`:
#' \itemize{
#'   \item{"all"}{List with two components: 
#'     \describe{
#'       \item{runs}{List of individual runs, each containing an object with the following structure:
#'         \itemize{
#'           \item{\code{plateID}: Character vector of plate identifier}
#'           \item{\code{ExecutionDetails}: Run metadata including command line, execution time, instrument details, lot information}
#'           \item{\code{RunSummary}: Data frame of read statistics for the run (total reads, parseable reads, matches, etc.)}
#'           \item{\code{targets}: Data frame of target metadata (target names, UniProt IDs, QC flags, etc.)}
#'           \item{\code{samples}: Data frame of sample metadata (sample names, matrix type, sample plate well positions, etc.)}
#'           \item{\code{Data_raw}: Raw count data (matrix, targets × samples)}
#'           \item{\code{attributes}: List of quality attributes including QCS and SN matrices}
#'           \item{\code{IC, IPC, SC, NC, Bridge}: Character vector of control target or sample names}
#'           \item{\code{IC_normed}: Internal control normalized data (matrix, targets × samples)}
#'           \item{\code{normed_untransformedReverse, normed}: Normalized data in multiple formats (matrix, targets × samples), before (`normed_untransformedReverse`) and after (`normed`) reverse-curve transformation is applied to any relevant targets}
#'           \item{\code{AQ}: Absolute quantification data if available, including:
#'             \itemize{
#'               \item{\code{Data_AQ_aM}: Absolute quantification in aM units (matrix, targets × samples)}
#'               \item{\code{Data_AQ}: Absolute quantification in pg/mL units (matrix, targets × samples)}
#'               \item{\code{targetAQ_param}: Data frame of curve parameters per target}
#'               \item{\code{withinDR}: Logical matrix indicating targets within dynamic range (targets × samples)}
#'             }
#'           }
#'           \item{\code{NPQ}: NULISA Protein Quantification (NPQ) data (log2 scale) (matrix, targets × samples)}
#'           \item{\code{lod}: Limit of detection information including:
#'             \itemize{
#'               \item{\code{LOD}: LOD values per target, on the unlogged normalized count scale}
#'               \item{\code{aboveLOD}: Logical matrix indicating detection above LOD (targets × samples)}
#'               \item{\code{blank_outlier_table}: NC outlier analysis results}
#'               \item{\code{LODNPQ}: LOD in NPQ units}
#'               \item{\code{LOD_pgmL, LOD_aM}: LOD in different units if AQ data available}
#'             }
#'           }
#'           \item{\code{detectability}: List of group-specific and overall detectability}
#'           \item{\code{qcTarget, qcSample, qcPlate}: Quality control flags for targets, samples, and plates}
#'           \item{\code{quantifiability}: List of quantifiability metrics by sample group and overall, if AQ data available}
#'           \item{\code{qcSamplebyTarget}: Sample-by-target QC matrices (read threshold, LOD, dynamic range)}
#'           \item{\code{AbsAssay, advancedQC}: Assay type and advanced target QC flags}
#'         }
#'       }
#'       \item{merged}{Merged dataset from all runs containing:
#'         \itemize{
#'           \item{\code{plateID}: Character vector of plate identifiers}
#'           \item{\code{fileNames}: Character vector of input file names}
#'           \item{\code{covariateNames}: Character vector of sample covariate names}
#'           \item{\code{ExecutionDetails}: Per-plate metadata (command line, run time, instrument, assay, lot info)}
#'           \item{\code{RunSummary}: Data frame of read statistics (total reads, parseable, matches, etc.)}
#'           \item{\code{IC}: Identifier of internal control}
#'           \item{\code{targets}: Data frame of target metadata including UniProt IDs, QC flags, LOD/ULOQ/LLOQ, detectability}
#'           \item{\code{samples}: Data frame of sample metadata including names, well positions, sample matrices, other covariates}
#'           \item{\code{qcSample}: Data frame of sample-level QC flags (e.g., IC_Median) with values and status}
#'           \item{\code{qcTarget}: Data frame of target-level QC flags (e.g., concentration accuracy, thresholds)}
#'           \item{\code{qcPlate}: Data frame of plate-level QC flags (e.g., SC CV, WARN targets, thresholds)}
#'           \item{\code{aqParams}: Data frame of curve parameters and quantification metrics 
#'           (LLOQ, ULOQ, LOD) (If AQ data available)}
#'           \item{\code{inconsistent_targets}: Placeholder for targets inconsistent across plates/runs (NULL if none)}
#'           \item{\code{detectability}: Data frame of detectability by target and sample matrix}
#'           \item{\code{quantifiability}: Data frame of quantifiability by target and sample matrix}
#'           \item{\code{Data_raw}, \code{Data_rawlog2}: Raw counts and log2-transformed counts (matrix, targets × samples)}
#'           \item{\code{Data_IC}, \code{Data_IClog2}: Internal control–normalized data (linear and log2)}
#'           \item{\code{Data_Reverse}, \code{Data_Reverselog2}: Reverse-transformed IC-IPC normalized data (linear and log2)}
#'           \item{\code{Data_AQ_aM}, \code{Data_AQlog2_aM}: Absolute quantitation in attomolar units (linear and log2), 
#'           if AQ data available}
#'           \item{\code{Data_AQ_pgmL}, \code{Data_AQlog2_pgmL}: Absolute quantitation in pg/mL units (linear and log2), 
#'           if AQ data available}
#'           \item{\code{aboveLOD}: Logical matrix indicating values above LOD (targets × samples)}
#'           \item{\code{unit}: Character string of concentration units (e.g., "pg/mL")}
#'           \item{\code{Data_NPQ}: Matrix of NULISA Protein Quantification (NPQ) values (log2)}
#'           \item{\code{Data_NPQ_long}: Data frame (long format) of NPQ data with sample and target annotations}
#'           \item{\code{Data_AQ_long}: Data frame (long format) of absolute quantification data with concentrations and LOD/LLOQ/ULOQ, 
#'           if AQ data available}
#'         }
#'       }
#'     }
#'   }
#'   \item{"run"}{List of individual runs (as described in the 'runs' component above)}
#'   \item{"merged"}{Merged dataset from all runs (as described in the 'merged' component above)}
#' }
#' 
#' 
#' @examples
#' \dontrun{
#' result <- importNULISAseq(
#'   files = c("run1.xml", "run2.xml"),
#'   plateName = c("Experiment1", "Experiment2"),
#'   excludeSamples = list("Experiment1" = c("SampleA"), "Experiment2" = c("SampleB")),
#'   return = "all"
#' )
#' }
#'
#' @seealso
#' \code{\link{loadNULISAseq}} for loading individual NULISAseq runs
#'
#' @importFrom dplyr select filter mutate arrange group_by summarise rename distinct pull left_join
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tibble column_to_rownames rownames_to_column
#' @export
importNULISAseq <- function(files,
                            plateName = NULL,
                            return_type = c("all", "run", "merged"),
                            IC = NULL,
                            IPC = NULL,
                            SC = NULL,
                            NC = NULL,
                            Bridge = NULL,
                            Calibrator = NULL,
                            sample_group_covar = 'SAMPLE_MATRIX',
                            excludeSamples = NULL,
                            excludeTargets = NULL,
                            include_qc = TRUE,
                            verbose = TRUE,
                            ...) {
  
  # Get NAS_compatible from option
  NAS_compatible = getOption("NAS_compatible", FALSE)
  
  # Validate return parameter
  return_type <- match.arg(return_type)
  
  # Input validation
  if (length(files) == 0) {
    stop("No files provided")
  }
  
  if (!all(file.exists(files))) {
    missing_files <- files[!file.exists(files)]
    stop(paste("Files not found:", paste(missing_files, collapse = ", ")))
  }
  
  # Validate plateName
  user_provided_plateName <- !is.null(plateName)
  if (user_provided_plateName) {
    if (length(plateName) != length(files)) {
      stop(sprintf("plateName length (%d) does not match files length (%d)",
                   length(plateName), length(files)))
    } else if (any(duplicated(plateName))) {
      stop(sprintf("plateName contains duplicate values: %s",
                   paste(unique(plateName[duplicated(plateName)]), collapse = ", ")))
    }
  }
  
  # Check IC only accept vector
  if(!is.null(IC) && !is.vector(IC)) {
    stop("IC must be the same across plates/runs. IC must be NULL or a vector")
  }
  
  # Process all parameters that can be named lists
  param_names_for_processing <- if (user_provided_plateName) plateName else seq_along(files)
  
  IC <- process_named_param(IC, param_names_for_processing) # should be a vector, same IC across plates
  IPC <- process_named_param(IPC, param_names_for_processing)
  SC <- process_named_param(SC, param_names_for_processing)
  NC <- process_named_param(NC, param_names_for_processing)
  Bridge <- process_named_param(Bridge, param_names_for_processing)
  Calibrator <- process_named_param(Calibrator, param_names_for_processing)
  excludeSamples <- process_named_param(excludeSamples, param_names_for_processing)
  excludeTargets <- process_named_param(excludeTargets, param_names_for_processing)
  
  if (verbose) {
    logger::log_info("Processing {length(files)} NULISAseq files.")
    logger::log_info("{paste('\n', seq_along(files), basename(files), collapse = '')}")
  }
  
  # Helper to select from list parameters
  pick_param <- function(param, i) {
    if (is.null(param)) return(NULL)
    return(param[[i]])  
  }
  
  # Load data with error handling
  runs <- list()
  for (i in seq_along(files)) {
    if (verbose) message(sprintf("Loading file %d/%d: %s", i, length(files), basename(files[i])))
    plateID_to_pass <- if (user_provided_plateName) plateName[i] else NULL
    tryCatch({
      runs[[i]] <- loadNULISAseq(
        file = files[i],
        IC = pick_param(IC, i),
        IPC = pick_param(IPC, i),
        SC = pick_param(SC, i),
        NC = pick_param(NC, i),
        Bridge = pick_param(Bridge, i),
        Calibrator = pick_param(Calibrator, i), 
        sample_group_covar = sample_group_covar,
        plateID = plateID_to_pass, # Pass user-provided plateID or NULL
        excludeSamples = pick_param(excludeSamples, i),
        excludeTargets = pick_param(excludeTargets, i),
        ...
      )
    }, error = function(e) {
      stop("Error loading file ", files[i], ": ", e$message)
    })
  }
  
  PlateNameMessage <- NULL
  
  # Check for UUID plateID 
  for (i in seq_along(runs)) {
    if (any(uuid::UUIDvalidate(runs[[i]]$samples$plateID))) {
      if (!is.null(runs[[i]]$ExecutionDetails) && 
          !is.null(runs[[i]]$ExecutionDetails$ExperimentName)) {
  
        if (verbose) message(sprintf("UUID plateID found for run %d, replacing with ExecutionDetails$ExperimentName", i))
        runs[[i]]$samples$plateID <- runs[[i]]$ExecutionDetails$ExperimentName
        runs[[i]]$plateID <- runs[[i]]$ExecutionDetails$ExperimentName
      } else {
        if (verbose) message(sprintf("UUID plateID found for run %d, but ExecutionDetails$ExperimentName is missing. Skipping replacement.", i))
      }
    }
  }
  
  # Check and assign plateName if not user-provided
  if (!user_provided_plateName) {
    # Get plate names from run objects, prioritizing the assigned plateID/AUTO_PLATE
    plateNames <- sapply(runs, function(x) {
      plate_name_i <- x$plateID # Use plateID after potential UUID replacement
      if (is.null(plate_name_i) && 
          !is.null(x$samples$AUTO_PLATE)) {
        plate_name_i <- unique(x$samples$AUTO_PLATE)[1]
      }
      if(is.null(plate_name_i)) plate_name_i <- NA
      return(plate_name_i)
    })
    
    if(any(is.na(plateNames)) || any(duplicated(plateNames))) {
      # Assign Plate_XX or Plate_XXX
      if(length(files) > 99) {
        plateNames <- paste0('Plate_', formatC(1:length(files) , width=3, format='d', flag='0'))
      } else {
        plateNames <- paste0('Plate_', formatC(1:length(files) , width=2, format='d', flag='0'))
      }
      PlateNameMessage <- 'NOTE: Plate ID variable had missing, UUID, or duplicate entries so plate IDs were automatically assigned (starting with Plate_01/Plate_001).'
    }
    
    # Update runs with the final determined plateName
    for (i in seq_along(runs)) {
      runs[[i]]$plateID <- plateNames[i]
      runs[[i]]$samples$plateID <- plateNames[i]
    }
    
    plateName <- plateNames
  }
  
  # Reorder runs if auto-assigned Plate_XX/Plate_XXX names were used
  plate_pattern <- "^Plate_\\d+$"
  
  if (all(grepl(plate_pattern, plateName))) {
    sorted_order <- order(plateName)
    runs <- runs[sorted_order]
    plateName <- plateName[sorted_order]
  }
  
  names(runs) <- plateName
  
  if (!is.null(PlateNameMessage) && verbose) {
    message(PlateNameMessage)
  }
  
  # If only runs are requested, return early
  if (return_type == "run") {
    if (verbose) message("Returning loaded runs only")
    return(runs)
  }
  
  if (verbose) message(sprintf("Processing %d loaded runs ", length(runs)))
  
  # Process each run with error handling
  processed_runs <- list()
  for (i in seq_along(runs)) {
    tryCatch({
      if (verbose) message(sprintf("Processing run %d/%d: %s", i, length(runs), basename(files[i])))
      processed_runs[[i]] <- process_loadNULISAseq(runs[[i]])
      if (verbose) message(sprintf("Successfully processed %s", basename(files[i])))
    }, error = function(e) {
      stop(sprintf("Error processing run %d (%s): %s", i, basename(files[i]), e$message))
    })
  }
  
  # Merge
  merged_data <- NULL
  original_file_names <- sapply(runs, function(x) x$xmlFile)
  named_list <- setNames(as.list(original_file_names), plateName)
  
  tryCatch({
    merged_data <- mergeNULISAseq(dataList = processed_runs, fileNameList = named_list, sample_group_covar = sample_group_covar)
    merged_data$Data_NPQ_long <- format_wide_to_long(merged_data, AQ = FALSE, include_qc = include_qc)
    merged_data$Data_NPQ <- merged_data$Data_IClog2

    if(any(grepl("^Data_AQ", names(merged_data)))) {
      merged_data$Data_AQ_long <- format_wide_to_long(merged_data, AQ = TRUE, include_qc = include_qc)
      
      names(merged_data)[names(merged_data) == "Data_AQ"] <- "Data_AQ_aM"
      names(merged_data)[names(merged_data) == "Data_AQlog2"] <- "Data_AQlog2_aM" 
    }
    if (verbose) {
      message("Data successfully merged")
      message(sprintf("Final merged dataset contains %d runs", length(processed_runs)))
    }
  }, error = function(e) {
    stop("Error during data merging: %s", conditionMessage(e))
  })
  
  if (NAS_compatible) {
    if (verbose) message("Applying output modifications for NAS")
    # NAS compatible mode: Remove long format data
    if ("Data_NPQ_long" %in% names(merged_data)) {
      merged_data$Data_NPQ_long <- NULL
    }
    if ("Data_AQ_long" %in% names(merged_data)) {
      merged_data$Data_AQ_long <- NULL
    }
  } else {
    # Default mode: Remove covariateName, RunSummary, NULL object and unused data
    # also removed IC from data, except from Data_raw
    if ("covariateNames" %in% names(merged_data)) {
      merged_data$covariateNames <- NULL
    }
    if ("RunSummary" %in% names(merged_data)) {
      merged_data$RunSummary <- NULL
    }
    data_to_drop <- c("Data_rawlog2", "Data_IC", "Data_IClog2", "Data_Reverse", "Data_Reverselog2")
    objects_to_remove <- intersect(data_to_drop, names(merged_data))
    if (length(objects_to_remove) > 0) {
      merged_data[objects_to_remove] <- NULL
    }
    merged_data <- Filter(Negate(is.null), merged_data)
    
    # Remove IC target from wide data matrices and long data frames
    if (!is.null(merged_data$IC) && length(merged_data$IC) > 0) {
      ic_targets <- merged_data$IC
      
      # Remove from wide matrices
      wide_matrices <- c("Data_NPQ", "Data_AQ_aM", "Data_AQlog2_aM", "aboveLOD")
      for (mat_name in wide_matrices) {
        if (mat_name %in% names(merged_data)) {
          merged_data[[mat_name]] <- remove_ic_from_matrix(merged_data[[mat_name]], ic_targets)
        }
      }
      
      # Remove from long data frames
      long_dfs <- c("Data_NPQ_long", "Data_AQ_long")
      for (df_name in long_dfs) {
        if (df_name %in% names(merged_data)) {
          merged_data[[df_name]] <- remove_ic_from_long(merged_data[[df_name]], ic_targets)
        }
      }
    }
  }
  
  # Return based on parameter
  if (return_type == "merged") {
    return(merged_data)
  } else if (return_type == "all") {
    return(list(runs = runs, merged = merged_data))
  }
}

##' Merge two processed NULISAseq data objects
#'
#' This function merges two NULISAseq data objects that have already been processed (e.g., by mergeNULISAseq), combining their metadata, sample, target, QC, and data matrices into a single unified object for downstream analysis. It handles duplicate plate IDs, file names, and sample names, and ensures compatibility of assay types and AQ units.
#'
#' @param existing_data A processed NULISAseq data object (as returned by mergeNULISAseq) representing the existing dataset.
#' @param new_data A processed NULISAseq data object (as returned by mergeNULISAseq) representing the new dataset to be merged.
#'
#' @return A named list containing:
#'   \item{plateID}{Character vector of plate IDs}
#'   \item{fileNames}{Character vector of file names}
#'   \item{ExecutionDetails}{List of execution details for each plate}
#'   \item{RunSummary}{Data frame summarizing each run}
#'   \item{IC}{Character, the internal control used}
#'   \item{targets}{Data frame of target information}
#'   \item{samples}{Data frame of sample information}
#'   \item{qcSample}{Data frame of sample-level QC metrics}
#'   \item{qcTarget}{Data frame of target-level QC metrics (if present)}
#'   \item{qcPlate}{Data frame of plate-level QC metrics}
#'   \item{aqParams}{Data frame of AQ target parameters (if present)}
#'   \item{excluded}{Vector of excluded targets (if any)}
#'   \item{Data_IC, Data_IClog2, Data_raw, Data_rawlog2, aboveLOD, Data_AQ, Data_AQlog2, Data_Reverse, Data_Reverselog2, Data_AQ_pgmL, Data_AQlog2_pgmL}{Matrices of merged data (if present)}
#'   \item{unit}{Character, AQ unit (if present)}
#'
#' @details
#' - Checks for duplicate plate IDs, file names, and sample names, and renames as needed.
#' - Ensures consistent internal control (IC) and AQ units across datasets.
#' - Combines all available data matrices, keeping only common targets.
#' - Handles missing or NULL QC and AQ parameter data frames gracefully.
#'
#' @examples
#' # merged <- mergeProcessedNULISAseq(existing_data, new_data)
#'
#' @keywords internal
mergeProcessedNULISAseq <- function(existing_data, new_data){
  # Check if there are any duplicate plateID in the existing and new data objects
  dups <- check_duplicates(existing_data, new_data, type = "plateID")
  dup_filename <- check_duplicates(existing_data, new_data, type = "fileNames")
  dup_sampleName <- check_duplicates(existing_data, new_data, type = "sampleName")
  
  if(dups & dup_filename){
    logger::log_error("Identical FileNames and plateIDs found, please re-name newer datafiles")
    stop("Identical FileNames and plateIDs found, please re-name newer datafiles")
  }
  
  # if duplicates are present, then rename the new object with its fileNames as plateID
  # changes are also applied to their sampleNames in the object
  if(dups | dup_sampleName){
    if(dups){
      logger::log_info("Found duplicate plateID, fixing dups...")
    } else{
      logger::log_info("Found duplicate sampleName, fixing dups...")
    }
    new_data <- fix_duplicate_plates(new_object = new_data, existingSamples = existing_data$samples$sampleName)
  }
  
  plateID <- c(existing_data$plateID, new_data$plateID)
  fileNames <- c(existing_data$fileNames, new_data$fileNames)
  
  unit <- unique(c(existing_data$unit, new_data$unit))
  if(length(unit) > 1){
    logger::log_error("Non unique AQ units")
    stop("Non unique AQ units")
  }
  
  dataList <- list(existing_data, new_data)
  # Combine Execution Details and check for compatibility
  ExecutionDetails <- c(existing_data[["ExecutionDetails"]], new_data[["ExecutionDetails"]])
  check_assay_type(ExecutionDetails)
  
  #Combine Run Summary object
  RunSummary <- dataList %>%
    purrr::map(., ~.x[["RunSummary"]]) %>%
    bind_rows() %>%
    tibble::as_tibble() %>%
    dplyr::relocate("plateID")
  
  # Combine the targets dataframe to keep only common targets between the two objects
  targets <- combine_targets(dataList, addPlateID = FALSE)
  if("list" == class(targets)){
    excluded <- targets$excluded
    targets <- targets$targets
  }
  
  # Generate a vector of the shared columns
  samples <- lapply(dataList, function(x) x$samples)
  shared_sample_columns <- Reduce(union, lapply(samples, colnames))
  # merge the two object together, keeping only the shared columns
  # TODO: Test if this will remove any covariate information, if it is missing in either objects
  samples <- dataList %>%
    purrr::map(., ~.x[['samples']]) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(plateID = forcats::fct(as.character(plateID))) %>%
    dplyr::select(dplyr::any_of(shared_sample_columns)) %>%
    tibble::as_tibble() 
  
  # Combine the qcSample objects
  qcSample <- dataList %>%
    purrr::map(., ~.x[["qcSample"]]) %>%
    dplyr::bind_rows() %>%
    tibble::as_tibble() %>%
    mutate(plateID = forcats::fct(as.character(plateID))) %>%
    dplyr::relocate(plateID)
  
  # Combine the qcSample objects
  qcTargets_exists <- "qcTarget" %in% unique(unlist(lapply(dataList, names)))
  if(qcTargets_exists){
    qcTarget <- tryCatch({
      dataList %>%
        purrr::map(., ~.x[["qcTarget"]]) %>%
        bind_rows() %>%
        tibble::as_tibble() %>%
        mutate(plateID = forcats::fct(as.character(plateID))) %>%
        dplyr::relocate(plateID)
    }, error = function(e){
      logger::log_info("qcTarget object found, but is NULL")
      NULL
    })
  } else {
    qcTarget <- NULL
  }
  
  # merge AQ_params if it exists
  target_params_exists <- "aqParams" %in% unique(unlist(lapply(dataList, names)))
  if(target_params_exists){
    logger::log_info("AQ Target params dataframe found")
    aqParams <- tryCatch({
      dataList %>%
        purrr::map(., ~.x[["aqParams"]]) %>%
        bind_rows() %>%
        tibble::as_tibble() %>%
        mutate(plateID = forcats::fct(as.character(plateID))) %>%
        dplyr::relocate(plateID)
    }, error = function(e){
      logger::log_info("AQ target params object found, but is NULL")
      NULL
    })
  } else{
    logger::log_info("AQ Target params dataframe not found")
    aqParams <- NULL 
  }
  
  # Combine the qcPlates objects
  qcPlate <- dataList %>%
    purrr::map(., ~.x[["qcPlate"]]) %>%
    bind_rows() %>%
    tibble::as_tibble() %>%
    dplyr::relocate(plateID)
  
  if(!is.factor(qcPlate$plateID)) {
    qcPlate$plateID <- forcats::fct(as.character(qcPlate$plateID))
  }
  
  # Initiate an empty list for all the matrix like objects
  dataMatrix <- list()
  for (i in c("Data_IC","Data_IClog2","Data_raw","Data_rawlog2","aboveLOD","Data_AQ",
              "Data_AQlog2", "Data_Reverse","Data_Reverselog2", "Data_AQ_pgmL", "Data_AQlog2_pgmL")) {
    logger::log_info("Merging ProcessedNULISAseq Data Matrix -- ",i)
    tryCatch({
      # Generate a list object for all the new/existing data 
      dataObj <- lapply(dataList, function(x){
        x[[i]] %>%
          tibble::as_tibble(rownames = "targetName")
      })
      # merge the list objects together and keep only the common targets and preserve the matrix format
      dataMatrix[[i]] <- dataObj %>%
        purrr::reduce(dplyr::full_join, by = "targetName") %>%
        filter(targetName %in% targets$targetName) %>%
        tibble::column_to_rownames("targetName") %>%
        as.matrix()
    }, error = function(e){
      logger::log_warn("Data Matrix -- ", i, " not found... skipped")
    })
  }
  
  IC <- unique(existing_data$IC, new_data$IC)
  
  if(length(IC) > 1){
    logger::log_error("Non-Unique IC targets")
    stop("Non-Unique IC targets")
  }
  
  # Return the data object
  c(
    list(
      plateID=plateID,
      fileNames=fileNames,
      ExecutionDetails=ExecutionDetails,
      RunSummary=RunSummary,
      IC = IC,
      targets=targets,
      samples=samples,
      qcSample=qcSample,
      qcTarget=qcTarget,
      qcPlate=qcPlate,
      aqParams=aqParams,
      excluded = excluded
    ),
    dataMatrix, # This is a named list of the matrices
    unit = unit
  )
}



##' Add plateID column to a data frame within a list
#'
#' This helper function adds or updates a `plateID` column in a specified data frame element within a list, ensuring the column is a factor. Used internally for merging NULISAseq data objects.
#'
#' @param x A list containing data frames as elements.
#' @param element Character. The name of the element in `x` to which the plateID should be added.
#' @param plateID Character or factor. The plate ID(s) to assign to the data frame.
#'
#' @return The data frame with an added or updated `plateID` column (as a factor).
#'
#' @keywords internal
add_plateID <- function(x, element, plateID) {
  data <- x[[element]]
  if(is.factor(plateID)){
    data$plateID <- plateID
  } else {
    data$plateID <- forcats::fct(plateID)
  }
  data
}

replace_empty_numeric <- function(x) {
  if(identical(x, numeric(0))) 0 else x
}

##' Combine target information across multiple NULISAseq data objects
#'
#' This helper function merges the target data frames from a list of processed NULISAseq data objects, returning a unified targets data frame and a vector of any excluded targets (those not present in all objects). Optionally, it adds a plateID column to each target data frame before merging.
#'
#' @param dataList A list of processed NULISAseq data objects (as returned by mergeNULISAseq), each containing a `targets` data frame.
#' @param plateID A character or factor vector of plate IDs for the data objects. Only required if `addPlateID = TRUE`.
#' @param addPlateID Logical. If TRUE (default), adds a `plateID` column to each targets data frame before merging.
#'
#' @return A list with two elements:
#'   \item{targets}{A unified data frame of all targets present in the input objects, with columns from the original targets data frames.}
#'   \item{excluded}{A character vector of target names that were excluded because they were not present in all input objects, or NULL if none.}
#'
#' @details
#' - All targets present in all input objects are retained in the merged targets data frame.
#' - If some targets are missing from one or more objects, their names are returned in the `excluded` vector and a WARN is logged.
#' - Used internally for merging and aligning data matrices across plates.
#'
#' @keywords internal
combine_targets <- function(dataList, plateID, addPlateID = TRUE){
  
  target_ids <- lapply(dataList, function(x){
    x$targets$targetName
  })
  
  target_intersect <- Reduce(intersect, target_ids)
  # print a WARN if some targets do not match
  n_targets <- sapply(target_ids, length)
  
  # MODIFY TO INCLUDE NUMBER OF EXCLUDED TARGETS IN WARNING
  excluded <- NULL
  if (sum(n_targets > length(target_intersect)) > 0){
    logger::log_warn("Some plates have targets that do not match those on other plates. These targets will be excluded.")
    target_ids <- unique(unname(do.call("c", target_ids)))
    excluded <- target_ids[!target_ids %in% target_intersect]
  }
  # save target data.frames in a list
  if(addPlateID){
    targets <- mapply(add_plateID, dataList, "targets", plateID, SIMPLIFY = F)
  } else {
    targets <- lapply(dataList, "[[", "targets")
  }
  #create a dataframe of the targets
  targets <- targets %>% 
    bind_rows() %>%
    unique() %>%
    mutate(targetType = tolower(targetType))
  #filter(targetName %in% target_intersect)
  
  list(
    "targets" = targets,
    "excluded" = excluded
  )
}

##' Check compatibility of assay types across NULISAseq data objects
#'
#' This helper function checks that all provided NULISAseq data objects have the same assay type, using the `Assay` field in their ExecutionDetails. If multiple assay types are found, an error is raised.
#'
#' @param ExecutionDetails A list of ExecutionDetails objects (one per plate), each containing an `Assay` field.
#' @param fileNames Optional character vector of file names for each assay, used for error reporting.
#'
#' @return Logical TRUE if all assay types are compatible (identical). Raises an error if not.
#'
#' @details
#' - Used internally to ensure only compatible NULISAseq data objects are merged.
#' - If incompatible, an informative error message is generated listing the file names and assay types.
#'
#' @keywords internal
check_assay_type <- function(ExecutionDetails, fileNames = NA){
  
  assays <- unname(do.call("c",lapply(ExecutionDetails, function(x){
    val <- x$Assay
    ifelse(is.null(val), NA, val)
  })))
  
  assay <- unique(assays)
  
  if(length(assay) > 1){
    if(any(is.na(fileNames))){
      fileNames <- names(ExecutionDetails)
    }
    
    msg <- paste("Found multiple assay types ", 
                 paste(
                   sprintf("%s - %s", fileNames, assays),
                   collapse = ", "
                 ))
    FALSE
    stop(msg)
  } else {
    message(format(Sys.time()),": INFO Assays Compatible")
    TRUE
  }
}

##' Check for duplicate plate IDs, file names, or sample names between two NULISAseq data objects
#'
#' This helper function checks whether there are any duplicate plate IDs, file names, or sample names between two processed NULISAseq data objects. Used to prevent accidental merging of non-unique data.
#'
#' @param existing_data A processed NULISAseq data object (as returned by mergeNULISAseq) representing the existing dataset.
#' @param new_data A processed NULISAseq data object (as returned by mergeNULISAseq) representing the new dataset to be checked.
#' @param type Character. The type of identifier to check for duplicates: one of `"plateID"`, `"fileNames"`, or `"sampleName"`.
#'
#' @return Logical TRUE if any duplicates are found, FALSE otherwise.
#'
#' @details
#' - Used internally to ensure unique merging of NULISAseq data objects.
#' - If `type = "sampleName"`, checks the `sampleName` column in the `samples` data frame of each object.
#'
#' @keywords internal
check_duplicates <- function(existing_data, new_data, type = "plateID"){
  if(type != "sampleName"){
    old <- existing_data[[type]]
    new <- new_data[[type]]
  } else {
    old <- existing_data[["samples"]][["sampleName"]]
    new <- new_data[["samples"]][["sampleName"]]
  }
  
  any(old %in% new)
}

##' Fix duplicate plate IDs and sample names in a NULISAseq data object
#'
#' This helper function reconstructs a processed NULISAseq data object by renaming duplicate plate IDs with their respective file names and updating sample names to ensure uniqueness. It is used internally when merging data objects to prevent identifier collisions.
#'
#' @param new_object A processed NULISAseq data object (as returned by mergeNULISAseq) that may contain duplicate plate IDs or sample names.
#' @param existingSamples Character vector of sample names already present in the existing dataset, used to identify duplicates.
#'
#' @return The input data object with plate IDs replaced by their respective file names and sample names updated for uniqueness. All relevant data frames and matrices are updated accordingly.
#'
#' @details
#' - Updates the `plateID` field and all relevant data frames and matrices in the object.
#' - Appends the plate ID to duplicate sample names to ensure uniqueness.
#' - Used internally during merging to avoid identifier conflicts.
#'
#' @keywords internal
fix_duplicate_plates <- function(new_object, existingSamples){
  
  plateIDs <- new_object$plateID
  file_name <- new_object$fileNames
  names_df <- tibble(
    plateID = plateIDs,
    fileName = file_name
  )
  
  objects_to_fix <- names(new_object)
  samplesToRename <- existingSamples[which(existingSamples %in% new_object$samples$sampleName)]
  
  for (i in objects_to_fix){
    message(format(Sys.time()),": INFO Renaming Plates in new_data, object - ", i)
    if(i == "plateID"){
      new_object[[i]] <- names_df$fileName
    }
    
    if(i == "ExecutionDetails"){
      names(new_object[[i]]) <- names_df$fileName
    }
    
    if("data.frame" %in% class(new_object[[i]])){
      
      new_object[[i]] <- new_object[[i]] %>%
        left_join(names_df, by = "plateID") %>%
        mutate(plateID = fileName) %>%
        select(-fileName) %>%
        dplyr::relocate("plateID")
      
      if(i %in% c("samples","qcSample")){
        new_object[[i]] <- new_object[[i]] %>%
          mutate(sampleName = dplyr::if_else(sampleName %in% samplesToRename, paste(sampleName, plateID, sep = "_"), sampleName))
      }
    } 
    
    if(all(class(new_object[[i]]) %in% c("matrix","array"))){
      sampleName_df <- new_object[["samples"]][,c("sampleID","sampleName")]
      new_object[[i]] <- new_object[[i]] %>%
        tibble::as_tibble(rownames = "targetName") %>%
        data.table::setnames(., 
                             old = sampleName_df$sampleID, 
                             new = sampleName_df$sampleName,
                             skip_absent = T) %>%
        column_to_rownames("targetName") %>%
        as.matrix()
    } 
  }
  new_object[["samples"]] <- new_object[["samples"]] %>%
    mutate(sampleID = sampleName)
  message(format(Sys.time()),": INFO Duplicate plates fixed in the new_data object")
  new_object
}

##' Clean and standardize a NULISAseq samples data frame
#'
#' This helper function processes a samples data frame from NULISAseq data, converting columns to appropriate types, cleaning covariate levels, renaming columns for consistency, and ensuring required identifiers are present. Used internally for data preparation and merging.
#'
#' @param df A data frame of sample metadata, typically from NULISAseq data.
#' @param numericCovariates Named logical vector indicating which columns should be treated as numeric (TRUE) or factor/character (FALSE).
#'
#' @return A cleaned and standardized data frame with:
#'   - Columns converted to appropriate types (numeric, factor, character)
#'   - Covariate levels cleaned of punctuation
#'   - Standardized column names (e.g., wellRow, wellCol, plateGUID, plateID)
#'   - A `sampleID` column ensuring unique sample identifiers
#'
#' @details
#' - Handles both current and legacy XML column names.
#' - Used internally by process_loadNULISAseq and related functions.
#'
#' @keywords internal
process_samples_df <- function(df, numericCovariates){
  # format columns based on provided format types in the processed data
  curr_type <- stack(sapply(df, class))
  req_type <- stack(numericCovariates)
  not_req <- c("sampleType", "sampleName", "sampleID", "sampleName_long")
  types <- curr_type %>%
    left_join(req_type %>%
                mutate(ptype = dplyr::if_else(values, "numeric", "factor"),
                       ptype = dplyr::if_else(ind %in% not_req, "character", ptype)) %>%
                select(-values), by = "ind") %>%
    mutate(ptype = dplyr::if_else(is.na(ptype) & values == "character", 
                           "factor", ptype),
           ptype = dplyr::if_else(ind == "plateID", "factor", ptype),
           ptype = substr(ptype,1,1)) %>%
    pull(ptype) %>%
    paste(collapse = "")
  
  # sanitize the samples dataframe of any punctuations marks
  colsToSanitize <-setdiff(names(df), protected_columns())
  logger::log_info("Columns to sanitize: ", paste(colsToSanitize, collapse = ", "))
  
  # convert the columns to the required types
  df <- readr::type_convert(df, col_types = types)
  
  # clean the covariate levels
  df <- df %>%
    mutate(across(all_of(colsToSanitize), clean_covariate_levels))
  
  # rename some columns 
  names(df)[names(df) == "AUTO_WELLROW"] <- "wellRow"
  names(df)[names(df) == "AUTO_WELLCOL"] <- "wellCol"
  names(df)[names(df) == "AUTO_PLATE"] <- "plateGUID"
  
  # work with older versions of xml
  if(!("plateID" %in% names(df))){
    names(df)[names(df) == "Annot1"] <- "plateID"
  }
  if(!("wellRow" %in% names(df))){
    names(df)[names(df) == "Annot2"] <- "wellRow"
  }
  if(!("wellCol" %in% names(df))){
    names(df)[names(df) == "Annot3"] <- "wellCol"
  }
  
  # edit samples dataframe 
  # add sample id (now it"s just the sampleName column)
  # This will come from NULISAseqR package now
  if(!"sampleID" %in% colnames(df)){
    df <- df %>%
      mutate(sampleID = dplyr::if_else(sampleType != "Sample", 
                                paste(plateID, sampleName, sep = "_"),
                                sampleName))
  }
  
  df
}

#' function to get the protected columns in the covariate data
#' This function returns a vector of column names that are protected and should not be cleaned
#' the column are created by ACC or upstream process and cannot be user edited

##' List of protected columns in NULISAseq sample data
#'
#' This helper function returns a character vector of column names that should not be sanitized or altered during data cleaning and processing. These columns are essential identifiers or metadata fields in NULISAseq sample data.
#'
#' @return Character vector of protected column names.
#'
#' @details
#' - Used internally to prevent renaming or cleaning of key columns during sample data processing.
#' - Includes plate, sample, and experiment identifiers, as well as legacy and alternative naming conventions.
#' - Used Primarily in NAS for internal control of the samples dataframe
#'
#' @keywords internal
protected_columns <- function(){
  c('plateGUID','sampleType','sampleBarcode','wellRow','wellCol','matching','non-matching','AUTO_WELLPOSITION','CONDITION_1','EXPERIMENT_ID','OPERATOR','SAMPLE_MATRIX','SEQ_INSTRUMENT_ID','SEQ_RUN_DATE','sampleID','plateID','sampleName','sampleName_original','AUTO_WELLCOL','AUTO_PLATE','AUTO_WELLROW','EXPERIMENT.ID','SAMPLE.MATRIX','SEQ.INSTRUMENT.ID','SEQ.RUN.DATE')
}

##' Clean and standardize covariate levels in a vector
#'
#' This helper function cleans the levels of a covariate vector by applying `janitor::make_clean_names` to non-numeric values, replacing empty strings with `NA`, and correcting cases where cleaned levels are incorrectly prefixed with "x". Numeric vectors are returned unchanged; all other vectors are returned as factors with cleaned, standardized levels.
#'
#' @param x A vector of covariate values (character, factor, or numeric).
#'
#' @return
#' If `x` is numeric, returns `x` unchanged. Otherwise, returns a factor with cleaned levels:
#'   - Empty strings are replaced with `NA`.
#'   - Non-missing values are cleaned using `janitor::make_clean_names` (with `allow_dupes = TRUE` and `case = 'lower_camel'`).
#'   - If any cleaned level starts with "x" followed by a digit, but the original value did not, the leading "x" is removed.
#'
#' @details
#' - Used internally to ensure covariate levels are consistent and safe for use as column names or factor levels.
#' - Prevents issues with R variable naming and downstream data processing.
#'
#' @examples
#' clean_covariate_levels(c("A-1", "B 2", "", NA, "3rd"))
#'
#' @keywords internal
#' This function is useful for standardizing categorical covariate levels for analysis or modeling,
#' especially when original levels may contain spaces, special characters, or begin with digits.
clean_covariate_levels <- function(x){
  val <- trimws(as.character(x), which = 'both')
  
  if(!is.numeric(x)){
    x <- trimws(x, which = 'both')
    val[val==''] <- NA
    val <- ifelse(
      is.na(val), 
      val,
      janitor::make_clean_names(val, allow_dupes = TRUE, case = 'lower_camel')
    )
    # check if the cleaned levels begin with `x`
    if(any(grepl("^x[[:digit:]]", as.character(val)))){
      # get the index of the new levels which begin with `x` and followed by digit
      idx <- which(grepl("^x[[:digit:]]", as.character(val)))
      
      # check if any of the original levels begin with the `x` followed by digit pattern
      # if not remove the beginning x from the renamed levels using the index from the original vector
      if(!any(grepl("^x[[:digit:]]", as.character(x[idx])))){
        orig_idx <- which(grepl("^[[:digit:]]", as.character(x)))
        val[orig_idx] <- gsub("^x", "", val[orig_idx])
        as.factor(val)
      }else{
        as.factor(val)
      }
    } else{
      as.factor(val)
    }
  }else{
    x
  }
}

##' Test if all (finite, non-missing) values in a numeric vector are (nearly) constant
#'
#' This helper function checks whether all finite, non-missing values in a numeric vector are equal within a specified tolerance. Used to detect zero-variance columns or features.
#'
#' @param x Numeric vector to test.
#' @param tol Numeric tolerance for equality (default: square root of machine epsilon).
#'
#' @return Logical TRUE if all (finite, non-missing) values are equal within tolerance, FALSE otherwise.
#'
#' @details
#' - Ignores NA, NaN, and infinite values.
#' - Returns TRUE for vectors of length 1.
#' - Used internally for data quality checks and filtering.
#' - Primarily used in NAS to determine viability of columns as covariates in the samples dataframe
#'
#' @examples
#' zero_range(c(1, 1, 1))        # TRUE
#' zero_range(c(1, 1.0000001))   # TRUE (within tolerance)
#' zero_range(c(1, 2, 1))        # FALSE
#'
#' @keywords internal
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- x[!is.na(x) & !is.infinite(x) & !is.nan(x)]
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

#' This helper function applies `janitor::make_clean_names` to a character vector of column names, except for those in the protected columns list. Used to ensure column names are syntactically valid and consistent for downstream processing.
#'
#' @param x Character vector of column names to clean.
#' @param ... Additional arguments PASS to `janitor::make_clean_names` (e.g., `case = 'lower_camel'`).
#'
#' @return Character vector of cleaned column names, with protected columns left unchanged.
#'
#' @details
#' - Protected columns (see `protected_columns()`) are not altered.
#' - Used internally to standardize column names in NULISAseq sample data frames.
#'
#' @examples
#' clean_covariate_names(c("Sample Name", "plateID", "Well Row"), case = "lower_camel")
#'
#' @keywords internal
clean_covariate_names <- function(x, ...){
  dont_clean <- protected_columns()
  to_clean_idx <- which(!x %in% dont_clean)
  logger::log_info("Columns to clean - {paste(x[to_clean_idx], collapse = ', ')}")
  x[to_clean_idx] <- janitor::make_clean_names(x[to_clean_idx], ...)
  logger::log_info("New column names - {paste(x[to_clean_idx], collapse = ', ')}")
  x
}

#' A function to rename the column names of a matrix given a dataframe of old and new
#' names
#' @param mat A named (rownames/column names) matrix
#' @param names_df A two column dataframe with sampleName and sampleID as the old and new column name
#' respectively
#' 
#' @return A matrix with the updated columns names
#' @importFrom data.table setnames
#' @keywords internal
rename_cols <- function(mat, names_df){
  mat %>%
    tibble::as_tibble(rownames = 'rn') %>%
    data.table::setnames(., old = names_df$sampleName, new = names_df$sampleID,
                         skip_absent = TRUE) %>%
    column_to_rownames('rn') %>%
    as.matrix()
}

#' A function that converts a matrix or data frame to long format by pivoting
#' sample columns into rows.
#'
#' @param data_matrix A matrix or data frame with IDs/targets as row names and
#' samples as columns.
#' @param id_col Character string specifying the name for the ID column
#' created from the row names (Default is "Target").
#' @param sample_name_col Character string specifying the name for the new
#' column that will contain the original sample column names (Default is "SampleName").
#' @param data_col Character string specifying the name for the values column
#' in the output (e.g., "Expression" or "Value"). This is a required argument.
#'
#' @return A data frame in long format with columns: ID, Sample Name, and Value.
#'
#' @keywords internal
convert_to_long <- function(data_matrix, id_col = "Target", sample_name_col = "SampleName", data_col) {
  data_matrix %>%
    tibble::as_tibble(rownames = id_col) %>% 
    tidyr::pivot_longer(
      cols = -dplyr::all_of(id_col),
      names_to = sample_name_col,
      values_to = data_col
    )
}


#' Prepare Sample QC Data for Long Format
#'
#' Transforms qcSample data from long to wide format suitable for joining
#' with the long format data output.
#'
#' @param qc_sample Data frame with columns: sampleName, plateID, flagName, val, status
#'   The status column can be logical (TRUE/FALSE), numeric (0/non-zero), or character.
#'   Character values are case-insensitive. Common passing values: "OK", "PASS", "PASSED", "FALSE", "0".
#'   Common warning values: "TRUE", "1", "FAIL", "FAILED", "FAILURE".
#'
#' @return A tibble with columns SampleName, PlateID, Sample_QC_Status ("PASS" or "WARN"), and
#'   Sample_QC_<metric> and Sample_QC_<metric>_Status ("PASS" or "WARN") for each QC metric
#'
#' @details Column naming convention: Sample_QC_<MetricName> and Sample_QC_<MetricName>_Status.
#'   Status values are "PASS" or "WARN" (uppercase).
#'   Overall Sample_QC_Status is "PASS" only if all individual metrics pass; otherwise "WARN".
#'
#' @keywords internal
prepare_sample_qc_for_long <- function(qc_sample) {
  
  required_cols <- c("sampleName", "plateID", "flagName", "val", "status")
  missing_cols <- setdiff(required_cols, names(qc_sample))
  if (length(missing_cols) > 0) {
    stop("qcSample missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (is.character(qc_sample$status)) {
    logger::log_info("Converting character status to logical (TRUE/non-zero = warning, FALSE/zero = passed)")
    # Convert to uppercase for case-insensitive matching
    # TRUE values (warning): "TRUE", "1", "FAIL", "FAILED", "FAILURE", or anything not "OK"/"PASS"
    status_upper <- toupper(trimws(qc_sample$status))
    qc_sample$status <- status_upper %in% c("TRUE", "1", "FAIL", "FAILED", "FAILURE") |
                        (!status_upper %in% c("OK", "PASS", "PASSED", "FALSE", "0", ""))
  } else if (is.numeric(qc_sample$status)) {
    logger::log_info("Converting numeric status to logical (non-zero = warning, zero = passed)")
    qc_sample$status <- as.logical(qc_sample$status)
  } else if (!is.logical(qc_sample$status)) {
    stop("qcSample 'status' column must be logical, numeric, or character, got: ", class(qc_sample$status)[1])
  }
  
  key_cols <- c("sampleName", "plateID", "flagName")
  dupes <- qc_sample %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(key_cols))) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup()
  
  if (nrow(dupes) > 0) {
    logger::log_warn(sprintf("Found %d duplicate QC sample entries, using first value", nrow(dupes)))
    qc_sample <- qc_sample %>%
      dplyr::distinct(dplyr::across(dplyr::all_of(key_cols)), .keep_all = TRUE)
  }
  
  qc_sample_clean <- qc_sample %>%
    dplyr::mutate(flagName = sub("^[Ss]ample_", "", flagName))
  
  sample_qc_wide <- qc_sample_clean %>%
    dplyr::select(sampleName, plateID, flagName, val, status) %>%
    tidyr::pivot_wider(
      id_cols = c(sampleName, plateID),
      names_from = flagName,
      values_from = c(val, status),
      names_glue = "Sample_QC_{flagName}_{.value}"
    )
  
  status_cols <- names(sample_qc_wide)[grepl("_status$", names(sample_qc_wide))]
  
  sample_qc_wide <- add_overall_qc_status(sample_qc_wide, status_cols, "Sample_QC_Status")
  sample_qc_wide <- convert_status_to_labels(sample_qc_wide, status_cols)
  
  names(sample_qc_wide) <- gsub("_val$", "", names(sample_qc_wide))
  names(sample_qc_wide) <- gsub("_status$", "_Status", names(sample_qc_wide))
  
  sample_qc_wide <- sample_qc_wide %>%
    dplyr::rename(SampleName = sampleName, PlateID = plateID)
  
  col_order <- c("SampleName", "PlateID", "Sample_QC_Status",
                 setdiff(names(sample_qc_wide), c("SampleName", "PlateID", "Sample_QC_Status")))
  sample_qc_wide %>% dplyr::select(dplyr::all_of(col_order))
}


#' Prepare Target QC Data for Long Format
#'
#' Transforms qcTarget data from long to wide format suitable for joining
#' with the long format data output.
#'
#' @param qc_target Data frame with columns: target, plateID, flagName, val, status (optional)
#'   The status column can be logical (TRUE/FALSE), numeric (0/non-zero), or character.
#'   Character values are case-insensitive. Common passing values: "OK", "PASS", "PASSED", "FALSE", "0".
#'   Common warning values: "TRUE", "1", "FAIL", "FAILED", "FAILURE".
#'
#' @return A tibble with columns Target, PlateID, Target_QC_Status ("PASS" or "WARN"), and
#'   Target_QC_<metric> and Target_QC_<metric>_Status ("PASS" or "WARN") for each QC metric.
#'   If status column is missing, Target_QC_Status will be NA_character_.
#'
#' @details Column naming convention: Target_QC_<MetricName> and Target_QC_<MetricName>_Status.
#'   Status values are "PASS" or "WARN" (uppercase).
#'
#' @keywords internal
prepare_target_qc_for_long <- function(qc_target) {
  
  required_cols <- c("target", "plateID", "flagName", "val")
  missing_cols <- setdiff(required_cols, names(qc_target))
  if (length(missing_cols) > 0) {
    stop("qcTarget missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  has_status <- "status" %in% names(qc_target)
  
  if (has_status) {
    if (is.character(qc_target$status)) {
      logger::log_info("Converting character status to logical (TRUE/non-zero = warning, FALSE/zero = passed)")
      # Convert to uppercase for case-insensitive matching
      status_upper <- toupper(trimws(qc_target$status))
      qc_target$status <- status_upper %in% c("TRUE", "1", "FAIL", "FAILED", "FAILURE") |
                          (!status_upper %in% c("OK", "PASS", "PASSED", "FALSE", "0", ""))
    } else if (is.numeric(qc_target$status)) {
      logger::log_info("Converting numeric status to logical (non-zero = warning, zero = passed)")
      qc_target$status <- as.logical(qc_target$status)
    } else if (!is.logical(qc_target$status)) {
      stop("qcTarget 'status' column must be logical, numeric, or character, got: ", class(qc_target$status)[1])
    }
  }
  
  key_cols <- c("target", "plateID", "flagName")
  dupes <- qc_target %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(key_cols))) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup()
  
  if (nrow(dupes) > 0) {
    logger::log_warn(sprintf("Found %d duplicate QC target entries, using first value", nrow(dupes)))
    qc_target <- qc_target %>%
      dplyr::distinct(dplyr::across(dplyr::all_of(key_cols)), .keep_all = TRUE)
  }
  
  qc_target_clean <- qc_target %>%
    dplyr::mutate(flagName = sub("^[Tt]arget_", "", flagName))
  
  if (has_status) {
    target_qc_wide <- qc_target_clean %>%
      dplyr::select(target, plateID, flagName, val, status) %>%
      tidyr::pivot_wider(
        id_cols = c(target, plateID),
        names_from = flagName,
        values_from = c(val, status),
        names_glue = "Target_QC_{flagName}_{.value}"
      )
    
    names(target_qc_wide) <- gsub("_val$", "", names(target_qc_wide))
    names(target_qc_wide) <- gsub("_status$", "_Status", names(target_qc_wide))
    
    status_cols <- names(target_qc_wide)[grepl("_Status$", names(target_qc_wide))]
    target_qc_wide <- add_overall_qc_status(target_qc_wide, status_cols, "Target_QC_Status")
    target_qc_wide <- convert_status_to_labels(target_qc_wide, status_cols)
  } else {
    target_qc_wide <- qc_target_clean %>%
      dplyr::select(target, plateID, flagName, val) %>%
      tidyr::pivot_wider(
        id_cols = c(target, plateID),
        names_from = flagName,
        values_from = val,
        names_glue = "Target_QC_{flagName}"
      )
    target_qc_wide$Target_QC_Status <- NA_character_
  }
  
  target_qc_wide <- target_qc_wide %>%
    dplyr::rename(Target = target, PlateID = plateID)
  
  col_order <- c("Target", "PlateID", "Target_QC_Status",
                 setdiff(names(target_qc_wide), c("Target", "PlateID", "Target_QC_Status")))
  target_qc_wide %>% dplyr::select(dplyr::all_of(col_order))
}


#' Calculate Overall QC Status from Individual Status Columns
#'
#' @param data Data frame containing status columns
#' @param status_cols Character vector of status column names to aggregate
#' @param status_col_name Name for the output status column
#'
#' @return Data frame with added overall status column
#' @keywords internal
add_overall_qc_status <- function(data, status_cols, status_col_name) {
  if (length(status_cols) > 0) {
    status_matrix <- dplyr::select(data, dplyr::all_of(status_cols))
    
    data[[status_col_name]] <- dplyr::if_else(
      rowSums(status_matrix, na.rm = TRUE) == 0,
      "PASS",
      "WARN"
    )
  } else {
    data[[status_col_name]] <- NA_character_
  }
  data
}


#' Convert Boolean Status Columns to Character Labels
#'
#' @param data Data frame containing boolean status columns
#' @param status_cols Character vector of status column names to convert
#'
#' @return Data frame with status columns converted to character
#' @keywords internal
convert_status_to_labels <- function(data, status_cols) {
  if (length(status_cols) == 0) return(data)
  
  for (col in status_cols) {
    if (col %in% names(data)) {
      data[[col]] <- dplyr::if_else(
        as.logical(data[[col]]),
        "WARN",
        "PASS"
      )
    }
  }
  data
}

#' Filter Target QC columns based on data mode (RQ vs AQ)
#'
#' Removes AQ-specific Target QC metrics when processing RQ data. This ensures that
#' RQ long format output only contains Target QC metrics appropriate for RQ data,
#' avoiding confusion from AQ-specific metrics like concentration accuracy and CV
#' that are calculated on AQ values and not applicable to RQ data.
#'
#' @param data Data frame containing Target QC columns
#' @param AQ Logical indicating whether this is AQ data (TRUE) or RQ data (FALSE)
#'
#' @return Data frame with appropriate Target QC columns for the data mode
#'
#' @details
#' Uses QCTargetCriteria() to dynamically determine which Target QC metrics are
#' AQ-specific by comparing criteria for AQ=TRUE vs AQ=FALSE. AQ-specific metrics
#' are automatically excluded from RQ data output.
#'
#' @keywords internal
filter_target_qc_by_mode <- function(data, AQ = FALSE) {
  # If AQ mode, keep all Target QC columns
  if (AQ) {
    return(data)
  }

  # For RQ mode, identify AQ-specific metrics by calling QCTargetCriteria
  rq_criteria <- QCTargetCriteria(AQ = FALSE, advancedQC = TRUE)
  aq_criteria <- QCTargetCriteria(AQ = TRUE, advancedQC = TRUE)

  # Validate QCTargetCriteria returned expected structure
  if (is.null(rq_criteria) || is.null(rq_criteria$thresholds)) {
    logger::log_warn("QCTargetCriteria returned NULL for RQ mode, no filtering applied")
    return(data)
  }
  if (is.null(aq_criteria) || is.null(aq_criteria$thresholds)) {
    logger::log_warn("QCTargetCriteria returned NULL for AQ mode, no filtering applied")
    return(data)
  }

  # Find metrics that exist in AQ but not in RQ (these are AQ-specific)
  rq_metrics <- names(rq_criteria$thresholds)
  aq_metrics <- names(aq_criteria$thresholds)
  aq_specific_metrics <- setdiff(aq_metrics, rq_metrics)

  if (length(aq_specific_metrics) == 0) {
    logger::log_debug("No AQ-specific metrics to filter for RQ data")
    return(data)
  }

  # Pre-filter to Target QC columns only for better performance
  all_cols <- names(data)
  target_qc_cols <- all_cols[grepl("^Target_QC_", all_cols)]

  if (length(target_qc_cols) == 0) {
    return(data)  # Early return if no Target QC columns exist
  }

  # Convert metric names to column name patterns
  # Target_Conc_Accuracy -> Target_QC_Conc_Accuracy, Target_QC_Conc_Accuracy_Status
  cols_to_remove <- character(0)

  for (metric in aq_specific_metrics) {
    # Remove "Target_" prefix if present to match QC column naming
    metric_suffix <- sub("^Target_", "", metric)

    # Defensive check for empty suffix
    if (nchar(metric_suffix) == 0) {
      logger::log_warn("Empty metric suffix after removing Target_ prefix from: {metric}")
      next
    }

    # Escape special regex characters to prevent regex injection
    metric_suffix_escaped <- gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", metric_suffix)

    # Match both the value column and _Status column
    pattern <- paste0("^Target_QC_", metric_suffix_escaped, "(_Status)?$")
    matching_cols <- target_qc_cols[grepl(pattern, target_qc_cols)]
    cols_to_remove <- c(cols_to_remove, matching_cols)
  }

  # Remove AQ-specific columns if any were found
  if (length(cols_to_remove) > 0) {
    data <- data %>% dplyr::select(-dplyr::any_of(cols_to_remove))
    logger::log_info(
      "Filtered out AQ-specific Target QC columns for RQ data: {paste(cols_to_remove, collapse = ', ')}"
    )
  }

  return(data)
}

#' A helper function that converts NULISAseq data from wide format (matrices) to long format (tidy data)
#' by joining sample metadata, target metadata, and execution details.
#'
#' @param merged A list containing NULISAseq data with required elements:
#'   \itemize{
#'     \item{\code{targets} - Target metadata data frame}
#'     \item{\code{samples} - Sample metadata data frame}
#'     \item{\code{ExecutionDetails} - Execution details list}
#'     \item{\code{Data_Reverselog2} or \code{Data_NPQ} - NPQ data matrix for RQ}
#'     \item{\code{Data_raw} - Raw count matrix for RQ}
#'     \item{\code{Data_AQ} or \code{Data_AQ_aM} - Absolute quantification matrix in aM (optional)}
#'     \item{\code{Data_AQ_pgmL} - Absolute quantification matrix in pg/mL (optional)}
#'     \item{\code{Data_AQlog2} or \code{Data_AQlog2_aM} - log2 Absolute quantification matrix in aM (optional)}
#'     \item{\code{Data_AQlog2_pgmL} - log2 Absolute quantification matrix in pg/mL (optional)}
#'   }
#' @param AQ Logical indicating whether to process AQ data.
#'   Defaults to \code{FALSE}.
#' @param exclude_sample_cols Character vector of column names to exclude from
#'   sample metadata. Defaults to \code{"plateGUID"}.
#' @param include_qc Logical indicating whether to include QC columns in the output.
#'   Defaults to \code{TRUE}. When TRUE, adds Sample_QC_* and Target_QC_* columns if
#'   qcSample and qcTarget data are available. Target QC columns are automatically
#'   filtered based on the \code{AQ} parameter: RQ data excludes AQ-specific metrics
#'   like concentration accuracy and CV calculated on AQ values.
#'
#' @return A tibble in long format with the following columns (depending on mode):
#'   \itemize{
#'     \item{Panel, PanelLotNumber, PlateID, InstrumentID, InstrumentBay}
#'     \item{SampleName, SampleType}
#'     \item{Target, AlamarTargetID, UniProtID, ProteinName}
#'     \item{For AQ : Conc_pgmL, log2Conc_pgmL, LOD_pgmL, LLOQ_pgmL, ULOQ_pgmL, Conc_aM, log2Conc_aM, LOD_aM, LLOQ_aM, ULOQ_aM, withinDR}
#'     \item{For RQ: UnnormalizedCount, NPQ, LOD}
#'     \item{If include_qc = TRUE: Sample_QC_Status, Sample_QC_*, Target_QC_Status, Target_QC_* (filtered by mode)}
#'   }
#' @keywords internal
format_wide_to_long <- function(merged, AQ = FALSE, exclude_sample_cols = "plateGUID", include_qc = TRUE) {
  
  # Input validation
  if (!inherits(merged, "list")) {
    stop("merged must be a list")
  }
  
  required_elements <- c("targets", "samples", "ExecutionDetails")
  missing_elements <- setdiff(required_elements, names(merged))
  if (length(missing_elements) > 0) {
    stop("Missing required elements in merged: ", paste(missing_elements, collapse = ", "))
  }
  
  # Check if AQ data exists when AQ = TRUE
  if (AQ) {
    has_aq_data <- any(grepl("^Data_AQ", names(merged)))
    if (!has_aq_data) {
      stop("No AQ data found in merged object")
    }
    
    # Only process AQ data
    long_data <- NULL
    
    # Check which AQ data column is available
    AQ_aM_data_used <- if("Data_AQ" %in% names(merged)) {"Data_AQ"} 
    else if("Data_AQ_aM" %in% names(merged)) {"Data_AQ_aM"} else {
      stop("Neither 'Data_AQ' nor 'Data_AQ_aM' was found in the merged data. 
           Please ensure that AQ data is present in the input.")
    }
    
    AQ_aM_log2_data_used <- if("Data_AQlog2" %in% names(merged)) {"Data_AQlog2"} 
    else if("Data_AQlog2_aM" %in% names(merged)) {"Data_AQlog2_aM"} else {
      stop("Neither 'Data_AQlog2' nor 'Data_AQlog2_aM' was found in the merged data. 
           Please ensure that AQ log2 data is present in the input.")
    }
    
    data_AQ_aM  <- convert_to_long(data_matrix = merged[[AQ_aM_data_used]], data_col = "Conc_aM")
    data_AQ_pgmL <- convert_to_long(data_matrix = merged$Data_AQ_pgmL, data_col = "Conc_pgmL")
    data_AQ_aM_log2  <- convert_to_long(data_matrix = merged[[AQ_aM_log2_data_used]], data_col = "log2Conc_aM")
    data_AQ_pgmL_log2 <- convert_to_long(data_matrix = merged$Data_AQlog2_pgmL, data_col = "log2Conc_pgmL")

    long_data <- data_AQ_aM %>%
      left_join(data_AQ_pgmL, by = c("Target", "SampleName")) %>%
      left_join(data_AQ_aM_log2, by = c("Target", "SampleName")) %>%
      left_join(data_AQ_pgmL_log2, by = c("Target", "SampleName"))

    # Add withinDR column if available
    if ("withinDR" %in% names(merged) && !is.null(merged$withinDR)) {
      data_withinDR <- convert_to_long(data_matrix = merged$withinDR, data_col = "withinDR")
      long_data <- long_data %>%
        left_join(data_withinDR, by = c("Target", "SampleName"))
    }
    
  } else {
    # Check which NPQ data column is available
    NPQ_data_used <- if("Data_Reverselog2" %in% names(merged)) {
      "Data_Reverselog2"
    } else if("Data_NPQ" %in% names(merged)) {
      "Data_NPQ" 
    } else {
      stop("Neither 'Data_Reverselog2' nor 'Data_NPQ' found in the merged data")
    }
    
    data_long <- convert_to_long(data_matrix = merged[[NPQ_data_used]], data_col ="NPQ")
    data_raw_long <- convert_to_long(data_matrix = merged$Data_raw, data_col = "UnnormalizedCount") 
    
    long_data <- data_raw_long %>% left_join(data_long, by = c("Target", "SampleName"))
  }
  
  # Prepare target metadata based on mode
  target_base_cols <- c("PlateID" = "plateID", "Target" = "targetName", 
                        "AlamarTargetID", "UniProtID", "ProteinName")
  
  if (AQ) {
    target_metadata <- merged$targets %>%
      select(any_of(c(target_base_cols, 
                      "LOD_pgmL", "LLOQ_pgmL" = "targetLLOQ_pg_ml", 
                      "ULOQ_pgmL" = "targetULOQ_pg_ml", "LOD_aM", 
                      "LLOQ_aM" = "targetLLOQ_aM", "ULOQ_aM" = "targetULOQ_aM"))) %>%
      distinct()
  } else {
    target_metadata <- merged$targets %>%
      select(any_of(c(target_base_cols, "LOD" = "logged_LOD"))) %>%
      distinct()
  }
  
  # Prepare sample metadata
  sample_metadata <- merged$samples %>%
    rename(
      PlateID = plateID, 
      SampleName = sampleName,
      SampleType = sampleType
    ) %>%
    select(
      PlateID, SampleName, SampleType, 
      everything(), 
      -any_of("plateGUID"),
      -where(~all(is.na(.)))
    )
  
  # Prepare execution metadata
  execution_df <- tryCatch({
    exec_data <- merged$ExecutionDetails %>%
      tibble::enframe(name = "PlateID", value = "Details") %>%
      tidyr::unnest_wider(Details)
    
    # Safely select columns that exist
    available_cols <- names(exec_data)
    select_cols <- c("PlateID")
    
    if ("Assay" %in% available_cols) select_cols <- c(select_cols, "Panel" = "Assay")
    if ("TargetKitLotNumber" %in% available_cols) select_cols <- c(select_cols, "PanelLotNumber" = "TargetKitLotNumber")
    if ("InstrumentSerial" %in% available_cols) select_cols <- c(select_cols, "InstrumentID" = "InstrumentSerial")
    if ("InstrumentBay" %in% available_cols) select_cols <- c(select_cols, "InstrumentBay")
    
    exec_data %>% select(any_of(select_cols))
  }, error = function(e) {
    warning("Error processing ExecutionDetails: ", conditionMessage(e))
    NULL
  })
  
  # Join all data
  long_data <- long_data %>%
    left_join(sample_metadata, by = "SampleName") %>% 
    left_join(target_metadata, by = c("Target", "PlateID")) 
  
  if (!is.null(execution_df)) {
    long_data <- long_data %>% left_join(execution_df, by = "PlateID")
  }
  
  # Define column order based on mode
  base_cols <- c("Panel", "PanelLotNumber", "PlateID", "InstrumentID", 
                 "InstrumentBay", "SampleName", "SampleType")
  target_cols <- c("Target", "AlamarTargetID", "UniProtID", "ProteinName")
  
  if (AQ) {
    value_cols <- c("Conc_pgmL", "log2Conc_pgmL", "LOD_pgmL", "LLOQ_pgmL", "ULOQ_pgmL",
                    "Conc_aM", "log2Conc_aM", "LOD_aM", "LLOQ_aM", "ULOQ_aM", "withinDR")
  } else {
    value_cols <- c("LOD", "UnnormalizedCount", "NPQ")
  }
  
  # columns which are integers
  integer_cols <- c("UnnormalizedCount", "matching", "non-matching", "wellCol")
  
  # Reorder columns
  final_data <- long_data %>%
    mutate(across(where(is.character), as.factor)) %>% # recode characters columns as factors
    mutate(across(any_of(integer_cols), as.integer)) %>% # recode integerCols as integers from dbl
    dplyr::relocate(
      any_of(base_cols),    
      any_of(target_cols),  
      any_of(value_cols),
      .before = 1
    ) 
  
  
  # Add QC columns if requested and available
  if (include_qc) {
    # Add Sample QC if available
    if ("qcSample" %in% names(merged) && !is.null(merged$qcSample)) {
      logger::log_info("Adding Sample QC columns")
      tryCatch({
        sample_qc_wide <- prepare_sample_qc_for_long(merged$qcSample)
        final_data <- final_data %>%
          dplyr::left_join(sample_qc_wide, by = c("SampleName", "PlateID"))
        logger::log_info("Sample QC columns added successfully")
      }, error = function(e) {
        logger::log_warn("Failed to add Sample QC columns: ", conditionMessage(e))
      })
    }
    
    # Add Target QC if available
    if ("qcTarget" %in% names(merged) && !is.null(merged$qcTarget)) {
      logger::log_info("Adding Target QC columns")
      tryCatch({
        target_qc_wide <- prepare_target_qc_for_long(merged$qcTarget)
        final_data <- final_data %>%
          dplyr::left_join(target_qc_wide, by = c("Target", "PlateID"))
        logger::log_info("Target QC columns added successfully")
      }, error = function(e) {
        logger::log_warn("Failed to add Target QC columns: ", conditionMessage(e))
      })

      # Filter Target QC columns based on data mode (removes AQ-specific metrics from RQ data)
      # Done outside tryCatch so filtering errors are visible
      if ("Target_QC_Status" %in% names(final_data)) {
        final_data <- filter_target_qc_by_mode(final_data, AQ = AQ)
      }
    }

    # Reorder columns to place QC columns appropriately
    all_cols <- names(final_data)
    
    # Sample QC columns: status first, then details
    sample_qc_status <- "Sample_QC_Status"
    sample_qc_details <- all_cols[grepl("^Sample_QC_", all_cols) & all_cols != sample_qc_status]
    
    # Target QC columns: status first, then details
    target_qc_status <- "Target_QC_Status"
    target_qc_details <- all_cols[grepl("^Target_QC_", all_cols) & all_cols != target_qc_status]
    
    # Remaining columns (covariates, etc.)
    priority_cols <- c(base_cols, target_cols, value_cols,
                       sample_qc_status, sample_qc_details,
                       target_qc_status, target_qc_details)
    remaining_cols <- setdiff(all_cols, priority_cols)
    
    # Build final column order (only include columns that exist)
    final_order <- c(
      intersect(base_cols, all_cols),
      intersect(target_cols, all_cols),
      intersect(value_cols, all_cols),
      intersect(sample_qc_status, all_cols),
      intersect(sample_qc_details, all_cols),
      intersect(target_qc_status, all_cols),
      intersect(target_qc_details, all_cols),
      remaining_cols
    )
    
    final_data <- final_data %>% dplyr::select(dplyr::all_of(final_order))
  }
  
  return(final_data)
}

#' An helper function that processes parameters that can be specified
#' as named lists (per plate), vectors (applied to all plates), or NULL
#' 
#' @param param The parameter to process. Can be:
#'   - `NULL`: No action for any plate
#'   - A vector: Applied to all plates
#'   - A named list: Specific values for specific plates, names must match plate_names, 
#'   length of plate_names must match length of list
#' @param param Character string of the parameter 
#' @param plate_names Character vector of plate names to match against
#'
#' @return A named list with one element per plate, where:
#'   - Names correspond to `plate_names`
#'   - Values are either the plate-specific parameter or `NULL`
#'   
#' @details
#' This function handles the following input patterns:
#' \itemize{
#'   \item{\code{param = NULL}: Returns `NULL`}
#'   \item{\code{param = c("value1", "value2")}: Returns a named list where all plates get the vector}
#'   \item{\code{param = list("Plate_01" = "value1", "Plate_02" = "value2")}: Returns a named list with plate-specific values}
#' }
#'
#' @keywords internal
process_named_param <- function(param, plate_names) {
  param_name <- deparse(substitute(param))
  
  if (is.null(param)) {
    return(NULL)
  }
  
  # Handle vectors - apply to all plates
  if (is.vector(param) && !is.list(param) && length(param) >= 1) {
    v_list <- rep(list(param), length(plate_names))
    names(v_list) <- plate_names
    return(v_list)
  } 
  
  # Handle named lists with strict validation
  if (is.list(param)) {
    if (is.null(names(param))) {
      stop(sprintf("'%s' must be a named list. Unnamed lists are not accepted.", param_name))
    }
    
    # Validate exact match
    if (!identical(sort(names(param)), sort(plate_names))) {
      invalid_names <- setdiff(names(param), plate_names)
      missing_names <- setdiff(plate_names, names(param))
      
      error_parts <- character()
      if (length(invalid_names) > 0) {
        error_parts <- c(error_parts, sprintf("invalid names: %s", paste(invalid_names, collapse = ", ")))
      }
      if (length(missing_names) > 0) {
        error_parts <- c(error_parts, sprintf("missing names: %s", paste(missing_names, collapse = ", ")))
      }
      
      stop(sprintf("'%s' names do not match plate_names. %s", 
                   param_name, paste(error_parts, collapse = "; ")))
    }
    
    if (length(param) != length(plate_names)) {
      stop(sprintf("'%s' length (%d) must match plate_names length (%d)", 
                   param_name, length(param), length(plate_names)))
    }
    
    # Reorder to match plate_names order
    result <- lapply(plate_names, function(plate) param[[plate]])
    names(result) <- plate_names
    return(result)
  }
  
  stop(sprintf("'%s' must be NULL, a vector, or a named list. Got: %s", 
               param_name, class(param)[1]))
}

#' Remove internal control targets from data matrices
#'
#' Helper function that removes specified internal control (IC) targets from 
#' wide-format data matrices by filtering out rows matching the IC target names.
#' 
#' @param mat A matrix or data frame with targets as rownames
#' @param ic_targets Character vector of internal control target names to remove
#'
#' @return The input matrix with IC target rows removed, or NULL if input is NULL.
#'   Maintains matrix structure even when only one row remains (drop = FALSE).
#'   
#' @details
#' This function is used internally to clean data matrices by removing internal
#' control targets that should not be included in downstream analysis.
#'
#' @keywords internal
remove_ic_from_matrix <- function(mat, ic_targets) {
  if (is.null(mat)) return(mat)
  targets_to_keep <- setdiff(rownames(mat), ic_targets)
  return(mat[targets_to_keep, , drop = FALSE])
}

#' Remove internal control targets from long-format data frames
#'
#' Helper function that removes specified internal control (IC) targets from 
#' long-format data frames by filtering out rows where the Target column 
#' matches the IC target names.
#' 
#' @param df A data frame in long format with a 'Target' column
#' @param ic_targets Character vector of internal control target names to remove
#'
#' @return The input data frame with IC target rows removed, or NULL if input is NULL.
#'   
#' @details
#' This function handles the Target column as either a factor or character vector,
#' converting factors to character for reliable comparison. Used internally to 
#' clean long-format data frames during data processing.
#'
#' @keywords internal
remove_ic_from_long <- function(df, ic_targets) {
  if (is.null(df)) return(df)
  # Handle factor columns by converting to character for comparison
  return(df[!as.character(df$Target) %in% ic_targets, ])
}
