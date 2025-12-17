#' Write NULISAseq data in long format Excel file
#'
#' @description

#' Reads in one or more NULISAseq XML files. Output is a long-format Excel file
#' where each row corresponds to a particular sample-target combination.
#' 
#' @details 
#' #' For RQ data, Excel file tabs include RQ data, target detectability 
#' and sample information.
#' For AQ data, tabs include RQ data, AQ data, target detectability, 
#' target quantifiability, and sample information.
#' 
#' The function takes as input an optional target info file containing target metadata, 
#' including Alamar target IDs and protein names, and an optional sample info 
#' file containing sample metadata that
#' is desired to appear in the output file.
#' 
#' Sample QC is currently based on the internal control within 40% of
#' median criterion. If a sample's IC count falls outside this threshold, 
#' it will be assigned a "WARN" status. Otherwise the sample will be assigned a 
#' "PASS" status.  
#' 
#' Target QC (if AQ data is available) is currently based on the target concentration accuracy
#' within 30% of median criterion. If a target's concentration accuracy falls outside this threshold, 
#' it will be assigned a "WARN" status. Otherwise the target will be assigned a 
#' "PASS" status.  
#' 
#' The function outputs long format data where each row corresponds to a particular
#' sample-target combination. Total number of rows is number of samples times
#' number of targets. By default the sample control (SC/AQSC) wells are included 
#' in the output, and IPC and NC wells are omitted.
#' 
#' 
#' 
#'
#' @param xml_files Vector of file names (character strings). Files should 
#' be in order of the desired plateID variable, unless that is otherwise
#' defined.
#' @param dataDir Data directory where xml_files, target_info_file, and sample_info_file 
#' reside (character string). 
#' @param target_info_file Optional. Path and filename for the target info CSV file 
#' (character string). Must include columns for TargetName (matches the targets in 
#' XML files), AlamarTargetID, UniProtID, and ProteinName. Only targets in the 
#' target_info_file will be output in the CSV file. If target_info_file is not provided, 
#' the function will use information from Barcode A in the xml files.
#' @param output_filename Filename for output xlsx file. 
#' @param Panel Name of multi-plex panel. Default is '200-plex Inflammation v1'
#' @param PanelLotNumber The panel lot number. Can be a single value or vector applied
#' to all plates, or a named list with plate-specific values where names correspond to
#' plate IDs.
#' @param plateIDs A vector of plate IDs. If `NULL`, default is to number plates from 01 to 
#' total number of plates based on the order of xml_files. Passed to `readNULISAseq()` function
#' and output as a column in the data file. 
#' @param sample_info_file Optional. Path and filename for the sample annotation CSV file. 
#' Must include columns "plateID" and "sampleName" which correspond to the 
#' plateID (matching the plateIDs input to this function) and sampleName in the `readNULISAseq()` 
#' `samples` data.frame.
#' @param sample_info_file_variables Subset of column names in 
#' sample_info_file that will be included in data file output. Other columns 
#' will be excluded. Otherwise, if `NULL` (default), all columns will be included.
#' @param sample_group_covar Optional column name in the Barcode B file 
#' and samples data matrix output by readNULISAseq that represents subgroups
#' for which detectability will be calculated separately. 
#' Default is 'SAMPLE_MATRIX', Function will check first to
#' be sure that the variable is present in the column names of the samples matrix.
#' Can be set to NULL to not use this feature.
#' @param ICs Vector of string(s) or a named list of character string vectors . 
#' Internal control names. Default is "mCherry". First IC in vector will be used in 
#' intra-plate IC-normalization (usually mCherry). 
#' ICs will be omitted from output file by default. Can be a single value or vector 
#' applied to all plates, or a named list with plate-specific values where names 
#' correspond to plate IDs.
#' @param IPC_string Optional vector of string(s) or a named list of character string vectors 
#' that identifies the sample names of IPC wells 
#'  (e.g. 'IPC'). This will override the default behavior which is 
#'  to use the sample 'type' column from Barcode B information. Can be a single value or vector 
#' applied to all plates, or a named list with plate-specific values where names 
#' correspond to plate IDs.
#' @param SC_string Optional vector of string(s) or a named list of character string vectors 
#' that identifies the sample names of SC wells 
#'  (e.g. 'SC'). This will override the default behavior which is 
#'  to use the sample 'type' column from Barcode B information. Can be a single value or vector 
#' applied to all plates, or a named list with plate-specific values where names 
#' correspond to plate IDs.
#' @param Bridge_string Optional vector of string(s) or a named list of character string vectors 
#' that identifies the sample names of Bridge wells 
#'  (e.g. 'Bridge'). This will override the default behavior which is 
#'  to use the sample 'type' column from Barcode B information. 
#'  (Bridge normalization not currently implemented.) Can be a single value or vector 
#' applied to all plates, or a named list with plate-specific values where names 
#' correspond to plate IDs.
#' @param Calibrator_string Optional vector of character string(s) or a named list 
#' of character string vectors that identifies the sample names of Calibrator wells 
#'  (e.g. 'Calibrator'). This will override the default behavior which is 
#'  to use the sample 'type' column from Barcode B information. 
#'  (Calibrator not currently implemented. AQ uses IPC as calibrator.)
#'  Can be a single value or vector applied to all plates, or a named 
#'  list with plate-specific values where names correspond to plate IDs.
#' @param NC_string Optional vector of character string(s) or a named list of character string vectors 
#' that identifies the sample names of NC wells 
#'  (e.g. 'NC'). This will override the default behavior which is 
#'  to use the sample 'type' column from Barcode B information. Can be a single value or vector 
#' applied to all plates, or a named list with plate-specific values where names 
#' correspond to plate IDs.
#' @param include_IPC Logical. Should IPC samples be included in output? Default is `FALSE`.
#' @param include_SC Logical. Should SC samples be included in output? Default is `TRUE`.
#' @param include_Bridge Logical. Should Bridge samples be included in output? Default is `TRUE`.
#' @param include_Calibrator Logical. Should Calibrator samples be included in output? Default is `FALSE`.
#' @param include_NC Logical. Should NC samples be included in output? Default is `FALSE`.
#' @param include_unnorm_counts Logical. Should unnormalized counts be included 
#' as an additional column in output? Default is `FALSE`.
#' @param include_IC_counts Logical. Should IC counts be included in the output? 
#' Default is `FALSE`. This is only useful when 
#' \code{include_unnorm_counts=TRUE}.
#' @param excludeSamples Optional vector of string(s) or a named list of character string vectors 
#' that give sample names to be excluded from the output file. 
#' Can be a single value or vector applied to all plates, or a named list with plate-specific 
#' values where names correspond to plate IDs. If no sample is to be excluded from a plate, 
#' use `NULL` in the list for that plate.
#' @param excludeTargets Optional vector of string(s) or a named list of character string vectors 
#' that give target names to be excluded from the output file. Can be a single value or vector 
#' applied to all plates, or a named list with plate-specific values where names correspond 
#' to plate IDs. If no target is to be excluded from a plate, use `NULL` in the list for that plate.
#' @param interPlateNorm_method Default is "IPC" for inter-plate control normalization. 
#' Use "IN" for IPC normalization followed by intensity normalization. 
#' @param IN_samples Optional argument passed to `interPlateNorm` function. 
#' A list of column names or indices specifying which subset of samples to use 
#' for intensity normalization step for each plate in data_list. 
#' By default, when this is set to `NULL`, all samples are used for IN. 
#' @param replaceNA Logical. Passed to `readNULISAseq()` function.
#' If `TRUE` (default), will replace any missing counts with 
#' zero for generating NPQ. (For AQ data, 
#' see \code{replace_zeros_with_NA}, below.)
#' @param TAP If `TRUE` (default), uses TAP detectability criteria in sample QC 
#' which includes more matrix types than non-TAP criteria. However, 
#' this function currently only flags samples based on IC median deviation.
#' @param output_TAP_AQ Logical. Default is `FALSE`. If xml files include AQ 
#' parameters and this is set to TRUE, the function will output the Excel 
#' file format used for TAP projects. File will have 5 tabs including RQ data, 
#' AQ data, detectability, quantifiability, and sample information.
#' @param replace_cal_blank_zeros Logical. Default is `FALSE`. If `TRUE`, AQ 
#' targets with a zero mean NC value will have the a_yint parameter replaced with 
#' the AQMC value, rather than be set to zero. 
#' @param replace_zeros_with_NA Logical. Default is `TRUE`. When `TRUE`, will replace
#' and zero AQ data values with NA.
#' @param metadata An optional named list of metadata, such as software package 
#' version numbers, that will be added as the last sheet in the Excel file. 
#' The list names will be the first column, and the values will form the second 
#' column. If `NULL` (default) then no metadata sheet is added.
#' @param verbose Logical. Should function output step completion info.
#' Default is `TRUE`.
#'
#' @return Outputs an Excel file.
#'
#' @export
#'
writeNULISAseq <- function (xml_files, dataDir, target_info_file = NULL, output_filename, 
                            Panel = "200-plex Inflammation v1", PanelLotNumber = "", 
                            plateIDs = NULL, sample_info_file = NULL, sample_info_file_variables = NULL, 
                            sample_group_covar = "SAMPLE_MATRIX", ICs = "mCherry", IPC_string = NULL, 
                            SC_string = NULL, Bridge_string = NULL, Calibrator_string = NULL, 
                            NC_string = NULL, include_IPC = FALSE, include_SC = TRUE, 
                            include_Bridge = TRUE, include_Calibrator = FALSE, include_NC = FALSE, 
                            include_unnorm_counts = FALSE, include_IC_counts = FALSE, 
                            excludeSamples = NULL, excludeTargets = NULL, interPlateNorm_method = "IPC", 
                            IN_samples = NULL, replaceNA = TRUE, TAP = TRUE, output_TAP_AQ = FALSE, 
                            replace_cal_blank_zeros = FALSE, replace_zeros_with_NA = TRUE, 
                            metadata = NULL, verbose = TRUE) {
  ############## Process xmls ##############
  all_data <- importNULISAseq(files = file.path(dataDir, xml_files), 
                              plateName = plateIDs,
                              sample_group_covar = sample_group_covar, 
                              IC = ICs,
                              IPC = IPC_string,
                              SC = SC_string,
                              NC = NC_string,
                              Bridge = Bridge_string,
                              Calibrator = Calibrator_string,
                              replaceNA = replaceNA, 
                              excludeSamples = excludeSamples, 
                              excludeTargets = excludeTargets,
                              TAP = TAP,
                              replace_cal_blank_zeros = replace_cal_blank_zeros,
                              replace_zeros_with_NA = replace_zeros_with_NA,  
                              security = FALSE)
  
  ############## Process Target ##############
  run_targets_cols <- Reduce(intersect, lapply(all_data$runs, function(x) colnames(x$targets)))
  run_targets_cols <- run_targets_cols[run_targets_cols %in% 
                                         c("targetName", "AlamarTargetID", "MW", "ProteinName")]
  run_targets <- do.call(rbind, lapply(all_data$runs, function(x) x$targets[, run_targets_cols]))
  if ("MW" %in% colnames(run_targets)) {
    run_targets$MW <- round(as.numeric(run_targets$MW))
  }
  run_targets <- unique(run_targets)
  if (nrow(run_targets) != length(unique(run_targets$targetName))) {
    stop("Inconsistent target metadata was found across runs. Please check targetName, AlamarTargetID, MW and ProteinName in the Barcode A information.")
  }
  if (!is.null(target_info_file)) {
    target_info <- read.csv(file.path(dataDir, target_info_file))
    run_target_cols <- unique(c("targetName", run_targets_cols[!(run_targets_cols %in% 
                                                                   colnames(target_info))]))
    run_targets <- run_targets[, run_target_cols]
    all_targets <- merge(target_info[target_info$TargetName %in% 
                                       run_targets$targetName, ], run_targets, by.x = "TargetName", 
                         by.y = "targetName", all.x = TRUE, all.y = FALSE, 
                         sort = FALSE)
  } else {
    all_targets <- run_targets
  }
  colnames(all_targets) <- toupper(colnames(all_targets))
  all_targets <- all_targets[order(all_targets$TARGETNAME), ]
  
  if (include_IC_counts == TRUE) {
    missing_ICs <- ICs[!ICs %in% all_targets$TARGETNAME]
    if (length(missing_ICs) > 0) {
      # Create new rows with proper NA types
      new_rows <- data.frame(TARGETNAME = missing_ICs)
      all_targets <- dplyr::bind_rows(all_targets, new_rows)
    }
    
    # Reorder
    all_targets <- dplyr::arrange(all_targets, !TARGETNAME %in% ICs, TARGETNAME)
  } else {
    all_targets <- all_targets[!all_targets$TARGETNAME %in% ICs, ]
  }
  
  all_targets <- all_targets %>%
    dplyr::rename(dplyr::any_of(c(
      "AlamarTargetID" = "ALAMARTARGETID",
      "UniProtID" = "UNIPROTID",
      "ProteinName" = "PROTEINNAME"
    ))) %>%
    dplyr::mutate(
      AlamarTargetID = if (!"AlamarTargetID" %in% colnames(.)) NA else AlamarTargetID,
      UniProtID = if (!"UniProtID" %in% colnames(.)) NA else UniProtID,
      ProteinName = if (!"ProteinName" %in% colnames(.)) NA else ProteinName
    )
  ############## Process Sample ##############
  sample_data <- all_data$merged$samples
  
  if (!is.null(sample_info_file)) {
    sample_info <- read.csv(file.path(dataDir, sample_info_file))
    if (!is.null(sample_info_file_variables)) {
      sample_info <- sample_info[, c("plateID", "sampleName", 
                                     sample_info_file_variables)]
    }
    else {
      sample_info_file_variables <- colnames(sample_info)
      sample_info_file_variables <- sample_info_file_variables[!(sample_info_file_variables %in% 
                                                                   c("plateID", "sampleName"))]
    }
    sample_data <- merge(sample_data, sample_info, by = c("plateID", 
                                                          "sampleName"), all.x = TRUE, all.y = FALSE, sort = FALSE)
  }
  
  sampleQC_IC_Median <- all_data$merged$qcSample %>% 
    dplyr::filter(flagName == "IC_Median") %>% 
    dplyr::mutate(SampleQC = ifelse(status == TRUE, "WARN", "PASS")) %>% 
    dplyr::select(plateID, sampleName, SampleQC)
  sample_data <- sample_data %>%
    dplyr::left_join(sampleQC_IC_Median, by = c("plateID", "sampleName"))
  
  if (verbose == TRUE) 
    cat("Sample QC completed.\n")
  
  # Define sample types to include based on parameters
  sample_types <- "Sample"
  
  if (include_IPC) {
    sample_types <- c(sample_types, "IPC")
  }
  if (include_SC) {
    sample_types <- c(sample_types, "SC")
  }
  if (include_NC) {
    sample_types <- c(sample_types, "NC")
  }
  if (include_Calibrator) {
    sample_types <- c(sample_types, "Calibrator")
  }
  
  SampleInfo <- all_data$merged$samples %>%
    dplyr::filter(sampleType %in% sample_types) %>% 
    dplyr::left_join(sampleQC_IC_Median, by = c("plateID", "sampleName")) %>% 
    dplyr::select(SampleName = sampleName, PlateID = plateID, SampleType = sampleType, sample_info_file_variables, SampleQC) %>% 
    dplyr::left_join(all_data$merged$Data_NPQ_long %>% 
                       dplyr::select(PlateID, SampleName, SampleType, Target, LOD, UnnormalizedCount, NPQ),
                     by = c("PlateID", "SampleName", "SampleType")) %>% 
    dplyr::arrange(PlateID, Target)
  
  SampleInfo_all <- all_data$merged$samples %>% 
    dplyr::select("plateID", "wellRow", "wellCol", "sampleName") %>% 
    dplyr::mutate(WellPosition = paste0(wellRow, "_", sprintf("%02s", wellCol))) %>%
    dplyr::select(-wellRow, -wellCol) %>% 
    dplyr::select(PlateID = plateID, WellPosition, SampleName = sampleName)
  
  ############## Process Data ##############
  if (interPlateNorm_method == "IN") {
    intensityNorm_data <- interPlateNorm(data_list = lapply(all_data$runs, function(x) x$NPQ), 
                                         IPC = FALSE, IN = TRUE, 
                                         IPC_wells = lapply(runs, function(x) x$IPC), NC_wells = lapply(runs, function(x) x$NC), 
                                         IN_samples = IN_samples, dataScale = "log", scaleFactor = 1, 
                                         transformReverse = NULL)
    Inten_Data <- lapply(names(intensityNorm_data$log2_interNormData), function(name) {
      x <- intensityNorm_data$log2_interNormData[[name]]
      long <- convert_to_long(x, "NPQ_intensity_norm") %>% 
        dplyr::mutate(PlateID = name) 
      return(long)
      
    })
    Inten_Data <- do.call(rbind, Inten_Data)
    
    SampleInfo <- SampleInfo %>%
      dplyr::left_join(Inten_Data, by = c("PlateID", "SampleName", "Target")) %>% 
      dplyr::mutate(LOD_intensity_norm = NPQ_intensity_norm - NPQ + LOD)
    
    if (verbose == TRUE) 
      cat("Intensity normalization completed.\n")
  }
  if (interPlateNorm_method != "IPC" & interPlateNorm_method != 
      "IN") {
    warning("No valid inter-plate normalization method was specified.")
  }
  
  targets <- all_targets$TARGETNAME
  
  base_cols <- c("Panel", "PanelLotNumber", "PlateID", "SampleName", "SampleType", 
                 sample_info_file_variables, "Target", "AlamarTargetID", "UniProtID", 
                 "ProteinName", "SampleQC")
  
  if (include_unnorm_counts) {
    target_cols <- c(base_cols, "LOD", "UnnormalizedCount", "NPQ")
  } else {
    target_cols <- c(base_cols, "LOD", "NPQ")
  }
  # Add PanelLotNumber by plate if named list or vector provided
  processed_panel_lot <- process_named_param(PanelLotNumber, as.character(unique(SampleInfo$PlateID)))
  
  allData <- SampleInfo %>% 
    dplyr::filter(Target %in% targets) %>% 
    dplyr::left_join(all_targets, by = c("Target" = "TARGETNAME")) %>% 
    dplyr::mutate(
      Panel = Panel,  
      PanelLotNumber = if (!is.null(processed_panel_lot)) {
        # Create a named vector for easy lookup
        lot_vector <- setNames(unlist(processed_panel_lot), names(processed_panel_lot))
        lot_vector[as.character(PlateID)]  # Lookup by PlateID
      } else {
        NA_character_
      },
      SampleType = factor(SampleType, 
                          levels = c("Sample", "IPC", "SC", "Bridge", "Calibrator", "NC"))) %>% 
    dplyr::arrange(PlateID, SampleType) %>%
    dplyr::select(dplyr::all_of(target_cols))
  
  ############## Detectability ##############
  detect_table <- all_data$merged$detectability
  
  ############## Write RQ Output ##############
  anyMissingAQ <- any(sapply(all_data$runs, function(x) is.null(x$AQ)))
  if (output_TAP_AQ == FALSE | (output_TAP_AQ == TRUE & anyMissingAQ)) {
    
    # Convert NaN to NA in all data frames before creating the output list
    allData <- allData %>%
      dplyr::mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .)))
    
    detect_table <- detect_table %>%
      dplyr::mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .)))
    
    # Create workbook
    wb <- openxlsx::createWorkbook()
    
    # Add sheets in order
    add_sheet_to_wb(wb, "NPQ Values", allData, na_as_string = FALSE)
    add_sheet_to_wb(wb, "Target Detectability", detect_table, na_as_string = TRUE)
    add_sheet_to_wb(wb, "Sample Information", SampleInfo_all, na_as_string = TRUE)
    
    if (!is.null(metadata)) {
      metadata_sheet <- data.frame(names(metadata), unlist(metadata))
      add_sheet_to_wb(wb, "Analysis Metadata", metadata_sheet, colNames = FALSE, na_as_string = TRUE)
    }
    
    openxlsx::saveWorkbook(wb, file = file.path(dataDir, output_filename), 
                           overwrite = TRUE)
    
    if (verbose == TRUE) 
      cat("writeNULISAseq completed. RQ data file was generated.")
  }
  
  ############## AQ ##############
  if (output_TAP_AQ == TRUE & !anyMissingAQ) {
    sample_cols <- c("SampleName", "PlateID", "SampleType", sample_info_file_variables, "SampleQC") 
    
    target_Conc_Accuracy <- all_data$merged$qcTarget %>% 
      dplyr::filter(flagName == "Target_Conc_Accuracy") %>% 
      dplyr::mutate(TargetQC = ifelse(status == TRUE, "WARN", "PASS")) %>% 
      dplyr::select(plateID, target, TargetQC)
    
    all_targets_AQ <- all_targets %>% 
      dplyr::filter(TARGETNAME %in% unique(all_data$merged$Data_AQ_long$Target)) %>% 
      dplyr::select(-dplyr::any_of("MW")) %>% 
      dplyr::left_join(target_Conc_Accuracy, by = c("TARGETNAME" = "target")) %>% 
      dplyr:: arrange(plateID)
    
    SampleInfo_AQ <- SampleInfo %>%
      dplyr::select(dplyr::any_of(sample_cols)) %>%
      dplyr::distinct() %>% 
      dplyr::left_join(all_data$merged$Data_AQ_long %>% 
                         dplyr::select(PlateID, SampleName, SampleType, Target, 
                                       dplyr::contains("LOD") | dplyr::contains("Conc") | dplyr::contains("LLOQ") | dplyr::contains("ULOQ")),
                       by = c("PlateID", "SampleName", "SampleType")) 
    
    sample_type_levels <- c('Sample', 'SC', 'AQSC', 'IPC', 'CAL', 'NC')
    
    allData_AQ <- SampleInfo_AQ %>% 
      dplyr::filter(Target %in% all_targets_AQ$TARGETNAME) %>% 
      dplyr::left_join(all_targets_AQ, by = c("Target" = "TARGETNAME", "PlateID" = "plateID")) %>% 
      dplyr::mutate(Panel = Panel,
                    PanelLotNumber = if (!is.null(processed_panel_lot)) {
                      # Create a named vector for easy lookup
                      lot_vector <- setNames(unlist(processed_panel_lot), names(processed_panel_lot))
                      lot_vector[as.character(PlateID)]  # Lookup by PlateID
                    } else {
                      NA_character_
                    },
                    SampleType = dplyr::case_when(
                      SampleType == "IPC" ~ "CAL", 
                      SampleType == "SC" ~ "AQSC",
                      TRUE ~ SampleType),
                    SampleType = factor(SampleType, levels = sample_type_levels),
                    is_IC = Target %in% ICs) %>% 
      # Sort properly: non-IC targets first (by Target, PlateID, SampleType, SampleName), then IC targets
      dplyr::arrange(is_IC, Target, PlateID, SampleType, SampleName) %>%
      dplyr::select(-is_IC) %>%  
      dplyr::select(dplyr::any_of(c("Panel", "PanelLotNumber", 
                                    "PlateID", "SampleName", "SampleType", sample_info_file_variables, 
                                    "Target", "AlamarTargetID", "UniProtID", "ProteinName", 
                                    "SampleQC", "TargetQC", "Conc_pgmL", "LOD_pgmL", 
                                    "LLOQ_pgmL", "ULOQ_pgmL", "Conc_aM", "LOD_aM", 
                                    "LLOQ_aM", "ULOQ_aM"))) %>%
      dplyr::rename(
        "Conc (pg/mL)" = "Conc_pgmL",
        "LOD (pg/mL)" = "LOD_pgmL",
        "LLOQ (pg/mL)" = "LLOQ_pgmL",
        "ULOQ (pg/mL)" = "ULOQ_pgmL",
        "Conc (aM)" = "Conc_aM",
        "LOD (aM)" = "LOD_aM",
        "LLOQ (aM)" = "LLOQ_aM",
        "ULOQ (aM)" = "ULOQ_aM"
      )
    
    allData <- allData %>% 
      dplyr::mutate(
        SampleType = dplyr::case_when(
          SampleType == "IPC" ~ "CAL", 
          SampleType == "SC" ~ "AQSC",
          TRUE ~ SampleType),
        SampleType = factor(SampleType, levels = sample_type_levels),
        is_IC = Target %in% ICs
      ) %>%
      dplyr::arrange(is_IC, Target, PlateID, SampleType, SampleName) %>%
      dplyr::select(-is_IC)
    
    ############## Quantifiability ##############
    quant_table <- all_data$merged$quantifiability
    
    ############## Write AQ Output ##############
    header1 <- "* All reported concentrations are functional concentrations (i.e., the concentration present in the sample)."
    header2 <- "* Please note that LOD or sample AQ values will be blank if they fall above or below the asymptotes of the AQ calibration curve."
    AQ_header <- rbind(header1, header2)
    
    # Convert NaN to NA in all data frames before creating the output list
    allData <- allData %>%
      dplyr::mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .)))
    
    allData_AQ <- allData_AQ %>%
      dplyr::mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .)))
    
    detect_table <- detect_table %>%
      dplyr::mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .)))
    
    quant_table <- quant_table %>%
      dplyr::mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .)))
    
    # Create workbook
    wb <- openxlsx::createWorkbook()
    
    # Add RQ NPQ Values sheet with empty NAs
    add_sheet_to_wb(wb, "RQ NPQ Values", allData, na_as_string = FALSE)
    
    # Add AQ sheet with header at row 1 and data at row 3 (both with empty NAs)
    openxlsx::addWorksheet(wb, "AQ (pg per mL and aM)", gridLines = TRUE)
    openxlsx::writeData(wb, sheet = "AQ (pg per mL and aM)", 
                        x = AQ_header, startCol = 1, startRow = 1, 
                        colNames = FALSE, keepNA = FALSE)
    openxlsx::writeData(wb, sheet = "AQ (pg per mL and aM)", 
                        x = allData_AQ, startCol = 1, startRow = 3, 
                        colNames = TRUE, keepNA = FALSE)
    
    # Add remaining sheets with NA as string
    add_sheet_to_wb(wb, "Target Detectability", detect_table, na_as_string = TRUE)
    add_sheet_to_wb(wb, "Target Quantifiability", quant_table, na_as_string = TRUE)
    add_sheet_to_wb(wb, "Sample Information", SampleInfo_all, na_as_string = TRUE)
    
    if (!is.null(metadata)) {
      metadata_sheet <- data.frame(names(metadata), unlist(metadata))
      add_sheet_to_wb(wb, "Analysis Metadata", metadata_sheet, colNames = FALSE, na_as_string = TRUE)
    }
    
    openxlsx::saveWorkbook(wb, file = file.path(dataDir, output_filename), 
                           overwrite = TRUE)
    
    if (verbose == TRUE) 
      cat("writeNULISAseq completed. AQ data file was generated.")
  }
}

#' Add a Sheet to an Excel Workbook with Configurable NA Handling
#'
#' This helper function adds a new worksheet to an openxlsx workbook object and 
#' writes data to it with customizable handling of NA values. It provides a 
#' convenient wrapper around openxlsx::addWorksheet() and openxlsx::writeData() 
#' with standardized formatting.
#'
#' @param wb An openxlsx workbook object created with openxlsx::createWorkbook()
#' @param sheet_name Character string specifying the name of the worksheet to create
#' @param data A data frame or matrix containing the data to write to the worksheet
#' @param colNames Logical indicating whether to write column names. Default is TRUE
#' @param na_as_string Logical indicating how to handle NA values:
#'   - If TRUE, NA values are written as the string "NA" in Excel cells
#'   - If FALSE, NA values are written as empty cells in Excel
#'   Default is FALSE (empty cells)
#' @param startRow Integer specifying the row number to start writing data. Default is 1
#' @param startCol Integer specifying the column number to start writing data. Default is 1
#'
#' @return None. The function modifies the workbook object in place by adding a 
#'   new worksheet with the specified data.
#'
#' @details
#' The function creates a worksheet with gridLines enabled by default. When 
#' na_as_string = FALSE, the function uses keepNA = FALSE to ensure NA values 
#' appear as empty cells rather than displaying any text. This is particularly 
#' useful for detectability and quantifiability table where empty cells are shown as "NA" text.
#'
#' @keywords internal
#'
#' @seealso \code{\link[openxlsx]{addWorksheet}}, \code{\link[openxlsx]{writeData}}
add_sheet_to_wb <- function(wb, sheet_name, data, colNames = TRUE, na_as_string = FALSE, 
                            startRow = 1, startCol = 1) {
  openxlsx::addWorksheet(wb, sheet_name, gridLines = TRUE)
  if (na_as_string) {
    openxlsx::writeData(wb, sheet = sheet_name, x = data, 
                        startCol = startCol, startRow = startRow, 
                        colNames = colNames, na.string = "NA",
                        keepNA = TRUE)
  } else {
    openxlsx::writeData(wb, sheet = sheet_name, x = data, 
                        startCol = startCol, startRow = startRow, 
                        colNames = colNames, keepNA = FALSE)
  }
}

