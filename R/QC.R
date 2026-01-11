# Sample QC criteria
MIN_FRAC_DETECTABILITY_TAP <- c(plasma=0.9, serum=0.9, csf=0.7, urine=0.65, cell_culture=0.3, nhp_plasma=0.55, nhp_serum=0.55, nhp_csf=0.35, dried_blood_spot=0.75, control=0.9, other=0.0 )  # Minimim fraction (Target_Detectability): # Targets with reads above LOD
MIN_FRAC_DETECTABILITY <- c(plasma=0.9, serum=0.9, csf=0.7, other=0.0 )  # Minimim fraction (Target_Detectability): # Targets with reads above LOD
MIN_IC_READS_PER_SAMPLE <- 1000    # Minimum number (ICReads) of IC reads within a sample
MIN_NUM_READS_PER_SAMPLE <- 500000 # Minimum number (NumReads) of reads within a sample
#MAX_QCS <- 3
#MIN_SN <- 4
MIN_IC_MEDIAN <- "-0.4,0.4"# +/- of sample IC read count about the median
#' QCSampleCriteria
#'
#' Sample QC Criteria
#' @param TAP Whether the run is a TAP experiment
#' @return criteria
#' @examples
#' # QCSampleCriteria()
#'
#' @export
QCSampleCriteria <- function(TAP=TRUE){
  retVal <- NULL
  detect_criteria <- if(TAP) MIN_FRAC_DETECTABILITY_TAP else MIN_FRAC_DETECTABILITY
  detect_criteria_name <- if(TAP) "MIN_FRAC_DETECTABILITY_TAP" else "MIN_FRAC_DETECTABILITY"
  retVal$thresholds <-c(Detectability=detect_criteria,
                        ICReads=MIN_IC_READS_PER_SAMPLE, 
                        NumReads=MIN_NUM_READS_PER_SAMPLE, 
                        IC_Median=MIN_IC_MEDIAN)#,
#                        QCS=MAX_QCS,
#                        SN=MIN_SN)
  retVal$operators <-c(Detectability="<",
                       ICReads="<", 
                       NumReads="<", 
                       IC_Median="<,>")#,
#                       QCS=">",
#                       SN="<")
  retVal$format <- c(Detectability="percentage",
                       ICReads="integer", 
                       NumReads="integer", 
                       IC_Median="percentage")#,
#                       QCS="float",
#                       SN="float")
  retVal$thresholdNames<-c(Detectability=detect_criteria_name, 
                        ICReads="MIN_IC_READS_PER_SAMPLE", 
                        NumReads="MIN_NUM_READS_PER_SAMPLE", 
                        IC_Median="MIN_IC_MEDIAN")#,
#                        QCS="MAX_QCS",
#                        SN="MIN_SN")
  retVal$properNames <- c(Detectability="Detectability",
                          ICReads="IC Reads",
                          NumReads="Reads",
                          IC_Median="IC Median"
                         )
  retVal$explanations <- c(Detectability=paste0("Percentage of targets with reads above the Limit of Detection (LOD) (Minimum threshold = plasma-", 100*MIN_FRAC_DETECTABILITY[["plasma"]], "%, serum-", 100*MIN_FRAC_DETECTABILITY[["serum"]], "%, csf-", 100*MIN_FRAC_DETECTABILITY[["csf"]], "%, other-", 100*MIN_FRAC_DETECTABILITY[["other"]], "%)"),
                          ICReads=paste0("Number of Internal Control (IC) reads within a sample (Minimum threshold = ", format(MIN_IC_READS_PER_SAMPLE, big.mark = ",", scientific = FALSE), ")"),
                          NumReads=paste0("Number of reads within a sample (Minimum threshold = ", format(MIN_NUM_READS_PER_SAMPLE, big.mark = ",", scientific = FALSE), ")"),
                          IC_Median=paste0("Sample IC reads relative to the median (Within +/-", 100*as.numeric(strsplit(MIN_IC_MEDIAN, ",")[[1]][2]), "% of the plate median)")
                         )
  if(TAP){
  retVal$explanations <- c(retVal$explanations, Detectability=paste0("Percentage of targets with reads above the Limit of Detection (LOD) (Minimum threshold = plasma-", 100*MIN_FRAC_DETECTABILITY_TAP[["plasma"]], "%, serum-", 100*MIN_FRAC_DETECTABILITY_TAP[["serum"]], "%, csf-", 100*MIN_FRAC_DETECTABILITY_TAP[["csf"]], "%, urine-", 100*MIN_FRAC_DETECTABILITY_TAP[["urine"]], "%, cell_culture-", 100*MIN_FRAC_DETECTABILITY_TAP[["cell_culture"]], "%, nhp_plasma-", 100*MIN_FRAC_DETECTABILITY_TAP[["nhp_plasma"]], "%, nhp_csf-", 100*MIN_FRAC_DETECTABILITY_TAP[["nhp_csf"]], "%, dried_blood_spot-", 100*MIN_FRAC_DETECTABILITY_TAP[["dried_blood_spot"]], "%, control-", 100*MIN_FRAC_DETECTABILITY_TAP[["control"]], "%, other-", 100*MIN_FRAC_DETECTABILITY[["other"]],"%)"))
  }
  return(retVal)
}

# Plate QC criteria
MAX_SC_CV <- 0.25                # (SCRead_CV) CV of IC reads across all samples
MAX_MEDIAN_SC_TARGET_CV <- 0.10  # (SCTarget_CV) CV of total read count for each IPC sample
MAX_IC_CV <- 0.25                # (ICRead_CV) CV of IC reads across all samples
MAX_IPC_CV <- 0.25               # (IPCRead_CV) CV of total read count for each IPC sample
MAX_MEDIAN_IPC_TARGET_CV <- 0.1  # (IPCTarget_CV) Median of CVs of all IPC targets (Performed on normalized data) 
DETECTABILITY_FRAC <- 0.90       # (Detectability) Target-wise Detectability fraction (target is detected if >50% of samples > LOD)
MIN_READS <- 1e8                 # (MinReads) Minimum number of reads
MAX_FAILED_TARGET_PERC <- 0.1    # (Failed_Targets)
MAX_FAILED_SC <- 0               # (Failed_SC)
MAX_FAILED_IPC <- 0              # (Failed_IPC)

#' QCPlateCriteria
#'
#' Plate QC Criteria
#' @param AQ whether the experiment is an AQ experiment 
#' @return criteria
#' @examples
#' # QCPlateCriteria()
#'
#' @export
QCPlateCriteria <- function(AQ=FALSE){
  retVal <- NULL
 
  retVal$thresholds <- c(ICRead_CV=as.numeric(MAX_IC_CV),
                         IPCRead_CV=as.numeric(MAX_IPC_CV), 
                         IPCTarget_CV=as.numeric(MAX_MEDIAN_IPC_TARGET_CV), 
                         Detectability=as.numeric(DETECTABILITY_FRAC), 
                         MinReads=as.numeric(MIN_READS))
  retVal$operators <- c(ICRead_CV=">", IPCRead_CV=">", IPCTarget_CV=">", Detectability="<", MinReads="<")
  retVal$format <- c(ICRead_CV="percentage", IPCRead_CV="percentage", IPCTarget_CV="percentage", Detectability="percentage", MinReads="integer")
  retVal$thresholdNames<-c(ICRead_CV="MAX_IC_CV",
                         IPCRead_CV="MAX_IPC_CV", 
                         IPCTarget_CV="MAX_MEDIAN_IPC_TARGET_CV", 
                         Detectability="DETECTABILITY_FRAC", 
                         MinReads="MIN_READS")

  retVal$properNames <- c(ICRead_CV="IC CV",
                          IPCRead_CV="IPC CV",
                          IPCTarget_CV="IPC Target CV",
                          Detectability="Run Detectability",
                          MinReads="Reads")

  retVal$explanations <- c(ICRead_CV=paste0("Coefficient of variation of internal control Parseable Matching reads across all wells (Maximum threshold = ", 100*as.numeric(MAX_IC_CV), "%)"),
                          IPCRead_CV=paste0("Coefficient of variation of Parseable Matching reads across all IPCs (Maximum threshold = ", 100*as.numeric(MAX_IPC_CV), "%)"),
                          IPCTarget_CV=paste0("Median coefficient of variation of all IPC targets (Maximum threshold = ", 100*as.numeric(MAX_MEDIAN_IPC_TARGET_CV), "%)"),
                          Detectability=paste0("Percentage of targets that are detectable (target is considered detectable if > 50% of samples are above LOD) (Minimum threshold = ",100*as.numeric(DETECTABILITY_FRAC), "%)"),
                          MinReads=paste0("Minimum number of Parseable Matching reads for the run (Minimum threshold = ", format(MIN_READS, big.mark = ",", scientific = FALSE), ")"))
  if(AQ){
    retVal$thresholds <- c(retVal$thresholds, 
                           SCRead_CV=as.numeric(MAX_SC_CV),
                           SCTarget_CV=as.numeric(MAX_MEDIAN_SC_TARGET_CV),
                           Failed_Targets=as.numeric(MAX_FAILED_TARGET_PERC),
                           Failed_SC = as.numeric(MAX_FAILED_SC),
                           Failed_IPC = as.numeric(MAX_FAILED_IPC)
                           )
    retVal$operators <- c(retVal$operators, SCRead_CV=">", SCTarget_CV=">", Failed_Targets=">", Failed_SC=">", Failed_IPC=">")
    retVal$format <- c(retVal$format, SCRead_CV="percentage", SCTarget_CV="percentage", Failed_Targets="percentage", Failed_SC="integer", Failed_IPC="integer")
    retVal$thresholdNames <- c(retVal$thresholdNames, 
                               SCRead_CV="MAX_SC_CV",
                               SCTarget_CV="MAX_MEDIAN_SC_TARGET_CV",
                               Failed_Targets="MAX_FAILED_TARGET_PERC",
                               Failed_SC="MAX_FAILED_SC",
                               Failed_IPC="MAX_FAILED_IPC")
    retVal$properNames <- c(ICRead_CV="IC CV",
                            IPCRead_CV="CAL CV",
                            IPCTarget_CV="CAL Target CV",
                            Detectability="Run Detectability",
                            MinReads="Reads",
                            SCRead_CV="AQSC CV",
                            SCTarget_CV="AQSC Target CV",
                            Failed_Targets="Target Warning",
                            Failed_SC="AQSC Sample Warning",
                            Failed_IPC="CAL Sample Warning"
                            )
    retVal$explanations <- c(ICRead_CV=paste0("Coefficient of variation of internal control Parseable Matching reads across all wells (Maximum threshold = ", 100*as.numeric(MAX_IC_CV), "%)"),
                            IPCRead_CV=paste0("Coefficient of variation of Parseable Matching reads across all IPCs (Maximum threshold = ", 100*as.numeric(MAX_IPC_CV), "%)"),
                            IPCTarget_CV=paste0("Median coefficient of variation of all IPC targets (Maximum threshold = ", 100*as.numeric(MAX_MEDIAN_IPC_TARGET_CV), "%)"),
                            Detectability=paste0("Percentage of targets that are detectable (target is considered detectable if > 50% of samples are above LOD) (Minimum threshold = ",100*as.numeric(DETECTABILITY_FRAC), "%)"),
                            MinReads=paste0("Minimum number of Parseable Matching reads for the run (Minimum threshold = ", format(MIN_READS, big.mark = ",", scientific = FALSE), ")"),
                            SCRead_CV=paste0("Coefficient of variation of Parseable Matching reads across all AQSCs (Maximum threshold = ", 100*MAX_SC_CV, "%)"),
                            SCTarget_CV=paste0("Median of coefficient of variation of normalized reads of all AQSC targets (Maximum threshold = ", 100*MAX_MEDIAN_SC_TARGET_CV, "%)"),
                            Failed_Targets=paste0("Percentage of Targets with Target QC warnings (Maximum threshold = ", 100*MAX_FAILED_TARGET_PERC, "%)"),
                            Failed_SC=paste0("Number of AQSC samples with Sample QC warnings (Maximum threshold = ", MAX_FAILED_SC,")"),
                            Failed_IPC=paste0("Number of CAL samples with Sample QC warnings (Maximum threshold = ", MAX_FAILED_IPC, ")")
                           )
  }
  return(retVal)
}

# Target QC criteria
TARGET_CONC_ACCURACY <- "-0.3,0.3"
TARGET_CONC_CV <- 0.3
MIN_TARGET_READS <- 200 # Minimum number of reads for a Target
PERC_MIN_TARGET_READS <- 0.5
TARGET_DETECTABILITY <- 0.5
TARGET_CONC_CV_RQ <- 0.3
TARGET_IPC_MIN_READS <- 200
#' QCTargetCriteria
#'
#' Target QC Criteria
#' @return criteria
#' @examples
#' # QCTargetCriteria()
#'
#' @export
QCTargetCriteria <- function(AQ=FALSE, advancedQC=FALSE){
  retVal <- NULL
  if(advancedQC){
    AQSCname = if(AQ) "AQSC" else "SC"
    retVal$thresholds <- c(Target_Min_Reads=as.numeric(PERC_MIN_TARGET_READS), 
                           Target_Conc_CV_RQ=as.numeric(TARGET_CONC_CV_RQ),
                           Target_Detectability=as.numeric(TARGET_DETECTABILITY)#,
  #                         Target_IPC_Min_Reads=as.numeric(TARGET_IPC_MIN_READS)
                           )
    retVal$operators <- c(Target_Min_Reads="<",
                          Target_Conc_CV_RQ=">",
                          Target_Detectability="<"#,
  #                        Target_IPC_Min_Reads="<"
                          )
    retVal$format <- c(Target_Min_Reads="percentage",
                       Target_Conc_CV_RQ="percentage",
                       Target_Detectability="percentage"#,
  #                     Target_IPC_Min_Reads="integer"
                       )
    retVal$thresholdNames <- c(Target_Min_Reads="TARGET_MIN_READS",
                               Target_Conc_CV_RQ="TARGET_CONC_CV_RQ",
                               Target_Detectability="TARGET_DETECTABILITY"#,
  #                             Target_IPC_Min_Reads="TARGET_IPC_MIN_READS"
                               )
    retVal$properNames <- c(Target_Min_Reads="Target Signal",
                            Target_Conc_CV_RQ="CV Normalized Reads",
                            Target_Detectability="Target Detectability"#,
  #                          Target_IPC_Min_Reads="IPC Reads"
                            )
    retVal$explanations <- c(Target_Min_Reads=paste0("Percentage of samples with raw reads of at least ", MIN_TARGET_READS, " (Minimum threshold = ", PERC_MIN_TARGET_READS*100, "%)"),
                             Target_Conc_CV_RQ=paste0("Target-specific ", AQSCname, " intra-plate coefficient of variation of normalized reads (Maximum threshold = ", TARGET_CONC_CV_RQ*100, "%)"),
                             Target_Detectability=paste0("Percentage of samples with signal above the limit of detection (Minimum threshold = ", TARGET_DETECTABILITY*100, "%)")#,
  #                           Target_IPC_Min_Reads=paste0("Median of IPC raw reads (Minimum threshold = ", TARGET_IPC_MIN_READS, ")")
                             )
  }
  if(AQ){
    retVal$thresholds <- c(retVal$thresholds,
                           Target_Conc_Accuracy=TARGET_CONC_ACCURACY,
                           Target_Conc_CV=as.numeric(TARGET_CONC_CV))
    retVal$operators <- c(retVal$operators,
                          Target_Conc_Accuracy="<,>",
                          Target_Conc_CV=">")
    retVal$format <- c(retVal$format,
                       Target_Conc_Accuracy="percentage",
                       Target_Conc_CV="percentage")
    retVal$thresholdNames <- c(retVal$thresholdNames, Target_Conc_Accuracy="TARGET_CONC_ACCURACY",
                             Target_Conc_CV="TARGET_CONC_CV")
    retVal$properNames <- c(retVal$properNames, 
                            Target_Conc_Accuracy="Accuracy",
                            Target_Conc_CV="CV")
    retVal$explanations <- c(retVal$explanations, 
                             Target_Conc_Accuracy=paste0("Experimental vs Theoretical average AQSC concentration accuracy (Within +/-", 100*as.numeric(strsplit(TARGET_CONC_ACCURACY, ",")[[1]][2]), "% of theoretical, calculated on AQ targets)"),
                             Target_Conc_CV=paste0("Target-specific AQSC coefficient of variation (Maximum threshold = ", TARGET_CONC_CV*100, "%, calculated on AQ targets)"))
  }
  return(retVal)
}

#' Write Processed XML from QC table
#'
#' Writes NULISAseq QC flags to XML
#' @param input QC table 
#' @param QC QC XML node to add to 
#' @param type whether the table is sample, target or plate
#' @param combined 
#' @return QC XML
#' @examples
#' # QC2XML(inputtable)
#'
#' @export
QC2XML <- function(input, QCNode, type="plate", combined=F){
  for (i in 1:nrow(input)){
    if(type=="sample"){
      if(combined){
        addChildren(QCNode, newXMLNode("QCFlag", input$val[i],
                                           attrs=c(
                                                  name=input$flagName[i],
                                                  set=input$status[i],
                                                  method=input$normMethod[i],
                                                  rep=input$sampleBarcode[i],
                                                  format=input$QCformat[i]
                                                  )
                )
        )
      }else{
        addChildren(QCNode, newXMLNode("QCFlag", input$val[i],
                                           attrs=c(
                                                  name=input$flagName[i],
                                                  set=input$status[i],
                                                  method=input$normMethod[i],
                                                  format=input$QCformat[i]
                                                  )
                )
        )
      }
    }else if(type=="plate" || type=="target") {
      addChildren(QCNode, newXMLNode("QCFlag", input$val[i],
                                           attrs=c(
                                                  name=input$flagName[i],
                                                  set=input$status[i],
                                                  method=input$normMethod[i],
                                                  format=input$QCformat[i]
                                                  )
                )
      )
    }
  }
  return(QCNode)
}

#' Write Sample QC table
#'
#' Writes NULISAseq QC sample flags 
#' @param raw unnormalized data
#' @param aboveLOD aboveLOD from the lod function
#' @param samples samples data frame
#' @param targets targets data frame
#' @param QCS QCS values
#' @param SN SN values
#' @param ICs name of IC samples
#' @param NCs names of NC samples
#' @param IPCs names of IPC samples
#' @param SCs names of SC samples
#' @param TAP whether this is a TAP run
#' @return QC table
#' @examples
#' # QCFlagSample(inputtable)
#'
#' @export

QCFlagSample <- function(raw, aboveLOD, samples, targets, QCS=NULL, SN=NULL, 
                         well_order=NULL, ICs=NULL, IPCs=NULL, NCs=NULL, SCs=NULL, TAP=TRUE){
  columns <- c("sampleName", "flagName", "normMethod", "status", "val", "text", 
               "sampleBarcode", "sampleType", "QCthreshold", "QCoperator", "QCformat")

  criteria <- QCSampleCriteria(TAP)

  QCFlagList <- vector("list", length(columns))
  if(is.null(well_order)){
    #well_order <- 1:length(raw[1,])
    well_order <-colnames(raw)
  }else{
    well_order <- rev(well_order)
  }
  ICs  <- if(!is.null(ICs))   ICs else targets$targetName[which(tolower(targets$targetType) == "control")]
  IPCs <- if(!is.null(IPCs)) IPCs else samples$sampleName[which(samples$sampleType == "IPC")]
  NCs  <- if(!is.null(NCs))   NCs else samples$sampleName[which(samples$sampleType == "NC")]
  SCs  <- if(!is.null(SCs))   SCs else samples$sampleName[which(samples$sampleType == "SC")]

  # QCS
#  qcTargets <- c("pTau-217", "CRH", "pTau-231", "IL4", "pTDP43-409", "IL33")
#  qcTargetBarcodes <- targets$targetBarcode[which(targets$targetName %in% qcTargets)]


#  op <- criteria$operators[which(criteria$thresholdNames=="MAX_QCS")]
#  format <- criteria$format[which(criteria$thresholdNames=="MAX_QCS")]

#  for(j in 1:ncol(raw)){
#    i <- well_order[j]
#    type <- "Sample"
#    if(i %in% NCs){
#      type <- "NC"
#    } else if (i %in% IPCs){
#      type <- "IPC"
#    } else if (i %in% SCs){
#      type <- "SC"
#    }

#    set = FALSE
#    value = NA
#    if(length(qcTargetBarcodes) > 0 & !is.null(QCS) & !is.null(SN)){
#      value = mean( log2(QCS[qcTargetBarcodes, samples$sampleBarcode[i]]), na.rm=T)
#      if(length(value) > 0){
#        if (value > MAX_QCS){
#          set = TRUE
#        }
#      }
#    }
#    QCFlagReturn <- rbind(QCFlagReturn, c(colnames(raw)[i], "QCS", "IC", set, value, "", samples$sampleBarcode[i], type,
#                                        as.character(paste(MAX_QCS, collapse=',')), as.character(op), format))
#  }

  # SN
#  qcTargets <- c("pTau-217")
#  qcTargetBarcodes <- targets$targetBarcode[which(targets$targetName %in% qcTargets)]

#  op <- criteria$operators[which(criteria$thresholdNames=="MIN_SN")]
#  format <- criteria$format[which(criteria$thresholdNames=="MIN_SN")]

#  for(j in 1:ncol(raw)){
#    i <- well_order[j]
#    type <- "Sample"
#    if(i %in% NCs){
#      type <- "NC"
#    } else if (i %in% IPCs){
#      type <- "IPC"
#    } else if (i %in% SCs){
#      type <- "SC"
#    }

#    set = FALSE
#    value = NA
#    if(!is.null(QCS) & !is.null(SN)){
#      value1 = QCS[qcTargetBarcodes, samples$sampleBarcode[i]]
#      value2 = SN[qcTargetBarcodes, samples$sampleBarcode[i]] / QCS[qcTargetBarcodes, samples$sampleBarcode[i]]
#      if (length(qcTargetBarcodes) > 0 & length(value1) > 0 & length(value2) > 0){
#        value = value2
#        if( value1 > MAX_QCS & value2 < MIN_SN){ 
#          set = TRUE
#        }
#      }
#    }
 #   QCFlagReturn <- rbind(QCFlagReturn, c(colnames(raw)[i], "SN", "IC", set, value, "", samples$sampleBarcode[i], type,
 #                                       as.character(paste(MIN_SN, collapse=',')), as.character(op), format))
#  }
#  SN[qcTargetBarcodes, samples$sampleBarcode[i]]

  # Median IC between -40% and 40% of median  
  mCherry_median <- median(raw[ICs[1],], na.rm=T)
  min_ic_median <- unname(unlist(strsplit(MIN_IC_MEDIAN, ",")))
#  medianMin30 <- mCherry_median - mCherry_median * abs(as.numeric(min_ic_median[0]))
#  medianMax30 <- mCherry_median + mCherry_median * abs(as.numeric(min_ic_median[1]))
  medVals <- (raw[ICs[1], ] - mCherry_median ) / mCherry_median 
  op <- criteria$operators[which(criteria$thresholdNames=="MIN_IC_MEDIAN")]
  format <- criteria$format[which(criteria$thresholdNames=="MIN_IC_MEDIAN")]
  for(j in 1:length(medVals)){
    i <- well_order[j]
    set <- if(is.na(medVals[i])) "T" else evalCriterion("IC_Median", medVals[i], op, MIN_IC_MEDIAN)
    type <- "Sample"
    type <- if(i %in% NCs) "NC" 
           else if(i %in% IPCs) "IPC"
           else if(i %in% SCs) "SC"
           else type
    QCFlagList[[j]] <- c(i, "IC_Median", "raw", set, medVals[i], "", samples$sampleBarcode[which(samples$sampleName==i)], type, 
                                          as.character(paste(MIN_IC_MEDIAN, collapse=',')), as.character(op), format)
  }
  curIndex <- length(medVals)
  # Minimim fraction (Target_Detectability): # Targets with reads above LOD
  # Exclude ICs and noDetectability targets from detectability calculation
  noDetectTargets <- if("noDetectability" %in% colnames(targets)) targets$targetName[which(targets$noDetectability == TRUE)] else character(0)
  aboveLOD <- aboveLOD[!rownames(aboveLOD) %in% c(ICs, noDetectTargets), ] # remove controls and noDetectability targets
  perc_tar <- colSums(aboveLOD == TRUE, na.rm=TRUE)/ colSums(!is.na(aboveLOD))
  perc_tar[NCs] <- NA
  min_frac_detect <- if(TAP) "MIN_FRAC_DETECTABILITY_TAP" else "MIN_FRAC_DETECTABILITY"
  op <- criteria$operators[which(criteria$thresholdNames==min_frac_detect)]
  format <- criteria$format[which(criteria$thresholdNames==min_frac_detect)]

  for (j in 1:length(perc_tar)){
    i <- well_order[j]
    #samplematrix <- if(tolower(samples$SAMPLE_MATRIX[i]) == "ipc") "plasma" else tolower(samples$SAMPLE_MATRIX[i])
    samplematrix <- tolower(samples$SAMPLE_MATRIX[which(samples$sampleName == i)])
    sampleType <- samples$sampleType[which(samples$sampleName ==i)]
    samplematrix <- ifelse(sampleType %in% c('SC', 'IPC'), 'plasma', samplematrix)
    set <- if(is.na(perc_tar[i])) NA else evalCriterion("Detectability", perc_tar[i], op, criteria$thresholds[paste0('Detectability.', samplematrix)])
    type <- "Sample"
    if(i %in% NCs){
      type <- "NC"
    } else if (i %in% IPCs){
      type <- "IPC"
    } else if (i %in% SCs){
      type <- "SC"
    }
    val <- samples$sampleBarcode[which(samples$sampleName==i)]
    if(length(val) == 0 ) val <- "NA"
    QCFlagList[[curIndex + j]] <- c(i, "Detectability", "IPC", set, perc_tar[i], "", 
                                          val, 
                                          type, 
                                          as.character(criteria$thresholds[paste0('Detectability.', samplematrix)]), 
                                          as.character(op), format)
  }

  curIndex <- curIndex + length(perc_tar)
  # Minimum number (ICReads) of IC reads within a sample
  ICvals <- raw[ICs[1], ]
  ICvals[is.na(ICvals)] <- 0
  op <- criteria$operators[which(criteria$thresholdNames=="MIN_IC_READS_PER_SAMPLE")]
  format <- criteria$format[which(criteria$thresholdNames=="MIN_IC_READS_PER_SAMPLE")]
  for (j in 1:length(ICvals)){
    i <- well_order[j]
    set <- evalCriterion("ICReads", ICvals[i], op, MIN_IC_READS_PER_SAMPLE)
    type <-"Sample"
    if(i %in% NCs){
      type <- "NC"
    } else if (i %in% IPCs){
      type <- "IPC"
    } else if (i %in% SCs){
      type <- "SC"
    }
    QCFlagList[[curIndex + j]] <- c(i, "ICReads", "raw", set, ICvals[i], "", samples$sampleBarcode[which(samples$sampleName==i)], type, 
                                          as.character(MIN_IC_READS_PER_SAMPLE), as.character(op), format)
  }

  # Minimum number (NumReads) of reads within a sample
  raw2 <- raw
  #barNames <- samples$sampleBarcode
  val <- colSums(raw2, na.rm=T)
  op <- criteria$operators[which(criteria$thresholdNames=="MIN_NUM_READS_PER_SAMPLE")]
  format <- criteria$format[which(criteria$thresholdNames=="MIN_NUM_READS_PER_SAMPLE")]
  curIndex <- curIndex + length(ICvals)
  for(j in 1:length(val)){
    i <- well_order[j]
    set <- NULL
    if(i %in% NCs){
      set <- NA  
      val[i] <- NA
    }else{
      set <- evalCriterion("NumReads", val[i], op, MIN_NUM_READS_PER_SAMPLE)
    }
    type <- "Sample"
    if(i %in% NCs){
      type <- "NC"
    } else if (i %in% IPCs){
      type <- "IPC"
    } else if (i %in% SCs){
      type <- "SC"
    }
    QCFlagList[[curIndex + j]] <- c(i, "NumReads", "raw", set, val[i], "", samples$sampleBarcode[which(samples$sampleName ==i)], type, 
                                          as.character(MIN_NUM_READS_PER_SAMPLE), as.character(op), format)
  }
  # Convert to matrix/data.frame at the end
  QCFlagReturn <- data.frame(do.call(rbind, QCFlagList), stringsAsFactors=FALSE)
  colnames(QCFlagReturn) <- columns
  QCFlagReturn$val <- as.numeric(QCFlagReturn$val)
  QCFlagReturn$QCthreshold <- as.character(QCFlagReturn$QCthreshold)
  return(QCFlagReturn)
}

# Return TRUE if criteria is violated (fail)
evalCriterion <- function(name, value, operator, threshold){
  # Convert to character only once
  op_char <- as.character(operator)
  th_char <- as.character(threshold)

  # Split strings
  operators <- strsplit(op_char, ",", fixed = TRUE)[[1]]
  thresholds <- as.numeric(strsplit(th_char, ",", fixed = TRUE)[[1]])

  # Length check
  if (length(operators) != length(thresholds)) {
    stop("Error: Specification of QC criteria is incorrect, unequal number of operators and thresholds for ", name)
  }

  # Direct comparison without eval/parse
  for (i in seq_along(operators)) {
    op <- operators[i]
    th <- thresholds[i]

    result <- switch(op,
      ">" = value > th,
      "<" = value < th,
      ">=" = value >= th,
      "<=" = value <= th,
      "==" = value == th,
      "!=" = value != th,
      stop("Unsupported operator: ", op)
    )

    if (is.na(value) || (!is.na(result) && result)) return(TRUE)
  }

  return(FALSE)
}

#' Write Target QC table
#'
#' Writes NULISAseq QC target flags 
#' @param AQdata Absolute quantification data
#' @param raw Raw data
#' @param IPCnormed IPC-normalized data
#' @param detectability detectability
#' @param aboveLOD aboveLOD
#' @param withinDR withinDR
#' @param absRun whether the assay is absolute quantification
#' @param targets targets data
#' @param samples samples data
#' @param SCparams SC concentrations 
#' @param ICs names of IC samples
#' @param NCs names of NC samples
#' @param IPCs names of IPC samples
#' @param SCs names of SC samples
#' @return QC table
#' @examples
#' # QCFlagTarget(inputtable)
#'
#' @export
QCFlagTarget <- function(AQdata, raw, IPCnormed, detectability, aboveLOD, withinDR, absRun, targets, samples, SCparams,
                        ICs = NULL, IPCs = NULL, NCs = NULL, SCs = NULL, advancedQC=FALSE) {
  # Load QC criteria
  criteria <- QCTargetCriteria(absRun, advancedQC)
  if(length(criteria) == 0){
    return(NULL)
  }
  # Determine SCs based on sample types if not provided
  SCs <- if (!is.null(SCs)) SCs else samples$sampleName[which(samples$sampleType == "SC")]
  IPCs <- if (!is.null(IPCs)) IPCs else samples$sampleName[which(samples$sampleType == "IPC")]
  # Determine Samples based on sample types if not provided
  Samples <- samples$sampleName[which(tolower(samples$sampleType) == "sample")]
  Targets <- targets$targetName[which(tolower(targets$targetType) == "target")]
  # For detectability metric only, also exclude noDetectability targets
  Targets_detect <- Targets
  if("noDetectability" %in% colnames(targets)){
    noDetectTargets <- targets$targetName[which(targets$noDetectability == TRUE)]
    Targets_detect <- Targets[!Targets %in% noDetectTargets]
  }
  AQtargets <- rownames(AQdata)

# Define column names for QCFlagReturn
  columns <- c("target", "flagName", "normMethod", "status", "val", "QCthreshold", "QCoperator", "QCformat")

  # Create copy of IPCnormed and set values to NA where aboveLOD is FALSE
  IPCnormed_filtered <- IPCnormed
  IPCnormed_filtered[!aboveLOD] <- NA
  meanSC_RQ <- apply(IPCnormed_filtered[Targets, SCs], 1, mean, na.rm = TRUE)
  sdevSC_RQ <- apply(IPCnormed_filtered[Targets, SCs], 1, sd, na.rm = TRUE)

  # Define evaluation items with their respective calculations
  eval_info <- list(
    list(name = "TARGET_MIN_READS", method = "raw", calc = apply(raw[Targets, Samples, drop=FALSE], 1, function(row){
      mean(row > MIN_TARGET_READS)
      })
    ),
    list(name = "TARGET_DETECTABILITY", method = "IPC", calc = detectability[Targets_detect]/100),
    list(name = "TARGET_CONC_CV_RQ", method = "IPC", calc = sdevSC_RQ / meanSC_RQ)#,
#    list(name = "TARGET_IPC_MIN_READS", method = "IPC", calc = apply(raw[Targets, IPCs], 1, median))
  )

  if(absRun){
    # Calculate mean and standard deviation for SC samples across AQdata
    AQdata_filtered <- AQdata
    if(!is.null(withinDR)){
      AQdata_filtered[!withinDR] <- NA
    }
    meanSC <- apply(AQdata_filtered[, SCs], 1, mean, na.rm = TRUE)
    sdevSC <- apply(AQdata_filtered[, SCs], 1, sd, na.rm = TRUE)
    SCparamsNull <- rep(NA, length(meanSC))
    names(SCparamsNull) <- names(meanSC)
    eval_info <- append(eval_info, list(
      list(name = "TARGET_CONC_ACCURACY", method = "AQ_aM", calc = if(is.null(SCparams)) SCparamsNull else (meanSC - SCparams) / SCparams),
      list(name = "TARGET_CONC_CV", method = "AQ_aM", calc = sdevSC / meanSC)
    ))
  }

  # Process all evaluation items into a single list with parameters
  eval_list <- lapply(seq_along(criteria$thresholdNames), function(i) {
    name <- names(criteria$properNames)[i]
    result_name <- criteria$thresholdNames[i]
    result_info <- eval_info[[which(sapply(eval_info, function(x) x$name == result_name))]]
    list(
      name = name,
      result = result_info$calc,
      method = result_info$method,
      op = criteria$operators[i],
      threshold = criteria$thresholds[i],
      format = criteria$format[i]
    )
  })

  # Pre-calculate total number of rows needed
  total_rows <- sum(sapply(eval_list, function(item) {
    if(item$name == "Target_Detectability") {
      sum(!is.na(item$result))
    } else {
      length(item$result)
    }
  }))

  QCFlagList <- vector("list", total_rows)
  row_idx <- 1
  # Loop through AQdata and evaluate criteria for each row
  for (j in 1:length(eval_list)) {
    eval_item <- eval_list[[j]]
    for (i in 1:length(eval_item$result)) {
      if(is.na(eval_item$result[i]) & eval_item$name =="Target_Detectability"){ # don't consider NA for Target_Detectability to be TRUE
        next;
      }
      # Evaluate the criterion for the current row
      set <- evalCriterion(eval_item$name, eval_item$result[i], eval_item$op, eval_item$threshold)
      if (is.na(set)) {
        set <- TRUE
      }
      # Add to list instead of rbind
      QCFlagList[[row_idx]] <- c(
        names(eval_item$result[i]),
        eval_item$name,
        eval_item$method,
        set,
        eval_item$result[i],
        eval_item$threshold,
        eval_item$op,
        eval_item$format
      )
      row_idx <- row_idx + 1
    }
  }
  # Remove any unused pre-allocated slots (in case of skipped NA values)
  if(row_idx <= total_rows) {
    QCFlagList <- QCFlagList[1:(row_idx-1)]
  }

  # Set column names for QCFlagReturn
  QCFlagReturn <- data.frame(do.call(rbind, QCFlagList), stringsAsFactors=FALSE)
  colnames(QCFlagReturn) <- columns

  # Convert column types as needed
  QCFlagReturn$val <- as.numeric(QCFlagReturn$val)
  QCFlagReturn$status <- as.logical(QCFlagReturn$status)
  return(QCFlagReturn)
}


#' Write Plate QC table
#'
#' Writes NULISAseq QC plate flags 
#' @param raw unnormalized data
#' @param normed normalized data
#' @param aboveLOD aboveLOD from lod function
#' @param ICs names of IC samples
#' @param NCs names of NC samples
#' @param IPCs names of IPC samples
#' @param SCs names of SC samples
#' @param AQ Whether this assay is AQ
#' @param AQ_QC AQ QC metrics, these are necessary to calculate a plate QC metric for AQ
#' @param Sample_QC Sample QC metrics, one of these is necessary to calculate a plate QC metric for AQ
#' @return QC table
#' @examples
#' # QCFlagPlate(inputtable)
#'
#' @export
QCFlagPlate <- function(raw, normed, aboveLOD, targets, samples, 
                        ICs=NULL, IPCs=NULL, NCs=NULL, SCs=NULL, AQ=TRUE, AQ_QC=NULL, Sample_QC=NULL){
  criteria <- QCPlateCriteria(AQ)
  ICs  <- if(!is.null(ICs))   ICs else targets$targetName[which(tolower(targets$targetType) == "control")]
  IPCs <- if(!is.null(IPCs)) IPCs else samples$sampleName[which(samples$sampleType == "IPC")]
  NCs  <- if(!is.null(NCs))   NCs else samples$sampleName[which(samples$sampleType == "NC")]
  SCs  <- if(!is.null(SCs))   SCs else samples$sampleName[which(samples$sampleType == "SC")]
  columns <- c("flagName", "normMethod", "status", "val", "QCthreshold", "QCoperator", "QCformat")
  
  max_rows <- if(AQ) 10 else 5
  QCFlagList <- vector("list", max_rows)
  row_idx <- 1
  # Calculate Plate-wide QC vals
  if(AQ){
    ## MAX_SC_CV
    SCvals <- raw[, SCs]
    SCvals[is.na(SCvals)] <- 0
    SCvals2 <- colSums(SCvals, na.rm=T)
    SC_CV <- sd(SCvals2, na.rm=T) / mean(SCvals2, na.rm=T)
    op <- criteria$operators[which(criteria$thresholdNames=="MAX_SC_CV")]
    format <- criteria$format[which(criteria$thresholdNames=="MAX_SC_CV")]
    set <- evalCriterion("SCRead_CV", SC_CV, op, MAX_SC_CV)
    QCFlagList[[row_idx]] <- c("SCRead_CV", "raw", set, SC_CV, MAX_SC_CV, op, format)
    row_idx <- row_idx + 1

    ## MAX_SC_Target_CV
    SCnormvals <- normed[, SCs]
    SCnormvals[is.na(SCvals)] <- 0
    median_SC_targetCV <- median(apply(SCnormvals, 1, sd) / rowMeans(SCnormvals, na.rm=T), na.rm=T) 
    op <- criteria$operators[which(criteria$thresholdNames=="MAX_MEDIAN_SC_TARGET_CV")]
    format <- criteria$format[which(criteria$thresholdNames=="MAX_MEDIAN_SC_TARGET_CV")]
    set <- evalCriterion("SCTarget_CV", median_SC_targetCV, op, MAX_MEDIAN_SC_TARGET_CV)
    QCFlagList[[row_idx]] <- c("SCTarget_CV", "IPC", set, median_SC_targetCV, MAX_MEDIAN_SC_TARGET_CV, op, format)
    row_idx <- row_idx + 1

    ## Failed Assays <10% of total AQ Targets 
    op <- criteria$operators[which(criteria$thresholdNames == "MAX_FAILED_TARGET_PERC")]
    format <- criteria$format[which(criteria$thresholdNames=="MAX_FAILED_TARGET_PERC")]
    val <- length(unique(AQ_QC[which(AQ_QC$status == "TRUE"),]$target)) / length(unique(AQ_QC$target))
    set <- evalCriterion("Failed_Targets", val, op, MAX_FAILED_TARGET_PERC) 
    QCFlagList[[row_idx]] <- c("Failed_Targets", "IPC", set, val, MAX_FAILED_TARGET_PERC, op, format)
    row_idx <- row_idx + 1

    ## Flagged sample control from sample QC
    op <- criteria$operators[which(criteria$thresholdNames == "MAX_FAILED_SC")]
    format <- criteria$format[which(criteria$thresholdNames=="MAX_FAILED_SC")]
    inds <- (which(Sample_QC$sampleName %in% SCs & Sample_QC$status == TRUE))
    val = length(unique(Sample_QC[inds,]$sampleName))
    set <- evalCriterion("Failed_SC", val, op, MAX_FAILED_SC) 
    QCFlagList[[row_idx]] <- c("Failed_SC", "IPC", set, val, MAX_FAILED_SC, op, format)
    row_idx <- row_idx + 1
    
    ## Flagged IPC from sample QC
    op <- criteria$operators[which(criteria$thresholdNames == "MAX_FAILED_IPC")]
    format <- criteria$format[which(criteria$thresholdNames=="MAX_FAILED_IPC")]
    inds <- (which(Sample_QC$sampleName %in% IPCs & Sample_QC$status == TRUE))
    val = length(unique(Sample_QC[inds,]$sampleName))
    set <- evalCriterion("Failed_IPC", val, op, MAX_FAILED_IPC) 
    QCFlagList[[row_idx]] <- c("Failed_IPC", "IPC", set, val, MAX_FAILED_IPC, op, format)
    row_idx <- row_idx + 1
  }

  ## MAX_IC_CV (V)
  ICvals <- raw[ICs,]
  ICvals[is.na(ICvals)] <- 0
  IC_CV <- sd(ICvals, na.rm=T) / mean(ICvals, na.rm=T)
  op <- criteria$operators[which(criteria$thresholdNames=="MAX_IC_CV")]
  format <- criteria$format[which(criteria$thresholdNames=="MAX_IC_CV")]
  set <- evalCriterion("ICRead_CV", IC_CV, op, MAX_IC_CV)
  QCFlagList[[row_idx]] <- c("ICRead_CV", "raw", set, IC_CV, MAX_IC_CV, op, format)
  row_idx <- row_idx + 1

  ## MAX_IPC_CV (I)
  IPCvals <- raw[, IPCs]
  IPCvals[is.na(IPCvals)] <- 0
  IPCvals2 <- colSums(IPCvals, na.rm=T)
  IPC_CV <- sd(IPCvals2, na.rm=T) / mean(IPCvals2, na.rm=T)
  op <- criteria$operators[which(criteria$thresholdNames=="MAX_IPC_CV")]
  format <- criteria$format[which(criteria$thresholdNames=="MAX_IPC_CV")]
  set <- evalCriterion("IPCRead_CV", IPC_CV, op, MAX_IPC_CV)
  QCFlagList[[row_idx]] <- c("IPCRead_CV", "raw", set, IPC_CV, MAX_IPC_CV, op, format)
  row_idx <- row_idx + 1

  ## MAX_MEDIAN_IPC_TARGET_CV (P)
  IPCnormvals <- normed[, IPCs]
  IPCnormvals[is.na(IPCvals)] <- 0
  median_IPC_targetCV <- median(apply(IPCnormvals, 1, function(x) sd(x, na.rm=TRUE)) / rowMeans(IPCnormvals, na.rm=T), na.rm=T) 
  op <- criteria$operators[which(criteria$thresholdNames=="MAX_MEDIAN_IPC_TARGET_CV")]
  format <- criteria$format[which(criteria$thresholdNames=="MAX_MEDIAN_IPC_TARGET_CV")]
  set <- evalCriterion("IPCTarget_CV", median_IPC_targetCV, op, MAX_MEDIAN_IPC_TARGET_CV)
  QCFlagList[[row_idx]] <- c("IPCTarget_CV", "IPC", set, median_IPC_targetCV, MAX_MEDIAN_IPC_TARGET_CV, op, format)
  row_idx <- row_idx + 1

  ## Detectability fraction (D)
  # Exclude ICs and noDetectability targets from detectability calculation
  noDetectTargets <- if("noDetectability" %in% colnames(targets)) targets$targetName[which(targets$noDetectability == TRUE)] else character(0)
  aboveLOD <- aboveLOD[!rownames(aboveLOD) %in% c(ICs, noDetectTargets), !colnames(aboveLOD) %in% c(IPCs, NCs, SCs)] # remove controls and noDetectability targets
  if(!is.null(dim(aboveLOD))){
    perc_tar <- rowSums(aboveLOD == TRUE, na.rm=TRUE)/ rowSums(!is.na(aboveLOD))
    perc_all <- sum(perc_tar > 0.5, na.rm=TRUE)/ nrow(aboveLOD)
  } else{
    perc_all <- mean(aboveLOD, na.rm=TRUE)
  }
  op <- criteria$operators[which(criteria$thresholdNames=="DETECTABILITY_FRAC")]
  format <- criteria$format[which(criteria$thresholdNames=="DETECTABILITY_FRAC")]
  set <- evalCriterion("Detectability", perc_all, op, DETECTABILITY_FRAC)
  QCFlagList[[row_idx]] <- c("Detectability", "IPC", set, perc_all, DETECTABILITY_FRAC, op, format)
  row_idx <- row_idx + 1

  ## Min number of reads (R)
  nReads <- sum(raw, na.rm=T)
  op <- criteria$operators[which(criteria$thresholdNames=="MIN_READS")]
  format <- criteria$format[which(criteria$thresholdNames=="MIN_READS")]
  set <- evalCriterion("MinReads", nReads, op, MIN_READS)
  QCFlagList[[row_idx]] <- c("MinReads", "raw", set, nReads, MIN_READS, op, format)

  # Trim list if needed (e.g. when AQ=FALSE)
  if(row_idx <= max_rows){
    QCFlagList <- QCFlagList[1:row_idx]
  }

  # Convert to data.frame
  QCFlagReturn <- as.data.frame(do.call(rbind, QCFlagList), stringsAsFactors = FALSE)
  colnames(QCFlagReturn) <- columns
  QCFlagReturn$val <- as.numeric(QCFlagReturn$val)
  QCFlagReturn$QCthreshold <- as.numeric(QCFlagReturn$QCthreshold)
  return(QCFlagReturn)
}
