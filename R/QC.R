# Sample QC criteria
MIN_FRAC_DETECTABILITY <- 0.8  # Minimim fraction (Target_Detectability): # Targets with reads above LOD
MIN_IC_READS_PER_SAMPLE <- 1000    # Minimum number (ICReads) of IC reads within a sample
MIN_NUM_READS_PER_SAMPLE <- 500000 # Minimum number (NumReads) of reads within a sample
MIN_IC_MEDIAN <- "-0.3,0.3"# +/- of sample IC read count about the median

#' QCSampleCriteria
#'
#' Sample QC Criteria
#' @return criteria
#' @examples
#' # QCSampleCriteria()
#'
#' @export
QCSampleCriteria <- function(){
  retVal <- NULL
  retVal$thresholds <-c(Detectability=MIN_FRAC_DETECTABILITY,
                        ICReads=MIN_IC_READS_PER_SAMPLE, 
                        NumReads=MIN_NUM_READS_PER_SAMPLE, 
                        IC_Median=MIN_IC_MEDIAN)
  retVal$operators <-c(Detectability="<",
                       ICReads="<", 
                       NumReads="<", 
                       IC_Median="<,>")
  retVal$format <- c(Detectability="percentage",
                       ICReads="integer", 
                       NumReads="integer", 
                       IC_Median="percentage")
  retVal$thresholdNames<-c(Detectability="MIN_FRAC_DETECTABILITY", 
                        ICReads="MIN_IC_READS_PER_SAMPLE", 
                        NumReads="MIN_NUM_READS_PER_SAMPLE", 
                        IC_Median="MIN_IC_MEDIAN")
  return(retVal)
}

# Plate QC criteria
MAX_IC_CV <- 0.5                # (ICRead_CV) CV of IC reads across all samples
MAX_IPC_CV <- 0.5               # (IPCRead_CV) CV of total read count for each IPC sample
MAX_MEDIAN_IPC_TARGET_CV <- 0.2 # (IPCTarget_CV) Median of CVs of all IPC targets (Performed on normalized data) 
DETECTABILITY_FRAC <- 0.75      # (Detectability) Target-wise Detectability fraction (target is detected if >50% of samples > LOD)
MIN_READS <- 5e5                # (MinReads) Minimum number of reads

#' QCPlateCriteria
#'
#' Plate QC Criteria
#' @return criteria
#' @examples
#' # QCPlateCriteria()
#'
#' @export
QCPlateCriteria <- function(){
  retVal <- NULL
  retVal$thresholds <- c(ICRead_CV=as.numeric(MAX_IC_CV),
                         IPCRead_CV=as.numeric(MAX_IPC_CV), 
                         IPCTarget_CV=as.numeric(MAX_MEDIAN_IPC_TARGET_CV), 
                         Detectability=as.numeric(DETECTABILITY_FRAC), 
                         MinReads=as.numeric(MIN_READS))
  retVal$operators <- c(ICRead_CV=">",IPCRead_CV=">", IPCTarget_CV=">", Detectability="<", MinReads="<")
  retVal$format <- c(ICRead_CV="percentage",IPCRead_CV="percentage", IPCTarget_CV="percentage", Detectability="percentage", MinReads="integer")
  retVal$thresholdNames<-c(IC_Read_CV="MAX_IC_CV",
                         IPCRead_CV="MAX_IPC_CV", 
                         IPCTarget_CV="MAX_MEDIAN_IPC_TARGET_CV", 
                         Detectability="DETECTABILITY_FRAC", 
                         MinReads="MIN_READS")
  return(retVal)
}


#' Write Processed XML from QC table
#'
#' Writes NULISAseq QC flags to XML
#' @param input QC table 
#' @param sample whether the table is sample or plate
#' @return QC XML
#' @examples
#' # QC2XML(inputtable)
#'
#' @export
QC2XML <- function(input, QCNode, sample=F, combined=F){
  for (i in 1:nrow(input)){
    if(sample){
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
  }
  return(QCNode)
}

#' Write Sample QC table
#'
#' Writes NULISAseq QC sample flags 
#' @param raw unnormalized data
#' @param normed normalized data
#' @param ICs indices of IC samples
#' @param NCs indices of NC samples
#' @param IPCs indices of IPC samples
#' @param SCs indices of SC samples
#' @return QC table
#' @examples
#' # QCFlagSample(inputtable)
#'
#' @export

QCFlagSample <- function(raw, normed, samples, targets, 
                         well_order=NULL, ICs=NULL, IPCs=NULL, NCs=NULL, SCs=NULL){
  columns <- c("sampleName", "flagName", "normMethod", "status", "val", "text", 
               "sampleBarcode", "sampleType", "QCthreshold", "QCoperator", "QCformat")

  criteria <- QCSampleCriteria()

  QCFlagReturn <- data.frame(matrix(nrow=0, ncol=length(columns)))
  if(is.null(well_order)){
    well_order <- 1:length(raw[1,])
  }else{
    well_order <- rev(well_order)
  }
  ICs  <- if(!is.null(ICs))   ICs else which(targets$targetType == "Control")
  IPCs <- if(!is.null(IPCs)) IPCs else which(samples$sampleType == "IPC")
  NCs  <- if(!is.null(NCs))   NCs else which(samples$sampleType == "NC")
  SCs  <- if(!is.null(SCs))   SCs else which(samples$sampleType == "SC")

  # Median IC between -30% and 30% of median  
  mCherry_median <- median(raw[ICs[1],], na.rm=T)
  min_ic_median <- unname(unlist(strsplit(MIN_IC_MEDIAN, ",")))
  medianMin30 <- mCherry_median - mCherry_median * abs(as.numeric(min_ic_median[0]))
  medianMax30 <- mCherry_median + mCherry_median * abs(as.numeric(min_ic_median[1]))
  medVals <- (raw[ICs[1], ] - mCherry_median ) / mCherry_median 
  op <- criteria$operators[which(criteria$thresholdNames=="MIN_IC_MEDIAN")]
  format <- criteria$format[which(criteria$thresholdNames=="MIN_IC_MEDIAN")]
  for(j in 1:length(medVals)){
    i <- well_order[j]
    set <- if(is.na(medVals[i])) "T" else evalCriterion("IC_Median", medVals[i], op, MIN_IC_MEDIAN)
    type <- "Sample"
    if(i %in% NCs){
      type <- "NC"
    } else if (i %in% IPCs){
      type <- "IPC"
    } else if (i %in% SCs){
      type <- "SC"
    }
    QCFlagReturn <- rbind(QCFlagReturn, c(names(medVals)[i], "IC_Median", "raw", set, medVals[i], "", samples$sampleBarcode[i], type, 
                                          as.character(paste(MIN_IC_MEDIAN, collapse=',')), as.character(op), format))
  }

  # Minimim fraction (Target_Detectability): # Targets with reads above LOD
  normed2 <- normed
  
  lod <- lod(data_matrix=normed2, blanks=NCs, min_count=0)
  lod$aboveLOD[which(is.na(lod$aboveLOD))] <- FALSE
  perc_tar <- colSums(lod$aboveLOD == TRUE)/ nrow(lod$aboveLOD)
  perc_tar[NCs] <-NA
  op <- criteria$operators[which(criteria$thresholdNames=="MIN_FRAC_DETECTABILITY")]
  format <- criteria$format[which(criteria$thresholdNames=="MIN_FRAC_DETECTABILITY")]
  for (j in 1:length(perc_tar)){
    i <- well_order[j]
    set <- if(is.na(perc_tar[i])) NA else evalCriterion("Detectability", perc_tar[i], op,  MIN_FRAC_DETECTABILITY)
    type <- "Sample"
    if(i %in% NCs){
      type <- "NC"
    } else if (i %in% IPCs){
      type <- "IPC"
    } else if (i %in% SCs){
      type <- "SC"
    }
    QCFlagReturn <- rbind(QCFlagReturn, c(names(perc_tar)[i], "Detectability", "IC", set, perc_tar[i], "", samples$sampleBarcode[i], type, 
                                          as.character(MIN_FRAC_DETECTABILITY), as.character(op), format))
  }

  
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
    QCFlagReturn <- rbind(QCFlagReturn, c(colnames(raw)[i], "ICReads", "raw", set, ICvals[i], "", samples$sampleBarcode[i], type, 
                                          as.character(MIN_IC_READS_PER_SAMPLE), as.character(op), format))
  }

  # Minimum number (NumReads) of reads within a sample
  raw2 <- raw
  barNames <- samples$sampleBarcode
  val <- colSums(raw2, na.rm=T)
  op <- criteria$operators[which(criteria$thresholdNames=="MIN_NUM_READS_PER_SAMPLE")]
  format <- criteria$format[which(criteria$thresholdNames=="MIN_NUM_READS_PER_SAMPLE")]
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
    QCFlagReturn <- rbind(QCFlagReturn, c(colnames(raw)[i], "NumReads", "raw", set, val[i], "", barNames[i], type, 
                                          as.character(MIN_NUM_READS_PER_SAMPLE), as.character(op), format))
  }
  colnames(QCFlagReturn) <- columns
  QCFlagReturn$val <- as.numeric(QCFlagReturn$val)
  QCFlagReturn$QCthreshold <- as.character(QCFlagReturn$QCthreshold)
  return(QCFlagReturn)
}

# Return TRUE if criteria is violated (fail)
evalCriterion <- function(name, value, operator, threshold){
  operators <- unlist(strsplit(unname(as.character(operator)), ","))
  thresholds <- unlist(strsplit(unname(as.character(threshold)), ","))
  if(length(operators) != length(thresholds)){
    stop(paste("Error: Specification of QC criteria is incorrect, unequal number of operators and thresholds for ", name)) 
  }else{
    input <- data.frame(rep(value, length(operators)), operators, thresholds)
    flags <- apply(input, 1, function(x) eval(parse(text=paste(x[1], x[2], x[3]))))
    return (any(flags))
  }
}

#' Write Plate QC table
#'
#' Writes NULISAseq QC plate flags 
#' @param raw unnormalized data
#' @param normednormalized data
#' @param ICs indices of IC samples
#' @param NCs indices of NC samples
#' @param IPCs indices of IPC samples
#' @param SCs indices of SC samples
#' @return QC table
#' @examples
#' # QCFlagPlate(inputtable)
#'
#' @export
QCFlagPlate <- function(raw, normed, targets, samples, 
                        ICs=NULL, IPCs=NULL, NCs=NULL, SCs=NULL){
  criteria <- QCPlateCriteria()
  ICs  <- if(!is.null(ICs))   ICs else which(targets$targetType == "Control")
  IPCs <- if(!is.null(IPCs)) IPCs else which(samples$sampleType == "IPC")
  NCs  <- if(!is.null(NCs))   NCs else which(samples$sampleType == "NC")
  SCs  <- if(!is.null(SCs))   SCs else which(samples$sampleType == "SC")
  columns <- c("flagName", "normMethod", "status", "val", "QCthreshold", "QCoperator", "QCformat")
  QCFlagReturn <- data.frame(matrix(nrow=0, ncol=length(columns)))

  # Calculate Plate-wide QC vals
  ## MAX_IC_CV (V)
  ICvals <- raw[ICs,]
  ICvals[is.na(ICvals)] <- 0
  IC_CV <- sd(ICvals, na.rm=T) / mean(ICvals, na.rm=T)
  op <- criteria$operators[which(criteria$thresholdNames=="MAX_IC_CV")]
  format <- criteria$format[which(criteria$thresholdNames=="MAX_IC_CV")]
  set <- evalCriterion("ICRead_CV", IC_CV, op, MAX_IC_CV)
  QCFlagReturn <- rbind(QCFlagReturn, c("ICRead_CV", "raw", set, IC_CV, MAX_IC_CV, op, format))

  ## MAX_IPC_CV (I)
  IPCvals <- raw[, IPCs]
  IPCvals[is.na(IPCvals)] <- 0
  IPCvals2 <- colMeans(IPCvals, na.rm=T)
  IPC_CV <- sd(IPCvals2, na.rm=T) / mean(IPCvals2, na.rm=T)
  op <- criteria$operators[which(criteria$thresholdNames=="MAX_IPC_CV")]
  format <- criteria$format[which(criteria$thresholdNames=="MAX_IPC_CV")]
  set <- evalCriterion("IPCRead_CV", IPC_CV, op, MAX_IPC_CV)
  QCFlagReturn <- rbind(QCFlagReturn, c("IPCRead_CV", "raw", set, IPC_CV, MAX_IPC_CV, op, format))

  ## MAX_MEDIAN_IPC_TARGET_CV (P)
  IPCnormvals <- normed[, IPCs]
  IPCnormvals[is.na(IPCvals)] <- 0
  median_IPC_targetCV <- median(apply(IPCnormvals, 1, sd) / rowMeans(IPCnormvals, na.rm=T), na.rm=T) 
  op <- criteria$operators[which(criteria$thresholdNames=="MAX_MEDIAN_IPC_TARGET_CV")]
  format <- criteria$format[which(criteria$thresholdNames=="MAX_MEDIAN_IPC_TARGET_CV")]
  set <- evalCriterion("IPCTarget_CV", median_IPC_targetCV, op, MAX_MEDIAN_IPC_TARGET_CV)
  QCFlagReturn <- rbind(QCFlagReturn, c("IPCTarget_CV", "IC", set, median_IPC_targetCV, MAX_MEDIAN_IPC_TARGET_CV, op, format))

  ## Detectability fraction (D)
  normed2 <- normed
  normed2 <- normed2[-ICs,]
  lod <- lod(data_matrix=normed2, blanks=NCs, min_count=0)
  lod$aboveLOD[which(is.na(lod$aboveLOD))] <- FALSE
  perc_tar <- rowSums(lod$aboveLOD == TRUE)/ ncol(lod$aboveLOD)
  perc_all <- length(which(perc_tar > 0.5))/ nrow(lod$aboveLOD)
  op <- criteria$operators[which(criteria$thresholdNames=="DETECTABILITY_FRAC")]
  format <- criteria$format[which(criteria$thresholdNames=="DETECTABILITY_FRAC")]
  set <- evalCriterion("Detectability", perc_all, op, DETECTABILITY_FRAC)
  QCFlagReturn <- rbind(QCFlagReturn, c("Detectability", "IC", set, perc_all, DETECTABILITY_FRAC, op, format))

  ## Min number of reads (R)
  nReads <- sum(raw, na.rm=T)
  op <- criteria$operators[which(criteria$thresholdNames=="MIN_READS")]
  format <- criteria$format[which(criteria$thresholdNames=="MIN_READS")]
  set <- evalCriterion("MinReads", nReads, op, MIN_READS)
  QCFlagReturn <- rbind(QCFlagReturn, c("MinReads", "raw", set, nReads, MIN_READS, op, format))
  colnames(QCFlagReturn) <- columns
  QCFlagReturn$val <- as.numeric(QCFlagReturn$val)
  QCFlagReturn$QCthreshold <- as.numeric(QCFlagReturn$QCthreshold)
  return(QCFlagReturn)
}
