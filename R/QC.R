# Sample QC criteria
MIN_FRAC_DETECTABILITY <- 0.8  # Minimim fraction (Target_Detectability): # Targets with reads above LOD
MIN_IC_READS_PER_SAMPLE <- 1000    # Minimum number (ICReads) of IC reads within a sample
MIN_NUM_READS_PER_SAMPLE <- 500000 # Minimum number (NumReads) of reads within a sample
MIN_IC_MEDIAN <- 0.3 # +/- of sample IC read count about the median

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
  retVal$thresholds <-c(Detectability=as.numeric(MIN_FRAC_DETECTABILITY),
                        ICReads=as.numeric(MIN_IC_READS_PER_SAMPLE), 
                        NumReads=as.numeric(MIN_NUM_READS_PER_SAMPLE), 
                        IC_Median=paste0(-MIN_IC_MEDIAN, ",", MIN_IC_MEDIAN))
  retVal$operators <-c(Detectability=">",
                       ICReads=">", 
                       NumReads=">", 
                       IC_Median=">,<")
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
  retVal$operator <- c(ICRead_CV=">",IPCRead_CV=">", IPCTarget_CV=">", Detectability=">", MinReads=">")
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
                                                  rep=input$sampleBarcode[i]
                                                  )
                )
        )
      }else{
        addChildren(QCNode, newXMLNode("QCFlag", input$val[i],
                                           attrs=c(
                                                  name=input$flagName[i],
                                                  set=input$status[i],
                                                  method=input$normMethod[i]
                                                  )
                )
        )
      }
    }else{
      addChildren(QCNode, newXMLNode("QCFlag", input$val[i],
                                           attrs=c(
                                                  name=input$flagName[i],
                                                  set=input$status[i],
                                                  method=input$normMethod[i]
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
#' @param normednormalized data
#' @param ICs indices of IC samples
#' @param NCs indices of NC samples
#' @param IPCs indices of IPC samples
#' @return QC table
#' @examples
#' # QCFlagSample(inputtable)
#'
#' @export
QCFlagSample <- function(raw, normed, samples, targets, well_order=NULL, ICs=NULL, IPCs=NULL, NCs=NULL){
  columns <- c("sampleName", "flagName", "normMethod", "status", "val", "text", "sampleBarcode", "sampleType", "QCthreshold")
  QCFlagReturn <- data.frame(matrix(nrow=0, ncol=length(columns)))
  if(is.null(well_order)){
    well_order <- 1:length(raw[1,])
  }else{
    well_order <- rev(well_order)
  }
  ICs  <- if(!is.null(ICs))   ICs else which(targets$targetType == "Control")
  IPCs <- if(!is.null(IPCs)) IPCs else which(samples$sampleType == "IPC")
  NCs  <- if(!is.null(NCs))   NCs else which(samples$sampleType == "NC")
  # Median IC between -30% and 30% of median  
  mCherry_median <- median(raw[ICs,], na.rm=T)
  medianMin30 <- mCherry_median - mCherry_median * MIN_IC_MEDIAN
  medianMax30 <- mCherry_median + mCherry_median * MIN_IC_MEDIAN
  medVals <- (raw[ICs, ] - mCherry_median ) / mCherry_median 
  for(j in 1:length(medVals)){
    i <- well_order[j]
    set <- NULL
    if(is.na(medVals[i])){
      set <- "T"
    }else{
      set <- if(abs(medVals[i]) > MIN_IC_MEDIAN) "T" else "F"
    }
    type <- "Sample"
    if(i %in% NCs){
      type <- "NC"
    }else if (i %in% IPCs){
      type <- "IPC"
    }
    QCFlagReturn <- rbind(QCFlagReturn, c(names(medVals)[i], "IC_Median", "raw", set, medVals[i], "", samples$sampleBarcode[i], type, MIN_IC_MEDIAN))

  }
  # Minimim fraction (Target_Detectability): # Targets with reads above LOD
  normed2 <- normed
  
  lod <- lod(data_matrix=normed2, blanks=NCs, min_count=0)
  lod$aboveLOD[which(is.na(lod$aboveLOD))] <- FALSE
  perc_tar <- colSums(lod$aboveLOD == TRUE)/ nrow(lod$aboveLOD)
  perc_tar[NCs] <-NA
  for (j in 1:length(perc_tar)){
    i <- well_order[j]
    set <- NULL
    if(is.na(perc_tar[i])){
      set <- NA
    }else{
      set <- if(perc_tar[i] < MIN_FRAC_DETECTABILITY) "T" else "F"
    }
    type <- "Sample"
    if(i %in% NCs){
      type <- "NC"
    }else if (i %in% IPCs){
      type <- "IPC"
    }
    QCFlagReturn <- rbind(QCFlagReturn, c(names(perc_tar)[i], "Detectability", "IC", set, perc_tar[i], "", samples$sampleBarcode[i], type, MIN_FRAC_DETECTABILITY))
  }

  
  # Minimum number (ICReads) of IC reads within a sample
  ICvals <- raw[ICs, ]
  ICvals[is.na(ICvals)] <- 0
  for (j in 1:length(ICvals)){
    i <- well_order[j]
    set <- if(ICvals[i] < MIN_IC_READS_PER_SAMPLE) "T" else "F"
    type <-"Sample"
    if(i %in% NCs){
      type <- "NC"
    }else if (i %in% IPCs){
      type <- "IPC"
    }
    QCFlagReturn <- rbind(QCFlagReturn, c(colnames(raw)[i], "ICReads", "raw", set, ICvals[i], "", samples$sampleBarcode[i], type, MIN_IC_READS_PER_SAMPLE))
  }

  # Minimum number (NumReads) of reads within a sample
  raw2 <- raw
  barNames <- samples$sampleBarcode
  val <- colSums(raw2, na.rm=T)
  for(j in 1:length(val)){
    i <- well_order[j]
    set <- NULL
    if(i %in% NCs){
      set <- NA  
      val[i] <- NA
    }else{
      set <- if(val[i] < MIN_NUM_READS_PER_SAMPLE) "T" else "F"
    }
    type <- "Sample"
    if(i %in% NCs){
      type <- "NC"
    }else if (i %in% IPCs){
      type <- "IPC"
    }
    QCFlagReturn <- rbind(QCFlagReturn, c(colnames(raw)[i], "NumReads", "raw", set, val[i], "", barNames[i], type, MIN_NUM_READS_PER_SAMPLE))
  }
  colnames(QCFlagReturn) <- columns
  QCFlagReturn$val <- as.numeric(QCFlagReturn$val)
  QCFlagReturn$QCthreshold <- as.numeric(QCFlagReturn$QCthreshold)
  return(QCFlagReturn)
}

#' Write Plate QC table
#'
#' Writes NULISAseq QC plate flags 
#' @param raw unnormalized data
#' @param normednormalized data
#' @param ICs indices of IC samples
#' @param NCs indices of NC samples
#' @param IPCs indices of IPC samples
#' @return QC table
#' @examples
#' # QCFlagPlate(inputtable)
#'
#' @export
QCFlagPlate <- function(raw, normed, targets, samples, ICs=NULL, IPCs=NULL, NCs=NULL){
  ICs  <- if(!is.null(ICs))   ICs else which(targets$targetType == "Control")
  IPCs <- if(!is.null(IPCs)) IPCs else which(samples$sampleType == "IPC")
  NCs  <- if(!is.null(NCs))   NCs else which(samples$sampleType == "NC")
  columns <- c("flagName", "normMethod", "status", "val", "QCthreshold", "QCoperator")
  QCFlagReturn <- data.frame(matrix(nrow=0, ncol=length(columns)))

  # Calculate Plate-wide QC vals
  ## MAX_IC_CV (V)
  ICvals <- raw[ICs,]
  ICvals[is.na(ICvals)] <- 0
  IC_CV <- sd(ICvals, na.rm=T) / mean(ICvals, na.rm=T)
  set <- if(IC_CV > MAX_IC_CV) "T" else "F"
  QCFlagReturn <- rbind(QCFlagReturn, c("ICRead_CV", "raw", set, IC_CV, MAX_IC_CV, ">"))

  ## MAX_IPC_CV (I)
  IPCvals <- raw[, IPCs]
  IPCvals[is.na(IPCvals)] <- 0
  IPCvals2 <- colMeans(IPCvals, na.rm=T)
  IPC_CV <- sd(IPCvals2, na.rm=T) / mean(IPCvals2, na.rm=T)
  set <- if(IPC_CV > MAX_IPC_CV) "T" else "F"
  QCFlagReturn <- rbind(QCFlagReturn, c("IPCRead_CV", "raw", set, IPC_CV, MAX_IPC_CV, ">"))

  ## MAX_MEDIAN_IPC_TARGET_CV (P)
  IPCnormvals <- normed[, IPCs]
  IPCnormvals[is.na(IPCvals)] <- 0
  median_IPC_targetCV <- median(apply(IPCnormvals, 1, sd) / rowMeans(IPCnormvals, na.rm=T), na.rm=T) 
  set <- if(median_IPC_targetCV > MAX_MEDIAN_IPC_TARGET_CV) "T" else "F"
  QCFlagReturn <- rbind(QCFlagReturn, c("IPCTarget_CV", "IC", set, median_IPC_targetCV, MAX_MEDIAN_IPC_TARGET_CV, ">"))

  ## Detectability fraction (D)
  normed2 <- normed
  normed2 <- normed2[-ICs,]
  lod <- lod(data_matrix=normed2, blanks=NCs, min_count=0)
  lod$aboveLOD[which(is.na(lod$aboveLOD))] <- FALSE
  perc_tar <- rowSums(lod$aboveLOD == TRUE)/ ncol(lod$aboveLOD)
  perc_all <- length(which(perc_tar > 0.5))/ nrow(lod$aboveLOD)
  set <- if(perc_all < DETECTABILITY_FRAC) "T" else "F"
  QCFlagReturn <- rbind(QCFlagReturn, c("Detectability", "IC", set, perc_all, DETECTABILITY_FRAC, ">"))

  ## Min number of reads (R)
  nReads <- sum(raw, na.rm=T)
  set <- if(nReads < MIN_READS) "T" else "F"
  QCFlagReturn <- rbind(QCFlagReturn, c("MinReads", "raw", set, nReads, MIN_READS, ">"))
  colnames(QCFlagReturn) <- columns
  QCFlagReturn$val <- as.numeric(QCFlagReturn$val)
  QCFlagReturn$QCthreshold <- as.numeric(QCFlagReturn$QCthreshold)
  return(QCFlagReturn)
}
