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
QC2XML <- function(input, QCNode, sample=F){
  for (i in 1:nrow(input)){
    if(sample){
      addChildren(QCNode, newXMLNode("QCFlag", input$val[i],
                                           attrs=c(
                                                  name=input$flagName[i],
                                                  set=input$status[i],
                                                  method=input$normMethod[i]
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
QCFlagSample <- function(raw, normed, ICs, NCs, IPCs, samples){
  columns <- c("sampleName", "flagName", "normMethod", "status", "val", "text", "sampleBarcode")
  QCFlagReturn <- data.frame(matrix(nrow=0, ncol=length(columns)))
  # Sample QC criteria
  MIN_FRAC_TARGETS_ABOVE_LOD <- 0.8  # Minimim fraction (Target_Detectability): # Targets with reads above LOD
  MIN_IC_READS_PER_SAMPLE <- 1000    # Minimum number (ICReads) of IC reads within a sample
  MIN_NUM_READS_PER_SAMPLE <- 500000 # Minimum number (NumReads) of reads within a sample

  # Minimim fraction (Target_Detectability): # Targets with reads above LOD
  normed2 <- normed
  normed2 <- normed2[-ICs, ]
  
  lod <- lod(data_matrix=normed2, blanks=NCs, min_count=0)
  lod$aboveLOD[which(is.na(lod$aboveLOD))] <- FALSE
  perc_tar <- colSums(lod$aboveLOD == TRUE)/ nrow(lod$aboveLOD)
  perc_tar <- perc_tar[-NCs]
  for (i in 1:length(perc_tar)){
    set <- if(perc_tar[i] < MIN_FRAC_TARGETS_ABOVE_LOD) "T" else "F"
    QCFlagReturn <- rbind(QCFlagReturn, c(names(perc_tar)[i], "TARGETS_ABOVE_LOD", "IC", set, perc_tar[i], "", samples$sampleBarcode[i]))
  }

  
  # Minimum number (ICReads) of IC reads within a sample
  ICvals <- raw[ICs, ]
  ICvals[is.na(ICvals)] <- 0
  for (i in 1:length(ICvals)){
    set <- if(ICvals[i] < MIN_IC_READS_PER_SAMPLE) "T" else "F"
    QCFlagReturn <- rbind(QCFlagReturn, c(colnames(raw)[i], "ICReads", "raw", set, ICvals[i], "", samples$sampleBarcode[i]))
  }

  # Minimum number (NumReads) of reads within a sample
  raw2 <- raw
  raw2 <- raw2[,-NCs]
  barNames <- samples$sampleBarcode
  barNames <- barNames[-NCs] 
  val <- colSums(raw2, na.rm=T)
  for(i in 1:length(val)){
    set <- if(val[i] < MIN_NUM_READS_PER_SAMPLE) "T" else "F"
    QCFlagReturn <- rbind(QCFlagReturn, c(colnames(raw)[i], "NumReads", "raw", set, val[i], "", barNames[i]))
  }
  colnames(QCFlagReturn) <- columns
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
QCFlagPlate <- function(raw, normed, ICs, NCs, IPCs){
  columns <- c("flagName", "normMethod", "status", "val", "text")
  QCFlagReturn <- data.frame(matrix(nrow=0, ncol=length(columns)))

  # Plate QC criteria
  MAX_IC_CV <- 0.5                # (ICRead_CV) CV of IC reads across all samples
  MAX_IPC_CV <- 0.5               # (IPCRead_CV) CV of total read count for each IPC sample
  MAX_MEDIAN_IPC_TARGET_CV <- 0.2 # (IPCTarget_CV) Median of CVs of all IPC targets (Performed on normalized data) 
  DETECTABILITY_FRAC <- 0.75      # (Detectability) Target-wise Detectability fraction (target is detected if >50% of samples > LOD)
  MIN_READS <- 5e5                # (MinReads) Minimum number of reads

  # Calculate Plate-wide QC vals
  ## MAX_IC_CV (V)
  ICvals <- raw[ICs,]
  ICvals[is.na(ICvals)] <- 0
  IC_CV <- sd(ICvals, na.rm=T) / mean(ICvals, na.rm=T)
  set <- if(IC_CV > MAX_IC_CV) "T" else "F"
  QCFlagReturn <- rbind(QCFlagReturn, c("ICRead_CV", "raw", set, IC_CV, ""))

  ## MAX_IPC_CV (I)
  IPCvals <- raw[, IPCs]
  IPCvals[is.na(IPCvals)] <- 0
  IPCvals2 <- colMeans(IPCvals, na.rm=T)
  IPC_CV <- sd(IPCvals2, na.rm=T) / mean(IPCvals2, na.rm=T)
  set <- if(IPC_CV > MAX_IPC_CV) "T" else "F"
  QCFlagReturn <- rbind(QCFlagReturn, c("IPCRead_CV", "raw", set, IPC_CV, ""))

  ## MAX_MEDIAN_IPC_TARGET_CV (P)
  IPCnormvals <- normed[, IPCs]
  IPCnormvals[is.na(IPCvals)] <- 0
  median_IPC_targetCV <- median(apply(IPCnormvals, 1, sd) / rowMeans(IPCnormvals, na.rm=T), na.rm=T) 
  set <- if(median_IPC_targetCV > MAX_MEDIAN_IPC_TARGET_CV) "T" else "F"
  QCFlagReturn <- rbind(QCFlagReturn, c("IPCTarget_CV", "IC", set, median_IPC_targetCV, ""))

  ## Detectability fraction (D)
  normed2 <- normed
  normed2 <- normed2[-ICs,]
  lod <- lod(data_matrix=normed2, blanks=NCs, min_count=0)
  lod$aboveLOD[which(is.na(lod$aboveLOD))] <- FALSE
  perc_tar <- rowSums(lod$aboveLOD == TRUE)/ ncol(lod$aboveLOD)
  perc_all <- length(which(perc_tar > 0.5))/ nrow(lod$aboveLOD)
  set <- if(perc_all < DETECTABILITY_FRAC) "T" else "F"
  QCFlagReturn <- rbind(QCFlagReturn, c("Detectability", "IC", set, perc_all, ""))

  ## Min number of reads (R)
  nReads <- sum(raw, na.rm=T)
  set <- if(nReads < MIN_READS) "T" else "F"
  QCFlagReturn <- rbind(QCFlagReturn, c("MinReads", "raw", set, nReads, ""))

  colnames(QCFlagReturn) <- columns
  return(QCFlagReturn)
}