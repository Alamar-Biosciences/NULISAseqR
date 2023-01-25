base_XML <- function(ExecutionDetails, bcodeA, bcodeB, RunSummary){
  children <- c()
  for(i in 1:length(ExecutionDetails)){
    val <- if(length(ExecutionDetails[[i]]) == 1) ExecutionDetails[i] else ExecutionDetails[[i]][[1]]
    children <- c(children, newXMLNode(names(ExecutionDetails)[i], val))
  }
  return( newXMLNode("NULISAseq", 
                      .children=c(
                        newXMLNode("ExecutionDetails", 
                          .children=children
                        ),
                        newXMLNode("plateID"),
                        newXMLNode("RunSummary",
                          .children=c(
                            newXMLNode("Barcodes", 
                              .children=c(
                                bcodeA,
                                bcodeB
                              )
                            ),
                            newXMLNode("Tags",
                              .children=c(
                                newXMLNode("TagA"),
                                newXMLNode("TagB")
                              )
                            ),
                            newXMLNode("TotalReads", RunSummary$TotalReads),
                            newXMLNode("Parseable", RunSummary$Parseable),
                            newXMLNode("ParseableMatch", RunSummary$ParseableMatch),
                            newXMLNode("Unparseable", RunSummary$Unparseable)
                          )
                        )
                      )
                    )
  )
}

bcodeA_XML <- function(targets){
  # Process BarcodeAs
  bcodeA <- newXMLNode("BarcodeA")
  for(i in 1:length(targets$targetBarcode)){
    addChildren(bcodeA, newXMLNode("Barcode", targets$targetName[i], attrs=c(name=targets$targetBarcode[i]))) 
  }
  return(bcodeA)
}
bcodeB_XML <- function(samples, barcodeB=NULL){
  # Process BarcodeBs
  inds <- c()
  sampColI <- names(samples)
  for( i in 1:length(sampColI)){
    if (sampColI[i] != "sampleBarcode" && sampColI[i] != "sampleName" && 
      sampColI[i] !="matching" && sampColI[i] != "non-matching"){
      inds <- c(inds, i)
    }
  }
  barcodeB_file <- NULL
  if(!is.null(barcodeB)){
    barcodeB_file <- read.table(barcodeB, sep='\t', comment.char='&', header=T, na.strings=c())
  }
  bcodeB <- newXMLNode("BarcodeB")
  for(i in 1:length(samples$sampleBarcode)){
    attrVals <- attrNames <- c()
    if(is.null(barcodeB_file)){
      for (j in 1:length(inds)){
        attrNames <- c(attrNames, sampColI[inds[j]])
        attrVals <- c(attrVals, eval(parse(text=paste0("(samples$", sampColI[inds[j]],")[i]"))))
      }
      attrVals <- c(attrVals, samples$sampleBarcode[i])
      attrNames <- c(attrNames, "name")
    }else{
      info <- barcodeB_file[which(samples$sampleBarcode[i] == barcodeB_file[,1]),]

      attrVals <- c(info[2])
      attrNames <- c("name")
      if (length(info)>2){
        for( k in 3:length(info)){
          val <- colnames(barcodeB_file)[k]
          attrVals <- c(attrVals, info[k])
          attrNames <- c(attrNames, names(info)[k])
        }
      }
    }
    names(attrVals) <- attrNames
    addChildren(bcodeB, newXMLNode("Barcode", rownames(info), attrs=attrVals)) 
  }
  return(bcodeB)
}

NCBkgdLevels_XML <- function(Data, ICs, NCs){
  NCmeansRad <- rowMeans(Data[, NCs], na.rm=T)
  methodRad <- newXMLNode("Method", attrs=c(name="raw"))

  normedData <- intraPlateNorm(data_matrix=Data, method="IC", IC=ICs, NC_wells=NCs)

  NCmeansIC <-rowMeans(normedData$normData[, NCs], na.rm=T)
  methodIC <- newXMLNode("Method", attrs=c(name="IC"))
  for(i in 1:length(NCmeansRad)){
    addChildren(methodRad, newXMLNode("Target", NCmeansRad[i], attrs=c(name=targets$targetBarcode[which(names(NCmeansRad[i]) == targets$targetName)])))
    addChildren(methodIC,  newXMLNode("Target", NCmeansIC[i], attrs=c(name=targets$targetBarcode[which(names(NCmeansIC[i]) == targets$targetName)])))
  }
  return( newXMLNode("NCBkgdLevels",
                      .children=c(
                        methodRad,
                        methodIC
                      )))
}

QCFlagSample <- function(raw, normed, ICs, NCs, IPCs){
  # Sample QC criteria
  MIN_FRAC_TARGETS_ABOVE_BKGD <- 0.8 # Minimim fraction (B): # Targets with reads above background / TotalReads
  MIN_FRAC_COGNATE_RATIO <- 0.9      # Minimum fraction (C): Cognate reads / total reads
  MIN_IC_READS_PER_SAMPLE <- 500.0   # Minimum number (I) of IC reads within a sample

  QCFlagSamples <- c()
  return(QCFlagSamples)
}

QCFlagPlate <- function(raw, normed, ICs, NCs, IPCs){
  QCFlagPlateXML <- newXMLNode("PlateQC")

  # Plate QC criteria
  MAX_IC_CV <- 0.5                # (V) CV of number of IC reads across all samples
  MAX_IPC_CV <- 0.5               # (I) CV of total read count for each IPC sample
#  MAX_NUM_NC_READ_FRACTION <- 0.1 # (N) Fraction of: NC_reads / total_cognate_reads
#  MAX_UNPARSEABLE_FRACTION <- 0.5 # (U) Fraction of Non-ParseablerReads / Total Reads
  MAX_MEDIAN_IPC_TARGET_CV <- 0.2 # (P) Median of CVs of all IPC targets (Performed on normalized data) 
  DETECTABILITY_FRAC <- 0.75      # (D) Detectability fraction (target is detected if >50% of samples > LOD)

  # Calculate Plate-wide QC vals
  ## MAX_IC_CV (V)
  ICvals <- raw[ICs,]
  ICvals[is.na(ICvals)] <- 0
  IC_CV <- sd(ICvals) / mean(ICvals)
  addChildren(QCFlagPlateXML, newXMLNode("QCFlag", IC_CV,
                                         attrs=c(
                                                name="V",
                                                set=if(IC_CV > MAX_IC_CV) "T" else "F",
                                                method="raw"
                                                )
              )
  )

  ## MAX_IPC_CV (I)
  IPCvals <- raw[, IPCs]
  IPCvals[is.na(IPCvals)] <- 0
  IPCvals2 <- colMeans(IPCvals)
  IPC_CV <- sd(IPCvals2) / mean(IPCvals2)
  addChildren(QCFlagPlateXML, newXMLNode("QCFlag", IPC_CV,
                                         attrs=c(
                                                name="I",
                                                set=if(IPC_CV > MAX_IPC_CV) "T" else "F",
                                                method="raw"
                                                )
              )
  )

  ## MAX_MEDIAN_IPC_TARGET_CV (P)
  IPCnormvals <- normed[, IPCs]
  IPCnormvals[is.na(IPCvals)] <- 0
  median_IPC_targetCV <- median(apply(IPCnormvals, 1, sd) / rowMeans(IPCnormvals)) 
  addChildren(QCFlagPlateXML, newXMLNode("QCFlag", median_IPC_targetCV,
                                         attrs=c(
                                                name="P",
                                                set=if(median_IPC_targetCV > MAX_MEDIAN_IPC_TARGET_CV) "T" else "F",
                                                method="IC"
                                                )
              )
  )

  return(QCFlagPlateXML)
}

#* Write Processed NULISAseq XML
#*
#* Writes NULISAseq XML file with unnormalized / normalized data and QC flags
#*
#* @param xml_file:[file] Character string. Path and name of the file.
#* @IPC IPC
#* @NC NC
#* @IC IC
#* @barcodeB barcodeB
#* @return 
#*
#* @examples
#* writeXML('filename.xml')
#*
#* @export
#* @post /processXML
processXML <- function(in_xml_file, IPC, NC, IC, barcodeB=NULL, out_XML=NULL){
  c(plateID, ExecutionDetails, RunSummary, targets, samples, Data) %<-%  
    readNULISAseq(in_xml_file,
                  plateID="",
                  file_type='xml_no_mismatches')
#  Data[is.na(Data)] <- 0
#  IPC <- c('IPC') #"InterPlateControl"
#  NC <- c('NC') #"NegativeControl"
#  IC <- c("mCherry")
#  Bridge <- c('Donor') #NULL
  IPCs <- which(grepl(paste(IPC, collapse="|"), colnames(Data)))
  NCs <- which(grepl(paste(NC, collapse="|"), colnames(Data)))
  ICs <- which(grepl(paste(IC, collapse="|"), rownames(Data)))
  bcodeA <- bcodeA_XML(targets)
  bcodeB <- bcodeB_XML(samples, barcodeB)

  #NC Bkgd Levels
  NCBkgdLevels <- NCBkgdLevels_XML(Data, ICs, NCs)
  normedData <- intraPlateNorm(data_matrix=Data, method="IC", IC=ICs)
  base <- base_XML(ExecutionDetails, bcodeA, bcodeB, RunSummary)  
  data <- newXMLNode("Data")
  addChildren(base, addChildren(data, NCBkgdLevels))
  qcPlate <- QCFlagPlate(Data, normedData$normData, IC=ICs, NC=NCs, IPC=IPCs)
  addChildren(data, qcPlate)

  qcSample <- QCFlagSample(Data, normedData$normData, IC=ICs, NC=NCs, IPC=IPCs)
  uniqSampleNames <- unique(samples$sampleName)
  for( i in 1:length(uniqSampleNames)){
    ind <- which(samples$sampleName == uniqSampleNames[i])
    sampleNode <- newXMLNode("Sample", attrs=c(name=uniqSampleNames[i], replicates=length(ind)))
    for (j in 1:length(ind)){
      rep <- newXMLNode("Replicate", attrs=c(
                                             name=samples$sampleBarcode[ind[j]],
                                             signal=samples$matching[ind[j]], 
                                             bkgd=samples$"non-matching"[ind[j]]) 
                        )
      for (k in 1:2){
        name <- if(k == 1) "raw" else "IC"
        method <- newXMLNode("Method", attrs=c(name=name))
        vals <- if (k == 1) Data[, ind[j]] else normedData$normData[, ind[j]]
        lod <- if (k == 1) lod(data_matrix=Data, blanks=NCs, min_count=0) else lod(data_matrix=normedData$normData, blanks=NCs, min_count=0)
        for (m in 1:length(vals)){
          name <- targets$targetBarcode[which(names(vals[m]) == targets$targetName)]
          aboveLODval <- if(!lod$aboveLOD[m, ind[j]] || is.na(lod$aboveLOD[m, ind[j]])) "N" else "Y"
          addChildren(method, newXMLNode("Target", 
                                          attrs=c(
                                            name=name,
                                            aboveBkgd=aboveLODval
                                          ),
                                          vals[m]
                                        )
          )
        }
        addChildren(rep, method)
      }
      addChildren(sampleNode, rep)

      #CalcQC
      methodIC  <- newXMLNode("Method", attrs=c(name="IC")) 
      methodRaw <- newXMLNode("Method", attrs=c(name="raw")) 
      
    }
    addChildren(data, sampleNode)
  }
  if(!is.null(out_XML)){
    write_xml(base, out_XML)
    return(NULL)
  }
  return(base)
}
