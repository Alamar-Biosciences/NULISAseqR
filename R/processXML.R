base_XML <- function(ExecutionDetails, bcodeA, bcodeB, RunSummary){
  children <- c()
  for(i in 1:length(ExecutionDetails)){
    val <- if(length(ExecutionDetails[[i]]) == 1) ExecutionDetails[i] else ExecutionDetails[[i]][[1]]
    children <- c(children, newXMLNode(names(ExecutionDetails)[i], val))
  }
  norm_ver <- "?"
  norm_ver_file <- ".git/refs/heads/main"
  if(file.exists(norm_ver_file)){
    norm_ver <- paste(readLines(norm_ver_file), collapse="\n")
  }else{
    if (file.exists("main")){
    norm_ver <- paste(readLines("main"), collapse="\n")
    }
  }

  children <- c(children,
                newXMLNode("Normalize_CommitID", norm_ver),
                newXMLNode("Normalize_Date", toString(Sys.time()))
                )
  return( newXMLNode("NULISAseq", 
                      .children=c(
                        newXMLNode("ExecutionDetails", 
                          .children=children
                        ),
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
  if(!is.null(barcodeB) && barcodeB != ""){
    if(file.exists(barcodeB)){
      barcodeB_file <- read.table(barcodeB, sep='\t', comment.char='&', header=T, na.strings=c())
    }else{
      barcodeB_file <- read.table(text=barcodeB, sep='\t', comment.char='&', header=T, na.strings=c())
    }
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
      names(attrVals) <- attrNames
      addChildren(bcodeB, newXMLNode("Barcode", samples$sampleName[i], attrs=attrVals)) 
    }else{
      info <- barcodeB_file[which(samples$sampleBarcode[i] == barcodeB_file[,1]),]

      attrVals <- c(info[2])
      attrNames <- c("name")
      if (length(info) > 2){
        for( k in 3:length(info)){
          val <- colnames(barcodeB_file)[k]
          attrVals <- c(attrVals, info[k])
          attrNames <- c(attrNames, names(info)[k])
        }
      }
      names(attrVals) <- attrNames
      addChildren(bcodeB, newXMLNode("Barcode", rownames(info), attrs=attrVals)) 
    }
  }
  return(bcodeB)
}

NCBkgdLevels_XML <- function(Data, normedData, targets, NCs){
  NCmeansRad <- rowMeans(Data[, NCs], na.rm=T)
  methodRad <- newXMLNode("Method", attrs=c(name="raw"))


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

QCThresholds_XML <- function(){
  QCthresh <- newXMLNode("QCThresholds")
  QCsample <- QCSampleCriteria()
  QCplate <- QCPlateCriteria()
  sample <- newXMLNode("Sample")
  plate <- newXMLNode("Plate")
  detectability_ind <- 0
  detectabilities <- 1
  for (i in 1:length(QCsample$thresholds)){
    if (startsWith(names(QCsample$thresholds[i]), "Detectability.") && detectability_ind == 0){
      detectability_ind <- i
      j = i
    } else if(startsWith(names(QCsample$thresholds[i]), "Detectability.")){
      j = detectability_ind
      detectabilities <- detectabilities + 1
    } else if (detectability_ind == 0){
      j = 1
    }else{
      j = i - detectabilities 
    }
    addChildren(sample, newXMLNode("Threshold", attrs=c(name=names(QCsample$thresholds[i]), value=QCsample$thresholds[[i]], operator=QCsample$operators[[j]])))
  }
  for (i in 1:length(QCplate$thresholds)){
    addChildren(plate, newXMLNode("Threshold", attrs=c(name=names(QCplate$thresholds[i]), value=QCplate$thresholds[[i]], operator=QCplate$operators[[j]])))
  }
  addChildren(QCthresh, sample)
  addChildren(QCthresh, plate)
}

#' Write Processed NULISAseq XML
#'
#' Writes NULISAseq XML file with unnormalized / normalized data and QC flags
#'
#' @param xml_file:[file] Character string. Path and name of the file.
#' @param IPC Name to search for Interprocess control (IPC) samples
#' @param NC Name to search for Negative control (NC) samples
#' @param IC Name to search for Internal Control (IC) targets
#' @param barcodeB optional BarcodeB file
#' @param out_XML optional output XML to write to
#'
#' @examples
#' # writeXML('filename.xml')
#'
#' @export
processXML <- function(in_xml, IPCs=c("InterProcessControl"), NCs=c("NegativeControl"), ICs=c("mCherry"), barcodeB="", out_XML=""){
  val <- readNULISAseq(in_xml, IPC=IPCs, NC=NCs, IC=ICs, file_type='xml_no_mismatches')
  Data <- val$Data
  RunSummary <- val$RunSummary
  targets <- val$targets
  samples <- val$samples
  ExecutionDetails <- val$ExecutionDetails
  plateID <- val$plateID
  bcodeA <- bcodeA_XML(targets)
  bcodeB <- bcodeB_XML(samples, barcodeB)

  #NC Bkgd Levels
  normedData <- intraPlateNorm(data_matrix=Data, method="IC", IC=val$IC)
  NCBkgdLevels <- NCBkgdLevels_XML(Data, normedData, targets, val$NC)
  base <- base_XML(ExecutionDetails, bcodeA, bcodeB, RunSummary)  
  data <- newXMLNode("Data")
  addChildren(base, QCThresholds_XML())

  addChildren(base, addChildren(data, NCBkgdLevels))
  qcPlate <- QCFlagPlate(Data, normedData$normData, targets, samples)
  qcSample <- QCFlagSample(Data, normedData$normData, samples, targets)
  plateNode <- newXMLNode("PlateQC")
  QC2XML(qcPlate, plateNode, sample=F)
  addChildren(base, plateNode)
  uniqSampleNames <- unique(samples$sampleName)
  for( i in 1:length(uniqSampleNames)){
    ind <- which(samples$sampleName == uniqSampleNames[i])
    sampleNode <- newXMLNode("Sample", attrs=c(name=uniqSampleNames[i], replicates=length(ind)))
#    combinedNode <- newXMLNode("Combined", attrs=c(name=uniqSampleNames[i], replicates=length(ind)))

    for (j in 1:length(ind)){
      rep <- newXMLNode("Replicate", attrs=c(
                                             name=samples$sampleBarcode[ind[j]],
                                             signal=samples$matching[ind[j]], 
                                             bkgd=samples$"non-matching"[ind[j]]) 
                        )
      indQC <- which(qcSample$sampleBarcode == samples$sampleBarcode[ind[j]])
      
      for (q in 1:length(indQC)){
        QC2XML(qcSample[indQC[q],], rep, sample=T)
#        QC2XML(qcSample[indQC[q],], combinedNode, sample=T, combined=T)
      }
      for (k in 2:2){
        name <- if(k == 1) "raw" else "IC"
        method <- newXMLNode("Method", attrs=c(name=name))
        vals <- if (k == 1) Data[, ind[j]] else normedData$normData[, ind[j]]
        lod <- if (k == 1) lod(data_matrix=Data, blanks=val$NC, min_count=0) else lod(data_matrix=normedData$normData, blanks=val$NC, min_count=0)
        for (m in 1:length(vals)){
          name <- targets$targetBarcode[which(names(vals[m]) == targets$targetName)]
          aboveLODval <- if(!lod$aboveLOD[m, ind[j]] || is.na(lod$aboveLOD[m, ind[j]])) "N" else "Y"
          addChildren(method, newXMLNode("Target", 
                                          attrs=c(
                                            name=name,
                                            aboveBkgd=aboveLODval
                                          ),
                                          if(is.na(vals[m])) 0 else vals[m]
                                        )
          )
        }
        addChildren(rep, method)
      }
      addChildren(sampleNode, rep)
    }
    addChildren(data, sampleNode)

    ## Handle combining replicates here
#    dataCombined <- if(length(ind) > 1) cbind(rowMeans(Data[, ind],na.rm=T), Data[, val$NC]) else cbind(Data[, ind], Data[, val$NC])
#    normCombined <- if(length(ind) > 1) cbind(rowMeans(normedData$normData[, ind]), normedData$normData[, val$NC]) else cbind(normedData$normData[, ind], normedData$normData[, val$NC])
#    lodData <- lod(data_matrix=dataCombined, blanks=val$NC, min_count=0)
#    lodNorm <- lod(data_matrix=normCombined, blanks=val$NC, min_count=0)
#    combRaw <- newXMLNode("Method", attrs=c(name="raw"))
#    combIC  <- newXMLNode("Method", attrs=c(name="IC"))
#    addChildren(combinedNode, combRaw)
#    addChildren(combinedNode, combIC) 
#    for (w in 1:nrow(dataCombined)){
#      valLoD <- if (lodData$aboveLOD[w] && !is.na(lodData$aboveLOD[w])) "Y" else "N"
#      barcode <- targets$targetBarcode[which(targets$targetName == rownames(dataCombined)[w])]
#      addChildren(combRaw, newXMLNode("Target",
#                                      attrs=c(
#                                              name=barcode, #rownames(dataCombined)[w],
#                                              aboveBkgd=valLoD
#                                              ),
#                                      if(is.na(dataCombined[w, 1])) 0 else dataCombined[w, 1]
#                                      )
#      )
#      valLoD <- if (lodNorm$aboveLOD[w] && !is.na(lodNorm$aboveLOD[w])) "Y" else "N"
#      barcode <- targets$targetBarcode[which(targets$targetName == rownames(normCombined)[w])]
#      addChildren(combIC, newXMLNode("Target",
#                                      attrs=c(
#                                              name=barcode,#rownames(normCombined)[w],
#                                              aboveBkgd=valLoD
#                                              ),
#                                      if(is.na(normCombined[w, 1])) 0 else normCombined[w, 1]
#                                      )
#      )
#    }
#    addChildren(sampleNode, combinedNode)
  }
  if(!is.null(out_XML) && out_XML != ""){
    cat(saveXML(base, indent=TRUE, prefix='<?xml version="1.0" encoding="UTF-8"?>\n'), file=out_XML)
    return(base)
  }
  return(saveXML(base, indent=F, prefix='<?xml version="1.0" encoding="UTF-8"?>\n'))
}
