base_XML <- function(ExecutionDetails, bcodeA, bcodeB, RunSummary){
  return( newXMLNode("NULISAseq", 
                      .children=c(
                        newXMLNode("ExecutionDetails", 
                          .children=c(
                            newXMLNode("CommandLine", ExecutionDetails$CommandLine),
                            newXMLNode("ExecutionTime", ExecutionDetails$ExecutionTime[[1]],
                              attrs=c(units=ExecutionDetails$ExecutionTime$unit)),
                            newXMLNode("CommitID", ExecutionDetails$CommitID),
                            newXMLNode("Date", ExecutionDetails$Date),
                            newXMLNode("Normalize_CommitID"),
                            newXMLNode("NormalizeDate")
                          )
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
bcodeB_XML <- function(samples){
  # Process BarcodeBs
  inds <- c()
  sampColI <- names(samples)
  for( i in 1:length(sampColI)){
    if (sampColI[i] != "sampleBarcode" && sampColI[i] != "sampleName" && 
      sampColI[i] !="matching" && sampColI[i] != "non-matching"){
      inds <- c(inds, i)
    }
  }
  bcodeB <- newXMLNode("BarcodeB")
  for(i in 1:length(samples$sampleBarcode)){
    attrVals <- attrNames <- c()
    for (j in 1:length(inds)){
      attrNames <- c(attrNames, sampColI[inds[j]])
      attrVals <- c(attrVals, eval(parse(text=paste0("(samples$", sampColI[inds[j]],")[i]"))))
    }
    attrVals <- c(attrVals, samples$sampleBarcode[i])
    attrNames <- c(attrNames, "name")
    names(attrVals) <- attrNames
    addChildren(bcodeB, newXMLNode("Barcode", samples$sampleName[i], attrs=attrVals)) 
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
#' Write Processed NULISAseq XML
#'
#' Writes NULISAseq XML file with unnormalized / normalized data and QC flags
#'
#' @param xml_file Character string. Path and name of the file.
#' @param plateID Character string that will be added to the beginning of
#' column names before the sample name. This is helpful for 
#' identifying the plate each sample came from 
#' after interplate normalization. If no plate ID is given, the function
#' will use the date and time in the execution details (this is 
#' very long so it is recommended to provide a plate ID!).
#' @param file_type Character string. Type of input file, as output from Galaxy. Options include
#' xml_full_output, xml_no_mismatches (default) (both from NULISAseq tool),
#' or xml_normalization (from NULISAseq Normalization tool).
#'
#' @return NULL
#'
#' @examples
#' writeXML('filename.xml')
#'
#' @export
#'
writeProcessedXML <- function(in_xml_file, out_xml_file, IPCs, NCs, ICs){
  c(plateID, ExecutionDetails, RunSummary, targets, samples, Data) %<-%  
    readNULISAseq(in_xml_file,
                  plateID="",
                  file_type='xml_no_mismatches')
#  Data[is.na(Data)] <- 0
  IPC <- c('IPC') #"InterPlateControl"
  NC <- c('NC') #"NegativeControl"
  IC <- c("mCherry")
  Bridge <- c('Donor') #NULL
 IPCs <- which(grepl(paste(IPC, collapse="|"), colnames(Data)))
 NCs <- which(grepl(paste(NC, collapse="|"), colnames(Data)))
 ICs <- which(grepl(paste(IC, collapse="|"), rownames(Data)))
  bcodeA <- bcodeA_XML(targets)
  bcodeB <- bcodeB_XML(samples)

  #NC Bkgd Levels
  NCBkgdLevels <- NCBkgdLevels_XML(Data, ICs, NCs)
  base <- base_XML(ExecutionDetails, bcodeA, bcodeB, RunSummary)  
  data <- newXMLNode("Data")
  addChildren(base, addChildren(data, NCBkgdLevels))

  normedData <- intraPlateNorm(data_matrix=Data, method="IC", IC=ICs)
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
#      for (k in 1:length()){

#      }
      addChildren(sampleNode, rep)

      #CalcQC
      methodIC  <- newXMLNode("Method", attrs=c(name="IC")) 
      methodRaw <- newXMLNode("Method", attrs=c(name="raw")) 
      
    }
    addChildren(data, sampleNode)
  }
  print(base)
}
