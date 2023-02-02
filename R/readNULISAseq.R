#' Read NULISAseq XML
#'
#' Reads NULISAseq XML file, where XML file is output from the Alamar
#' Biosciences Galaxy NULISAseq (Beta) tool or the NULISAseq Normalization (Alpha) tool.
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
#' @return List of lists, data frames, and matrices.
#' Output will differ slightly depending on the input file type.
#'
#' @examples
#' # plate1 <- readNULISAseq('filename.xml')
#'
#' @export
#'
readNULISAseq <- function(xml_file, 
                          plateID=NULL, 
                          file_type='xml_no_mismatches'){
  # read in xml file
  xml <- xml2::read_xml(xml_file)
  
  ###########################
  # save Execution Details
  ###########################
  ExecutionDetails <- xml2::xml_find_first(xml, './/ExecutionDetails')
  ExecutionDetails <- xml2::as_list(ExecutionDetails)
  # if we need to save the Execution Time units
  ExecutionTimeUnits <- attributes(ExecutionDetails$ExecutionTime)
  ExecutionDetails <- lapply(ExecutionDetails, unlist)
  ExecutionDetails$ExecutionTime <- c(ExecutionDetails$ExecutionTime[1],
                                      ExecutionTimeUnits)
  ###########################
  # save Run Summary 
  ###########################
  RunSummary <- xml2::xml_find_first(xml, './/RunSummary')
  RunSummary <- xml2::as_list(RunSummary)
  # save target data in data frame
  targetBarcode <- unlist(lapply(RunSummary$Barcodes$BarcodeA, function(x) attributes(x)$name))
  targetName <- unlist(RunSummary$Barcodes$BarcodeA)
  # get other target annotations
  barcodeA_attrs <- names(attributes(RunSummary$Barcodes$BarcodeA[[1]]))
  barcodeA_attrs <- barcodeA_attrs[barcodeA_attrs!='name']
  if (length(barcodeA_attrs) > 0){
    target_metadata <- data.frame(matrix(nrow=length(targetBarcode),
                                           ncol=length(barcodeA_attrs)))
    colnames(target_metadata) <- barcodeA_attrs
    for (i in 1:length(barcodeA_attrs)){
      target_metadata[,i] <- unlist(lapply(RunSummary$Barcodes$BarcodeA, function(x) attributes(x)[barcodeA_attrs[i]]))
    }
    targets <- data.frame(targetBarcode=targetBarcode,
                          targetName=targetName,
                          target_metadata)
  } else if (length(barcodeA_attrs) == 0) {
    targets <- data.frame(targetBarcode=targetBarcode,
                          targetName=targetName)
  }
  # sort by target name
  targets <- targets[order(targets$targetName),]
  # save sample data in data frame
  sampleBarcode <- unlist(lapply(RunSummary$Barcodes$BarcodeB, function(x) attributes(x)$name))
  sampleName <- unlist(RunSummary$Barcodes$BarcodeB)
  # get other sample annotations
  barcodeB_attrs <- names(attributes(RunSummary$Barcodes$BarcodeB[[1]]))
  barcodeB_attrs <- barcodeB_attrs[barcodeB_attrs!='name']
  if (length(barcodeB_attrs) > 0){
    sample_metadata <- data.frame(matrix(nrow=length(sampleBarcode),
                                           ncol=length(barcodeB_attrs)))
    colnames(sample_metadata) <- barcodeB_attrs
    for (i in 1:length(barcodeB_attrs)){
      sample_metadata[,i] <- unlist(lapply(RunSummary$Barcodes$BarcodeB, function(x) attributes(x)[barcodeB_attrs[i]]))
    }
    samples <- data.frame(sampleBarcode=sampleBarcode, 
                          sampleName=sampleName, 
                          sample_metadata)
  } else if (length(barcodeB_attrs) == 0) {
    samples <- data.frame(sampleBarcode, sampleName)
  }
  # save remaining relevant RunSummary data
  BalancerNames <- unlist(lapply(RunSummary$Balancers, attr, 'name'))
  RunSummary <- lapply(RunSummary[c('TotalReads', 
                                    'Parseable', 
                                    'ParseableMatch', 
                                    'Unparseable',
                                    'Balancers')], unlist)
  names(RunSummary$Balancers) <- BalancerNames
  
  ###########################
  # save Data section
  ###########################
  SampleData <- xml2::xml_find_all(xml, './/Sample')
  # save matching / non-matching
  matchNonMatch <- do.call(rbind, xml2::xml_attrs(SampleData))
  matchNonMatch <- as.data.frame(matchNonMatch)
  # convert to numeric values
  matchNonMatch[, c("matching", "non-matching")] <- sapply(matchNonMatch[, c("matching", "non-matching")], as.numeric)
  # add matching / non-matching samples data frame
  samples <- merge(samples, matchNonMatch, 
                   by.x='sampleBarcode', by.y='barcode', all=TRUE)
  # sort samples by sample name
  samples <- samples[order(samples$sampleName),]
  # loop over sample replicates and save data
  Data <- targets[,c('targetBarcode', 'targetName')]
  sampleBarcode <- c()
  for(i in 1:length(SampleData)){
    # get data for sample i
    sampleBarcode[i] <- xml2::xml_attr(SampleData[[i]], 'barcode')
    targetBarcode <- xml2::xml_attr(xml2::xml_children(SampleData[[i]]), 'target')
    Data_i <- cbind(targetBarcode, 
                    xml2::xml_double(xml2::xml_children(SampleData[[i]])))
    colnames(Data_i)[2] <- sampleBarcode[i]
    # merge with rest of data
    Data <- merge(Data, Data_i, 
                  all.x=TRUE, by='targetBarcode')
  }
  # sort rows by target name
  Data <- Data[order(Data$targetName),]
  # sort columns to match sample data frame order
  sampleOrder <- match(samples$sampleBarcode, colnames(Data)[3:ncol(Data)])
  Data <- Data[,c(1,2,sampleOrder+2)]
  # convert Data to numeric matrix
  DataMatrix <- as.matrix(Data[,3:ncol(Data)])
  DataMatrix <- apply(DataMatrix, 2, as.numeric)
  # add row names
  rownames(DataMatrix) <- Data$targetName
  # add column names
  # if no plate ID given, use date and time
  if (is.null(plateID)){
    # replace space with underscore
    plateID <- gsub(" ", "_", ExecutionDetails$Date, fixed = TRUE)
    # remove day of week
    plateID <- substr(plateID, 5, nchar(plateID))
  }
  colnames(DataMatrix) <- if(plateID != "") paste(plateID, samples$sampleName, sep='_') else samples$sampleName
  
  ###########################
  # return the output
  ###########################
  return(list(
    plateID=plateID,
    ExecutionDetails=ExecutionDetails,
    RunSummary=RunSummary,
    targets=targets,
    samples=samples,
    Data=DataMatrix
  ))
}
