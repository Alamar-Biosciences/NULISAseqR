#' Read NULISAseq XML
#'
#' Reads NULISAseq XML file, where XML file is output from the Alamar
#' Biosciences Galaxy NULISAseq (Beta) tool or the NULISAseq Normalization (Alpha) tool.
#'
#' @param file Character string. Path and name of the file.
#' @param plateID Character string that denotes plate ID.
#' @param file_type Character string. Type of input file.
#' Options include \code{xml_no_mismatches} (default, output from Galaxy NULISAseq tool),
#' \code{xml_full_output} (output from Galaxy NULISAseq tool, not currently implemented),
#' or \code{csv_long} (output from \code{writeNULISAseq()} function).
#' @param target_column_names A vector of column names for target-specific data
#' that will appear in the target data frame output. 
#' Only used for \code{csv_long} format, and only needed when target column names
#' are different from the default, which are \code{c('Target', 'AlamarTargetID', 
#' 'UniProtID', 'ProteinName')}.
#' @param sample_column_names A vector of column names for sample-specific data
#' that will appear in the sample data frame output. 
#' Only used for \code{csv_long} format, and only needed when sample column names
#' are different from the default. Default includes all column names in the csv except 
#' for the target column names, given above (either the default or whatever is specified), 
#' and the following: \code{c('Panel', 'PanelLotNumber', 'LOD', 'log2NormalizedCount')}. 
#' The sample-specific data
#' will include by default \code{c('PlateID', 'SampleName', 'SampleType', 'SampleQC')}.
#' In addition, by default, it will include any other variables in the dataset.
#' If \code{sample_column_names} is specified, however, it will only include
#' those variables in the sample data frame output. 
#' @param sample_group_covar string that represents the name of the covariate
#' that defines sample group information. Defaults to NULL. Typically, this covariate
#' defines whether each sample is plasma, serum, csf, etc.
#' This information can be used when reporting detectability (i.e., Report sample
#' type group specific detectability values, etc.).
#' @param IC string(s) that represents the names of internal control targets. 
#' Default is \code{'mCherry'}.Only used for xml file formats.
#' @param IPC string(s) that represent the inter-plate control wells. 
#' For example, \code{'IPC'} (default). Set to \code{NULL} if there are no IPCs.
#' Only used for xml file formats.
#' @param NC string(s) that represent the negative control wells. 
#' For example, \code{'NC'} (default). Set to \code{NULL} if there are no NCs.
#' Only used for xml file formats.
#' @param Bridge string(s) that represent the bridge sample wells. 
#' Set to \code{NULL} if there are no bridge samples (default).
#' Only used for xml file formats.
#' @param replaceNA Logical. If TRUE (default), will replace missing counts with 
#' zeros.
#'
#' @return List of lists, data frames, and matrices.
#' Output will differ slightly depending on the input file type.
#'
#'
#' @export
#'
readNULISAseq <- function(file, 
                          plateID=NULL, 
                          file_type='xml_no_mismatches',
                          target_column_names=NULL,
                          sample_column_names=NULL,
                          sample_group_covar=NULL,
                          IC='mCherry', 
                          IPC='IPC', SC='SC', NC='NC', Bridge=NULL,
                          replaceNA=TRUE){
  
  if(file_type == 'xml_no_mismatches'){

    # read in xml file
    xml <- xml2::read_xml(file)
    
    ###########################
    # save Execution Details
    ###########################
    ExecutionDetails <- xml2::xml_find_first(xml, './/ExecutionDetails')
    ExecutionDetails <- xml2::as_list(ExecutionDetails)
    # if we need to save the Execution Time units
    ExecutionTimeUnits <- attributes(ExecutionDetails$ExecutionTime)
    ExecutionDetails <- lapply(ExecutionDetails, unlist)
    ExecutionDetails$ExecutionTime <- c(as.numeric(ExecutionDetails$ExecutionTime[1]),
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
      targetMetadata <- data.frame(matrix(nrow=length(targetBarcode),
                                          ncol=length(barcodeA_attrs)))
      colnames(targetMetadata) <- barcodeA_attrs
      for (i in 1:length(barcodeA_attrs)){
        targetMetadata[,i] <- unlist(lapply(RunSummary$Barcodes$BarcodeA, function(x) attributes(x)[barcodeA_attrs[i]]))
      }
      targets <- data.frame(targetBarcode=targetBarcode,
                            targetName=targetName,
                            targetMetadata)
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
      sampleMetadata <- data.frame(matrix(nrow=length(sampleBarcode),
                                          ncol=length(barcodeB_attrs)))
      colnames(sampleMetadata) <- barcodeB_attrs
      for (i in 1:length(barcodeB_attrs)){
        sampleMetadata[,i] <- unlist(lapply(RunSummary$Barcodes$BarcodeB, function(x) attributes(x)[barcodeB_attrs[i]]))
      }
      samples <- data.frame(sampleBarcode=sampleBarcode, 
                            sampleName=sampleName, 
                            sampleMetadata)
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
    RunSummary <- lapply(RunSummary[c('TotalReads',
                                      'Parseable',
                                      'ParseableMatch',
                                      'Unparseable',
                                      'Balancers')], as.numeric)
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
    colnames(DataMatrix) <- samples$sampleName
    # if replaceNA=TRUE, replace NAs with zeros
    if (replaceNA==TRUE){
      DataMatrix[is.na(DataMatrix)] <- 0
    }
    
    val <- if(is.null(NC) && is.null(IPC) && is.null(SC) && is.null(Bridge)) NA else "Sample"
    sampleType <- rep(val, length(samples$sampleName))
    if(!is.null(samples$type)){
      sampleType <- unlist(lapply(samples$type, FUN=function(t) gsub(pattern="sample", replacement="Sample", x=t, fixed=T)))
      samples$type <- NULL
    } 
    
    # add well type information
    if(!is.null(NC)){     sampleType[grep(paste(NC, collapse="|"), samples$sampleName)]  <- "NC" }
    if(!is.null(SC)){     sampleType[grep(paste(SC, collapse="|"), samples$sampleName)] <- "SC" }
    if(!is.null(IPC)){    sampleType[grep(paste(IPC, collapse="|"), samples$sampleName)] <- "IPC" }
    if(!is.null(Bridge)){ sampleType[grep(paste(Bridge, collapse="|"), samples$sampleName)] <- "Bridge"}
    samples$sampleType <- sampleType
    
    # add IC target information
    val <- if(is.null(IC)) NA else "Target"
    targetType <- rep(val, length(targets$targetName))
    if(!is.null(IC)){ targetType[grep(paste(IC, collapse="|"), targets$targetName)] <- "Control"}
    targets$targetType <- targetType
    
    # save the special well type column names
    specialWellsTargets <- list()
    if(!is.null(IPC)) {
      specialWellsTargets[['IPC']] <- samples$sampleName[samples$sampleType=='IPC']
    }
    if(!is.null(NC)) { 
      specialWellsTargets[['NC']] <- samples$sampleName[samples$sampleType=='NC']
    }
    if(!is.null(SC)) {
      specialWellsTargets[['SC']] <- samples$sampleName[samples$sampleType=='SC']
    }
    if(!is.null(Bridge)) {
      specialWellsTargets[['Bridge']] <- samples$sampleName[samples$sampleType=='Bridge']
    }
    specialWellsTargets[['SampleNames']] <- samples$sampleName[samples$sampleType=='Sample']
    
    # add sample identity information
    # users can specify the covariate with sample_group_covar parameter
    if (length(sample_group_covar) == 1 && !is.null(sample_group_covar) && !all(is.na(samples[, sample_group_covar]) | samples[, sample_group_covar] == "NA")) {
      specialWellsTargets[[sample_group_covar]] <- samples[,sample_group_covar][samples$sampleType=='Sample']
    }
    
    # save IC row names
    if(!is.null(IC)) {
      specialWellsTargets[['IC']] <- targets$targetName[targets$targetType=='Control']
    }
    
    # if given, add plateID to samples 
    if( !is.null(plateID) ){
      samples$plateID <- plateID
    }else{
      if ( !is.null(samples$AUTO_PLATE) ){
        plateID <- unique(samples$AUTO_PLATE)[1]
        samples$plateID <- samples$AUTO_PLATE
      }
    }
    
    # Remove empty sample / target covariates
    rEmpty <- function(samples){
      return (if( sum(samples == "NA", na.rm=T) == length(samples) || # all "NA"
                  sum(is.na(samples), na.rm=T) == length(samples) || # all NA (e.g. NULL)
                  is.null(samples)) TRUE else FALSE)  # check if the covariate is NULL
    }
    inds <- which(lapply(samples, rEmpty) == TRUE)
    if(length(inds) > 0){
      samples[, inds] <- NULL
    }
    
    inds <- which(lapply(targets, rEmpty) == TRUE)
    if(length(inds) > 0){
      targets[, inds] <- NULL
    }
    
    ###########################
    # return the output
    ###########################
    return(c(list(
      plateID=plateID,
      ExecutionDetails=ExecutionDetails,
      RunSummary=RunSummary,
      targets=targets,
      samples=samples,
      Data=DataMatrix),
      specialWellsTargets))
    # end file type xml_no_mismatches
  } 
  
  if(file_type == 'csv_long'){
    Data <- read.csv(file)
    # make target and sample data frames
    Data_colnames <- colnames(Data)
    if(is.null(target_column_names)){
      target_column_names <- c('Target', 'AlamarTargetID','UniProtID', 'ProteinName')
    }
    if(is.null(sample_column_names)){
      sample_column_names <- Data_colnames[!(Data_colnames %in% c(target_column_names, 
                                                                  'Panel', 'PanelLotNumber',
                                                                  'LOD', 'log2NormalizedCount'))]
    }
    targets <- unique(Data[,target_column_names])
    samples <- unique(Data[,sample_column_names])
    # make an LOD data frame
    LOD <- unique(Data[,c('PlateID', 'Target', 'LOD')])
    # reformat the log2NormalizedCount data
    # reformat into wide, targets in columns
    Data <- reshape(Data[,c('SampleName', 'Target', 'log2NormalizedCount')],
                    direction= 'wide',
                    idvar = 'Target',
                    timevar='SampleName')
    
    rownames(Data) <- Data$Target
    Data <- Data[,2:ncol(Data)]
    colnames(Data) <- substr(colnames(Data), start=21, stop=nchar(colnames(Data)))
    Data <- as.matrix(Data)
    # set rownames 
    rownames(targets) <- NULL
    rownames(samples) <- NULL
    rownames(LOD) <- NULL
    
    # return output
    return(list(
      targets=targets,
      samples=samples,
      LOD=LOD,
      Data=Data
    ))
  } # end file type csv_long
}



#' Convenience function for use by SAM
#' Read NULISAseq XML, perform normalization, and QC 
#'
#' Reads NULISAseq XML file, where XML file is output from the Alamar
#' Biosciences Galaxy NULISAseq (Beta) tool or the NULISAseq Normalization (Alpha) tool.
#'
#' @param file Character string. Path and name of the file.
#' @param IPC String(s) that represent the inter-plate control wells. 
#' For example, \code{'IPC'} (default). Set to \code{NULL} if there are no IPCs.
#' Only used for xml file formats.
#' @param IC string(s) that represents the names of internal control targets. 
#' Default is \code{'mCherry'}.Only used for xml file formats.
#'
#' @return List of lists, data frames, and matrices.
#' Output will differ slightly depending on the input file type.
#'
#'
#' @export
#'
loadNULISAseq <- function(file, IPC, IC, ...){
  raw <- readNULISAseq(file, IPC=IPC, IC=IC, ...)
  raw$IC_normed <- intraPlateNorm(raw$Data, IC=IC[1], ...)
  raw$normed <- interPlateNorm(list(raw$IC_normed$normData), IPC_wells=list(raw$IPC), ...)
  raw$qcPlate <- QCFlagPlate(raw$Data, raw$IC_normed$normData, raw$targets, raw$samples)
  raw$qcSample <- QCFlagSample(raw$Data, raw$IC_normed$normData, raw$samples, raw$targets)
  return(raw)
}

#' Read Covariate file and add to NULISAseq object
#'
#' Reads NULISAseq covariate file and adds new covariates to list object. Merge based on plateID, row, col and sampleName. 
#' Fails if plateID is not set AND more than one entry in list
#'
#' @param txt_file Character string. File containing covariates (columns with name not already existing
#' @param NULISAseqRuns List of NULISAseq run objects 
#'
#' @return List of lists, data frames, and matrices.
#'
#' @export
#'
readCovariateFile <- function(txt_file, NULISAseqRuns){
  newInfo <- read.table(txt_file, header=T, sep="\t")
  for(i in 1:length(NULISAseqRuns)){
    NULISAseqRuns[[i]]$samples <- merge(runs[[i]]$samples, newInfo, nodups=T, all.x=T,
                                        by.x=c("AUTO_PLATE", "AUTO_WELLROW", "AUTO_WELLCOL", "sampleName", "sampleBarcode"), 
                                        by.y=c("AUTO_PLATE", "AUTO_WELLROW", "AUTO_WELLCOL", "sampleName", "sampleBarcode"))
  }
  return(NULISAseqRuns)
}
