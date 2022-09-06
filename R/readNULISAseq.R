#' Read NULISAseq XML
#'
#' Reads NULISAseq XML file, where XML file is output from the Alamar
#' Biosciences Galaxy NULISAseq (Beta) tool or the NULISAseq Normalization (Alpha) tool.
#'
#' @param xml_file Character string. Path and name of the file.
#' @param file_type Character string. Type of input file, as output from Galaxy. Options include
#' xml_full_output, xml_no_mismatches (default) (both from NULISAseq tool),
#' or xml_normalization (from NULISAseq Normalization tool).
#'
#' @return List of lists, data frames, and matrices.
#' Output will differ slightly depending on the input file type.
#'
#' @examples
#' plate1 <- readNULISAseq('filename.xml')
#'
#' @export
#'
readNULISAseq <- function(xml_file, file_type='xml_no_mismatches'){
  # read in xml file
  xml <- xml2::read_xml(xml_file)
  # NOTE: NEED TO ADD file_type CONDITIONAL SOMEWHERE
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
  # parse the Data section
  ###########################
  # background levels
  ###########################
  # get background levels
  NCBkgdLevels <- xml2::xml_find_first(xml, './/NCBkgdLevels')
  NCBkgdMethods <- xml2::xml_find_all(NCBkgdLevels, './/Method')
  # raw background
  NCBkgd_Raw <- NCBkgdMethods[1]
  bkgd_raw <- xml2::xml_double(xml2::xml_find_all(NCBkgd_Raw, './/Target'))
  targetBarcode <- xml2::xml_attr(xml2::xml_find_all(NCBkgd_Raw, './/Target'), 'name')
  NCBkgd_Raw <- data.frame(targetBarcode, bkgd_raw)
  # IC background
  NCBkgd_IC <- NCBkgdMethods[2]
  bkgd_IC <- xml2::xml_double(xml2::xml_find_all(NCBkgd_IC, './/Target'))
  targetBarcode <- xml2::xml_attr(xml2::xml_find_all(NCBkgd_IC, './/Target'), 'name')
  NCBkgd_IC <- data.frame(targetBarcode, bkgd_IC)
  # merge background levels with target data frame
  targets <- merge(targets, NCBkgd_Raw, all=TRUE)
  targets <- merge(targets, NCBkgd_IC, all=TRUE)
  # sort alphabetically by target full name 
  targets <- targets[order(targets$targetName),]
  
  ###########################
  # target counts
  ###########################
  SampleData <- xml2::xml_find_all(xml, './/Sample')
  SampleReplicates <- xml2::xml_find_all(SampleData, './/Replicate')
  # save signal / background
  SignalBkgd <- do.call(rbind, xml2::xml_attrs(SampleReplicates))
  SignalBkgd <- as.data.frame(SignalBkgd)
  # convert "signal", "bkgd" to numeric values
  SignalBkgd[, c("signal", "bkgd")] <- sapply(SignalBkgd[, c("signal", "bkgd")], as.numeric)
  # loop over sample replicates and save data
  Data_IC <- targets[,c('targetBarcode', 'targetName')]
  Data_TC <- targets[,c('targetBarcode', 'targetName')]
  Data_raw <- targets[,c('targetBarcode', 'targetName')]
  aboveBkgd_IC <- targets[,c('targetBarcode', 'targetName')]
  aboveBkgd_TC <- targets[,c('targetBarcode', 'targetName')]
  aboveBkgd_raw <- targets[,c('targetBarcode', 'targetName')]
  sampleBarcode <- c()
  # sample QC values (threshold):
  # B (0.8): Minimum fraction of targets with reads above background
  # C (0.9): Minimum fraction cognate reads / total reads
  # I (500.0): Minimum number of IC reads in a sample 
  SampleQC <- NULL
  for(i in 1:length(SampleReplicates)){
    sampleBarcode[i] <- xml2::xml_attr(SampleReplicates[[i]], 'name')
    QCFlags <- xml2::xml_find_all(SampleReplicates[[i]], './/QCFlag')
    QCFlag_names <- unlist(xml2::xml_attrs(QCFlags))
    QCFlag_values <- xml2::xml_double(QCFlags)
    SampleQC_i <- data.frame(sampleBarcode=xml2::xml_attr(SampleReplicates[[i]], 'name'),
                             B=NA, C=NA, I=NA)
    for (flag in c('B', 'C', 'I')){
      if (sum(QCFlag_names==flag) > 0){
        SampleQC_i[,flag] <- QCFlag_values[QCFlag_names==flag]
      }
    }
    SampleQC <- rbind(SampleQC, SampleQC_i)
    methods <- xml2::xml_find_all(SampleReplicates[[i]], './/Method')
    method_names <- xml2::xml_attr(methods, 'name')
    # IC data
    targetBarcode <- xml2::xml_attr(xml2::xml_children(methods[method_names=='IC']), 'name')
    IC_Data <- cbind(targetBarcode, 
                     xml2::xml_double(xml2::xml_children(methods[method_names=='IC'])))
    colnames(IC_Data)[2] <- sampleBarcode[i]
    IC_aboveBkgd <- cbind(targetBarcode, 
                          xml2::xml_attr(xml2::xml_children(methods[method_names=='IC']), 'aboveBkgd'))
    colnames(IC_aboveBkgd)[2] <- sampleBarcode[i]
    # merge with rest of data
    Data_IC <- merge(Data_IC, IC_Data, 
                     all.x=TRUE, by='targetBarcode')
    aboveBkgd_IC <- merge(aboveBkgd_IC, IC_aboveBkgd,
                          all.x=TRUE, by='targetBarcode')
    # TC data
    targetBarcode <- xml2::xml_attr(xml2::xml_children(methods[method_names=='TC']), 'name')
    TC_Data <- cbind(targetBarcode, 
                     xml2::xml_double(xml2::xml_children(methods[method_names=='TC'])))
    colnames(TC_Data)[2] <- sampleBarcode[i]
    TC_aboveBkgd <- cbind(targetBarcode, 
                          xml2::xml_attr(xml2::xml_children(methods[method_names=='TC']), 'aboveBkgd'))
    colnames(TC_aboveBkgd)[2] <- sampleBarcode[i]
    # merge with rest of data
    Data_TC <- merge(Data_TC, TC_Data, 
                     all.x=TRUE, by='targetBarcode')
    aboveBkgd_TC <- merge(aboveBkgd_TC, TC_aboveBkgd,
                          all.x=TRUE, by='targetBarcode')
    # raw data
    targetBarcode <- xml2::xml_attr(xml2::xml_children(methods[method_names=='raw']), 'name')
    raw_Data <- cbind(targetBarcode, 
                      xml2::xml_double(xml2::xml_children(methods[method_names=='raw'])))
    colnames(raw_Data)[2] <- sampleBarcode[i]
    raw_aboveBkgd <- cbind(targetBarcode, 
                           xml2::xml_attr(xml2::xml_children(methods[method_names=='raw']), 'aboveBkgd'))
    colnames(raw_aboveBkgd)[2] <- sampleBarcode[i]
    # merge with rest of data
    Data_raw <- merge(Data_raw, raw_Data, 
                      all.x=TRUE, by='targetBarcode')
    aboveBkgd_raw <- merge(aboveBkgd_raw, raw_aboveBkgd,
                           all.x=TRUE, by='targetBarcode')
  }
  # add signal / background and QC to samples data frame
  sample_covariates <- samples[,c('sampleBarcode', covariate_names)]
  samples <- merge(samples[,1:7], SignalBkgd, 
                   by.x='sampleBarcode', by.y='name', all=TRUE)
  colnames(SampleQC) <- c('sampleBarcode',
                          'QCFlag_B', 'QCFlag_C', 'QCFlag_I')
  samples <- merge(samples, SampleQC, by='sampleBarcode')
  samples <- merge(samples, sample_covariates, by='sampleBarcode')
  # sort samples by sample name and then sample type
  samples <- samples[order(samples$sampleName),]
  sampleTypeOrder <- c('Sample', 'InterPlateControl', 'NegativeControl')
  samples <- samples[order(match(samples$sampleType, sampleTypeOrder)),]
  # format aboveBkgd as 0/1 logical
  aboveBkgd_IC[,3:ncol(aboveBkgd_IC)] <- aboveBkgd_IC[,3:ncol(aboveBkgd_IC)] == 'Y'
  aboveBkgd_TC[,3:ncol(aboveBkgd_TC)] <- aboveBkgd_TC[,3:ncol(aboveBkgd_TC)] == 'Y'
  aboveBkgd_raw[,3:ncol(aboveBkgd_raw)] <- aboveBkgd_raw[,3:ncol(aboveBkgd_raw)] == 'Y'
  ###########################
  # combined replicate data
  # one unique column per subject
  ###########################
  SampleCombined <- xml2::xml_find_all(SampleData, './/Combined')
  # loop over samples and save data
  DataCombined_IC <- targets[,c('targetBarcode', 'targetName')]
  DataCombined_TC <- targets[,c('targetBarcode', 'targetName')]
  DataCombined_raw <- targets[,c('targetBarcode', 'targetName')]
  aboveBkgdCombined_IC <- targets[,c('targetBarcode', 'targetName')]
  aboveBkgdCombined_TC <- targets[,c('targetBarcode', 'targetName')]
  aboveBkgdCombined_raw <- targets[,c('targetBarcode', 'targetName')]
  sampleName <- c()
  replicates <- c()
  SampleCombinedQC <- NULL
  for(i in 1:length(SampleCombined)){
    sampleName[i] <- xml2::xml_attr(SampleCombined[[i]], 'name')
    replicates[i] <- xml2::xml_attr(SampleCombined[[i]], 'replicates')
    QCFlags <- xml2::xml_find_all(SampleCombined[[i]], './/QCFlag')
    QCFlag_names <- unlist(xml2::xml_attrs(QCFlags))
    QCFlag_values <- xml2::xml_double(QCFlags)
    SampleQC_i <- data.frame(sampleName=xml2::xml_attr(SampleCombined[[i]], 'name'),
                             B=0, C=0, I=0)
    # save number of replicates with the given flag
    for (flag in c('B', 'C', 'I')){
      SampleQC_i[,flag] <- sum(QCFlag_names==flag)
    }
    SampleCombinedQC <- rbind(SampleCombinedQC, SampleQC_i)
    methods <- xml2::xml_find_all(SampleCombined[[i]], './/Method')
    method_names <- xml2::xml_attr(methods, 'name')
    # IC data
    targetBarcode <- xml2::xml_attr(xml2::xml_children(methods[method_names=='IC']), 'name')
    IC_Data <- cbind(targetBarcode, 
                     xml2::xml_double(xml2::xml_children(methods[method_names=='IC'])))
    colnames(IC_Data)[2] <- sampleName[i]
    IC_aboveBkgd <- cbind(targetBarcode, 
                          xml2::xml_attr(xml2::xml_children(methods[method_names=='IC']), 'aboveBkgd'))
    colnames(IC_aboveBkgd)[2] <- sampleName[i]
    # merge with rest of data
    DataCombined_IC <- merge(DataCombined_IC, IC_Data, 
                             all.x=TRUE, by='targetBarcode')
    aboveBkgdCombined_IC <- merge(aboveBkgdCombined_IC, IC_aboveBkgd,
                                  all.x=TRUE, by='targetBarcode')
    # TC data
    targetBarcode <- xml2::xml_attr(xml2::xml_children(methods[method_names=='TC']), 'name')
    TC_Data <- cbind(targetBarcode, 
                     xml2::xml_double(xml2::xml_children(methods[method_names=='TC'])))
    colnames(TC_Data)[2] <- sampleName[i]
    TC_aboveBkgd <- cbind(targetBarcode, 
                          xml2::xml_attr(xml2::xml_children(methods[method_names=='TC']), 'aboveBkgd'))
    colnames(TC_aboveBkgd)[2] <- sampleName[i]
    # merge with rest of data
    DataCombined_TC <- merge(DataCombined_TC, TC_Data, 
                             all.x=TRUE, by='targetBarcode')
    aboveBkgdCombined_TC <- merge(aboveBkgdCombined_TC, TC_aboveBkgd,
                                  all.x=TRUE, by='targetBarcode')
    # raw data
    targetBarcode <- xml2::xml_attr(xml2::xml_children(methods[method_names=='raw']), 'name')
    raw_Data <- cbind(targetBarcode, 
                      xml2::xml_double(xml2::xml_children(methods[method_names=='raw'])))
    colnames(raw_Data)[2] <- sampleName[i]
    raw_aboveBkgd <- cbind(targetBarcode, 
                           xml2::xml_attr(xml2::xml_children(methods[method_names=='raw']), 'aboveBkgd'))
    colnames(raw_aboveBkgd)[2] <- sampleName[i]
    # merge with rest of data
    DataCombined_raw <- merge(DataCombined_raw, raw_Data, 
                              all.x=TRUE, by='targetBarcode')
    aboveBkgdCombined_raw <- merge(aboveBkgdCombined_raw, raw_aboveBkgd,
                                   all.x=TRUE, by='targetBarcode')
  }
  # create a samples combined data.frame
  samples_combined <- cbind(sampleName, replicates)
  colnames(SampleCombinedQC) <- c('sampleName',
                                  'QCFlag_B', 'QCFlag_C', 'QCFlag_I')
  samples_combined <- merge(samples_combined, SampleCombinedQC, by='sampleName')
  samples_combined <- merge(samples_combined, unique(samples[,c('plateID', 
                                                                'sampleName',
                                                                'sampleType',
                                                                covariate_names)]),
                            by='sampleName')
  samples_combined$sampleID <- paste(samples_combined$plateID[1], samples_combined$sampleName, sep='_')
  samples_combined <- samples_combined[,c('sampleID', 'plateID','sampleName','sampleType',
                                          'replicates','QCFlag_B','QCFlag_C','QCFlag_I',
                                          covariate_names)]
  # sort samples_combined by sample name and then sample type
  samples_combined <- samples_combined[order(samples_combined$sampleName),]
  samples_combined <- samples_combined[order(match(samples_combined$sampleType, sampleTypeOrder)),]
  # format aboveBkgd as 0/1 logical
  aboveBkgdCombined_IC[,3:ncol(aboveBkgdCombined_IC)] <- aboveBkgdCombined_IC[,3:ncol(aboveBkgdCombined_IC)] == 'Y'
  aboveBkgdCombined_TC[,3:ncol(aboveBkgdCombined_TC)] <- aboveBkgdCombined_TC[,3:ncol(aboveBkgdCombined_TC)] == 'Y'
  aboveBkgdCombined_raw[,3:ncol(aboveBkgdCombined_raw)] <- aboveBkgdCombined_raw[,3:ncol(aboveBkgdCombined_raw)] == 'Y'
  ###########################
  # convert numeric data frames to matrices
  # and add row / column names
  ########################### 
  # all sample replicate data
  convert_to_matrix <- function(dataFrame){
    # sort rows to match target data frame order
    dataFrame <- dataFrame[order(dataFrame$targetName),]
    # sort columns to match sample data frame order
    sampleOrder <- match(samples$sampleBarcode, colnames(dataFrame)[3:ncol(dataFrame)])
    dataFrame <- dataFrame[,c(1,2,sampleOrder+2)]
    dataMatrix <- as.matrix(dataFrame[,3:ncol(dataFrame)])
    dataMatrix <- apply(dataMatrix, 2, as.numeric)
    rownames(dataMatrix) <- dataFrame$targetBarcode
    sampleID_order <- match(colnames(dataFrame)[3:ncol(dataFrame)], samples$sampleBarcode)
    colnames(dataMatrix) <- samples$sampleID[sampleID_order]
    return(dataMatrix)
  }
  Data_IC_matrix <- convert_to_matrix(Data_IC)
  Data_TC_matrix <- convert_to_matrix(Data_TC)
  Data_raw_matrix <- convert_to_matrix(Data_raw)
  aboveBkgd_IC_matrix <- convert_to_matrix(aboveBkgd_IC)
  aboveBkgd_TC_matrix <- convert_to_matrix(aboveBkgd_TC)
  aboveBkgd_raw_matrix <- convert_to_matrix(aboveBkgd_raw)
  # combined sample replicate data
  convert_to_matrix_combined <- function(dataFrame){
    # sort rows to match target data frame order
    dataFrame <- dataFrame[order(dataFrame$ ),]
    # sort columns to match sample data frame order
    sampleOrder <- match(samples_combined$sampleName, colnames(dataFrame)[3:ncol(dataFrame)])
    dataFrame <- dataFrame[,c(1,2,sampleOrder+2)]
    dataMatrix <- as.matrix(dataFrame[,3:ncol(dataFrame)])
    dataMatrix <- apply(dataMatrix, 2, as.numeric)
    rownames(dataMatrix) <- dataFrame$targetBarcode
    colnames(dataMatrix) <- paste(plateID[1], colnames(dataFrame)[3:ncol(dataFrame)], sep='_')
    return(dataMatrix)
  }
  DataCombined_IC_matrix <- convert_to_matrix_combined(DataCombined_IC)
  DataCombined_TC_matrix <- convert_to_matrix_combined(DataCombined_TC)
  DataCombined_raw_matrix <- convert_to_matrix_combined(DataCombined_raw)
  aboveBkgdCombined_IC_matrix <- convert_to_matrix_combined(aboveBkgdCombined_IC)
  aboveBkgdCombined_TC_matrix <- convert_to_matrix_combined(aboveBkgdCombined_TC)
  aboveBkgdCombined_raw_matrix <- convert_to_matrix_combined(aboveBkgdCombined_raw)
  
  # log transform original (raw) counts
  Data_rawlog2_matrix <- log(Data_raw_matrix+1, base=2)
  DataCombined_rawlog2_matrix <- log(DataCombined_raw_matrix+1, base=2)
  
  ###########################
  # return the output
  ###########################
  return(list(
    plateID=plateID[1],
    ExecutionDetails=ExecutionDetails,
    RunSummary=RunSummary,
    targets=targets,
    samples=samples,
    samples_combined=samples_combined,
    Data_IC=Data_IC_matrix,
    Data_TC=Data_TC_matrix,
    Data_raw=Data_raw_matrix,
    Data_rawlog2=Data_rawlog2_matrix,
    aboveBkgd_IC=aboveBkgd_IC_matrix,
    aboveBkgd_TC=aboveBkgd_TC_matrix,
    aboveBkgd_raw=aboveBkgd_raw_matrix,
    DataCombined_IC=DataCombined_IC_matrix,
    DataCombined_TC=DataCombined_TC_matrix,
    DataCombined_raw=DataCombined_raw_matrix,
    DataCombined_rawlog2=DataCombined_rawlog2_matrix,
    aboveBkgdCombined_IC=aboveBkgdCombined_IC_matrix,
    aboveBkgdCombined_TC=aboveBkgdCombined_TC_matrix,
    aboveBkgdCombined_raw=aboveBkgdCombined_raw_matrix
  ))

}
