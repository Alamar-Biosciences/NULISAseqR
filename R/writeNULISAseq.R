#' Write NULISAseq normalized data CSV in long format
#'
#' @description
#' Reads in one or more NULISAseq XML files, where XML file is output from the Alamar
#' Biosciences Galaxy NULISAseq (Beta) tool. Output is a long-format CSV file
#' where each row corresponds to a particular sample-target combination.
#' 
#' @details 
#' This function performs intra-plate 
#' normalization (default is internal control (IC) mCherry), and when multiple 
#' XML files are input, also performs inter-plate normalization (default is inter-plate 
#' control (IPC)). Finally, data is log2-transformed (after adding small constant, 0.01). 
#' 
#' The function also takes as input a target info file containing target metadata, 
#' including Alamar target IDs and protein names, and an optional sample info 
#' file containing sample metadata that
#' is desired to appear in the output file.
#' 
#' Sample QC is currently based on the internal control (mCherry) within 30% of
#' median criterion. If a sample's mCherry count falls outside this threshold, 
#' it will be assigned a "WARN" status. Otherwise the sample will be assigned a 
#' "PASS" status.  
#' 
#' The function outputs a CSV file where each row corresponds to a particular
#' sample-target combination. Total number of rows is number of samples times
#' number of targets. By default the sample control (SC) and inter-plate control 
#' (IPC) wells are included in the output, and NC wells are omitted.
#' 
#' 
#' Columns of output data (in order) include:
#' \itemize{
#'  \item{"PanelLotNumber"}{}
#'  \item{"SampleID"}{}
#'  \item{"(optional sample metadata columns)"}{}
#'  \item{"Target"}{}
#'  \item{"AlamarTargetID"}{}
#'  \item{"UniProtID"}{}
#'  \item{"ProteinName"}{}
#'  \item{"SampleQC"}{}
#'  \item{"LOD"}{}
#'  \item{"log2NormalizedCount"}{}
#' }
#' 
#'
#' @param xml_files Vector of filenames (character strings). Files should 
#' be in order of the desired plateID variable, unless that is otherwise
#' defined.
#' @param dataDir Data directory where xml_files, target_info_file, and sample_info_file 
#' reside (character string). 
#' @param target_info_file Path and filename for the target info CSV file 
#' (character string). Must include columns for TargetName (matches the targets in 
#' XML files), AlamarTargetID, UniProtID, and ProteinName. Only targets in the 
#' target_info_file will be output in the CSV file.
#' @param sample_info_file Optional. Path and filename for the sample annotation CSV file. 
#' Must include columns "plateID" and "sampleName" which correspond to the 
#' plateID (matching the plateIDs input to this function) and sampleName in the readNULISAseq 
#' samples data.frame.
#' @param sample_info_file_variables Subset of column names in 
#' sample_info_file that will be included in data file output. Other columns 
#' will be excluded. Otherwise, if NULL (default), all columns will be included.
#' @param output_filename Filename for output CSV file. 
#' @param Panel Name of multi-plex panel. Default is '200-plex Inflammation v1'
#' @param PanelLotNumber The panel lot number.
#' @param plateIDs A vector of plate IDs. If NULL, default is to number plates from 01 to 
#' total number of plates based on the order of xml_files. Passed to readNULISAseq() function
#' and output as a column in the data file. 
#' @param ICs vector of string(s). Internal control names. Default is "mCherry". 
#' First IC in vector will be used in intra-plate IC-normalization (usually mCherry). 
#' ICs will be omitted from CSV output by default. 
#' If NULL, no intra-plate normalization is done.
#' @param IPC_string IPC_string (or IPC_wells) is required for IPC and IN normalization.
#' Vector of character string(s) that represents IPCs in the column names (e.g. 'IPC'). 
#' @param SC_string Optional (or SC_wells) Vector of character string(s) that represents SCs in the column 
#' names (e.g. 'SC'). SC wells are included in the data file output.
#' @param Bridge_string Optional (or Bridge_wells) Vector of character string(s) that represents 
#' Bridge samples in the column names (e.g. 'Bridge'). Bridge wells 
#' can be used for inter-plate normalization (not currently implemented).
#' @param NC_string Required (or NC wells). Vector of character string(s) that represents NCs in the 
#' column names (e.g. 'NC'). 
#' @param include_IPC Logical. Should IPC samples be included in output? Default is FALSE.
#' @param include_SC Logical. Should SC samples be included in output? Default is TRUE.
#' @param include_Bridge Logical. Should Bridge samples be included in output? Default is TRUE.
#' @param include_NC Logical. Should NC samples be included in output? Default is FALSE.
#' @param excludeSamples List of character string vectors that give sample names to be
#' excluded from the output file. List should be in order of xml files. If 
#' no sample is to be excluded from a plate, use NULL in the list for that plate.
#' @param intraPlateNorm_method intra-plate normalization method passed to 
#' intraPlateNorm() function. Default is 'IC' (internal control). Other option is
#' 'TC' (total count), but this method is not currently implemented.
#' @param intraPlateNorm_scaleFactor Optional scaling factor to apply after 
#' intra-plate normalization. Passed to intraPlateNorm() function.
#' @param interPlateNorm_method Default is "IPC" for inter-plate control normalization. 
#' Use "IN" for intensity normalization. If neither, no interplate normalization is done.
#' @param IPC_method Passed to interPlateNorm function. 'median' is the default. 
#' Other options include 'mean' (arithmetic mean) and 'geom_mean' (geometric mean). 
#' Determines how the counts are summarized across the IPC wells on a given plate.
#' @param IN_samples Optional argument. Passed to interPlateNorm function. 
#' A list of column names or indices specifying which subset of samples to use 
#' for intensity normalization step for each plate in data_list. Will over-ride 
#' the IPC_wells and NC_wells arguments for IN.
#' @param interPlateNorm_dataScale Passed to interPlateNorm function. 'count' is 
#' the default and interplate normalization is multiplicative. Use option 'log' 
#' for log-transformed data; normalization is additive on the log scale.
#' @param interPlateNorm_scaleFactor Passed to interPlateNorm function. Optional 
#' numeric value used to rescale all data after normalizing. Default is 1. 
#' This may be desirable to avoid normalized quantities between 0 and 1 
#' (which will be negative in the log scale). Only useful for count scale data.
#' @param replaceNA Logical. Passed to readNULISAseq() function.
#' If TRUE (default), will replace missing counts with 
#' zeros.
#' @param verbose Logical. Should function output step completion into.
#' Default is TRUE.
#'
#' @return writes CSV file
#'
#'
#' @export
#'
writeNULISAseq <- function(xml_files,
                           dataDir,
                           target_info_file,
                           sample_info_file=NULL,
                           sample_info_file_variables=NULL,
                           output_filename,
                           Panel='200-plex Inflammation v1',
                           PanelLotNumber='',
                           plateIDs=NULL,
                           ICs='mCherry',
                           IPC_string='IPC',
                           SC_string='SC',
                           Bridge_string='Bridge',
                           NC_string='NC',
                           include_IPC=FALSE,
                           include_SC=TRUE,
                           include_Bridge=TRUE,
                           include_NC=FALSE,
                           excludeSamples=NULL,
                           intraPlateNorm_method='IC',
                           intraPlateNorm_scaleFactor=1,
                           interPlateNorm_method='IPC',
                           IPC_method='median',
                           IN_samples=NULL,
                           interPlateNorm_dataScale='count',
                           interPlateNorm_scaleFactor=1,
                           replaceNA=TRUE,
                           verbose=TRUE){
  n_plates <- length(xml_files)
  if(is.null(plateIDs)) plateIDs <- paste0('Plate_', c(paste0('0', 1:9), 10:99)[1:length(xml_files)])
  # read in NULISAseq data
  runs <- vector(mode="list", length=n_plates)
  for (i in 1:n_plates){
    runs[[i]] <- readNULISAseq(file=file.path(dataDir, xml_files[i]),
                               plateID=plateIDs[i],
                               IPC=IPC_string,
                               SC=SC_string,
                               NC=NC_string,
                               Bridge=Bridge_string,
                               replaceNA=replaceNA,
                               IC=ICs)
    # exclude samples if listed
    if (!is.null(excludeSamples[[i]])){
      runs[[i]]$Data <- runs[[i]]$Data[,!(colnames(runs[[i]]$Data) %in% excludeSamples[[i]])]
      runs[[i]]$samples <- runs[[i]]$samples[!(runs[[i]]$samples$sampleName %in% excludeSamples[[i]]),]
      runs[[i]]$SampleNames <- runs[[i]]$SampleNames[!(runs[[i]]$SampleNames %in% excludeSamples[[i]])]
    }
  }
  names(runs) <- plateIDs
  if(verbose==TRUE) cat(paste0(n_plates, ' XML files were read.\n'))
  # read target info files
  target_info <- read.csv(file.path(dataDir, target_info_file))
  # merge target info file with targets
  all_targets <- merge(target_info, runs[[1]]$targets,
                       by.x='TargetName',
                       by.y='targetName',
                       all.x=TRUE, all.y=FALSE)
  # get all sample data 
  sample_data <- lapply(runs, function(x) x$samples)
  sample_data <- do.call(rbind, sample_data)
  # read sample metadata file if given
  if(!is.null(sample_info_file)){
    sample_info <- read.csv(file.path(dataDir, sample_info_file))
    # get subset of sample metadata to merge
    if (!is.null(sample_info_file_variables)){
      sample_info <- sample_info[,c('plateID', 'sampleName', sample_info_file_variables)]
    } else {
      sample_info_file_variables <- colnames(sample_info)
      sample_info_file_variables <- sample_info_file_variables[!(sample_info_file_variables %in% c('plateID', 'sampleName'))]
    }
    # merge sample annotations with sample data from readNULISAseq
    sample_data <- merge(sample_data, sample_info, 
                         by=c('plateID', 'sampleName'), 
                         all.x=TRUE, all.y=FALSE)
  }
  
  # do intra-plate normalization -- IC 
  if(is.null(ICs)){
    warning('Argument ICs is NULL. Intra-plate normalization was not done.\n')
    Data <- lapply(runs, function(x) x$Data)
  }
  if(!is.null(ICs) & intraPlateNorm_method=='IC'){
    intraPlateNorm_data <- lapply(runs, function(x){
      intraPlateNorm(data_matrix=x$Data,
                     method=intraPlateNorm_method,
                     IC=ICs[1],
                     scaleFactor=intraPlateNorm_scaleFactor)
    })
    if(verbose==TRUE) cat('Intra-plate IC normalization completed.\n')
    Data <- lapply(intraPlateNorm_data, function(x) x$normData)
  }
  
  # if more than one plate, do inter-plate normalization
  if(n_plates > 1){
    if(interPlateNorm_method=='IPC'){
      interPlateNorm_data <- interPlateNorm(data_list=lapply(intraPlateNorm_data, function(x) x$normData),
                                            IPC=TRUE, IN=FALSE,
                                            IPC_wells=lapply(runs, function(x) x$IPC),
                                            IPC_method=IPC_method,
                                            dataScale=interPlateNorm_dataScale,
                                            scaleFactor=interPlateNorm_scaleFactor)
      Data <- interPlateNorm_data$interNormData
      if(verbose==TRUE) cat('Inter-plate IPC normalization completed.\n')
    } else if(interPlateNorm_method=='IN'){
      interPlateNorm_data <- interPlateNorm(data_list=lapply(intraPlateNorm_data, function(x) x$normData),
                                            IPC=FALSE, IN=TRUE,
                                            IPC_wells=NULL,
                                            NC_wells=lapply(runs, function(x) x$NC),
                                            IN_samples=IN_samples,
                                            dataScale=interPlateNorm_dataScale,
                                            scaleFactor=interPlateNorm_scaleFactor)
      Data <- interPlateNorm_data$interNormData
      if(verbose==TRUE) cat('Inter-plate intensity normalization completed.\n')
    } else if(interPlateNorm_method!='IPC' & interPlateNorm_method!='IPC'){
      warning('No valid inter-plate normalization method was specified.')
    }
  }
  
  # calculate LODs 
  LODs <- vector(mode='list', length=n_plates)
  for(i in 1:n_plates){
    LODs[[i]] <- lod(data_matrix=Data[[i]],
                     blanks=runs[[i]]$NC,
                     min_count=0)$LOD
  }
  if(verbose==TRUE) cat('LOD calculation completed.\n')
  
  # sample QC
  # check if mCherry within +/- 30% of median
  sampleQC <- vector(mode='list', length=n_plates)
  sampleQC_IC_Median <- vector(mode='list', length=n_plates)
  for(i in 1:n_plates){
    sampleQC[[i]] <- QCFlagSample(raw=runs[[i]]$Data,
                                  normed=intraPlateNorm_data[[i]]$normData,
                                  samples=runs[[i]]$samples,
                                  targets=runs[[i]]$targets,
                                  well_order=NULL,
                                  ICs=runs[[i]]$IC,
                                  IPCs=runs[[i]]$IPC,
                                  NCs=runs[[i]]$NC)
    sampleQC_IC_Median[[i]] <- sampleQC[[i]]$sampleName[sampleQC[[i]]$flagName=='IC_Median' & sampleQC[[i]]$status=='T']
  }
  sample_data$SampleQC <- 'PASS'
  sample_data$SampleQC[sample_data$sampleName %in% unlist(sampleQC_IC_Median)] <- 'WARN'
  if(verbose==TRUE) cat('Sample QC completed.\n')
  
  
  # samples to include
  sampleNames_output <- lapply(runs, function(x) {
    samples <- x$SampleNames
    if(include_IPC==TRUE & !is.null(x$IPC)) samples <- c(samples, x$IPC)
    if(include_SC==TRUE & !is.null(x$SC)) samples <- c(samples, x$SC)
    if(include_Bridge==TRUE & !is.null(x$Bridge)) samples <- c(samples, x$Bridge)
    if(include_NC==TRUE & !is.null(x$NC)) samples <- c(samples, x$NC)
    plateID <- rep(x$plateID, length(samples))
    return(list(samples=samples, plateID=plateID))
  })
  
  # transform to long format for each target
  targets <- all_targets$TargetName
  # create the long data list
  long_data <- vector('list', length=length(targets))
  for (i in 1:length(targets)){
    target <- targets[i]
    SampleInfo <- data.frame(SampleName=unlist(lapply(sampleNames_output, function(x) x$samples)),
                             plateID=unlist(lapply(sampleNames_output, function(x) x$plateID)))
    SampleInfo <- merge(SampleInfo, 
                        sample_data[,c("plateID","sampleName","sampleType",sample_info_file_variables,'SampleQC')],
                        by.x=c('plateID', 'SampleName'), by.y=c('plateID', 'sampleName'),
                        all.x=TRUE, all.y=FALSE)
    colnames(SampleInfo)[colnames(SampleInfo)=="sampleType"] <- "SampleType"
    colnames(SampleInfo)[colnames(SampleInfo)=="plateID"] <- "PlateID"
    # get target specific data
    AlamarTargetID <- all_targets$AlamarTargetID[all_targets$TargetName==target]
    UniProtID <- all_targets$UniProtID[all_targets$TargetName==target]
    ProteinName <- all_targets$ProteinName[all_targets$TargetName==target]
    
    # get plate-specific and target-specific data
    target_data <- vector(mode='list', length=n_plates)
    for(j in 1:n_plates){
      plate_j_LOD <- log2(LODs[[j]][target] + 0.01)
      plate_j_log2NormalizedCount <- log2(Data[[j]][target,] + 0.01)
      data_length <- length(plate_j_log2NormalizedCount)
      target_data[[j]] <- data.frame(SampleName=names(plate_j_log2NormalizedCount), 
                                     Panel=rep(Panel, data_length),
                                     PanelLotNumber=rep(PanelLotNumber, data_length),
                                     PlateID=rep(plateIDs[j], data_length),
                                     Target=rep(target, data_length),
                                     AlamarTargetID=rep(AlamarTargetID, data_length),
                                     UniProtID=rep(UniProtID, data_length),
                                     ProteinName=rep(ProteinName, data_length),
                                     LOD=rep(plate_j_LOD, data_length),
                                     log2NormalizedCount=plate_j_log2NormalizedCount)
    }
    target_data <- do.call(rbind, target_data)
    
    
    target_data <- merge(SampleInfo, target_data, all.x=TRUE, all.y=FALSE,
                         by=c('PlateID', 'SampleName'))
    target_data$SampleType <- factor(target_data$SampleType, levels=c("Sample", "IPC", "SC", "Bridge", "NC"))
    target_data <- target_data[order(target_data$PlateID, target_data$SampleType, target_data$SampleName),]
    target_data <- target_data[,c("Panel","PanelLotNumber","PlateID",
                                  "SampleName","SampleType",
                                  sample_info_file_variables,"Target",
                                  "AlamarTargetID","UniProtID","ProteinName",
                                  "SampleQC","LOD","log2NormalizedCount")]
    
    long_data[[i]] <- target_data
  }
  
  # merge targets
  allData <- do.call(rbind, long_data)
  
  # save data
  write.csv(allData, 
            file.path(dataDir, output_filename),
            row.names=FALSE)
  
  if(verbose==TRUE) cat('writeNULISAseq completed.')
}

