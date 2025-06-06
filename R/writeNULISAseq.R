#' Write NULISAseq normalized data CSV in long format
#'
#' @description
#' Reads in one or more NULISAseq XML files. Output is a long-format CSV file
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
#'  \item{"NPQ"}{}
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
#' @param Calibrator_string Optional (or Calibrator_wells) Vector of character string(s) that represents 
#' Calibrator samples in the column names (e.g. 'Calibrator'). Calibrator wells 
#' can be used for absolute quantification (not currently implemented).
#' @param NC_string Required (or NC wells). Vector of character string(s) that represents NCs in the 
#' column names (e.g. 'NC'). 
#' @param include_IPC Logical. Should IPC samples be included in output? Default is FALSE.
#' @param include_SC Logical. Should SC samples be included in output? Default is TRUE.
#' @param include_Bridge Logical. Should Bridge samples be included in output? Default is TRUE.
#' @param include_Calibrator Logical. Should Calibrator samples be included in output? Default is TRUE.
#' @param include_NC Logical. Should NC samples be included in output? Default is FALSE.
#' @param include_unnorm_counts Logical. Should unnormalized counts be included 
#' as am additional column in output? Default is FALSE.
#' @param include_IC_counts Logical. Should IC counts be included in the output? 
#' Default is FALSE. This is probably only useful when 
#' \code{include_unnorm_counts=TRUE}.
#' @param excludeSamples List of character string vectors that give sample names to be
#' excluded from the output file. List should be in order of xml files. If 
#' no sample is to be excluded from a plate, use NULL in the list for that plate.
#' @param intraPlateNorm_method intra-plate normalization method passed to 
#' intraPlateNorm() function. Default is 'IC' (internal control). Other option is
#' 'TC' (total count), but this method is not currently implemented.
#' @param intraPlateNorm_scaleFactor Optional scaling factor to apply after 
#' intra-plate normalization. Passed to intraPlateNorm() function.
#' @param interPlateNorm_method Default is "IPC" for inter-plate control normalization. 
#' Use "IN" for IPC normalization followed by intensity normalization. 
#' If neither, no interplate normalization is done.
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
#' @param interPlateNorm_transformReverse_covariateName The name of the 
#' covariate in the Barcode A file that indicates whether or not to use
#' reverse curve transformation for each target. This column will have 
#' an R entry for reverse curve targets, and F for forward curve targets. 
#' Default is \code{"Curve_Quant"}. Please note that if there are multiple
#' runs, writeNULISAseq will only use the first run target data to determine
#' the reverse curve targets!
#' @param interPlateNorm_transformReverse_scaleFactor The scaling factor used in the 
#' reverse curve transformation. Default is 1e4. Reverse curve transformation is 
#' \code{transformReverse_scaleFactor / (IPC normalized count + 1)}. Then the log2
#' tranformation is applied to this value, after adding 1, to obtain NPQ.
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
                           Calibrator_string='Calibrator',
                           NC_string='NC',
                           include_IPC=FALSE,
                           include_SC=TRUE,
                           include_Bridge=TRUE,
                           include_Calibrator=TRUE,
                           include_NC=FALSE,
                           include_unnorm_counts=FALSE,
                           include_IC_counts=FALSE,
                           excludeSamples=NULL,
                           intraPlateNorm_method='IC',
                           intraPlateNorm_scaleFactor=1,
                           interPlateNorm_method='IPC',
                           IPC_method='median',
                           IN_samples=NULL,
                           interPlateNorm_dataScale='count',
                           interPlateNorm_scaleFactor=10^4,
                           interPlateNorm_transformReverse_covariateName='Curve_Quant',
                           interPlateNorm_transformReverse_scaleFactor=1e4,
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
                               Calibrator=Calibrator_string,
                               replaceNA=replaceNA,
                               IC=ICs)
    # remove AlamarTargetID, UniProtID, and ProteinName from the targets data.frame if present
    # so we don't have issue with merging target info file which has same column
    runs[[i]]$targets <- runs[[i]]$targets[,!(colnames(runs[[i]]$targets) %in% c('AlamarTargetID', 'UniProtID', 'ProteinName'))]
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
  run_targets_cols <- Reduce(intersect, lapply(runs, function(x) colnames(x$targets)))
  run_targets <- unique(do.call(rbind, lapply(runs, function(x) x$targets[,run_targets_cols])))
  # merge target info file with targets
  all_targets <- merge(target_info[target_info$TargetName %in% run_targets$targetName,], 
                       run_targets,
                       by.x='TargetName',
                       by.y='targetName',
                       all.x=TRUE, all.y=FALSE, sort=FALSE)
  all_targets <- all_targets[order(all_targets$TargetName),]
  if(include_IC_counts==TRUE){
    # assumes IC for run 1 is the IC used for all runs
    all_targets <- rbind(all_targets, c(ICs[1], rep(NA, ncol(all_targets) - 1)))
  }
  # get all sample data 
  sample_data <- lapply(runs, function(x) x$samples)
  # use only matching columns 
  matching_column <- Reduce(intersect, lapply(sample_data, colnames))
  sample_data <- lapply(sample_data, function(x) x[,matching_column])
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
                         all.x=TRUE, all.y=FALSE, sort=FALSE)
  }
  
  # save unnorm data if include_unnorm_counts==TRUE
  if(include_unnorm_counts==TRUE){
    unnorm_data <- lapply(runs, function(x) x$Data)
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
  
  # do IPC normalization
  if(interPlateNorm_method=='IPC'){
    transformReverse_targets <- runs[[1]]$targets$targetName[runs[[1]]$targets[,interPlateNorm_transformReverse_covariateName]=="R"]
    interPlateNorm_data <- interPlateNorm(data_list=lapply(intraPlateNorm_data, function(x) x$normData),
                                          IPC=TRUE, IN=FALSE,
                                          IPC_wells=lapply(runs, function(x) x$IPC),
                                          IPC_method=IPC_method,
                                          dataScale=interPlateNorm_dataScale,
                                          scaleFactor=interPlateNorm_scaleFactor,
                                          transformReverse=transformReverse_targets,
                                          transformReverse_scaleFactor=interPlateNorm_transformReverse_scaleFactor)
    Data <- interPlateNorm_data$interNormData
    if(verbose==TRUE) cat('Inter-plate IPC normalization completed.\n')
  } else if(interPlateNorm_method=='IN'){
    interPlateNorm_data <- interPlateNorm(data_list=lapply(intraPlateNorm_data, function(x) x$normData),
                                          IPC=TRUE, IN=TRUE,
                                          IPC_wells=lapply(runs, function(x) x$IPC),
                                          NC_wells=lapply(runs, function(x) x$NC),
                                          IN_samples=IN_samples,
                                          dataScale=interPlateNorm_dataScale,
                                          scaleFactor=interPlateNorm_scaleFactor,
                                          transformReverse=transformReverse_targets,
                                          transformReverse_scaleFactor=interPlateNorm_transformReverse_scaleFactor)
    Data <- interPlateNorm_data$interNormData
    if(verbose==TRUE) cat('Inter-plate IPC and intensity normalization completed.\n')
  } else if(interPlateNorm_method!='IPC' & interPlateNorm_method!='IN'){
    warning('No valid inter-plate normalization method was specified.')
  }
  
  # calculate LODs 
  LODs <- vector(mode='list', length=n_plates)
  for(i in 1:n_plates){
    # Determine which targets should NOT have outlier detection performed
    # Currently we require that at least one NC sample have at least 100 
    # raw reads for that target for outlier detection to be applied 
    # Therefore, if all NC samples have < 100 raw reads, no outlier detection
    # should be performed
    targetNoOutlierDetection <- names(which(apply(runs[[i]]$Data[, runs[[i]]$NC] < 100, 1, all)))
    LODs[[i]] <- lod(data_matrix=Data[[i]],
                     blanks=runs[[i]]$NC,
                     min_count=0,
                     targetNoOutlierDetection=targetNoOutlierDetection)$LOD
    # set LODs for reverse curve targets to NA
    if (length(transformReverse_targets) > 0){
      LODs[[i]][transformReverse_targets] <- NA
    }
  }
  if(verbose==TRUE) cat('LOD calculation completed.\n')
  
  # sample QC
  # check if mCherry within +/- 40% of median
  sampleQC <- vector(mode='list', length=n_plates)
  sampleQC_IC_Median <- vector(mode='list', length=n_plates)
  for(i in 1:n_plates){
    sampleQC[[i]] <- QCFlagSample(raw=runs[[i]]$Data,
                                  aboveLOD=LODs[[i]]$aboveLOD,
                                  samples=runs[[i]]$samples,
                                  targets=runs[[i]]$targets,
                                  well_order=NULL,
                                  ICs=which(runs[[i]]$targets$targetName %in% runs[[i]]$IC),
                                  IPCs=which(runs[[i]]$sample$sampleName %in% runs[[i]]$IPC),
                                  NCs=which(runs[[i]]$sample$sampleName %in% runs[[i]]$NC))
    sampleQC_IC_Median[[i]] <- sampleQC[[i]]$sampleName[sampleQC[[i]]$flagName=='IC_Median' & sampleQC[[i]]$status==TRUE]
  }
  sample_data$SampleQC <- 'PASS'
  for(i in 1:length(plateIDs)){
    sample_data$SampleQC[sample_data$plateID==plateIDs[i] & sample_data$sampleName %in% sampleQC_IC_Median[[i]]] <- 'WARN'
  }
  if(verbose==TRUE) cat('Sample QC completed.\n')
  
  # samples to include
  sampleNames_output <- lapply(runs, function(x) {
    samples <- x$SampleNames
    if(include_IPC==TRUE & !is.null(x$IPC)) samples <- c(samples, x$IPC)
    if(include_SC==TRUE & !is.null(x$SC)) samples <- c(samples, x$SC)
    if(include_Bridge==TRUE & !is.null(x$Bridge)) samples <- c(samples, x$Bridge)
    if(include_Calibrator==TRUE & !is.null(x$Calibrator)) samples <- c(samples, x$Calibrator)
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
                        all.x=TRUE, all.y=FALSE, sort=FALSE)
    colnames(SampleInfo)[colnames(SampleInfo)=="sampleType"] <- "SampleType"
    colnames(SampleInfo)[colnames(SampleInfo)=="plateID"] <- "PlateID"
    # get target specific data
    AlamarTargetID <- all_targets$AlamarTargetID[all_targets$TargetName==target]
    UniProtID <- all_targets$UniProtID[all_targets$TargetName==target]
    ProteinName <- all_targets$ProteinName[all_targets$TargetName==target]
    
    # get plate-specific and target-specific data
    target_data <- vector(mode='list', length=n_plates)
    for(j in 1:n_plates){
      plate_j_LOD <- log2(LODs[[j]][target] + 1)
      plate_j_NPQ <- log2(Data[[j]][target,] + 1)
      data_length <- length(plate_j_NPQ)
      
      if(include_unnorm_counts==FALSE){
        target_data[[j]] <- data.frame(SampleName=names(plate_j_NPQ), 
                                       Panel=rep(Panel, data_length),
                                       PanelLotNumber=rep(PanelLotNumber, data_length),
                                       PlateID=rep(plateIDs[j], data_length),
                                       Target=rep(target, data_length),
                                       AlamarTargetID=rep(AlamarTargetID, data_length),
                                       UniProtID=rep(UniProtID, data_length),
                                       ProteinName=rep(ProteinName, data_length),
                                       LOD=rep(plate_j_LOD, data_length),
                                       NPQ=plate_j_NPQ)
        target_data_colnames <- c("Panel","PanelLotNumber","PlateID",
                                  "SampleName","SampleType",
                                  sample_info_file_variables,"Target",
                                  "AlamarTargetID","UniProtID","ProteinName",
                                  "SampleQC","LOD","NPQ")
      } else if(include_unnorm_counts==TRUE){
        plate_j_UnnormalizedCount <- unnorm_data[[j]][target,]
        target_data[[j]] <- data.frame(SampleName=names(plate_j_NPQ), 
                                       Panel=rep(Panel, data_length),
                                       PanelLotNumber=rep(PanelLotNumber, data_length),
                                       PlateID=rep(plateIDs[j], data_length),
                                       Target=rep(target, data_length),
                                       AlamarTargetID=rep(AlamarTargetID, data_length),
                                       UniProtID=rep(UniProtID, data_length),
                                       ProteinName=rep(ProteinName, data_length),
                                       LOD=rep(plate_j_LOD, data_length),
                                       UnnormalizedCount=plate_j_UnnormalizedCount,
                                       NPQ=plate_j_NPQ)
        target_data_colnames <- c("Panel","PanelLotNumber","PlateID",
                                  "SampleName","SampleType",
                                  sample_info_file_variables,"Target",
                                  "AlamarTargetID","UniProtID","ProteinName",
                                  "SampleQC","LOD","UnnormalizedCount","NPQ")
      }
    }
    target_data <- do.call(rbind, target_data)
    
    
    target_data <- merge(SampleInfo, target_data, all.x=TRUE, all.y=FALSE,
                         by=c('PlateID', 'SampleName'), sort=FALSE)
    target_data$SampleType <- factor(target_data$SampleType, levels=c("Sample", "IPC", "SC", "Bridge", "Calibrator", "NC"))
    target_data <- target_data[order(target_data$PlateID, target_data$SampleType),]
    target_data <- target_data[,target_data_colnames]
    
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

