#' Write NULISAseq normalized data CSV in long format
#'
#' Reads in one or more NULISAseq XML files, where XML file is output from the Alamar
#' Biosciences Galaxy NULISAseq (Beta) tool. This function then performs intra-plate 
#' normalization (default is internal control (IC) mCherry), and when multiple 
#' XML files are input, also performs inter-plate normalization (default is inter-plate 
#' control (IPC)). Finally, data is log2-transformed. 
#' 
#' The function also takes as input a target info file containing target IDs 
#' and protein names, and an optional sample info file containing sample metadata that
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
#' Columns of output data (in order) include:
#' PanelLotNumber
#' PlateID
#' SampleID
#' (optional sample metadata columns)
#' Target
#' AlamarTargetID
#' UniProtID
#' ProteinName
#' SampleQC
#' LOD
#' log2NormalizedCount
#' 
#'
#' @param xml_files Vector of filenames (character strings). Files should 
#' be in order of the desired plateID variable, unless that is otherwise
#' defined.
#' @param dataDir Data directory where xml_files, target_info_file, and sample_info_file 
#' reside (character string). 
#' @param target_info_file Path and filename for the target info CSV file 
#' (character string). Must include columns for Target (matches the targets in 
#' XML files), AlamarTargetID, UniProtID, and ProteinName.
#' @param sample_info_file Path and filename for the sample annotation CSV file. 
#' One column should be "sampleName" which corresponds to the column names 
#' of the readNULISAseq Data matrix and the sampleName in the readNULISAseq 
#' samples data.frame.
#' @param sample_info_file_variables Optional subset of column names in 
#' sample_info_file that will be included in data file output. Other columns 
#' will be excluded. Otherwise, if NULL (default), all columns will be included.
#' @param output_filename Filename for output CSV file. 
#' @param IPC_string IPC_string (or IPC_wells) is required for IPC and IN normalization.
#' Vector of character string(s) that represents IPCs in the column names (e.g. 'IPC'). 
#' @param SC_string Optional (or SC_wells) Vector of character string(s) that represents SCs in the column 
#' names (e.g. 'SC'). SC wells are included in the data file output.
#' @param Bridge_string Optional (or Bridge_wells) Vector of character string(s) that represents 
#' Bridge samples in the column names (e.g. 'Bridge'). Bridge wells 
#' can be used for inter-plate normalization (not currently implemented).
#' @param NC_string Required (or NC wells). Vector of character string(s) that represents NCs in the 
#' column names (e.g. 'NC'). 
#' @param IPC_wells Required for IPC and IN normalization. (IPC_string overrides 
#' this if provided.) A list of vectors of column names for each plate that match the 
#' inter-plate control (IPC) column names in the Data matrix.
#' @param SC_wells Optional. (SC_string overrides this if provided.) 
#' A list of vectors of column names for each plate that match the 
#' sample control (SC) column names in the Data matrix. These samples will be 
#' labelled as SC in the data output file.
#' @param Bridge_wells Optional (Bridge_string overrides this if provided.) 
#' List of vectors of column names for Bridge samples. Bridge wells 
#' can be used for inter-plate normalization (not currently implemented).
#' @param NC_wells Required. (NC_string overrides this if provided.)  
#' A list of vectors of column names for each plate that match the 
#' negative control (NC) column names in the Data matrix.
#' @param omit_samples Optional. A list of vectors of columns names for each plate that match
#' samples to omit from the data output file (e.g. bridge wells).
#' @param omit_targets Optional. Vector of targets (such as ICs) to omit from data 
#' output file. Default is 'mCherry'.
#' @param Panel Name of multi-plex panel. Default is '200plex v1 Inflammation'
#' @param PanelLotNumber The panel lot number.
#' @param plateIDs A vector of plate IDs. Default is to number plates from 01 to 
#' total number of plates based on the order of xml_files. Passed to readNULISAseq() function
#' and output as a column in the data file. 
#' @param intraPlateNorm_method intra-plate normalization method passed to 
#' intraPlateNorm() function. Default is 'IC' (internal control). Other option is
#' 'TC' (total count), but this method is not currently implemented.
#' @param intraPlateNorm_IC String vector denoting the IC / ICs. Default is 'mCherry'.
#' Passed to intraPlateNorm() function.
#' @param intraPlateNorm_scaleFactor Optional scaling factor to apply after 
#' intra-plate normalization. Passed to intraPlateNorm() function.
#' @param interPlateNorm_method_IPC Logical TRUE/FALSE. Should IPC normalization
#' be used.
#' @param interPlateNorm_method_IN Logical TRUE/FALSE. Should IN normalization 
#' be used. If both of the above are false, no inter-plate normalization is done.
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
#' @param file_type Character string. Type of input file, as output from Galaxy. 
#' Passed to readNULISAseq() function. Options include
#' xml_no_mismatches (default), xml_full_output (not currently implemented) 
#' (both from NULISAseq tool),
#' or xml_normalization (from NULISAseq Normalization tool, 
#' not currently implemented).
#' @param replaceNA Logical. Passed to readNULISAseq() function.
#' If TRUE (default), will replace missing counts with 
#' zeros.
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
                           IPC_string=NULL,
                           SC_string=NULL,
                           Bridge_string=NULL,
                           NC_string=NULL,
                           IPC_wells=NULL,
                           SC_wells=NULL,
                           Bridge_wells=NULL,
                           NC_wells=NULL,
                           omit_samples=NULL,
                           omit_targets='mCherry',
                           Panel='200plex v1 Inflammation',
                           PanelLotNumber='',
                           plateIDs=paste0('Plate_', c(paste0('0', 1:9), 10:99)[1:length(xml_files)]),
                           intraPlateNorm_method='IC',
                           intraPlateNorm_IC='mCherry',
                           intraPlateNorm_scaleFactor=1,
                           interPlateNorm_method_IPC=TRUE,
                           interPlateNorm_method_IN=FALSE, 
                           IPC_method='median',
                           IN_samples=NULL,
                           interPlateNorm_dataScale='count',
                           interPlateNorm_scaleFactor=1,
                           file_type='xml_no_mismatches',
                           replaceNA=TRUE){
  n_plates <- length(xml_files)
  # read in NULISAseq data
  plate_data <- vector(mode="list", length=n_plates)
  for (i in 1:n_plates){
    plate_data[[i]] <- readNULISAseq(xml_file=paste0(dataDir, xml_files[i]),
                                     plateID=plateIDs[i])
  }
  names(plate_data) <- plateIDs
  # save the IPC wells as elements of the data list
  if (!is.null(IPC_string)){
    lapply(plate_data, function(x){
      x$IPC_wells <- x$samples$sampleName[grep(paste(IPC_string, collapse="|"), x$samples$sampleName)]
    })
  } 
  # save the SC wells as elements of the data list
  if (!is.null(SC_string)){
    plate_data$SC_wells <- lapply(plate_data, function(x){
      x$samples$sampleName[grep(paste(SC_string, collapse="|"), x$samples$sampleName)]
    })
  } 
  # save the Bridge wells as elements of the data list
  if (!is.null(Bridge_string)){
    plate_data$Bridge_wells <- lapply(plate_data, function(x){
      x$samples$sampleName[grep(paste(Bridge_string, collapse="|"), x$samples$sampleName)]
    })
  } 
  # save the NC wells as elements of the data list
  if (!is.null(NC_string)){
    plate_data$NC_wells <- lapply(plate_data, function(x){
      x$samples$sampleName[grep(paste(NC_string, collapse="|"), x$samples$sampleName)]
    })
  } 
  # read target and sample info files
  target_info <- read.csv(file.path(dataDir, target_info_file))
  sample_info <- read.csv(file.path(dataDir, sample_info_file))
  # get subset of sample metadata to merge
  if (!is.null(sample_info_file_variables)){
    sample_info <- sample_info[,c('sampleName', sample_info_file_variables)]
  }
  # merge sample annotations with sample data from readNULISAseq
  sample_data <- lapply(plate_data, function(x) x$samples)
  sample_data <- do.call(rbind, sample_data)
  sample_data <- merge(sample_data, sample_info, 
                       by='sampleName', 
                       all.x=TRUE, all.y=FALSE)
  # need to figure out how to match sample ID ... partial match to sampleName
  
  
  # do intra-plate normalization -- IC 
  if(intraPlateNorm_method=='IC'){
    intraPlateNorm_data <- lapply(plate_data, function(x){
      intraPlateNorm(data_matrix=x$Data,
                     method=intraPlateNorm_method,
                     IC=intraPlateNorm_IC,
                     scaleFactor=intraPlateNorm_scaleFactor)
    })
  }
  
  # do inter-plate normalization
  
  
  
  
  # do interplate normalization mCherry + IPC
  mCherryIPC <- interPlateNorm(data_list = list(plate1_mCherry$normData,
                                                plate2_mCherry$normData),
                               IPC=TRUE, IN=FALSE, 
                               IPC_wells=list(c("C_12_P01_IPC_rep01","E_12_P01_IPC_rep03","D_12_P01_IPC_rep02"),
                                              c("D_12_P02_IPC_rep02","C_12_P02_IPC_rep01","E_12_P02_IPC_rep03")))
  
  # merge target info file with targets
  all_targets <- merge(plate1$targets, target_panel_info,
                       by.x='targetName',
                       by.y='Primary.Gene.Name',
                       all.x=FALSE, all.y=TRUE)
  
}




# merge annotations with the other sample data
plate1_samples <- plate1$samples
plate1_samples$SampleID <- substr(plate1_samples$sampleName, start=10, stop=20)
plate1_samples$SampleID[grep('NC', plate1_samples$SampleID)] <- NA
plate1_samples$SampleID[grep('IPC', plate1_samples$SampleID)] <- 'IPC'

plate2_samples <- plate2$samples
plate2_samples$SampleID <- substr(plate2_samples$sampleName, start=10, stop=20)
plate2_samples$SampleID[grep('NC', plate2_samples$SampleID)] <- NA
plate2_samples$SampleID[grep('IPC', plate2_samples$SampleID)] <- 'IPC'





# calculate LODs 
plate1_lod <- lod(data_matrix=plate1_mCherry$normData, 
                  blanks=c('F_12_P01_NC_rep01', 
                           'G_12_P01_NC_rep02', 
                           'H_12_P01_NC_rep03'), 
                  min_count=0)

plate2_lod <- lod(data_matrix=plate2_mCherry$normData, 
                  blanks=c('F_12_P02_NC_rep01',
                           'G_12_P02_NC_rep02',
                           'H_12_P02_NC_rep03'), 
                  min_count=0)

# calculate detectability
detect_samples1 <- colnames(plate1_mCherry$normData)
detect_samples1 <- detect_samples1[-grep('IPC', detect_samples1)]
detect_samples1 <- detect_samples1[-grep('NC', detect_samples1)]
plate1_detect <- detectability(aboveLOD_matrix=plate1_lod$aboveLOD,
                               sample_subset=detect_samples1,
                               exclude_targets=c('mCherry'))

detect_samples2 <- colnames(plate2_mCherry$normData)
detect_samples2 <- detect_samples2[-grep('IPC', detect_samples2)]
detect_samples2 <- detect_samples2[-grep('NC', detect_samples2)]
plate2_detect <- detectability(aboveLOD_matrix=plate2_lod$aboveLOD,
                               sample_subset=detect_samples2,
                               exclude_targets=c('mCherry'))

# calculate intra-plate CVs mCherry
# use all 3 IPCs
IPCs <- rep(NA, 96)
IPCs[grep('IPC', colnames(plate1_mCherry$normData))] <- 'IPC'
plate1_mCherry_intraCV <- intraCV(data_matrix=plate1_mCherry$normData,
                                  samples=IPCs,
                                  aboveLOD=plate1_lod$aboveLOD,
                                  exclude_targets=c('mCherry'))

mean(plate1_mCherry_intraCV, na.rm=TRUE)
median(plate1_mCherry_intraCV, na.rm=TRUE)

IPCs <- rep(NA, 96)
IPCs[grep('IPC', colnames(plate2_mCherry$normData))] <- 'IPC'
plate2_mCherry_intraCV <- intraCV(data_matrix=plate2_mCherry$normData,
                                  samples=IPCs,
                                  aboveLOD=plate2_lod$aboveLOD,
                                  exclude_targets=c('mCherry'))

mean(plate2_mCherry_intraCV, na.rm=TRUE)
median(plate2_mCherry_intraCV, na.rm=TRUE)


# sample QC
# total sample counts >= 500*(number of targets)
table(colSums(plate1$Data) >= 500*203) 
colSums(plate1$Data)[colSums(plate1$Data) < 500*203]
# check if mCherry within +/- 30% of median
mCherry_median <- median(plate1$Data['mCherry',])
mCherry_lower <- 0.7*mCherry_median
mCherry_upper <- 1.3*mCherry_median
sum(plate1$Data['mCherry',] < mCherry_lower) # 2 sample failed
sum(plate1$Data['mCherry',] > mCherry_upper) # all pass
colnames(plate1$Data)[plate1$Data['mCherry',] < mCherry_lower] # "F_05_P01_S-001992264" "G_07_P01_S-001992276" failed

mCherry_median <- median(plate2$Data['mCherry',])
mCherry_lower <- 0.7*mCherry_median
mCherry_upper <- 2.3*mCherry_median
sum(plate2$Data['mCherry',] < mCherry_lower) # 1 sample failed
sum(plate2$Data['mCherry',] > mCherry_upper) # all pass
colnames(plate2$Data)[plate2$Data['mCherry',] < mCherry_lower] # "F_05_P02_S-001992348" failed

QC_warn_samples <- c("F_05_P01_S-001992264",
                     "G_07_P01_S-001992276",
                     "F_05_P02_S-001992348")
all_samples$SampleQC <- 'PASS'
all_samples$SampleQC[all_samples$sampleName %in% QC_warn_samples] <- 'WARN'

# transform to long format
# for each target
targets <- all_targets$targetName
sampleNames <- c(detect_samples1,
                 detect_samples2,
                 all_samples$sampleName[all_samples$SampleID=='IPC' & !is.na(all_samples$SampleID)])

long_data <- vector('list', length=length(targets))
for (i in 1:length(targets)){
  target <- targets[i]
  Panel <- '200plex v1 Inflammation'
  PanelLotNumber <- 'p20230117_20221215'
  SampleInfo <- data.frame(sampleName=c(colnames(plate1_mCherry$normData),
                                        colnames(plate2_mCherry$normData)))
  SampleInfo <- merge(SampleInfo, all_samples[,c("plateID","sampleName","SampleID","Condition",'SampleQC')])
  colnames(SampleInfo)[1] <- 'SampleName'
  # remove non-samples
  SampleInfo <- SampleInfo[SampleInfo$SampleName %in% sampleNames,]
  AlamarTargetID <- all_targets$AlamarTargetID[all_targets$targetName==target]
  UniProtID <- all_targets$UniProtID[all_targets$targetName==target]
  ProteinName <- all_targets$ProteinName[all_targets$targetName==target]
  plate1_LOD <- log2(plate1_lod$LOD[target] + 0.01)
  plate2_LOD <- log2(plate2_lod$LOD[target] + 0.01)
  plate1_log2NormalizedCount <- log2(mCherryIPC$interNormData[[1]][target,] + 0.01)
  plate2_log2NormalizedCount <- log2(mCherryIPC$interNormData[[2]][target,] + 0.01)
  plate1_target_data <- data.frame(SampleName=names(plate1_log2NormalizedCount), 
                                   Panel=Panel,
                                   PanelLotNumber=PanelLotNumber,
                                   PlateID='Plate_01',
                                   Target=target,
                                   AlamarTargetID=AlamarTargetID,
                                   UniProtID=UniProtID,
                                   ProteinName=ProteinName,
                                   LOD=plate1_LOD,
                                   log2NormalizedCount=plate1_log2NormalizedCount)
  plate2_target_data <- data.frame(SampleName=names(plate2_log2NormalizedCount), 
                                   Panel=Panel,
                                   PanelLotNumber=PanelLotNumber,
                                   PlateID='Plate_02',
                                   Target=target,
                                   AlamarTargetID=AlamarTargetID,
                                   UniProtID=UniProtID,
                                   ProteinName=ProteinName,
                                   LOD=plate2_LOD,
                                   log2NormalizedCount=plate2_log2NormalizedCount)
  target_data <- rbind(plate1_target_data, plate2_target_data)
  target_data <- merge(SampleInfo, target_data, all.x=TRUE, all.y=FALSE)
  target_data <- target_data[order(target_data$PlateID),]
  target_data <- target_data[,c("Panel","PanelLotNumber","PlateID",
                                "SampleName","SampleID","Condition","Target",
                                "AlamarTargetID","UniProtID","ProteinName",
                                "SampleQC","LOD","log2NormalizedCount")]
  
  # remove non-samples
  target_data <- target_data[target_data$SampleName %in% sampleNames,]
  long_data[[i]] <- target_data
}

allData <- do.call(rbind, long_data)

# save data
write.csv(allData, 
          paste0(dataDir, 'Harvard_Infection_Proneness_NULISAseq_200plex_v1_Inflammation_2022-02-09.csv'),
          row.names=FALSE)


