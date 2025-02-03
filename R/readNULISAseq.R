#' Rename Duplicates in a List with Incrementing Values
#'
#' This function renames duplicate entries in a list by adding an incrementing
#' value to each duplicate while leaving unique elements unchanged.
#'
#' @param input_list A vector or list containing elements to be renamed.
#'
#' @return A vector with duplicate entries renamed by adding incrementing
#' values, while unique elements remain unchanged.
#'
#' @examples
#' original_list <- c("a", "b", "a", "b", "c", "d", "a")
#' new_list <- renameDuplicateNames(original_list)
#' print(new_list)
#'
#' # Output: "a_1" "b_1" "a_2" "b_2" "c" "d" "a_3"
#'
#' @export
renameDuplicateNames <- function(input_list, order_indices=NULL){
  unique_names <- unique(input_list)
  counts <- table(input_list)
  new_names <- input_list
  
  if(is.null(order_indices)){
    order_indices <- 1:length(input_list)
  }
  for (name in unique_names) {
    if (counts[name] > 1) {
      indices <- which(input_list == name)
      order_indices2 <- order_indices[which(order_indices %in% indices)]
      for (i in 1:length(order_indices2)) {
        new_names[order_indices2[i]] <- paste0(name, "_", i)
      }
    }
  }
  
  return(new_names)
}

#' Calculate the number of SCs (Sample Controls) in the dataset
#'
#' This function calculates the number of SCs (Sample Controls) in a dataset based on 
#' unique replicate naming. It works by matching sample types and target combinations in 
#' the data, using a pattern to identify SCs or AQSCs (depending on the `AQ` parameter).
#'
#' @param Data A data frame that includes at least two columns: `SampleType` and `Target`.
#'             The `SampleType` represents the type of sample, and `Target` refers to the 
#'             target analyte being measured.
#' @param AQ A logical flag indicating whether the function should search for AQSCs 
#'           (if `TRUE`) or standard SCs (if `FALSE`). AQSC refers to 'Absolute Quantification Sample Control'.
#'
#' @return The function returns an integer indicating the number of SCs found in the data.
#'
#' @details The function identifies replicates based on the concatenation of `SampleType` 
#' and `Target`. It looks for sample names that match either "AQSC_" or "SC_" patterns, depending 
#' on the value of the `AQ` flag. The number of replicates for the identified SC type is 
#' then calculated and returned.
#'
#' @examples
#' # Example usage:
#' numSCs(data, AQ = TRUE)
#' numSCs(data, AQ = FALSE)
#'
numSCs <- function(Data, AQ){
  # Calculate # of SCs in a hacky way, as the replicates are named the same in ACC
  r <- paste0(Data$SampleType, "_", Data$Target)
  uR <- unique(r)
  exampleAQSC <- ifelse (AQ, which(grepl("^AQSC_", uR))[1], which(grepl("^SC_", uR))[1])
  num <- length(which(uR[exampleAQSC] == r))
  return(num)
}  

#' Rename Sample Controls (SC) and Absolute Quantification Sample Controls (AQSC)
#'
#' This function renames instances of "SC" and "AQSC" in the input column by appending
#' a counter to each occurrence. The counter cycles through values from 1 to `n`,
#' and resets back to 1 after reaching the specified limit.
#'
#' @param column A character vector or column of values to be renamed. It should contain
#'               "SC" and "AQSC" values that need to be renamed.
#' @param label A string specifying the first label to be renamed (default is "SC").
#' @param label2 A string specifying the second label to be renamed (default is "AQSC").
#' @param n An integer specifying the maximum value for the counter. The counter cycles
#'          from 1 to `n` for each label occurrence (default is 3).
#'
#' @return A character vector with the renamed values. The labels "SC" and "AQSC" will
#'         have a counter appended to them (e.g., "SC_1", "AQSC_1", etc.), while other
#'         values will remain unchanged.
#'
#' @details The function works by iterating over the input vector and checking each value.
#'          If the value is equal to `label` or `label2`, a counter is appended to it.
#'          The counter resets to 1 after reaching the value `n`. This is useful for ensuring
#'          unique naming of standard controls and assay quality standard controls in a dataset.
#'
#' @examples
#' # Example usage:
#' renameSC(c("SC", "AQSC", "SC", "Other", "AQSC"))
#' renameSC(c("SC", "AQSC", "SC", "Other", "AQSC"), n = 2)
#'
renameSC <- function(column, label="SC", label2="AQSC", n=3){
  counter <- 1
  sapply(column, function(val) {
    if (val == label | val == label2) {
      result <- paste0(val, "_", counter)
      counter <<- ifelse(counter < n, counter + 1, 1)  # Cycle through 1 to n
      return(result)
    }
    return(val)
  })
}


#' Read NULISAseq XML
#'
#' Reads NULISAseq XML file
#'
#' @param file Character string. Path and name of the file.
#' @param plateID Character string that denotes plate ID.
#' @param file_type Character string. Type of input file.
#' Options include \code{xml_no_mismatches}, \code{xlsx} (output from ACC),
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
#' and the following: \code{c('Panel', 'PanelLotNumber', 'LOD', 'NPQ')}. 
#' The sample-specific data
#' will include by default \code{c('PlateID', 'SampleName', 'SampleType', 'SampleQC')}.
#' In addition, by default, it will include any other variables in the dataset.
#' If \code{sample_column_names} is specified, however, it will only include
#' those variables in the sample data frame output. 
#' @param sample_group_covar string that represents the name of the covariate
#' that defines sample group information. Defaults to NULL. Typically, this covariate
#' defines whether each sample is plasma, serum, csf, etc.
#' This information can be used when reporting detectability. For instance report sample
#' group specific detectability values, etc.).
#' @param IC string(s) that represents the names of internal control targets. 
#' Default is \code{'mCherry'}.Only used for xml file formats.
#' @param IPC string(s) present in the sample names
#' that represent the inter-plate control wells. 
#' For example, \code{'IPC'}. Set to \code{NULL} (default) to use the 
#' type variable from the Barcode B file or if there are no IPCs.
#' Only used for xml file formats.
#' @param SC string(s) present in the sample names
#' that represent the sample control wells. 
#' For example, \code{'SC'}. Set to \code{NULL} (default) to use the 
#' type variable from the Barcode B file or if there are no IPCs.
#' Only used for xml file formats.
#' @param NC string(s) present in the sample names
#' that represent the negative control wells. 
#' For example, \code{'NC'}. Set to \code{NULL} (default) to use the 
#' type variable from the Barcode B file or if there are no NCs.
#' Only used for xml file formats.
#' @param Bridge string(s) present in the sample names
#' that represent the bridge sample wells. 
#' Set to \code{NULL} (default) to use the 
#' type variable from the Barcode B file or 
#' if there are no bridge samples (default).
#' Only used for xml file formats.
#' @param Calibrator string(s) present in the sample names
#' that represent the calibrator wells. 
#' Set to \code{NULL} (default) to use the 
#' type variable from the Barcode B file or 
#' if there are no calibrator samples (default).
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
                          IPC=NULL, SC=NULL, NC=NULL, Bridge=NULL, Calibrator=NULL,
                          replaceNA=TRUE, nonCogEncrypt=TRUE){
  
  if(!file.exists(file)){
    stop(paste0("Error: The file \'", file, "\' does not exist!\n"))
  }

  if(file_type == 'xml_no_mismatches'){
    
    # Fix rare problem where the XML cannot be loaded because an XML attribute value 
    # has been duplicated. This violates the XML standard and prevents xml2::read_xml
    # from sucessfully loading the file. This finds and removes duplicate attributes
    # read in xml file
    tryCatch({
      if(requireNamespace("NULISAseqAQ", quietly=T)){
        xml <- NULISAseqAQ::readXML_remove_duplicate_attributes(file)
      } else{
        xml <- xml2::read_xml(file)
      }
    },
      error = function(cond){
        stop("Error: Input file is not a valid NAS XML file\n")
      }
    )
    #writeLines(xml, 'temporaryOut.xml')
    #xml<-xml2::read_xml('temporaryOut.xml')
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
    if("Abs" %in% names(ExecutionDetails) & !is.null(ExecutionDetails$Abs)){
      ExecutionDetails$Abs <- data.frame(
        seq = unname(ExecutionDetails$Abs)[seq(1, length(ExecutionDetails$Abs), 2)],
        val = unname(ExecutionDetails$Abs)[seq(2, length(ExecutionDetails$Abs), 2)])
    }
    if("ACC" %in% names(ExecutionDetails)){
      stop("Error: Input file is not a valid NAS XML file!\n")
    }
    
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
      # target Barcode A metadata fill in missing values with NA
      for (i in 1:length(barcodeA_attrs)){
        for (j in 1:length(targetName)){
          value <- attributes(RunSummary$Barcodes$BarcodeA[[j]])[barcodeA_attrs[i]][[1]]
          if(is.null(value)) value <- NA
          targetMetadata[j,i] <- value
        }
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
    num_elements <- length(RunSummary$Barcodes$BarcodeB)
    
    # If AUTO_WELLCOL and AUTO_WELLROW are either null or blank, rewrite them based on AUTO_WELLPOSITION
    # AUTO_WELLCOL and AUTO_WELLROW are required for NAS and are blank in some lot definitions
    RunSummary$Barcodes$BarcodeB <- lapply(RunSummary$Barcodes$BarcodeB, function(barcode){
      if(!is.null(attr(barcode, "AUTO_WELLPOSITION")) && attr(barcode, "AUTO_WELLPOSITION") != "" &&
          (is.null(attr(barcode,"AUTO_WELLCOL")) || is.null(attr(barcode,"AUTO_WELLROW")) ||
          attr(barcode, "AUTO_WELLCOL") == "" || attr(barcode, "AUTO_WELLROW") == "")){
        well_pos <- attr(barcode, "AUTO_WELLPOSITION")
        attr(barcode, "AUTO_WELLROW") <- substr(well_pos, 1, 1)
        attr(barcode, "AUTO_WELLCOL") <- substring(well_pos, 2)
      }
      return(barcode)
    })    

    sorted_index <- seq(num_elements)
    if(!is.null(attr(RunSummary$Barcodes$BarcodeB[[1]], "AUTO_WELLCOL")) && !is.null(attr(RunSummary$Barcodes$BarcodeB[[1]], "AUTO_WELLROW"))){
      # Use lapply to iterate through each element and extract AUTO_WELLCOL
      auto_wellcol_values <- as.numeric(unlist(lapply(seq(num_elements), function(i) {
        attr(RunSummary$Barcodes$BarcodeB[[i]], "AUTO_WELLCOL")
      })))
      auto_wellrow_values <- unlist(lapply(seq(num_elements), function(i) {
        attr(RunSummary$Barcodes$BarcodeB[[i]], "AUTO_WELLROW")
      }))
      
      combined_values <- paste0(auto_wellrow_values, auto_wellcol_values)
      
      custom_sort <- function(vec) {
        # Extract letters and numbers separately
        letters_part <- gsub("\\d", "", vec)
        numbers_part <- as.integer(gsub("\\D", "", vec))
        
        # Sort based on letters_part and numbers_part
        sorted_indices <- order(letters_part, numbers_part)
        return(sorted_indices)
      }
      # Create a sorted index list
      sorted_index <- custom_sort(combined_values)
    }
    sampleName <- unlist(RunSummary$Barcodes$BarcodeB)
    sampleName <- renameDuplicateNames(sampleName, sorted_index) 
    sampleID <- if(!is.null(plateID)) paste0(plateID, ".", sampleName) else sampleName
    # get other sample annotations
    barcodeB_attrs <- unique(
      unname(
        do.call('c',
                lapply(RunSummary$Barcodes$BarcodeB, function(x) {
                  names(attributes(x))
                }))))
    
    barcodeB_attrs <- barcodeB_attrs[barcodeB_attrs!='name']
    if (length(barcodeB_attrs) > 0){
      sampleMetadata <- data.frame(matrix(nrow=length(sampleBarcode),
                                          ncol=length(barcodeB_attrs)))
      colnames(sampleMetadata) <- barcodeB_attrs
      # sample Barcode B metadata fill in missing values with NA
      for (i in 1:length(barcodeB_attrs)){
        for (j in 1:length(sampleName)){
          value <- attributes(RunSummary$Barcodes$BarcodeB[[j]])[barcodeB_attrs[i]][[1]]
          if(is.null(value)) value <- NA
          sampleMetadata[j,i] <- value
        }
      }
      samples <- data.frame(sampleBarcode=sampleBarcode, 
                            sampleName=sampleName, 
                            sampleID=sampleID,
                            sampleMetadata)
    } else if (length(barcodeB_attrs) == 0) {
      samples <- data.frame(sampleBarcode, sampleName)
    }

    # Sort samples dataframe by sorted_index
    # When AUTO_WELLCOL and AUTO_WELLROW are present, the samples dataframe is sorted by the well position
    samples <- samples[sorted_index,]

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


    ##############################
    # QC variable initialization section for SNR
    ##############################
    QCS <- matrix(0, nrow=nrow(targets), ncol=nrow(samples), dimnames=list(targets$targetBarcode, samples$sampleBarcode))
    SN <- matrix(0, nrow=nrow(targets), ncol=nrow(samples), dimnames=list(targets$targetBarcode, samples$sampleBarcode))

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
                     by.x='sampleBarcode', by.y='barcode', all=TRUE, sort=FALSE)
    # sort samples by sample name
    # loop over sample replicates and save data
    Data <- targets[,c('targetBarcode', 'targetName')]
    sampleBarcode <- character(length(SampleData))
    noncog <- vector("list", length(SampleData))
    key <- "awebke3j38djwj22dodl5kkk6" |>charToRaw() |>sodium::sha256() |> cyphr::key_sodium()
    A <- if(nonCogEncrypt) 'A' else 'barcodeA1'
    B <- if(nonCogEncrypt) 'B' else 'barcodeA2'
    for(i in 1:length(SampleData)){
      # get data for sample i
      sampleBarcode[i] <- xml2::xml_attr(SampleData[[i]], 'barcode')
      sampleDataChildren <- xml2::xml_children(SampleData[[i]])

      sampleDataChildrenAttr <- xml2::xml_has_attr(sampleDataChildren, 'target')

      offtarget_ind <- which(sampleDataChildrenAttr != TRUE)
      target_ind <- which(sampleDataChildrenAttr)  
      vals <- xml2::xml_attrs(sampleDataChildren)
      targetBarcode <- sapply(vals[target_ind], '[[', 'target')
      ## Assign the QC SNR vals
      if(length(targetBarcode) > 0){
        q<- as.numeric(xml2::xml_attr(sampleDataChildren, 'QCS', default="0"))
        if(any(q != 0)){
          QCS[targetBarcode, i] <- q
        }
        s <- as.numeric(xml2::xml_attr(sampleDataChildren, 'SN', default="0"))
        if(any(s != 0)){
          SN[targetBarcode, i] <- s
        }
      }
      A1 <- sapply(vals[offtarget_ind], '[[', A)#function(x) x[[A]])
      A2 <- sapply(vals[offtarget_ind], '[[', B)#function(x) x[[B]])
      uA1 <- unique(A1)
      uA2 <- unique(A2)

      mat <- matrix(0, nrow = length(uA1), ncol = length(uA2), dimnames=list(uA1, uA2))

      mat[cbind(match(A1, uA1), 
                match(A2, uA2))] <- if(nonCogEncrypt & requireNamespace("NULISAseqAQ", quietly=T)) as.numeric(NULISAseqAQ::decrypt(xml2::xml_text(sampleDataChildren[offtarget_ind]), key)) else as.numeric(xml2::xml_text(sampleDataChildren[offtarget_ind]))
      noncog[[sampleBarcode[i]]] = mat
      Data_i <- cbind(targetBarcode, 
                      xml2::xml_double(sampleDataChildren[target_ind]))
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
    
    # add sampleType information
    if(!is.null(samples$type) & is.null(IPC)){
      sampleType <- unlist(lapply(samples$type, FUN=function(t) gsub(pattern="sample", replacement="Sample", x=t, fixed=F, ignore.case=TRUE)))
      samples$type <- NULL
    } else { # try to infer sample type from sample names
      sampleType <- rep("Sample", length(samples$sampleName))
      sampleType[grep(paste(NC, collapse="|"), samples$sampleName)]  <- "NC"
      sampleType[grep(paste(SC, collapse="|"), samples$sampleName)] <- "SC"
      sampleType[grep(paste(IPC, collapse="|"), samples$sampleName)] <- "IPC" 
      if(!is.null(Bridge)) sampleType[grep(paste(Bridge, collapse="|"), samples$sampleName)] <- "Bridge"
      if(!is.null(Calibrator)) sampleType[grep(paste(Calibrator, collapse="|"), samples$sampleName)] <- "Calibrator"
    } 
    samples$sampleType <- sampleType
    
    # add IC target information
    val <- if(is.null(IC)) NA else "Target"
    targetType <- rep(val, length(targets$targetName))
    if(!is.null(IC)){ targetType[grep(paste(IC, collapse="|"), targets$targetName)] <- "Control"}
    targets$targetType <- targetType
    
    # save the special well type column names
    specialWellsTargets <- list()
    specialWellsTargets[['IPC']] <- samples$sampleName[samples$sampleType=='IPC']
    specialWellsTargets[['NC']] <- samples$sampleName[samples$sampleType=='NC']
    specialWellsTargets[['SC']] <- samples$sampleName[samples$sampleType=='SC']
    specialWellsTargets[['Bridge']] <- samples$sampleName[samples$sampleType=='Bridge']
    specialWellsTargets[['Calibrator']] <- samples$sampleName[samples$sampleType=='Calibrator']
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

    # Determine if covariates are numeric
    nonControl <- which(toupper(samples$sampleType) != "SAMPLE")
    numericCovariates <- sapply(samples[-nonControl,], function(lst) all(sapply(lst, function(x) is.na(x) || (is.character(x) && x == "NA") || grepl("^\\-?\\d+\\.?\\d*$", x))))
    if ("AUTO_PLATE" %in% names(numericCovariates)){
      numericCovariates['AUTO_PLATE'] = FALSE
    }
    if ("sampleBarcode" %in% names(numericCovariates)){
      numericCovariates['sampleBarcode'] = FALSE
    }
    if ("plateID" %in% names(numericCovariates)){
      numericCovariates['plateID'] = FALSE
    }
    
    # if no sample matrix given assign "PLASMA"
    if(is.null(samples$SAMPLE_MATRIX)){
      if(is.null(samples$Matrix_Type)){
        samples$SAMPLE_MATRIX <- "PLASMA"
      } else{
        samples$SAMPLE_MATRIX <- toupper(samples$Matrix_Type)
        samples$Matrix_Type <- NULL
      }
    }else{
      samples$SAMPLE_MATRIX <- toupper(samples$SAMPLE_MATRIX)
    }
    samples$SAMPLE_MATRIX[which(samples$sampleType != "Sample")] <- "CONTROL"
    samples$SAMPLE_MATRIX[samples$SAMPLE_MATRIX == ""] <- "PLASMA"

    # if no Curve_Quant attribute given, assign "F" for forward curve
    if(is.null(targets$Curve_Quant)){
      targets$Curve_Quant <- "F"
    }
    
    # if no AlamarTargetID given assign it a NA value
    if(is.null(targets$AlamarTargetID)){
      targets$AlamarTargetID <- NA
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
      Data=DataMatrix,
      noncog=noncog,
      numericCovariates=numericCovariates,
      QCS=QCS,
      SN=SN),
      specialWellsTargets))
    # end file type xml_no_mismatches
  } 
  
  if(file_type == 'csv_long' | file_type == 'xlsx'){
    if(file_type == 'csv_long'){
      Data <- read.csv(file)
    } else if(file_type == 'xlsx'){
      Data <- suppressWarnings(data.frame(readxl::read_excel(file, sheet='RQ NPQ Values')))
    } else {
      stop(paste0("Error: Invalid file_type \'", file_type, "\'!\n"))
    }

    # make target and sample data frames
    aM <- "Conc_aM"
    pgmL <- "Conc_pgmL"
    LOD <- NULL
    if(file_type == 'xlsx'){
      sheets <- readxl::excel_sheets(file)
      AQ <- ifelse(any(grepl("^AQ", sheets)), TRUE, FALSE)
      if(AQ){
        aqData <- suppressWarnings(data.frame(readxl::read_excel(file, sheet=sheets[which(grepl("^AQ", sheets))])))
        aqData$SampleName <- renameSC(aqData$SampleName, n=numSCs(aqData, AQ))
        AQdata <- list(Data_AQ=NULL, Data_AQ_aM=NULL)

        # Reshape the data so that they are matrices
        suppressWarnings(AQdata$Data_AQ <- reshape(aqData[,c('SampleName', 'Target', make.names(pgmL))],
                          direction = 'wide',
                          idvar = 'Target',
                          timevar = 'SampleName'))
        suppressWarnings(AQdata$Data_AQ_aM <- reshape(aqData[,c('SampleName', 'Target', make.names(aM))],
                          direction = 'wide',
                          idvar = 'Target',
                          timevar = 'SampleName'))
        
        # Replace "-" with NA
        AQdata$Data_AQ[AQdata$Data_AQ == "-"] <- NA
        AQdata$Data_AQ_aM[AQdata$Data_AQ_aM == "-"] <- NA
        
        rowNames <- AQdata$Data_AQ_aM$Target
        rowNames <- AQdata$Data_AQ$Target
        AQdata$Data_AQ <- as.data.frame(lapply(AQdata$Data_AQ, as.character), stringsAsFactors = FALSE)
        AQdata$Data_AQ_aM <- as.data.frame(lapply(AQdata$Data_AQ_aM, as.character), stringsAsFactors = FALSE)
        
        AQdata$Data_AQ <- suppressWarnings(apply(AQdata$Data_AQ, 2, function(x) as.numeric(x)))
        AQdata$Data_AQ_aM <- suppressWarnings(apply(AQdata$Data_AQ_aM, 2, function(x) as.numeric(x)))
        
        # Write row and column names
        rownames(AQdata$Data_AQ) <- rowNames
        rownames(AQdata$Data_AQ_aM) <- rowNames

        AQdata$Data_AQ <- (AQdata$Data_AQ[,-1])
        AQdata$Data_AQ_aM <- (AQdata$Data_AQ_aM[,-1])
        colnames(AQdata$Data_AQ) <- unique(aqData$SampleName)
        colnames(AQdata$Data_AQ) <- unique(aqData$SampleName)

        LOD <- list(LOD_aM=NULL, LOD_pgmL=NULL)
        lod_aM <- unique(cbind(aqData$Target, aqData[ ,make.names("LoD_aM")]))
        lod_pgmL <- unique(cbind(aqData$Target, aqData[ ,make.names("LoD_pgmL")]))

        LOD$LOD_aM <- suppressWarnings(as.numeric(lod_aM[,2]))
        LOD$LOD_pgmL <- suppressWarnings(as.numeric(lod_pgmL[,2]))

        names(LOD$LOD_aM) <- lod_aM[,1]
        names(LOD$LOD_pgmL) <- lod_pgmL[,1]
      }

      Data$SampleName <- renameSC(Data$SampleName, n=numSCs(Data, AQ))
    }

    Data_colnames <- colnames(Data)
    if(is.null(target_column_names)){
      target_column_names <- c('Target', 'AlamarTargetId','UniProtId', 'ProteinName')
    }
    if(is.null(sample_column_names)){
      sample_column_names <- Data_colnames[!(Data_colnames %in% c(target_column_names, 
                                                                  'Panel', 'PanelLotNumber',
                                                                  'LOD', 'NPQ'))]
    }
    targets <- unique(Data[,target_column_names])
    samples <- unique(Data[,sample_column_names])
    # make an LOD data frame
    LOD$info <- unique(Data[,c('PlateId', 'Target', 'LoD')])
    # reformat the NPQ data
    # reformat into wide, targets in columns
    Data <- reshape(Data[,c('SampleName', 'Target', 'NPQ')],
                    direction= 'wide',
                    idvar = 'Target',
                    timevar='SampleName')
    LOD$LODNPQ <- LOD$info$LoD
    names(LOD$LODNPQ) <- LOD$info$Target
    rownames(Data) <- Data$Target
    Data <- Data[,2:ncol(Data)]
    colnames(Data) <- substr(colnames(Data), start=5, stop=nchar(colnames(Data)))
    Data <- as.matrix(Data)
    # set rownames 
    rownames(targets) <- NULL
    rownames(samples) <- NULL
#    rownames(LOD) <- NULL
    
    # return output
    return(list(
      targets=targets,
      samples=samples,
      lod=LOD,
      Data=Data,
      AQ=AQdata
    ))
  } # end file type csv_long

}

#' Convenience function for use by NAS
#' Read NULISAseq XML, perform normalization, and QC 
#'
#' Reads NULISAseq XML file 
#'
#' @param file Character string. Path and name of the file.
#' @param IC string(s) that represents the names of internal control targets. 
#' Default is \code{'mCherry'}.Only used for xml file formats.
#' @param IPC string(s) present in the sample names
#' that represent the inter-plate control wells. 
#' For example, \code{'IPC'}. Set to \code{NULL} (default) to use the 
#' type variable from the Barcode B file or if there are no IPCs.
#' Only used for xml file formats.
#' @param SC string(s) present in the sample names
#' that represent the sample control wells. 
#' For example, \code{'SC'}. Set to \code{NULL} (default) to use the 
#' type variable from the Barcode B file or if there are no IPCs.
#' Only used for xml file formats.
#' @param NC string(s) present in the sample names
#' that represent the negative control wells. 
#' For example, \code{'NC'}. Set to \code{NULL} (default) to use the 
#' type variable from the Barcode B file or if there are no NCs.
#' Only used for xml file formats.
#' @param Bridge string(s) present in the sample names
#' that represent the bridge sample wells. 
#' Set to \code{NULL} (default) to use the 
#' type variable from the Barcode B file or 
#' if there are no bridge samples (default).
#' Only used for xml file formats.
#' @param sample_group_covar Optional column name in the Barcode B file 
#' and samples data matrix output by readNULISAseq that represents subgroups
#' for which detectability will be calculated separately, in addition to 
#' overall detectability. 
#' @param replace_cal_blank_zeros Logical TRUE / FALSE. Default is FALSE.
#' This parameter is passed to the applyAQ function.
#' If FALSE, the "a" parameter for targets with a zero blank calibrator (NC) mean
#' is set using the blank calibrator mean value (which would equal zero in these cases).
#' IF TRUE, the "a" parameter for targets with a zero blank calibrator mean
#' is set using the nonzero master curve "a" parameter
#' estimate instead.
#' @param replace_zeros_with_NA Logical TRUE / FALSE.
#' This parameter is passed to the applyAQ function.
#' When TRUE (default), any zero values
#' in the AQ data output will be replaced with NA. When FALSE, these values remain as zero.
#'
#' @return List of lists, data frames, and matrices.
#' Output will differ slightly depending on the input file type.
#'
#'
#' @export
#'
loadNULISAseq <- function(file,
                          IC='mCherry',
                          IPC=NULL,
                          SC=NULL,
                          NC=NULL,
                          TAP=TRUE,
                          Bridge=NULL,
                          sample_group_covar=NULL,
                          plateID=NULL,
                          scaleFactor=10^4,
                          transformReverse_scaleFactor=10^4,
                          nonCogEncrypt=FALSE,
                          replace_cal_blank_zeros = FALSE,
                          replace_zeros_with_NA = TRUE,
                          ...){
  raw <- readNULISAseq(file, IPC=IPC, IC=IC, SC=SC, NC=NC, 
                       nonCogEncrypt=nonCogEncrypt,...)
  raw$IC_normed <- intraPlateNorm(data_matrix=raw$Data, IC=IC[1])
  ind <-  which(raw$targets$Curve_Quant == "R")
  reverseCurve <- if (length(ind) > 0) raw$targets$targetName[ind] else NULL
  raw$normed_untransformedReverse <- interPlateNorm(data_list=list(raw$IC_normed$normData),
                                                    IPC_wells=list(raw$IPC),
                                                    scaleFactor=scaleFactor,
                                                    transformReverse_scaleFactor=transformReverse_scaleFactor)
  raw$normed <- interPlateNorm(list(raw$IC_normed$normData), 
                               transformReverse=reverseCurve, 
                               IPC_wells=list(raw$IPC),
                               scaleFactor=scaleFactor,
                               transformReverse_scaleFactor=transformReverse_scaleFactor)
  
  # if Execution Details has an Abs (Absolute quantification) section, add AQ results
  AbsAssay <- "Abs" %in% names(raw$ExecutionDetails) & !is.null(raw$ExecutionDetails$Abs)
  if (AbsAssay){
    rawNormed <- interPlateNorm(list(raw$IC_normed$normData),
                                transformReverse=reverseCurve,
                                IPC_wells=list(raw$IPC),
                                scaleFactor=scaleFactor,
                                transformReverse_scaleFactor=transformReverse_scaleFactor)
    if(requireNamespace("NULISAseqAQ", quietly=T)){
      tryCatch(
        {
          raw$AQ <- NULISAseqAQ::applyAQ(
              normDataIPC=rawNormed$interNormData[[1]],
              params=NULL, #raw$ExecutionDetails$Abs, 
              file=file,
              targetDataFrame=raw$targets, 
              IPC=raw$IPC, 
              NC=raw$NC,
              replace_cal_blank_zeros = replace_cal_blank_zeros,
              replace_zeros_with_NA = replace_zeros_with_NA
          )
        },
        error = function(cond){
          message("Could not perform absolute quantification\n")
          stop(conditionMessage(cond))
        }
      )
    }
  }
  
  # set assay name if not exists
  if(!AbsAssay & is.null(raw$ExecutionDetails$Assay)){
    raw$ExecutionDetails$Assay <- "NULISAseq"
  }

   
  # calculate NPQ from IPC-normalized data
  raw$NPQ <- log2(raw$normed$interNormData[[1]] + 1)
  
  # Determine which targets should NOT have outlier detection performed
  # Currently we require that at least one NC sample have at least 100 
  # raw reads for that target for outlier detection to be applied 
  # Therefore, if all NC samples have < 100 raw reads, no outlier detection
  # should be performed
  targetNoOutlierDetection = names(which(apply(raw$Data[, raw$NC] < 100, 1, all)))

  # calculate LODs on IPC-normalized data
  # Note that the reverse-curve is already applied since this went through interPlateNorm
  # This means that aboveLOD is wrong, so we need to recalculate this using normed_unstranformedReverse$interNormData[[1]]
  raw$lod <- lod(data_matrix=raw$normed$interNormData[[1]], blanks=raw$NC, min_count=0, targetNoOutlierDetection=targetNoOutlierDetection)
  
  # Recalculate lod without performing the reverse curve calculation we did in interPlateNorm so we can get the aboveLOD correct
  lodTemp <- lod(data_matrix=raw$normed_untransformedReverse$interNormData[[1]], blanks=raw$NC, min_count=0, targetNoOutlierDetection=targetNoOutlierDetection)

  # Create LOD values on NPQ data. Need to perform IPC normalization on IC-normalized 
  # LOD vals so that we can report these to users
  raw$lod$LODNPQ <- log2(raw$lod$LOD + 1)

  ## manually apply reverse curve correction since we didn't go through the interPlateNorm
  if(length(reverseCurve) > 0){
    raw$lod$LODNPQ[reverseCurve] <- log2((transformReverse_scaleFactor * scaleFactor) / (lodTemp$LOD[reverseCurve] + 1) + 1) # transformReverse_scaleFactor * scaleFactor 
    ## replace the reverse curve aboveLOD values with ones calculated from untransformed Reverse
    raw$lod$aboveLOD[reverseCurve, ] <- lodTemp$aboveLOD[reverseCurve, ]
  }

  # Perform QC at the target, plate and sample levels, also calculate the LODs in AQ units
  if(AbsAssay){
    if(requireNamespace("NULISAseqAQ", quietly=T)){ 
      raw$qcTarget <- QCFlagTarget(raw$AQ$Data_AQ_aM, raw$targets, raw$samples, raw$AQ$targetAQ_param$SC_conc)
      tryCatch(
        {
          s <- cbind(raw$lod$LOD, rawNormed$interNormData[[1]][, raw$IPC], rawNormed$interNormData[[1]][, raw$NC])
          rownames(s) <- names(raw$lod$LOD)
          temp <- NULISAseqAQ::applyAQ(
              normDataIPC=s,
              params=NULL, #raw$ExecutionDetails$Abs, 
              file=file,
              targetDataFrame=raw$targets, 
              IPC=raw$IPC, 
              NC=raw$NC,
              replace_cal_blank_zeros = replace_cal_blank_zeros,
              replace_zeros_with_NA = replace_zeros_with_NA
          )
          raw$lod$LOD_pgmL <- temp$Data_AQ[, 1]
          raw$lod$LOD_aM <- temp$Data_AQ_aM[, 1]
          ind <-names(which(s[,1] == 0));
          if(length(ind) > 0){
            raw$lod$LOD_aM[ind] <- NA
            raw$lod$LOD_pgmL[ind] <- NA
          }
        },
        error = function(cond){
          message("Could not perform absolute quantification\n")
          stop(conditionMessage(cond))
        }
      )
    }
  }
  raw$qcSample <- QCFlagSample(raw$Data, raw$lod$aboveLOD, raw$samples, raw$targets, QCS=raw$QCS, SN=raw$SN, TAP=TAP)
  raw$qcPlate <- QCFlagPlate(raw$Data, raw$IC_normed$normData, raw$lod$aboveLOD, raw$targets, raw$samples, AQ=AbsAssay, AQ_QC=raw$qcTarget, Sample_QC=raw$qcSample)
  # calculate Detectability
  sample_groups <- NULL
  if(!is.null(sample_group_covar)) sample_groups <- raw$samples[raw$samples$sampleType == "Sample", sample_group_covar]
  raw$detectability <- detectability(aboveLOD_matrix=raw$lod$aboveLOD, 
                                     sample_subset=raw$samples$sampleName[which(raw$samples$sampleType == "Sample")], 
                                     sample_groups=sample_groups, 
                                     exclude_targets=raw$IC)

  # do not return an NPQ value if the sample_matrix is "other" and target is reverse curve
  noQuantElems <- c("other")
  noQuantify <- which(toupper(raw$samples$SAMPLE_MATRIX) %in% toupper(noQuantElems))
  ind <-  which(raw$targets$Curve_Quant == "R")
  if(length(noQuantify) > 0 & length(ind) > 0){
    raw$normed$log2_interNormData[[1]][ind, noQuantify] <- NA
    raw$NPQ[ind, noQuantify] <- NA
  }
  return(raw)
}

#' Read Covariate file and add to NULISAseq object
#'
#' Reads NULISAseq covariate file and adds new covariates to list object. Merge based on plateID and sampleName. 
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
  newInfo <- read.csv(txt_file, header=T)
  for(i in 1:length(NULISAseqRuns)){
    NULISAseqRuns[[i]]$samples <- merge(NULISAseqRuns[[i]]$samples, newInfo, nodups=T, all.x=T,
                                        by.x=c("plateID", "sampleName"), 
                                        by.y=c("plateID", "sampleName"))
    # Determine if covariates are numeric
    NULISAseqRuns[[i]]$numericCovariates <- sapply(NULISAseqRuns[[i]]$samples, function(lst) all(sapply(na.omit(lst), function(x) suppressWarnings(!is.na(as.numeric(x))))))
  }
  return(NULISAseqRuns)
}
