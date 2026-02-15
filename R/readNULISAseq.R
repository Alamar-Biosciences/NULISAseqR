#' Identify Missing Target Entries in Sample Curve Quant
#'
#' This function scans a vector of semicolon-delimited target strings (`curve_quant`)
#' and identifies which entries contain missing targets, indicated by the presence
#' of `"<sample_matrix>-NA"` for each target listed in `sample_matrices`.
#'
#' @param curve_quant A character vector where each element contains semicolon-separated
#'   target identifiers, e.g., `"R;CSF-NA"`. The "R" identifier indicates a reverse
#'   curve target.
#' @param sample_matrices A character vector of sample_matrix values
#'
#' @return A two-column matrix where each row gives the (row, col) index of a match:
#'   - `row`: the index in `curve_quant`
#'   - `col`: the index in `sample_matrices` where a corresponding `"-NA"` entry was found.
#'   If no matches are found, returns an empty matrix with 0 rows and 2 columns.
#'   
#' @keywords internal
#'
#' @details
#' Reverse curve targets (indicated by "R" in `curve_quant`) receive special handling:
#' they are marked as NA for sample matrices other than "PLASMA", "SERUM", or "CONTROL".
#' This ensures that NPQ values are not reported for reverse curve targets in
#' incompatible sample matrices.
#'
#' @examples
#' curve_quant <- c("R;CSF-NA", "F", "R;PLASMA-NA")
#' sample_matrices <- c("CSF", "PLASMA", "SERUM")
#' calcSampleTargetNAs(curve_quant, sample_matrices)
#'
#' @export
calcSampleTargetNAs <- function(curve_quant, sample_matrices){
  matches <- list()

  # do not return an NPQ value for reverse curve targets unless sample_matrix is "PLASMA" or "SERUM" or "CONTROL" 
  for (i in seq_along(curve_quant)) {
    parts <- unlist(strsplit(toupper(curve_quant[i]), ";"))
    for (j in seq_along(sample_matrices)) {
      target <- paste0(toupper(sample_matrices[j]), "-NA")
      if (target %in% parts || ("R" %in% parts && !sample_matrices[j] %in% c("PLASMA", "SERUM", "CONTROL"))) {
        matches[[length(matches) + 1]] <- c(row = i, col = j)
      }
    }
  }
  # Combine into a matrix of indices
  if(length(matches) == 0) {
    # Return empty matrix with correct structure instead of NULL
    match_matrix <- matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("row", "col")))
  } else {
    match_matrix <- do.call(rbind, matches)
    colnames(match_matrix) <- c("row", "col")
  }
  return(match_matrix)
}

#' Parse QC XML Nodes into a Data Frame
#'
#' This function processes XML nodes representing QC data (from either `<PlateQC>` or `<SampleQC>`) and
#' converts them into a tidy data frame. It handles both flat `<QCFlag>` nodes and nested `<Sample>` nodes containing `<QCFlag>` children.
#'
#' @param nodes A list of XML nodes, typically retrieved with `xml2::xml_find_all()`, representing either `<QCFlag>` elements directly
#'              or `<Sample>` elements that contain `<QCFlag>` children.
#' @param rename Logical; if `TRUE`, column names will be renamed to more descriptive names (e.g., `name` → `flagName`, `set` → `status`, `value` → `val`).
#'
#' @return A data frame with one row per QC flag. For nested `<Sample>` nodes, a `sample` column will be included. Returns `NULL` if `nodes` is `NULL`.
#'
#' @details
#' - If `nodes` contains only `<QCFlag>` elements (as under `<PlateQC>`), the output will contain attributes and values for each flag.
#' - If `nodes` contains `<Sample>` elements (as under `<SampleQC>`), each `<QCFlag>` inside will be parsed and include an additional `sample` column.
#' - The `rename` parameter allows renaming of columns using standard recoding.
#'
#' @examples
#' # For Plate-level QC:
#' plate_nodes <- xml2::xml_find_all(xml_doc, ".//PlateQC/QCFlag")
#' readQCXMLNode(plate_nodes)
#'
#' # For Sample-level QC:
#' sample_nodes <- xml2::xml_find_all(xml_doc, ".//SampleQC/Sample")
#' readQCXMLNode(sample_nodes)
#'
#' @importFrom xml2 xml_name xml_attrs xml_attr xml_text xml_find_all
#' @importFrom dplyr bind_rows rename_with recode
#' @export

readQCXMLNode <- function(nodes, rename = TRUE, tag="QCFlag") {
  rows <- list()
  if(is.null(nodes) || length(nodes) == 0){
    return (NULL);
  }
  
  for (node in nodes) {
    node_name <- xml2::xml_name(node)
    
    if (node_name == "QCFlag") {
      # Case 1: Plate-level QCFlag
      attrs <- xml2::xml_attrs(node)
      name <- attrs[["name"]]
      
      row <- as.list(attrs)
      row$value <- xml2::xml_text(node)
      
      rows[[name]] <- row
      
    } else {
      # Case 2: Sample-level <Sample> or Target-level <Target> containing multiple <QCFlag>
      sample_name <- xml2::xml_attr(node, "name")
      qc_flags <- xml2::xml_find_all(node, paste0(".//", tag))
      
      for (qc_node in qc_flags) {
        attrs <- xml2::xml_attrs(qc_node)
        row <- as.list(attrs)
        if(node_name == "Target"){
          row$target <- sample_name
        } else if (node_name == "Sample"){
          row$sample <- sample_name
          row$value <- xml2::xml_text(qc_node)
        } 
        rows[[length(rows)+1]]<-row
      }
    }
  }
  
  df <- dplyr::bind_rows(rows)
  
  if (rename) {
    df <- dplyr::rename_with(df, ~ dplyr::recode(.,
                                                 aq_aM = paste0(tag, "_aM"),
                                                 aq_pgmL = paste0(tag, "_pg_ml"),
                                                 name = "flagName",
                                                 set = "status",
                                                 value = "val",
                                                 method = "normMethod",
                                                 format = "QCformat",
                                                 explanation = "explanations"
    ))
  }
  
  return(as.data.frame(df))
}

#' Read and Parse QC Attributes from XML Nodes
#'
#' Parses a list of XML `<Threshold>` nodes and extracts their attributes into
#' named vectors grouped by attribute. The output is a named list of named vectors,
#' where each vector corresponds to a different attribute and is indexed by the
#' `name` attribute from each node.
#'
#' If the XML node has a text value (e.g. `<Threshold>0.3</Threshold>`), it will
#' override the "value" attribute for that node.
#'
#' Optionally, the list element names can be renamed to more descriptive keys such as
#' `thresholds`, `operators`, `properNames`, and `explanations`.
#'
#' @param nodes An `xml_nodeset` or list of XML nodes (typically `<Threshold>` nodes).
#' @param rename Logical (default `TRUE`). If `TRUE`, renames selected attribute names:
#'   - `"value"` → `"thresholds"`
#'   - `"operator"` → `"operators"`
#'   - `"properName"` → `"properNames"`
#'   - `"explanation"` → `"explanations"`
#'
#' @return A named list where each element is a named character vector, with threshold
#'   names as names and attribute values as values.
#'
#' @examples
#' # xml <- xml2::read_xml("your_file.xml")
#' # nodes <- xml2::xml_find_all(xml, ".//Threshold")
#' # result <- readQCThresholdXMLNode(nodes)
#'
#' @export
readQCThresholdXMLNode <- function(nodes, rename=TRUE){
  if(is.null(nodes) || length(nodes) == 0){
    return(NULL)
  }
  threshs <- list()
  for (node in nodes){
    attrs <- xml2::xml_attrs(node)
    name <- attrs[["name"]]
    for (attr_name in names(attrs)){
      if(is.null(threshs[[attr_name]])){
        threshs[[attr_name]] <- c()
      }
      threshs[[attr_name]][name] <- attrs[[attr_name]]
    }
    if(xml2::xml_text(node) != ""){
      threshs[["value"]][name] <- xml2::xml_text(node)
    }
  }
  if(rename){
    names(threshs) <- dplyr::recode(names(threshs), value="thresholds", operator="operators", properName="properNames", explanation="explanations")
  }
  return(threshs)
}    

#' Read and Parse XML, Removing Duplicate Attributes
#'
#' This function reads an XML file as plain text, removes duplicate attributes within XML tags,
#' and parses the cleaned XML using the `xml2::read_xml` function. For each tag, only the first 
#' occurrence of each attribute is retained, and any subsequent duplicate attributes are removed.
#' The cleaned XML is then returned as an `xml_document` object.
#'
#' @param xml_file A character string representing the path to the XML file to be read and cleaned.
#'
#' @return An `xml_document` object representing the cleaned XML with duplicate attributes removed.
#'
#' @details
#' The function processes the XML file in three main steps:
#' 1. Reads the XML file as plain text and combines it into a single string.
#' 2. For each XML tag, the function identifies and removes any duplicate attributes, retaining 
#'    only the first occurrence of each attribute.
#' 3. The cleaned XML string is parsed and returned as an `xml_document` object using `xml2::read_xml`.
#'
#' The function is designed to handle the violation of XML specifications where tags may contain
#' duplicated attribute names.
#'
#' @examples
#' # Example usage:
#' # xml_doc <- readXML_remove_duplicate_attributes("path_to_file.xml")
#'
#' @import xml2
#' @export
readXML_remove_duplicate_attributes <- function(xml_file) {
  # Try reading the file using traditional method
  tryCatch({
    cleaned_xml_doc <- xml2::read_xml(xml_file)
    return(cleaned_xml_doc)
  },
  error = function(cond){
    # Read the XML file as plain text
    xml_data <- readLines(xml_file, warn = FALSE)
    
    # Concatenate lines into a single string
    xml_text <- paste(xml_data, collapse = "\n")
    
    # Function to clean a single tag's attributes
    clean_tag <- function(tag_text) {
      # Extract the tag name and attributes
      tag_name <- sub("^<([a-zA-Z0-9]+).*", "\\1", tag_text)
      attr_string <- sub("^<[a-zA-Z0-9]+\\s([^>]+)>", "\\1", tag_text)
      
      # Split the attribute string into individual attributes
      attrs <- unlist(strsplit(attr_string, "\\s+"))
      
      # Set to track seen attribute names
      seen_attrs <- c()
      cleaned_attrs <- c()
      
      # Loop through the attributes and remove duplicates
      for (attr in attrs) {
        # Extract the attribute name (before the "=" sign)
        attr_name <- sub("=.*", "", attr)
        
        # If the attribute hasn't been seen, keep it
        if (!(attr_name %in% seen_attrs)) {
          seen_attrs <- c(seen_attrs, attr_name)
          cleaned_attrs <- c(cleaned_attrs, attr)
        }
      }
      
      # Rebuild the cleaned tag with unique attributes
      cleaned_tag <- paste0("<", tag_name, " ", paste(cleaned_attrs, collapse = " "), ">")
      return(cleaned_tag)
    }
    
    # Use regular expressions to match tags and clean attributes
    tag_pattern <- "<[a-zA-Z0-9]+\\s[^>]+>"
    
    # Find all matching tags in the XML text
    matches <- gregexpr(tag_pattern, xml_text)
    matched_tags <- regmatches(xml_text, matches)[[1]]
    
    # Clean each tag by removing duplicate attributes
    cleaned_tags <- sapply(matched_tags, clean_tag, USE.NAMES = FALSE)
    
    # Replace the original tags with the cleaned versions
    cleaned_xml_text <- xml_text
    regmatches(cleaned_xml_text, matches) <- list(cleaned_tags)
    
    tempStore <- paste(cleaned_xml_text, collapse="\n")
    
    # Now that attributes are cleaned, use xml2::read_xml to parse the cleaned XML
    cleaned_xml_doc <- read_xml(tempStore)
    
    return(cleaned_xml_doc)
  })
}

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
#' @keywords internal
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
#' @keywords internal
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
#' Options include \code{xml_no_mismatches}, \code{xlsx} (output from ACC)
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
                          IC=NULL, 
                          IPC=NULL, SC=NULL, NC=NULL, Bridge=NULL, Calibrator=NULL,
                          replaceNA=TRUE){
  
  if(!file.exists(file)){
    stop(paste0("Error: The file \'", file, "\' does not exist!\n"))
  }
  
  if(file_type == 'xml_no_mismatches'){
    
    # Fix rare problem where the XML cannot be loaded because an XML attribute value 
    # has been duplicated. This violates the XML standard and prevents xml2::read_xml
    # from sucessfully loading the file. This finds and removes duplicate attributes
    # read in xml file
    tryCatch({
      xml <- readXML_remove_duplicate_attributes(file)
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
    ExecutionDetails$Assay <- trimws(ExecutionDetails$Assay)
    ExecutionDetails$Assay <- if(length(ExecutionDetails$Assay) == 0) NULL else ExecutionDetails$Assay
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
    # Force targetName to character to prevent downstream join type mismatches
    # XML parsing can preserve numeric types which causes errors when merging plates
    targetName <- as.character(unlist(RunSummary$Barcodes$BarcodeA))
    # get other target annotations
    # Collect all unique attribute names across ALL BarcodeA elements, not just the first
    # This ensures attributes like 'modifiers' are captured even if only some targets have them
    barcodeA_attrs <- unique(unlist(lapply(RunSummary$Barcodes$BarcodeA, function(x) names(attributes(x)))))
    barcodeA_attrs <- barcodeA_attrs[barcodeA_attrs != 'name']
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
    barcodeIndex <- unlist(lapply(seq(num_elements), function(i) {
      attr(RunSummary$Barcodes$BarcodeB[[i]], "name")
    }))
    sorted_indexRename <- sort(barcodeIndex, index.return=T)
    sampleName <- unlist(RunSummary$Barcodes$BarcodeB)
    sampleName <- renameDuplicateNames(sampleName, sorted_indexRename$ix) 
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
          sampleMetadata[j, i] <- value
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
    # Initialize attributes data structure
    ##############################
    # Get all possible attributes from first sample's data children
    SampleData <- xml2::xml_find_all(xml, './/Data//Sample')
    all_attrs <- unique(unlist(lapply(SampleData, function(sampleNode) {
      sampleDataChildren <- xml2::xml_children(sampleNode)
      # Filter out OffTarget nodes
      node_names <- xml2::xml_name(sampleDataChildren)
      keep_ind <- which(node_names != "OffTarget")
      sampleDataChildren <- sampleDataChildren[keep_ind]
      unlist(lapply(xml2::xml_attrs(sampleDataChildren), names))
    })))
    all_attrs <- setdiff(all_attrs, c('target')) # Remove 'target' as it's handled separately
    
    # Initialize list to store attribute matrices
    attrib <- list()
    for(at in all_attrs) {
      attrib[[at]] <- matrix(0, 
                               nrow=nrow(targets), 
                               ncol=nrow(samples), 
                               dimnames=list(targets$targetName, samples$sampleName))
    }    # save matching / non-matching
    
    ###########################
    # Save the QC sections if they exist
    ###########################
    qcXML <- NULL
    qcXML$PlateThresh  <- readQCThresholdXMLNode(xml2::xml_find_all(xml, './/QCThresholds/Plate/Threshold'))
    qcXML$TargetThresh <- readQCThresholdXMLNode(xml2::xml_find_all(xml, './/QCThresholds/Target/Threshold'))
    qcXML$SampleThresh <- readQCThresholdXMLNode(xml2::xml_find_all(xml, './/QCThresholds/Sample/Threshold'))
    
    qcXML$qcPlate  <- readQCXMLNode(xml2::xml_find_all(xml, './/PlateQC/QCFlag'), tag="QCFlag")
    qcXML$qcTarget <- readQCXMLNode(xml2::xml_find_all(xml, './/TargetQC/Target'), tag="QCFlag")
    qcXML$qcSample <- readQCXMLNode(xml2::xml_find_all(xml, './/SampleQC/Sample'), tag="QCFlag")
    
    if(!is.null(qcXML$qcPlate)){
      qcXML$qcPlate$row_id <- NULL # remove column that doesn't exist in qcPlate
      qcXML$qcPlate$QCoperator   <- qcXML$PlateThresh$operators[qcXML$qcPlate$flagName]
      qcXML$qcPlate$QCthreshold  <- qcXML$PlateThresh$thresholds[qcXML$qcPlate$flagName]
    }
    if(!is.null(qcXML$qcSample)){
      qcXML$qcSample <- dplyr::rename(qcXML$qcSample, sampleBarcode=dplyr::any_of("sample"))
      qcXML$qcSample$QCoperator  <- qcXML$SampleThresh$operators[qcXML$qcSample$flagName]
      qcXML$qcSample$QCthreshold <- qcXML$SampleThresh$thresholds[qcXML$qcSample$flagName]
      qcXML$qcSample$sampleName <- dplyr::left_join(data.frame(sampleBarcode=qcXML$qcSample$sampleBarcode), samples, by="sampleBarcode")$sampleName
    }
    
    if(!is.null(qcXML$qcTarget) && nrow(qcXML$qcTarget) > 0){
      qcXML$qcTarget <- dplyr::rename(qcXML$qcTarget, targetBarcode=dplyr::any_of("sample"))
      qcXML$qcTarget$QCoperator  <- qcXML$TargetThresh$operators[qcXML$qcTarget$flagName]
      qcXML$qcTarget$QCthreshold <- qcXML$TargetThresh$thresholds[qcXML$qcTarget$flagName]
      qcXML$qcTarget$target     <- dplyr::left_join(data.frame(target=qcXML$qcTarget$target), targets, by=c("target"="targetBarcode"))$targetName
      qcXML$target <- list()
      qcXML$target$lod  <- readQCXMLNode(xml2::xml_find_all(xml, './/TargetQC/Target'), tag="LOD")
      qcXML$target$ULOQ <- readQCXMLNode(xml2::xml_find_all(xml, './/TargetQC/Target'), tag="ULOQ")
      qcXML$target$LLOQ <- readQCXMLNode(xml2::xml_find_all(xml, './/TargetQC/Target'), tag="LLOQ")
      qcXML$target$SC_conc <- readQCXMLNode(xml2::xml_find_all(xml, './/TargetQC/Target'), tag="SC_conc")
      if(!is.null(qcXML$target$lod) && nrow(qcXML$target$lod) > 0){
        qcXML$target$lod  <- dplyr::left_join(data.frame(targets), qcXML$target$lod, by=c("targetBarcode"="target"))
      }
      if(!is.null(qcXML$target$LLOQ) && nrow(qcXML$target$LLOQ) > 0){
        qcXML$target$LLOQ <- dplyr::left_join(data.frame(targets), qcXML$target$LLOQ, by=c("targetBarcode"="target"))
      }
      if(!is.null(qcXML$target$ULOQ) && nrow(qcXML$target$ULOQ) > 0){
        qcXML$target$ULOQ <- dplyr::left_join(data.frame(targets), qcXML$target$ULOQ, by=c("targetBarcode"="target"))
      }
      if(!is.null(qcXML$target$SC_conc) && nrow(qcXML$target$SC_conc) > 0){
        qcXML$target$SC_conc <- dplyr::left_join(data.frame(targets), qcXML$target$SC_conc, by=c("targetBarcode"="target"))
      }
    }
    ###########################
    # save Data section
    ###########################
    
    QCS <- matrix(0, nrow=nrow(targets), ncol=nrow(samples), dimnames=list(targets$targetBarcode, samples$sampleBarcode))
    SN <- matrix(0, nrow=nrow(targets), ncol=nrow(samples), dimnames=list(targets$targetBarcode, samples$sampleBarcode))
    
    ###########################
    # save Data section
    ###########################
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
    A <- 'barcodeA1'
    B <- 'barcodeA2'
    for(i in 1:length(SampleData)){
      # get data for sample i
      sampleBarcode[i] <- xml2::xml_attr(SampleData[[i]], 'barcode')
      sampleDataChildren <- xml2::xml_children(SampleData[[i]])

      # Filter out OffTarget nodes - they should be ignored
      node_names <- xml2::xml_name(sampleDataChildren)
      keep_ind <- which(node_names != "OffTarget")
      sampleDataChildren <- sampleDataChildren[keep_ind]

      sampleDataChildrenAttr <- xml2::xml_has_attr(sampleDataChildren, 'target')

      offtarget_ind <- which(sampleDataChildrenAttr != TRUE)
      target_ind <- which(sampleDataChildrenAttr)  
      vals <- xml2::xml_attrs(sampleDataChildren)
      targetBarcode <- sapply(vals[target_ind], '[[', 'target')
      
      ## Assign all available attributes
      if(length(targetBarcode) > 0){
        for(at in all_attrs) {
          # XML may contain "nan" as string - suppress coercion warnings as NAs are expected
          attr_vals <- suppressWarnings(as.numeric(xml2::xml_attr(sampleDataChildren, at, default="0")))
          # Get target names corresponding to target barcodes
          targetNames <- targets$targetName[match(targetBarcode, targets$targetBarcode)]
          # Get sample name corresponding to current barcode
          sampleName <- samples$sampleName[match(sampleBarcode[i], samples$sampleBarcode)]
          attrib[[at]][targetNames, sampleName] <- attr_vals
        }
      }
      
      A1 <- sapply(vals[offtarget_ind], '[[', A)
      A2 <- sapply(vals[offtarget_ind], '[[', B)
      uA1 <- unique(A1)
      uA2 <- unique(A2)
      
      mat <- matrix(0, nrow = length(uA1), ncol = length(uA2), dimnames=list(uA1, uA2))
      
      mat[cbind(match(A1, uA1), 
                match(A2, uA2))] <- as.numeric(xml2::xml_text(sampleDataChildren[offtarget_ind]))
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
    
    specialWellsTargets <- list()
    
    if(!is.null(IC)){
      if(!all(IC %in% targets$targetName)){
        missing_IC <- IC[!IC %in% targets$targetName]
        stop("Error: The IC target name(s) provided do not match any target names in the XML file: ",
             paste(missing_IC, collapse = ", "), "\n")
      }
    }
    
    # add IC target information
    if(is.null(targets$type)){ # type was not defined in the barcodeA (older XMLs)
      targets$targetType <- rep("target", length(targets$targetName))
      IC_inUse <- ifelse(is.null(IC), 'mCherry', IC)
      targets$targetType[grep(paste(IC_inUse, collapse="|"), targets$targetName)] <- "control"
      targets$modifiers <- rep(NA, length(targets$targetName))
      targets$hide <- rep(FALSE, length(targets$targetType))
      targets$noDetectability <- rep(FALSE, length(targets$targetType))
      specialWellsTargets[['IC']] <- IC_inUse
      IC <- IC_inUse
    } else{ # type was defined for targets / barcodeA in the XML
      targets$targetType <- targets$type
      if(!is.null(IC)){ # Override what was in the XML 
        targets$targetType[grep(paste(IC, collapse="|"), targets$targetName)] <- "control"
      } 
      if(is.null(targets$modifiers)){ # its possible that type was defined, but modifiers was not
        targets$modifiers <- rep(NA, length(targets$targetName))
        targets$hide <- rep(FALSE, length(targets$targetName))
        targets$noDetectability <- rep(FALSE, length(targets$targetName))
      } else{
        targets$modifiers <- sapply(strsplit(as.character(targets$modifiers), ";"), `[`, 1)
        targets$hide <- grepl("hide", targets$modifiers)
        targets$noDetectability <- grepl("noDetectability", targets$modifiers)
      }
      IC <- targets$targetName[which(targets$targetType == "control" & targets$hide == FALSE)]
      specialWellsTargets[['IC']] <- targets$targetName[which(tolower(targets$targetType) == "control" & targets$hide == FALSE)]
    }
    
    # save the special well type column names
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
    
    # if given, add plateID to samples 
    if( !is.null(plateID) ){
      samples$plateID <- plateID
    } else{
      if ( !is.null(samples$AUTO_PLATE) ){
        plateID <- unique(samples$AUTO_PLATE)[1]
        samples$plateID <- samples$AUTO_PLATE
      } else {
        plateID <- "Plate_01"
        samples$plateID <- plateID
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
    inds <- inds[!names(inds) %in% "modifiers"]
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
      numericCovariates=numericCovariates,
      attributes=attrib,
      IC=IC,
      qcXML = qcXML),
      specialWellsTargets))
    # end file type xml_no_mismatches
  } 
  
  if(file_type == 'csv_long' | file_type == 'xlsx'){
    if(file_type == 'csv_long'){
      Data <- read.csv(file)
      AQdata <- NULL
    } else if(file_type == 'xlsx'){
      sheets <- readxl::excel_sheets(file)
      sheet_name <- sheets[grepl('^RQ NPQ values$', sheets, ignore.case = TRUE)]
      Data <- data.frame(readxl::read_excel(file, na="-", sheet=sheet_name))
    } else {
      stop(paste0("Error: Invalid file_type \'", file_type, "\'!\n"))
    }
    
    # make target and sample data frames
    if(file_type == 'xlsx'){
      AQ <- ifelse(any(grepl("^AQ", sheets)), TRUE, FALSE)
      AQdata <- NULL
      if(AQ){
        aqData <- data.frame(readxl::read_excel(file, na="-", sheet=sheets[which(grepl("^AQ", sheets))]))
        aqData$SampleName <- renameSC(aqData$SampleName, n=numSCs(aqData, AQ))
        AQdata <- list(Data_AQ=NULL, Data_AQ_aM=NULL)
        
        aM <- c("Conc_aM", "Conc (aM)")
        pgmL <- c("Conc_pgmL", "Conc (pg/mL)")
        LOD <- NULL
        aM_reshape <- intersect(c('SampleName', 'Target', make.names(aM)), colnames(aqData))
        pgmL_reshape <- intersect(c('SampleName', 'Target', make.names(pgmL)), colnames(aqData))
        # Reshape the data so that they are matrices
        suppressWarnings(AQdata$Data_AQ <- reshape(aqData[, pgmL_reshape],
                                                   direction = 'wide',
                                                   idvar = 'Target',
                                                   timevar = 'SampleName'))
        suppressWarnings(AQdata$Data_AQ_aM <- reshape(aqData[, aM_reshape],
                                                      direction = 'wide',
                                                      idvar = 'Target',
                                                      timevar = 'SampleName'))
        
        
        # Replace "-" with NA
        AQdata$Data_AQ[AQdata$Data_AQ == "-"] <- NA
        AQdata$Data_AQ_aM[AQdata$Data_AQ_aM == "-"] <- NA
        rowNames_AQ_aM <- AQdata$Data_AQ_aM$Target
        rowNames_AQ <-  AQdata$Data_AQ$Target
        AQdata$Data_AQ <- as.data.frame(lapply(AQdata$Data_AQ, as.character), stringsAsFactors = FALSE)
        AQdata$Data_AQ_aM <- as.data.frame(lapply(AQdata$Data_AQ_aM, as.character), stringsAsFactors = FALSE)
        
        AQdata$Data_AQ <- suppressWarnings(apply(AQdata$Data_AQ, 2, function(x) as.numeric(x)))
        AQdata$Data_AQ_aM <- suppressWarnings(apply(AQdata$Data_AQ_aM, 2, function(x) as.numeric(x)))
        
        # Write row and column names
        rownames(AQdata$Data_AQ) <- rowNames_AQ
        rownames(AQdata$Data_AQ_aM) <- rowNames_AQ_aM
        
        AQdata$Data_AQ <- (AQdata$Data_AQ[,-1])
        AQdata$Data_AQ_aM <- (AQdata$Data_AQ_aM[,-1])
        colnames(AQdata$Data_AQ) <- unique(aqData$SampleName)
        colnames(AQdata$Data_AQ_aM) <- unique(aqData$SampleName)
        
        LOD <- list(LOD_aM=NULL, LOD_pgmL=NULL)
        lod_aMnames <- intersect(c("LoD_aM", make.names("LOD (aM)")), colnames(aqData))
        lod_pgmLnames <- intersect(c("LoD_pgmL", make.names("LOD (pg/mL)")), colnames(aqData))
        
        lod_aM <- unique(cbind(aqData$Target, aqData[ , lod_aMnames]))
        lod_pgmL <- unique(cbind(aqData$Target, aqData[ , lod_pgmLnames]))
        
        LOD$LOD_aM <- suppressWarnings(as.numeric(lod_aM[,2]))
        LOD$LOD_pgmL <- suppressWarnings(as.numeric(lod_pgmL[,2]))
        
        names(LOD$LOD_aM) <- lod_aM[,1]
        names(LOD$LOD_pgmL) <- lod_pgmL[,1]
      }
      
      Data$SampleName <- renameSC(Data$SampleName, n=numSCs(Data, AQ))
    }
    
    Data_colnames <- colnames(Data)
    if(is.null(target_column_names)){
      target_column_names <- intersect(c('Target', 'AlamarTargetId','AlamarTargetID', 'UniProtId', 'UniProtID', 'ProteinName'), colnames(Data))
    }
    if(is.null(sample_column_names)){
      sample_column_names <- intersect(Data_colnames[!(Data_colnames %in% c(target_column_names, 
                                                                            'Panel', 'PanelLotNumber',
                                                                            'LoD', 'LOD', 'NPQ'))], colnames(Data))
    }
    targets <- unique(Data[,target_column_names])
    samples <- unique(Data[,sample_column_names])
    # make an LOD data frame
    LODnames <- intersect(c('PlateId', 'PlateID', 'Target', 'LOD', 'LoD'), colnames(Data))
    LOD$info <- unique(Data[, LODnames])
    # reformat the NPQ data
    # reformat into wide, targets in columns
    Data <- reshape(Data[,c('SampleName', 'Target', 'NPQ')],
                    direction= 'wide',
                    idvar = 'Target',
                    timevar='SampleName')
    if("LoD" %in% names(LOD$info)){
      LOD$LODNPQ <- LOD$info$LoD
    } else if ("LOD" %in% names(LOD$info)){
      LOD$LODNPQ <- LOD$info$LOD
    }
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
#' @param TAP If TRUE (default), uses TAP detectability criteria in sample QC 
#' which includes more matrix types than non-TAP criteria.
#' @param Bridge string(s) present in the sample names
#' that represent the bridge sample wells. 
#' Set to \code{NULL} (default) to use the 
#' type variable from the Barcode B file or 
#' if there are no bridge samples (default).
#' Only used for xml file formats.
#' @param sample_group_covar Optional column name in the Barcode B file 
#' and samples data matrix output by readNULISAseq that represents subgroups
#' for which detectability will be calculated separately, in addition to 
#' overall detectability. Default is 'SAMPLE_MATRIX', Function will check first to
#' be sure that the variable is present in the column names of the samples matrix.
#' Can be set to NULL to not use this feature.
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
#' @param security Logical. Default is TRUE. Should security checks be performed 
#' before generating AQ data.
#' @param excludeSamples A vector of sample names that will be excluded from all outputs.
#' @param excludeTargets A vector of target names that will be excluded from all outputs.
#' @param advancedQC Whether to use advancedQC metrics
#'
#' @return List of lists, data frames, and matrices.
#' Output will differ slightly depending on the input file type.
#'
#'
#' @export
#'
loadNULISAseq <- function(file,
                          IC=NULL,
                          IPC=NULL,
                          SC=NULL,
                          NC=NULL,
                          TAP=TRUE,
                          Bridge=NULL,
                          sample_group_covar='SAMPLE_MATRIX',
                          plateID=NULL,
                          scaleFactor=10^4,
                          transformReverse_scaleFactor=10^4,
                          replace_cal_blank_zeros = FALSE,
                          replace_zeros_with_NA = TRUE,
                          security=TRUE,
                          excludeSamples=NULL,
                          excludeTargets=NULL,
                          advancedQC = FALSE,
                          ...){
  raw <- readNULISAseq(file, IPC=IPC, IC=IC, SC=SC, NC=NC, 
                       plateID=plateID, ...)
  
  # xml file name without path 
  raw$xmlFile <- basename(file)
  
  if(is.null(IC)){ # if IC is not specified, get the ICs from the XML
    IC <- raw$IC 
  } else{ # if the IC is specified, use it instead of wha's in the XML
    raw$IC <- IC
  }
  
  if(length(raw$IPC) == 0 || length(raw$NC) == 0){
    stop("The input file does not contain IPC or NC sample types!\n")
  }
  if(length(raw$SC) == 0){
    warning("The input file does not contain the SC sample type!\n")
  }
  if(!is.null(excludeSamples)){
    raw$samples <- raw$samples[!(raw$samples$sampleName %in% excludeSamples),]
    raw$Data <- raw$Data[,!(colnames(raw$Data) %in% excludeSamples)]
    raw$IPC <- raw$IPC[!(raw$IPC %in% excludeSamples)]
    raw$SC <- raw$SC[!(raw$SC %in% excludeSamples)]
    raw$NC <- raw$NC[!(raw$NC %in% excludeSamples)]
    raw$Bridge <- raw$Bridge[!(raw$Bridge %in% excludeSamples)]
    raw$Calibrator <- raw$Calibrator[!(raw$Calibrator %in% excludeSamples)]
    raw$SampleNames <- raw$SampleNames[!(raw$SampleNames %in% excludeSamples)]
  }
  
  if(!is.null(excludeTargets)){
    raw$targets <- raw$targets[!(raw$targets$targetName %in% excludeTargets),]
    raw$Data <- raw$Data[!(rownames(raw$Data) %in% excludeTargets),]
  }

  # Filter hidden targets from data (hide=TRUE targets are excluded from all analyses)
  if("hide" %in% colnames(raw$targets)){
    hiddenTargets <- raw$targets$targetName[which(raw$targets$hide == TRUE)]
    if(length(hiddenTargets) > 0){
      raw$targets <- raw$targets[!raw$targets$targetName %in% hiddenTargets,]
      raw$Data <- raw$Data[!rownames(raw$Data) %in% hiddenTargets,]
      raw$IC <- raw$IC[!(raw$IC %in% hiddenTargets)]
    }
  }

  raw$IC_normed <- intraPlateNorm(data_matrix=raw$Data, IC=IC)
  ind <-  which(substr(raw$targets$Curve_Quant, 1, 1) == "R")
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
            replace_zeros_with_NA = replace_zeros_with_NA,
            security=security
          )
        },
        error = function(cond){
          message("Could not perform absolute quantification\n")
          stop(conditionMessage(cond))
        }
      )
    } else{ # Use AQ values if they were embedded (requires XML version >= 1.3.0)
      # Check XML version - fallback mode only works with v1.3.0+
      xml_version <- raw$ExecutionDetails$XMLversion
      fallback_supported <- !is.null(xml_version) &&
        as.numeric(sub("^([0-9]+\\.[0-9]+).*", "\\1", xml_version)) >= 1.3

      if(fallback_supported && "attributes" %in% names(raw)){
        if("AQ" %in% names(raw$attributes)){
          # Initialize Data_AQ with same dimensions as raw$Data to handle excludeSamples/excludeTargets
          raw$AQ$Data_AQ <- matrix(NA_real_, nrow = nrow(raw$Data), ncol = ncol(raw$Data),
                                   dimnames = dimnames(raw$Data))
          aq_attr <- raw$attributes$AQ
          if(!is.null(aq_attr)){
            aq_clean <- replace(aq_attr, is.nan(aq_attr), NA)
            # Align by row/column names if available, otherwise use direct assignment if dims match
            if(!is.null(rownames(aq_clean)) && !is.null(colnames(aq_clean))){
              common_rows <- intersect(rownames(raw$Data), rownames(aq_clean))
              common_cols <- intersect(colnames(raw$Data), colnames(aq_clean))
              if(length(common_rows) > 0 && length(common_cols) > 0){
                raw$AQ$Data_AQ[common_rows, common_cols] <- aq_clean[common_rows, common_cols]
              }
            } else if(identical(dim(aq_clean), dim(raw$Data))){
              raw$AQ$Data_AQ[,] <- aq_clean
            }
          }
          # Substitute NA when reads are zero AND concentration is 0
          zero_reads_zero_conc <- (raw$Data == 0) & (raw$AQ$Data_AQ == 0)
          raw$AQ$Data_AQ[zero_reads_zero_conc] <- NA
        }
        if("aM" %in% names(raw$attributes)){
          # Initialize Data_AQ_aM with same dimensions as raw$Data to handle excludeSamples/excludeTargets
          raw$AQ$Data_AQ_aM <- matrix(NA_real_, nrow = nrow(raw$Data), ncol = ncol(raw$Data),
                                      dimnames = dimnames(raw$Data))
          am_attr <- raw$attributes$aM
          if(!is.null(am_attr)){
            am_clean <- replace(am_attr, is.nan(am_attr), NA)
            # Align by row/column names if available, otherwise use direct assignment if dims match
            if(!is.null(rownames(am_clean)) && !is.null(colnames(am_clean))){
              common_rows <- intersect(rownames(raw$Data), rownames(am_clean))
              common_cols <- intersect(colnames(raw$Data), colnames(am_clean))
              if(length(common_rows) > 0 && length(common_cols) > 0){
                raw$AQ$Data_AQ_aM[common_rows, common_cols] <- am_clean[common_rows, common_cols]
              }
            } else if(identical(dim(am_clean), dim(raw$Data))){
              raw$AQ$Data_AQ_aM[,] <- am_clean
            }
          }
          # Substitute NA when reads are zero AND concentration is 0
          zero_reads_zero_conc <- (raw$Data == 0) & (raw$AQ$Data_AQ_aM == 0)
          raw$AQ$Data_AQ_aM[zero_reads_zero_conc] <- NA
        }

        # Filter to only include targets that have AQ parameters in ExecutionDetails$Abs
        # This ensures fallback mode returns the same targets as NULISAseqAQ would
        if("Abs" %in% names(raw$ExecutionDetails) && !is.null(raw$AQ$Data_AQ_aM)){
          # Validate raw$targets exists with required columns before using it
          if(is.null(raw$targets) || !is.data.frame(raw$targets) ||
             !all(c("targetName", "targetBarcode") %in% colnames(raw$targets))){
            warning("Cannot filter AQ targets: raw$targets missing or lacks required columns")
          } else {
            abs_data <- raw$ExecutionDetails$Abs
            # ExecutionDetails$Abs is a data.frame with 'seq' and 'val' columns
            # The 'seq' column contains target barcodes (and special entries like HashVal, ColumnNames)
            if(is.data.frame(abs_data) && "seq" %in% colnames(abs_data)){
              aq_barcodes <- abs_data$seq[!abs_data$seq %in% c("ColumnNames", "HashVal")]
            } else {
              # Fallback for older format: named vector
              seq_vals <- unname(abs_data)[seq(1, length(abs_data), 2)]
              aq_barcodes <- seq_vals[!seq_vals %in% c("ColumnNames", "HashVal")]
            }

            # Map barcodes to target names
            barcode_to_name <- setNames(raw$targets$targetName, raw$targets$targetBarcode)
            aq_target_names <- barcode_to_name[aq_barcodes]
            aq_target_names <- aq_target_names[!is.na(aq_target_names)]

            # Filter Data_AQ_aM to only include these targets
            keep_targets <- rownames(raw$AQ$Data_AQ_aM) %in% aq_target_names
            raw$AQ$Data_AQ_aM <- raw$AQ$Data_AQ_aM[keep_targets, , drop=FALSE]
            if(!is.null(raw$AQ$Data_AQ)){
              raw$AQ$Data_AQ <- raw$AQ$Data_AQ[keep_targets, , drop=FALSE]
            }
          }
        }
      } else {
        # Fallback not supported (XML version < 1.3) - warn and treat as RQ
        warning("This XML file contains AQ metadata but XML version < 1.3.0 and NULISAseqAQ is not installed. ",
                "AQ data will not be available. Processing as RQ (relative quantification) instead.")
        # Set AbsAssay to FALSE so downstream code treats this as RQ-only
        AbsAssay <- FALSE
      }
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
  
  # Get indices where sample_matrix is mentioned in curve_quant and has suffix -NA
  raw$match_matrix <- calcSampleTargetNAs(raw$targets$Curve_Quant, raw$samples$SAMPLE_MATRIX)
  
  # calculate LODs on IPC-normalized data
  # Note that the reverse-curve is already applied since this went through interPlateNorm
  # This means that aboveLOD is wrong, so we need to recalculate this using normed_unstranformedReverse$interNormData[[1]]
  raw$lod <- lod(data_matrix=raw$normed$interNormData[[1]], blanks=raw$NC, min_count=0, targetNoOutlierDetection=targetNoOutlierDetection, match_matrix=raw$match_matrix)
  
  # Recalculate lod without performing the reverse curve calculation we did in interPlateNorm so we can get the aboveLOD correct
  lodTemp <- lod(data_matrix=raw$normed_untransformedReverse$interNormData[[1]], blanks=raw$NC, min_count=0, targetNoOutlierDetection=targetNoOutlierDetection, match_matrix=raw$match_matrix)
  
  # Create LOD values on NPQ data. Need to perform IPC normalization on IC-normalized 
  # LOD vals so that we can report these to users
  raw$lod$LODNPQ <- log2(raw$lod$LOD + 1)
  raw$lod$untransformedReverse_LODNPQ <- log2(lodTemp$LOD + 1)
  
  ## manually apply reverse curve correction since we didn't go through the interPlateNorm
  if(length(reverseCurve) > 0){
    
    #raw$lod$LODNPQ[reverseCurve] <- log2((transformReverse_scaleFactor * scaleFactor) / (lodTemp$LOD[reverseCurve] + 1) + 1) # transformReverse_scaleFactor * scaleFactor 
    raw$lod$LODNPQ[reverseCurve] <- NA 
    raw$lod$LOD[reverseCurve] <- NA
    ## replace the reverse curve aboveLOD values with ones calculated from untransformed Reverse
    raw$lod$aboveLOD[reverseCurve, ] <- lodTemp$aboveLOD[reverseCurve, ]
  }
  # calculate Detectability
  sample_groups <- NULL
  if(!is.null(sample_group_covar) & sample_group_covar %in% colnames(raw$samples)) sample_groups <- raw$samples[raw$samples$sampleType == "Sample", sample_group_covar]
  raw$detectability <- detectability(aboveLOD_matrix=raw$lod$aboveLOD, 
                                     sample_subset=raw$samples$sampleName[which(raw$samples$sampleType == "Sample")], 
                                     sample_groups=sample_groups, 
                                     exclude_targets=raw$IC)
  if(AbsAssay){
    if(requireNamespace("NULISAseqAQ", quietly=T)){ 
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
            replace_zeros_with_NA = replace_zeros_with_NA,
            security=security
          )
          raw$lod$LOD_pgmL <- temp$Data_AQ[, 1]
          raw$lod$LOD_aM <- temp$Data_AQ_aM[, 1]
          # set zero normCount LODs to NA in the AQ LOD, they are likely already NA but maybe there are exceptions
          ind <- intersect(names(which(s[, 1] == 0 & !is.na(s[, 1]))), names(raw$lod$LOD_pgmL))
          if(length(ind) > 0){
            raw$lod$LOD_aM[ind] <- NA
            raw$lod$LOD_pgmL[ind] <- NA
          }
          
          # add LODs to targetAQ_param
          if(!is.null(raw$lod$LOD_pgmL)){
            AQ_LODs <- data.frame(target=names(raw$lod$LOD_aM),
                                  LOD_aM=raw$lod$LOD_aM, 
                                  LOD_pg_ml=raw$lod$LOD_pgmL)
          } else {
            AQ_LODs <- data.frame(target=names(raw$lod$LOD_aM),
                                  LOD_aM=raw$lod$LOD_aM)
          }
          
          raw$AQ$targetAQ_param <- merge(raw$AQ$targetAQ_param, AQ_LODs,
                                         by.x='targetName', by.y='target',
                                         all.x=TRUE)

          # Check if LLOQ/ULOQ were calculated - if not, will need to use qcXML later
          lloq_uloq_calculated <- "LLOQ" %in% colnames(raw$AQ$targetAQ_param) &&
                                  "ULOQ" %in% colnames(raw$AQ$targetAQ_param)

          if(lloq_uloq_calculated){
            # set run LLOQ = LOD if the LLOQ is below run LOD
            raw$AQ$targetAQ_param$LLOQ[raw$AQ$targetAQ_param$LLOQ < raw$AQ$targetAQ_param$LOD_aM & !is.na(raw$AQ$targetAQ_param$LOD_aM)] <- raw$AQ$targetAQ_param$LOD_aM[raw$AQ$targetAQ_param$LLOQ < raw$AQ$targetAQ_param$LOD_aM & !is.na(raw$AQ$targetAQ_param$LOD_aM)]
            if(!is.null(raw$AQ$targetAQ_param$LLOQ_pg_ml)){
              raw$AQ$targetAQ_param$LLOQ_pg_ml[raw$AQ$targetAQ_param$LLOQ_pg_ml < raw$AQ$targetAQ_param$LOD_pg_ml & !is.na(raw$AQ$targetAQ_param$LOD_pg_ml)] <- raw$AQ$targetAQ_param$LOD_pg_ml[raw$AQ$targetAQ_param$LLOQ_pg_ml < raw$AQ$targetAQ_param$LOD_pg_ml & !is.na(raw$AQ$targetAQ_param$LOD_pg_ml)]
            }
            raw$AQ$withinDR <- ((raw$AQ$Data_AQ_aM >= raw$AQ$targetAQ_param$LLOQ) & (raw$AQ$Data_AQ_aM <= raw$AQ$targetAQ_param$ULOQ))
          }
        },
        error = function(cond){
          message("Could not perform absolute quantification\n")
          stop(conditionMessage(cond))
        }
      )
    } else{ # Use AQ values embedded in XML if available
      if("attributes" %in% names(raw)){
        if("dr" %in% names(raw$attributes)){
          raw$AQ$withinDR <- matrix(as.logical(raw$attributes$dr), 
                                    nrow=nrow(raw$attributes$dr),
                                    ncol=ncol(raw$attributes$dr),
                                    dimnames=dimnames(raw$attributes$dr))
        }
      }
    }
  }

  # Fallback: Use qcXML values only if LLOQ/ULOQ were not calculated
  # This happens when NULISAseqAQ is not available OR when applyAQ didn't generate them
  if(AbsAssay && (is.null(raw$AQ$targetAQ_param) || !"LLOQ" %in% colnames(raw$AQ$targetAQ_param))){
    if("qcXML" %in% names(raw) && !is.null(raw$qcXML$target)){
      message("LLOQ/ULOQ not calculated - using values from XML")

      # Create targetAQ_param from raw$targets if it doesn't exist
      if(is.null(raw$AQ$targetAQ_param)){
        if(is.null(raw$targets) || !is.data.frame(raw$targets)){
          warning("Cannot create targetAQ_param: raw$targets is missing or not a data.frame")
        } else if(!"targetName" %in% colnames(raw$targets)){
          warning("Cannot create targetAQ_param: raw$targets missing required column 'targetName'")
        } else {
          raw$AQ$targetAQ_param <- raw$targets
        }
      }

      # Filter targetAQ_param to only include quantifiable targets (those in Data_AQ_aM)
      # Guard against NULL targetAQ_param (can happen if raw$targets was missing/invalid)
      if(!is.null(raw$AQ$targetAQ_param) && !is.null(raw$AQ$Data_AQ_aM)){
        quantifiable_targets <- rownames(raw$AQ$Data_AQ_aM)
        raw$AQ$targetAQ_param <- raw$AQ$targetAQ_param[
          raw$AQ$targetAQ_param$targetName %in% quantifiable_targets, , drop=FALSE]
      }

      # Merge LLOQ and ULOQ data from qcXML into targetAQ_param
      # Skip if targetAQ_param is NULL
      if(is.null(raw$AQ$targetAQ_param)){
        warning("Skipping LLOQ/ULOQ merge: targetAQ_param could not be created")
      } else {
      if(!is.null(raw$qcXML$target$LLOQ) &&
         nrow(raw$qcXML$target$LLOQ) > 0 &&
         all(c("targetName", "LLOQ_aM") %in% colnames(raw$qcXML$target$LLOQ))){
        # Remove existing LLOQ columns to avoid .x/.y suffixes from merge
        raw$AQ$targetAQ_param <- raw$AQ$targetAQ_param[,
          !colnames(raw$AQ$targetAQ_param) %in% c("LLOQ", "LLOQ_aM", "LLOQ_pg_ml"), drop = FALSE]

        # Explicitly create a clean dataframe with only needed columns
        lloq_cols <- intersect(c("targetName", "LLOQ_aM", "LLOQ_pg_ml"),
                               colnames(raw$qcXML$target$LLOQ))
        lloq_df <- data.frame(raw$qcXML$target$LLOQ[, lloq_cols, drop = FALSE])

        # Convert to numeric (qcXML values may be character from XML parsing)
        if("LLOQ_aM" %in% colnames(lloq_df)) lloq_df$LLOQ_aM <- as.numeric(lloq_df$LLOQ_aM)
        if("LLOQ_pg_ml" %in% colnames(lloq_df)) lloq_df$LLOQ_pg_ml <- as.numeric(lloq_df$LLOQ_pg_ml)

        # Substitute NA for zero LLOQ values
        if("LLOQ_aM" %in% colnames(lloq_df)) lloq_df$LLOQ_aM[lloq_df$LLOQ_aM == 0] <- NA
        if("LLOQ_pg_ml" %in% colnames(lloq_df)) lloq_df$LLOQ_pg_ml[lloq_df$LLOQ_pg_ml == 0] <- NA

        raw$AQ$targetAQ_param <- merge(raw$AQ$targetAQ_param,
                                       lloq_df,
                                       by="targetName", all.x=TRUE, sort=FALSE)
        raw$AQ$targetAQ_param$LLOQ <- raw$AQ$targetAQ_param$LLOQ_aM
      }

      if(!is.null(raw$qcXML$target$ULOQ) &&
         nrow(raw$qcXML$target$ULOQ) > 0 &&
         all(c("targetName", "ULOQ_aM") %in% colnames(raw$qcXML$target$ULOQ))){
        # Remove existing ULOQ columns to avoid .x/.y suffixes from merge
        raw$AQ$targetAQ_param <- raw$AQ$targetAQ_param[,
          !colnames(raw$AQ$targetAQ_param) %in% c("ULOQ", "ULOQ_aM", "ULOQ_pg_ml"), drop = FALSE]

        # Explicitly create a clean dataframe with only needed columns
        uloq_cols <- intersect(c("targetName", "ULOQ_aM", "ULOQ_pg_ml"),
                               colnames(raw$qcXML$target$ULOQ))
        uloq_df <- data.frame(raw$qcXML$target$ULOQ[, uloq_cols, drop = FALSE])

        # Convert to numeric (qcXML values may be character from XML parsing)
        if("ULOQ_aM" %in% colnames(uloq_df)) uloq_df$ULOQ_aM <- as.numeric(uloq_df$ULOQ_aM)
        if("ULOQ_pg_ml" %in% colnames(uloq_df)) uloq_df$ULOQ_pg_ml <- as.numeric(uloq_df$ULOQ_pg_ml)

        # Substitute NA for zero ULOQ values
        if("ULOQ_aM" %in% colnames(uloq_df)) uloq_df$ULOQ_aM[uloq_df$ULOQ_aM == 0] <- NA
        if("ULOQ_pg_ml" %in% colnames(uloq_df)) uloq_df$ULOQ_pg_ml[uloq_df$ULOQ_pg_ml == 0] <- NA

        raw$AQ$targetAQ_param <- merge(raw$AQ$targetAQ_param,
                                       uloq_df,
                                       by="targetName", all.x=TRUE, sort=FALSE)
        raw$AQ$targetAQ_param$ULOQ <- raw$AQ$targetAQ_param$ULOQ_aM
      }

      # Merge SC_conc data from qcXML into targetAQ_param
      # Expected XML format: <SC_conc aq_aM="X" aq_pgmL="Y"/>
      # Note: readQCXMLNode renames aq_aM -> SC_conc_aM, aq_pgmL -> SC_conc_pg_ml
      if(!is.null(raw$qcXML$target$SC_conc) &&
         nrow(raw$qcXML$target$SC_conc) > 0 &&
         all(c("targetName", "SC_conc_aM") %in% colnames(raw$qcXML$target$SC_conc))){
        # Remove existing SC_conc columns to avoid .x/.y suffixes from merge
        raw$AQ$targetAQ_param <- raw$AQ$targetAQ_param[,
          !colnames(raw$AQ$targetAQ_param) %in% c("SC_conc", "SC_conc_aM", "SC_conc_pg_ml"), drop = FALSE]

        # Explicitly create a clean dataframe with only needed columns
        sc_cols <- intersect(c("targetName", "SC_conc_aM", "SC_conc_pg_ml"),
                             colnames(raw$qcXML$target$SC_conc))
        sc_df <- data.frame(raw$qcXML$target$SC_conc[, sc_cols, drop = FALSE])

        # Convert to numeric (qcXML values may be character from XML parsing)
        if("SC_conc_aM" %in% colnames(sc_df)) suppressWarnings(sc_df$SC_conc_aM <- as.numeric(sc_df$SC_conc_aM))
        if("SC_conc_pg_ml" %in% colnames(sc_df)) suppressWarnings(sc_df$SC_conc_pg_ml <- as.numeric(sc_df$SC_conc_pg_ml))

        # Substitute NA for zero SC_conc values
        if("SC_conc_aM" %in% colnames(sc_df)) sc_df$SC_conc_aM[sc_df$SC_conc_aM == 0] <- NA
        if("SC_conc_pg_ml" %in% colnames(sc_df)) sc_df$SC_conc_pg_ml[sc_df$SC_conc_pg_ml == 0] <- NA

        raw$AQ$targetAQ_param <- merge(raw$AQ$targetAQ_param,
                                       sc_df,
                                       by="targetName", all.x=TRUE, sort=FALSE)
        raw$AQ$targetAQ_param$SC_conc <- raw$AQ$targetAQ_param$SC_conc_aM
      }

      if(!is.null(raw$qcXML$target$lod) &&
         nrow(raw$qcXML$target$lod) > 0 &&
         all(c("targetName", "LOD_aM") %in% colnames(raw$qcXML$target$lod))){
        # Explicitly create a clean dataframe with only needed columns
        lod_cols <- intersect(c("targetName", "LOD_aM", "LOD_pg_ml"),
                              colnames(raw$qcXML$target$lod))
        lod_df <- data.frame(raw$qcXML$target$lod[, lod_cols, drop = FALSE])

        # Convert to numeric (qcXML values may be character from XML parsing)
        if("LOD_aM" %in% colnames(lod_df)) suppressWarnings(lod_df$LOD_aM <- as.numeric(lod_df$LOD_aM))
        if("LOD_pg_ml" %in% colnames(lod_df)) suppressWarnings(lod_df$LOD_pg_ml <- as.numeric(lod_df$LOD_pg_ml))

        # Substitute NA for zero LOD values
        if("LOD_aM" %in% colnames(lod_df)) lod_df$LOD_aM[lod_df$LOD_aM == 0] <- NA
        if("LOD_pg_ml" %in% colnames(lod_df)) lod_df$LOD_pg_ml[lod_df$LOD_pg_ml == 0] <- NA

        raw$AQ$targetAQ_param <- merge(raw$AQ$targetAQ_param,
                                       lod_df,
                                       by="targetName", all.x=TRUE, sort=FALSE,
                                       suffixes=c("", "_qcXML"))
        
        # replace LLOQ with LOD when LOD > LLOQ: aM
        idx <- which(raw$AQ$targetAQ_param$LOD_aM > raw$AQ$targetAQ_param$LLOQ_aM)
        raw$AQ$targetAQ_param$LLOQ_aM[idx] <- raw$AQ$targetAQ_param$LOD_aM[idx]
        raw$AQ$targetAQ_param$LLOQ <- raw$AQ$targetAQ_param$LLOQ_aM
        # replace LLOQ with LOD when LOD > LLOQ: pg_ml
        idx <- which(raw$AQ$targetAQ_param$LOD_pg_ml > raw$AQ$targetAQ_param$LLOQ_pg_ml)
        raw$AQ$targetAQ_param$LLOQ_pg_ml[idx] <- raw$AQ$targetAQ_param$LOD_pg_ml[idx]

        # Also create raw$lod$LOD_aM as named vector (needed by quantifiability)
        lod_aM_vec <- lod_df$LOD_aM
        names(lod_aM_vec) <- lod_df$targetName
        # Warn if target names don't match between lod and qcXML
        if(!is.null(raw$lod) && !is.null(raw$lod$LOD) && !is.null(names(raw$lod$LOD))){
          matched_names <- intersect(names(raw$lod$LOD), names(lod_aM_vec))
          if(length(matched_names) < length(names(raw$lod$LOD))){
            warning("Some targets in raw$lod$LOD not found in qcXML LOD data")
          }
          raw$lod$LOD_aM <- lod_aM_vec[names(raw$lod$LOD)]
        } else {
          # If raw$lod$LOD doesn't exist, just use lod_aM_vec directly
          if(is.null(raw$lod)) raw$lod <- list()
          raw$lod$LOD_aM <- lod_aM_vec
        }

        # Create raw$lod$LOD_pgmL as named vector
        if("LOD_pg_ml" %in% colnames(lod_df)){
          lod_pgmL_vec <- lod_df$LOD_pg_ml
          names(lod_pgmL_vec) <- lod_df$targetName
          if(!is.null(raw$lod) && !is.null(raw$lod$LOD) && !is.null(names(raw$lod$LOD))){
            raw$lod$LOD_pgmL <- lod_pgmL_vec[names(raw$lod$LOD)]
          } else {
            raw$lod$LOD_pgmL <- lod_pgmL_vec
          }
        }
      }

      # Calculate withinDR if LLOQ/ULOQ and Data_AQ_aM are available
      if(!is.null(raw$AQ$Data_AQ_aM) &&
         !is.null(raw$AQ$targetAQ_param) &&
         "LLOQ" %in% colnames(raw$AQ$targetAQ_param) &&
         "ULOQ" %in% colnames(raw$AQ$targetAQ_param)){
        # Ensure proper row alignment between matrix and targetAQ_param
        target_order <- match(rownames(raw$AQ$Data_AQ_aM), raw$AQ$targetAQ_param$targetName)

        # Check for unmatched targets - NA in target_order means target not found
        unmatched <- is.na(target_order)
        if(any(unmatched)){
          unmatched_targets <- rownames(raw$AQ$Data_AQ_aM)[unmatched]
          warning(sprintf(
            "withinDR calculation: %d targets in Data_AQ_aM not found in targetAQ_param: %s",
            sum(unmatched),
            paste(head(unmatched_targets, 5), collapse=", ")
          ))
        }

        # Only calculate withinDR for matched targets to avoid NA propagation
        if(any(!unmatched)){
          raw$AQ$withinDR <- matrix(NA, nrow=nrow(raw$AQ$Data_AQ_aM), ncol=ncol(raw$AQ$Data_AQ_aM),
                                    dimnames=dimnames(raw$AQ$Data_AQ_aM))
          matched_idx <- which(!unmatched)
          raw$AQ$withinDR[matched_idx, ] <-
            ((raw$AQ$Data_AQ_aM[matched_idx, , drop=FALSE] >= raw$AQ$targetAQ_param$LLOQ[target_order[matched_idx]]) &
             (raw$AQ$Data_AQ_aM[matched_idx, , drop=FALSE] <= raw$AQ$targetAQ_param$ULOQ[target_order[matched_idx]]))
        }
      }
      } # End of else block for NULL targetAQ_param check
    }
  }

  # Use pre-calculated QC flags from XML only if NULISAseqAQ is not available
  if(AbsAssay && !requireNamespace("NULISAseqAQ", quietly=T)){
    if("qcXML" %in% names(raw)){
      if(!is.null(raw$qcXML$qcTarget)){
        raw$qcTarget <- raw$qcXML$qcTarget
      }
      if(!is.null(raw$qcXML$qcPlate)){
        raw$qcPlate <- raw$qcXML$qcPlate
      }
    }
  }

  # Perform QC at the target, plate and sample levels
  raw$qcTarget <- QCFlagTarget(AQdata=raw$AQ$Data_AQ_aM,
                               raw=raw$Data,
                               IPCnormed=raw$normed$interNormData[[1]],
                               detectability=raw$detectability$all$detectability,
                               aboveLOD=raw$lod$aboveLOD,
                               withinDR=raw$AQ$withinDR,
                               absRun=AbsAssay,
                               targets=raw$targets,
                               samples=raw$samples,
                               SCparams=raw$AQ$targetAQ_param$SC_conc, advancedQC=advancedQC)

  # QCS and SN are stored in raw$attributes if present in XML (will be NULL if not available)
  QCS <- raw$attributes$QCS
  SN <- raw$attributes$SN
  raw$qcSample <- QCFlagSample(raw$Data, raw$lod$aboveLOD, raw$samples, raw$targets, QCS=QCS, SN=SN, TAP=TAP)
  raw$qcPlate <- QCFlagPlate(raw$Data, raw$IC_normed$normData, raw$lod$aboveLOD, raw$targets, raw$samples, AQ=AbsAssay, AQ_QC=raw$qcTarget, Sample_QC=raw$qcSample)
  # calculate Detectability
  sample_groups <- NULL
  forDetectability <- raw$targets$targetName[which(raw$targets$noDetectability == FALSE)]
  if(!is.null(sample_group_covar) & sample_group_covar %in% colnames(raw$samples)) sample_groups <- raw$samples[raw$samples$sampleType == "Sample", sample_group_covar]
  noDetectability <- raw$targets$targetName[union(union(which(raw$targets$noDetectability == TRUE), 
                                                        which(raw$targets$hide == TRUE)),
                                                  which(raw$targets$targetType == "control"))]
  raw$detectability <- detectability(aboveLOD_matrix=lodTemp$aboveLOD[forDetectability, ],
                                     sample_subset=raw$samples$sampleName[which(raw$samples$sampleType == "Sample")], 
                                     sample_groups=sample_groups, 
                                     exclude_targets=noDetectability)
  
  # calculate Quantifiability
  if(AbsAssay){
    
    noQuantifiability <- raw$targets$targetName[which(raw$targets$hide == TRUE)]
    quant <- tryCatch({
      suppressWarnings(quantifiability(runs = raw, exclude_targets = noQuantifiability, sampleGroupCovar = sample_group_covar))
    }, error = function(cond) {
      message("Could not perform absolute quantification: ", conditionMessage(cond))
      return(NULL)
    })
    
    # restructure the output to match detectability
    raw$quantifiability <- list(
      sample_group = list(
        sampleNumber = as.list(quant$combined_quantifiability$n_samples[-1]),
        quantifiability = lapply(
          quant$combined_quantifiability$quant[-1],
          function(x) {
            vec <- as.numeric(x)
            names(vec) <- rownames(quant$combined_quantifiability$quant)
            return(vec)
          }
        )
      ),
      all = list(
        sampleNumber = as.integer(quant$combined_quantifiability$n_samples["overall"]),
        quantifiability = {
          vec <- as.numeric(quant$combined_quantifiability$quant$overall)
          names(vec) <- rownames(quant$combined_quantifiability$quant)
          vec
        }
      )
    )
  }
  # Do not report values for a specific target / sample_matrix combination if specified in barcodeA
  if(!is.null(raw$match_matrix) && nrow(raw$match_matrix) > 0 ) {
    raw$normed$log2_interNormData[[1]][cbind(raw$match_matrix[,"row"], raw$match_matrix[,"col"])] <- NA
    raw$NPQ[cbind(raw$match_matrix[,"row"], raw$match_matrix[,"col"])] <- NA

    # Helper to apply matrix-specific filtering to AQ data
    # Translates target indices via target names since AQ matrices may have different row ordering
    apply_aq_filter <- function(aq_matrix, targets, match_matrix) {
      if (is.null(aq_matrix)) return(aq_matrix)
      target_names_to_filter <- targets$targetName[match_matrix[,"row"]]
      aq_row_indices <- match(target_names_to_filter, rownames(aq_matrix))
      valid_idx <- !is.na(aq_row_indices)
      if (any(valid_idx)) {
        aq_matrix[cbind(aq_row_indices[valid_idx], match_matrix[valid_idx,"col"])] <- NA
      }
      aq_matrix
    }

    # Apply filtering to all AQ matrices
    raw$AQ$Data_AQ_aM <- apply_aq_filter(raw$AQ$Data_AQ_aM, raw$targets, raw$match_matrix)
    raw$AQ$Data_AQ <- apply_aq_filter(raw$AQ$Data_AQ, raw$targets, raw$match_matrix)
    raw$AQ$withinDR <- apply_aq_filter(raw$AQ$withinDR, raw$targets, raw$match_matrix)
  }
  
  # Add a new "qcSamplebyTarget" list as a convenience
  raw$qcSamplebyTarget$aboveReadThreshold = (raw$Data >= MIN_TARGET_READS)
  raw$qcSamplebyTarget$aboveLOD = raw$lod$aboveLOD
  if(AbsAssay){
    raw$qcSamplebyTarget$withinDR = raw$AQ$withinDR
  }
  raw$AbsAssay <- AbsAssay
  raw$advancedQC <- advancedQC
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
