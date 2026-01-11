#' Get XML Version from XML File
#'
#' This function reads an XML file and returns the XMLversion value.
#' If the XMLversion node does not exist, it returns "<1.3.0".
#'
#' @param infile The path to the XML file.
#'
#' @return A character string containing the XML version or "<1.3.0" if not found.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' version <- getXMLVersion("your_file.xml")
#' }
#'
#' @importFrom xml2 read_xml xml_find_first xml_text
#'
#' @export
getXMLVersion <- function(infile){
  # Read the XML file
  doc <- xml2::read_xml(infile)

  # Try to find the XMLversion node
  xmlVer <- xml2::xml_find_first(doc, "//XMLversion")

  # If the node doesn't exist, return "<1.3.0"
  if(length(xmlVer) == 0){
    return("<1.3.0")
  }

  # Return the text content of the XMLversion node
  return(xml2::xml_text(xmlVer))
}

#' Read Covariates from XML File
#'
#' This function reads covariate data from an XML file, extracts attributes
#' of Barcode elements under <BarcodeB>, and returns a data structure where
#' specified columns are excluded, and all remaining columns are converted to factors.
#'
#' @param infile The path to the XML file.
#' @param cols_to_exclude A character vector specifying the columns to exclude.
#'   Default is c("V1", "AUTO_PLATE", "AUTO_WELLCOL", "AUTO_WELLPOSITION", "AUTO_WELLROW",
#'   "EXPERIMENT_ID", "OPERATOR", "SEQ_INSTRUMENT_ID", "SEQ_RUN_DATE", "SMI_ID").
#'
#' @return A data structure with excluded columns and all remaining columns converted to factors.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' readCovariates("your_file.xml")
#' }
#'
#' @importFrom xml2 read_xml xml_find_all
#'
#' @export
readCovariates <- function(infile, cols_to_exclude = c("V1", "AUTO_PLATE", "AUTO_WELLCOL", "AUTO_WELLPOSITION", "AUTO_WELLROW", "EXPERIMENT_ID", "OPERATOR", "SEQ_INSTRUMENT_ID", "SEQ_RUN_DATE", "SMI_ID", "type")){
  # read in xml file
  xml_file <- xml2::read_xml(infile)
    
  # Extract attributes of Barcode elements under <BarcodeB>
  barcode_b_elements <- xml2::xml_find_all(xml_file, "//BarcodeB/Barcode")

  # Create a function to extract attributes from a Barcode element
  extract_attributes <- function(barcode_element) {
    attributes <- xml2::xml_attrs(barcode_element)
    text <- xml2::xml_text(barcode_element)
    return(c(text, attributes))
  }

  # Apply the function to each Barcode element under <BarcodeB>
  barcode_b_attributes <- lapply(barcode_b_elements, extract_attributes)

  # Convert the list to a data frame for easy manipulation
  barcode_b_df <- as.data.frame(do.call(rbind, barcode_b_attributes))

  if (!is.null(cols_to_exclude)) {
    barcode_b_df <- barcode_b_df[, !colnames(barcode_b_df) %in% cols_to_exclude, drop = FALSE]
  }

  barcode_b_df <- lapply(barcode_b_df, as.factor)
}

#' Removes BarcodeB nodes with a specified attribute and their associated Sample nodes.
#'
#' This function takes an XML document, along with the name and value of an attribute,
#' and removes the BarcodeB nodes that match the specified attribute. Additionally,
#' it identifies the associated Sample nodes with the same name attribute as the removed
#' BarcodeB nodes and removes them as well. These removals only occur if the sample/barcode
#' node in question is not a control (e.g. IPC, NC, SC, Bridge, etc.)
#'
#' @param doc An XML document.
#' @param attribute_name The name of the attribute to match in BarcodeB nodes.
#' @param attribute_value The value of the attribute to match in BarcodeB nodes. Can be multiple values.
#'
#' @return The modified XML document with BarcodeB and associated Sample nodes removed.
#'
#' @examples
#' doc <- read_xml("your_file.xml")
#' new_doc <- removeBarcodeAndSamples(doc, "name", c("774951", "23112"))
#'
#' @import xml2
#' @import glue
#' @export
removeBarcodeAndSamples <- function(doc, attribute_name, attribute_value) {

  # Make a copy of the input document
  doc_copy <- xml2::xml_unserialize(xml2::xml_serialize(doc, NULL))

  for (attr_val in attribute_value) {
    # Find BarcodeB nodes with the specified attribute
    xpath <- glue::glue("//BarcodeB/Barcode[@{attribute_name}='{attr_val}']")
    barcode_nodes <- xml2::xml_find_all(doc_copy, xpath)

    # Find Sample nodes with the same name as removed Barcodes
    for (barcode_node in barcode_nodes) {
      # Extract Barcode names before removing Barcode nodes
      barcodeName <- xml2::xml_attr(barcode_node, "name")

      # Find associated Sample nodes for the current Barcode node
      sample_xpath <- glue::glue("//Sample[@barcode='{barcodeName}']")
      sample_nodes <- xml2::xml_find_all(doc_copy, sample_xpath)
  
      # Check if any Sample node has the 'type' attribute with specific values
      # If found, retain the Sample nodes; otherwise, remove both Barcode and Sample nodes
      if (!any(xml2::xml_attr(barcode_node, "type") %in% c("IPC", "NC", "SC", "Bridge"))) {
        xml2::xml_remove(barcode_node)
        xml2::xml_remove(sample_nodes)
      }
    }
  }

  return(doc_copy)
}

#' Inserts Covariates into XML
#'
#' This function reads an XML file, modifies specific attributes within the <Barcode> elements 
#' under <BarcodeB>, and inserts corresponding covariate values from a provided data frame.
#'
#' @param input_XML Path to the input XML file.
#' @param covariates A data frame containing covariate information, including 'sampleName' 
#'                   and other covariate columns.
#'
#' @return The modified XML document.
#'
#' @examples
#' # Example usage:
#' # data <- loadNULISAseq('input.xml')
#' # data$samples$newAnnot <- c(1:length(1:nrow(data$samples)))
#' # insertCovariatesXML("input.xml", data$samples)
#'
#' @export
insertCovariatesXML <- function(input_XML, covariates){

  # Read the XML file
  doc <- xml2::read_xml(input_XML)

  if(is.null(covariates)){
    return (doc)
  }

  # Create wellRow and wellCol columns if they don't exist but AUTO_WELLPOSITION does
  if(!("wellRow" %in% colnames(covariates)) && !("wellCol" %in% colnames(covariates)) &&
     ("AUTO_WELLPOSITION" %in% colnames(covariates))) {
    covariates$wellRow <- substr(covariates$AUTO_WELLPOSITION, 1, 1)
    covariates$wellCol <- as.numeric(substring(covariates$AUTO_WELLPOSITION, 2))
  }

  # Attributes to skip during modification
  skip <- c("type", "TYPE", "name", "AUTO_PLATE", "AUTO_WELLCOL", "AUTO_WELLROW",
            "AUTO_WELLPOSITION", 'OPERATOR', 'SEQ_INSTRUMENT_ID', 'SEQ_RUN_DATE',
            'wellRow', 'wellCol')
  
  # Find all <Barcode> elements under <BarcodeB>
  barcode_b_elements <- xml2::xml_find_all(doc, "//BarcodeB/Barcode")

  # Loop through each <Barcode> element
  for (element in barcode_b_elements){
    
    # check for similar sampleNames by regex
    multi_check <- grep(xml2::xml_text(element), covariates$sampleName, fixed = TRUE)
    
    if(length(multi_check) == 1){
      
      # Find the index where the sampleName matches
      ind <- which(xml2::xml_text(element) == covariates$sampleName)
      if (length(ind) == 0){
        next;
      }
      # Loop through each covariate column
      # Set attribute with corresponding covariate value
      for (cov in colnames(covariates)){
        if (!(cov %in% skip)){
          value <- covariates[[cov]][ind]
          if (!is.null(value) && !is.na(value) && length(value) > 0) {
            xml2::xml_set_attr(element, cov, as.character(value))
          }
        }
      }
    } else if(length(multi_check) > 1){
      # get wellcol and wellrow info for current element
      if (!is.null(xml2::xml_attr(element, "AUTO_WELLPOSITION")) && xml2::xml_attr(element, "AUTO_WELLPOSITION") != "" &&
        (is.null(xml2::xml_attr(element, "AUTO_WELLCOL")) || is.null(xml2::xml_attr(element, "AUTO_WELLROW")) ||
         xml2::xml_attr(element, "AUTO_WELLCOL") == "" || xml2::xml_attr(element, "AUTO_WELLROW") == "")) {
        well_pos <- xml2::xml_attr(element, "AUTO_WELLPOSITION")
        req_attrs_val <- c(
          AUTO_WELLROW = substr(well_pos, 1, 1), 
          AUTO_WELLCOL = substring(well_pos, 2),
          wellRow = substr(well_pos, 1, 1), 
          wellCol = substring(well_pos, 2)
          )
      } else{
        req_attrs_val <- sapply(c('AUTO_WELLROW', 'AUTO_WELLCOL'), xml2::xml_attr, x = element)
      }
      #subset covariate information for the current element
      sub_covariates <- subset(covariates[multi_check, ], 
                               wellRow == req_attrs_val[1] & wellCol == sub('^0+', '', req_attrs_val[2]))
      
      # loop and apply to xml
      for (cov in colnames(covariates)) {
        if (!(cov %in% skip)) {
          value <- sub_covariates[[cov]]
          if (!is.null(value) && !is.na(value) && length(value) > 0) {
            xml2::xml_set_attr(element, cov, as.character(value))
          }
        }
      }
    } else{
      next;
    }
  }
  # Return the modified XML document
  return(doc)
}

#' Write Updated XML
#'
#' Reads and processes an XML file containing NULISAseq data. This function reads the input XML file
#' and returns the XML document with updated info (v1.3.0) from loadNULISAseq (e.g. AQ and QC info). 
#'
#' @param input_XML Path to the input XML file to be processed
#' @param rawReadCutoff The raw read value cutoff to apply to individual sample/target combinations
#' @param data output of loadNULISAseq (optional)
#' @return XML document object containing the NULISAseq data
#' @examples
#' # writeUpdatedXML("path/to/input.xml")
#'
#' @export
writeUpdatedXML <- function(input_XML, rawReadCutoff=200, data=NULL){
  # Read the XML file
  doc <- xml2::read_xml(input_XML)
  if(is.null(data)){
    data <- loadNULISAseq(input_XML)
  }
  root <- xml2::xml_root(doc)

  # Update the XML version #
  xmlVer <- xml_find_first(doc, "XMLversion")
  if(length(xmlVer) == 0){
    xmlVer <- xml_add_child(xml_find_first(doc, "ExecutionDetails"), "XMLversion")
  }
  xml_text(xmlVer) <- "1.3.0"

  # Create the QC thresholds
  QCthreshNode <- xml_add_child(root, "QCThresholds")
  
  createQCThresholdNode <- function(xmlNode, criteria){
    qcNames <- names(criteria[[1]])
    # Loop over each threshold name
    for (qcName in qcNames) {
      threshNode <- xml_add_child(xmlNode, "Threshold")
      xml_set_attr(threshNode, "name", gsub("\\.", "_", qcName))
      
      # Loop over each section in criteria and add it as an attribute
      for (section_name in names(criteria)) {
        if(qcName %in% names(criteria[[section_name]])){
          value <- criteria[[section_name]][[qcName]]
        } else{
          if (startsWith(qcName, "Detectability.")){
            newKey <- sub("\\..*", "", qcName)        
            value <- criteria[[section_name]][[newKey]]
          } else{
            warning(paste0('NOTE: QC [', qcName, ']', ' could not be found in [', section_name, ']'))
            next;
          }
        }
        # Optional: skip NULL or NA values
        if (!is.null(value) && !is.na(value)) {
          xml_set_attr(threshNode, section_name, as.character(value))
        }
      }
    }
  }

  # Write the QC thresholds into the XML
  PlateQCthreshNode <- xml_add_child(QCthreshNode, "PlateQC")
  PlateQCthresh <- QCPlateCriteria(AQ=data$AbsAssay)
  createQCThresholdNode(PlateQCthreshNode, PlateQCthresh)

  TargetQCthreshNode <- xml_add_child(QCthreshNode, "TargetQC")
  TargetQCthresh <- QCTargetCriteria(AQ=data$AbsAssay, advancedQC=data$advancedQC)
  createQCThresholdNode(TargetQCthreshNode, TargetQCthresh)

  SampleQCthreshNode <- xml_add_child(QCthreshNode, "SampleQC")
  SampleQCthresh <- QCSampleCriteria()
  createQCThresholdNode(SampleQCthreshNode, SampleQCthresh)

  # Write the Plate QC values
  PlateQCNode <- xml_add_child(root, "PlateQC")
  vals <- data$qcPlate %>% dplyr::select(val)
  plateFlags <- dplyr::rename(data$qcPlate, name="flagName", method="normMethod", set="status", format="QCformat") %>%
                  dplyr::select(-QCthreshold, -QCoperator, -val)
  for(i in 1:nrow(data$qcPlate)){
    QCFlag <- xml_add_child(PlateQCNode, "QCFlag")
    for(j in 1:ncol(plateFlags)){
      xml_set_attr(QCFlag, colnames(plateFlags)[j], plateFlags[i,j])
    } 
    xml_text(QCFlag) <- as.character(vals[[1]][i])
  }

  # Write the Sample QC values
  SampleQCNode <- xml_add_child(root, "SampleQC")
  vals <- data$qcSample %>% dplyr::select(val) 
  barcode <- data$qcSample %>% dplyr::select(sampleBarcode)
  sampleFlags <- dplyr::rename(data$qcSample, name="flagName", method="normMethod", set="status", format="QCformat") %>%
                  dplyr::select(-sampleName, -QCthreshold, -QCoperator, -val, -sampleType, -sampleBarcode, -text)
  uniqBarcodes <- unique(data$qcSample$sampleBarcode)
  for(i in 1:length(uniqBarcodes)){
    sampleQC <- xml_add_child(SampleQCNode, "Sample")
    xml_set_attr(sampleQC, "name", uniqBarcodes[i])
    qcflags <- sampleFlags[barcode == uniqBarcodes[i],]
    qcflagVals <- vals[barcode == uniqBarcodes[i]]
    for(j in 1:nrow(qcflags)){
      qcFlagNode <- xml_add_child(sampleQC, "QCFlag")
      xml_text(qcFlagNode) <- as.character(qcflagVals[j])
      for(k in 1:ncol(qcflags)){
        xml_set_attr(qcFlagNode, colnames(qcflags)[k], qcflags[j,k])
      }
    } 
  }

  # Write the Target QC values
  if(data$AbsAssay){
    TargetQCNode <- xml_add_child(root, "TargetQC")
    vals <- data$qcTarget %>% dplyr::select(val)
    targetFlags <- dplyr::rename(data$qcTarget, name="flagName", method="normMethod", set="status", format="QCformat") %>%
                  dplyr::select(-QCthreshold, -QCoperator)
    uniqTargets <- unique(data$qcTarget$target)
    for(i in 1:length(uniqTargets)){
      targetQC <- xml_add_child(TargetQCNode, "Target")
      targetBarcode <- data$targets %>% 
                          dplyr::filter(targetName == uniqTargets[i]) %>%
                          dplyr::pull(targetBarcode)
      xml_set_attr(targetQC, "name", targetBarcode)
      qcflags <- targetFlags %>% 
                          dplyr::filter(target == uniqTargets[i]) %>%
                          dplyr::select(-val, -target)
      qcflagVals <- targetFlags %>% 
                          dplyr::filter(target == uniqTargets[i]) %>% 
                          dplyr::pull(val)
      for(j in 1:nrow(qcflags)){
        qcFlagNode <- xml_add_child(targetQC, "QCFlag")
        xml_text(qcFlagNode) <- as.character(qcflagVals[j])
        for(k in 1:ncol(qcflags)){
          xml_set_attr(qcFlagNode, colnames(qcflags)[k], qcflags[j,k])
        }
      }

      # Function to create ULOQ/LLOQ/LOD nodes and set attributes      
      createQCnode <- function(parent_node, node_name, target_name, data_source, attr_map, use_named_vector = FALSE, name_column = "targetName") {
        node <- xml_add_child(parent_node, node_name)
        
        for (attr_name in names(attr_map)) {
          col_name <- attr_map[[attr_name]]
          
          value <- NA
          if (use_named_vector) {
            index <- which(names(data_source[[col_name]]) == target_name)
            if (length(index) != 0) {
              value <- data_source[[col_name]][index]
            }
          } else {
            index <- which(data_source[[name_column]] == target_name)
            if (length(index) != 0) {
              value <- data_source[[col_name]][index]
            }
          }
          if (!is.na(value)) {
            xml_set_attr(node, attr_name, as.character(value))
          }
        }
      }
      # Write LLOQ XML nodes
      createQCnode(
        parent_node = targetQC,
        node_name = "LLOQ",
        target_name = uniqTargets[i],
        data_source = data$AQ$targetAQ_param,
        attr_map = c("aq_aM" = "LLOQ", "aq_pgmL" = "LLOQ_pg_ml")
      )

      # Write ULOQ XML nodes
      createQCnode(
        parent_node = targetQC,
        node_name = "ULOQ",
        target_name = uniqTargets[i],
        data_source = data$AQ$targetAQ_param,
        attr_map = c("aq_aM" = "ULOQ", "aq_pgmL" = "ULOQ_pg_ml")
      )

      # Write LOD XML nodes — specify use_named_vector = TRUE
      createQCnode(
        parent_node = targetQC,
        node_name = "LOD",
        target_name = uniqTargets[i],
        data_source = data$lod,
        attr_map = c("rq" = "LODNPQ", "aq_aM" = "LOD_aM", "aq_pgmL" = "LOD_pgmL"),
        use_named_vector = TRUE
      )
    }
  }

  # Write the NPQ/AQ/etc. values into <Sample><ReadCount> XML nodes
  for(sample in xml2::xml_find_all(root, './/Data//Sample')){
    sBarcode <- xml2::xml_attr(sample, "barcode")
    sName <- data$samples %>% dplyr::filter(sampleBarcode == sBarcode) %>% dplyr::pull(sampleName)
    for(readcount in xml2::xml_find_all(sample, "ReadCount")){
      tBarcode <- xml2::xml_attr(readcount, "target")
      tName <- data$targets %>% dplyr::filter(targetBarcode == tBarcode) %>% dplyr::pull(targetName)
      if(tName %in% rownames(data$NPQ) && sName %in% colnames(data$NPQ)){
        xml_set_attr(readcount, "NPQ", as.character(data$NPQ[tName, sName]))
      }
      value <- as.integer(xml_text(readcount)) > rawReadCutoff
      xml_set_attr(readcount, "aR", as.character(as.integer(value)))

      if(tName %in% rownames(data$lod$aboveLOD) && sName %in% colnames(data$lod$aboveLOD)){
        xml_set_attr(readcount, "aB", as.character(as.integer(data$lod$aboveLOD)))
      }  
      if(data$AbsAssay){
        if(tName %in% rownames(data$AQ$Data_AQ) && sName %in% colnames(data$AQ$Data_AQ)){
          xml_set_attr(readcount, "AQ", as.character(data$AQ$Data_AQ[tName, sName]))
        }
        if(tName %in% rownames(data$AQ$Data_AQ_aM) && sName %in% colnames(data$AQ$Data_AQ_aM)){
          xml_set_attr(readcount, "aM", as.character(data$AQ$Data_AQ_aM[tName, sName]))
        }
        if(tName %in% rownames(data$AQ$withinDR) && sName %in% colnames(data$AQ$withinDR)){
          xml_set_attr(readcount, "dr", as.character(as.integer(data$AQ$withinDR[tName, sName])))
        }
      }
    }
  }
  return(as.character(doc))
}


#' Split XML file based on attribute and create a ZIP archive
#'
#' This endpoint reads an XML file, identifies specific elements based on an attribute name,
#' removes nodes with the specified attribute value, creates new XML files, and zips them together.
#'
#' @param xmlFile:file Character string vector. Path(s) and name(s) of the file(s).
#' @param attr_name2:[str] The attribute name used to identify nodes for splitting.
#' @export
splitXML <- function(xmlFile, attr_name2){
  doc <- xml2::read_xml(toString(xmlFile))
  barcode_b_elements <- xml2::xml_find_all(doc, "//BarcodeB/Barcode")
  attr_name <- trimws(strsplit(attr_name2, ":")[[1]][1])
  attr_name <- sub(" ", "_", attr_name)
  xpath <- paste0("//BarcodeB/Barcode[@", attr_name, "]")
  sample_nodes <- xml2::xml_find_all(doc, xpath)
  attr_vals <- unique(xml2::xml_attr(sample_nodes, attr_name))
  newFiles <- c()#vector(length = length(attr_vals))
  outFile <- paste0(Sys.Date(), ".xml")
  count <- 1
  for(i in 1:length(attr_vals)){
    attr_vals_toRemove <- attr_vals[-i]
    doc_removed <- removeBarcodeAndSamples(doc, attr_name, attr_vals_toRemove)

    # check the number of type=sample
    if(length(which(tolower(xml2::xml_attr(xml2::xml_find_all(doc_removed, "//Barcode[@type]"), "type")) == "sample")) > 0){
      # write XML
      newFiles[count] <- sub("\\.xml$", paste0("_", attr_name, "-", attr_vals[i], ".xml"), outFile)
      xml2::write_xml(doc_removed, newFiles[count])
      count <- count + 1
    }
  }
  # zip up the XMLs
  outFile <- paste0(outFile, ".zip")
  zip(outFile, files = normalizePath(newFiles), extras="-j")
  unlink(newFiles)
  bin <- readBin(outFile, "raw", n = file.info(outFile)$size)
  unlink(outFile)
  return(bin)
}
