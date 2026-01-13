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

# Internal helpers to make writeUpdatedXML faster and clearer
# Note: These are intentionally not exported and preserve existing behavior.

# Ensure the XMLversion node exists and is set to 1.3.0
#' Ensure XMLversion node
#'
#' Ensures the XML document contains an `ExecutionDetails/XMLversion` node
#' and sets its text to "1.3.0". If missing, creates the required nodes.
#'
#' @param doc An xml2 document.
#'
#' @return The same xml2 document with XMLversion set to 1.3.0.
#' @noRd
#' @keywords internal
#'
#' @importFrom xml2 xml_find_first xml_add_child xml_root xml_text
.ensure_xml_version <- function(doc){
  xmlVer <- xml2::xml_find_first(doc, "XMLversion")
  if(length(xmlVer) == 0){
    # Create XMLversion under ExecutionDetails if missing
    exec <- xml2::xml_find_first(doc, "ExecutionDetails")
    if(length(exec) == 0){
      exec <- xml2::xml_add_child(xml2::xml_root(doc), "ExecutionDetails")
    }
    xmlVer <- xml2::xml_add_child(exec, "XMLversion")
  }
  xml2::xml_text(xmlVer) <- "1.3.0"
  doc
}

# Write QC thresholds (PlateQC, TargetQC, SampleQC)
#' Add QC threshold nodes
#'
#' Writes `QCThresholds` and its child sections (`PlateQC`, `TargetQC`, `SampleQC`),
#' populating per-threshold attributes from computed criteria.
#'
#' @param root XML root node.
#' @param data Data list from `loadNULISAseq()` used to compute criteria.
#'
#' @return Invisibly returns NULL; mutates `root` by adding threshold nodes.
#' @noRd
#' @keywords internal
#'
#' @importFrom xml2 xml_add_child xml_set_attr
.add_qc_thresholds <- function(root, data){
  QCthreshNode <- xml2::xml_add_child(root, "QCThresholds")

  createQCThresholdNode <- function(xmlNode, criteria){
    qcNames <- names(criteria[[1]])
    for (qcName in qcNames) {
      threshNode <- xml2::xml_add_child(xmlNode, "Threshold")
      xml2::xml_set_attr(threshNode, "name", gsub("\\.", "_", qcName))
      for (section_name in names(criteria)) {
        if(qcName %in% names(criteria[[section_name]])){
          value <- criteria[[section_name]][[qcName]]
        } else{
          if (startsWith(qcName, "Detectability.")){
            newKey <- sub("\\..*", "", qcName)
            value <- criteria[[section_name]][[newKey]]
          } else{
            next
          }
        }
        if (!is.null(value) && !is.na(value)) {
          xml2::xml_set_attr(threshNode, section_name, as.character(value))
        }
      }
    }
  }

  PlateQCthreshNode <- xml2::xml_add_child(QCthreshNode, "PlateQC")
  PlateQCthresh <- QCPlateCriteria(AQ = data$AbsAssay)
  createQCThresholdNode(PlateQCthreshNode, PlateQCthresh)

  TargetQCthreshNode <- xml2::xml_add_child(QCthreshNode, "TargetQC")
  TargetQCthresh <- QCTargetCriteria(AQ = data$AbsAssay, advancedQC = data$advancedQC)
  createQCThresholdNode(TargetQCthreshNode, TargetQCthresh)

  SampleQCthreshNode <- xml2::xml_add_child(QCthreshNode, "SampleQC")
  SampleQCthresh <- QCSampleCriteria()
  createQCThresholdNode(SampleQCthreshNode, SampleQCthresh)
}

# Write Plate QC flags
#' Add Plate QC flags
#'
#' Adds `PlateQC/QCFlag` nodes, setting attributes from `data$qcPlate` and
#' the flag value as node text.
#'
#' @param root XML root node.
#' @param data Data list from `loadNULISAseq()`.
#'
#' @return Invisibly returns NULL; mutates `root` with plate QC flags.
#' @noRd
#' @keywords internal
#'
#' @importFrom xml2 xml_add_child xml_set_attrs xml_text
.add_plate_qc <- function(root, data){
  PlateQCNode <- xml2::xml_add_child(root, "PlateQC")
  vals <- data$qcPlate %>% dplyr::select(val)
  plateFlags <- dplyr::rename(data$qcPlate, name = "flagName", method = "normMethod", set = "status", format = "QCformat") %>%
    dplyr::select(-QCthreshold, -QCoperator, -val)
  for(i in seq_len(nrow(data$qcPlate))){
    QCFlag <- xml2::xml_add_child(PlateQCNode, "QCFlag")
    attrs <- setNames(as.character(plateFlags[i, ]), colnames(plateFlags))
    xml2::xml_set_attrs(QCFlag, attrs)
    xml2::xml_text(QCFlag) <- as.character(vals[[1]][i])
  }
}

# Write Sample QC flags grouped by barcode
#' Add Sample QC flags
#'
#' Adds `SampleQC/Sample` and child `QCFlag` nodes grouped by sample barcode,
#' setting attributes and flag text from `data$qcSample`.
#'
#' @param root XML root node.
#' @param data Data list from `loadNULISAseq()`.
#'
#' @return Invisibly returns NULL; mutates `root` with sample QC flags.
#' @noRd
#' @keywords internal
#'
#' @importFrom xml2 xml_add_child xml_set_attr xml_set_attrs xml_text
.add_sample_qc <- function(root, data){
  SampleQCNode <- xml2::xml_add_child(root, "SampleQC")
  vals <- data$qcSample %>% dplyr::select(val)
  barcode <- data$qcSample %>% dplyr::select(sampleBarcode)
  sampleFlags <- dplyr::rename(data$qcSample, name = "flagName", method = "normMethod", set = "status", format = "QCformat") %>%
    dplyr::select(-sampleName, -QCthreshold, -QCoperator, -val, -sampleType, -sampleBarcode, -text)
  uniqBarcodes <- unique(data$qcSample$sampleBarcode)
  for(bc in uniqBarcodes){
    sampleQC <- xml2::xml_add_child(SampleQCNode, "Sample")
    xml2::xml_set_attr(sampleQC, "name", bc)
    qcflags <- sampleFlags[barcode == bc, ]
    qcflagVals <- vals[barcode == bc]
    for(j in seq_len(nrow(qcflags))){
      qcFlagNode <- xml2::xml_add_child(sampleQC, "QCFlag")
      xml2::xml_text(qcFlagNode) <- as.character(qcflagVals[j])
      attrs <- setNames(as.character(qcflags[j, ]), colnames(qcflags))
      xml2::xml_set_attrs(qcFlagNode, attrs)
    }
  }
}

#' Create QC child node
#'
#' Creates a QC child node under the provided parent and sets attributes
#' based on a mapping from data source columns to attribute names.
#'
#' @param parent_node Parent XML node to attach the new child node to.
#' @param node_name Name of the child node to create (e.g., "LLOQ").
#' @param target_name Target identifier used to look up values in `data_source`.
#' @param data_source Data frame/list that contains values used for attributes.
#' @param attr_map Named character vector mapping attribute names to column keys.
#' @param use_named_vector Logical; when TRUE, look up by names in a vector; otherwise by `name_column`.
#' @param name_column Column name to use for non-vector lookups (default "targetName").
#'
#' @return Invisibly returns NULL; mutates `parent_node` by adding the child.
#' @noRd
#' @keywords internal
.createQCnode <- function(parent_node, node_name, target_name, data_source, attr_map,
                          use_named_vector = FALSE, name_column = "targetName") {
  node <- xml2::xml_add_child(parent_node, node_name)
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
      xml2::xml_set_attr(node, attr_name, as.character(value))
    }
  }
}

# Write Target QC nodes (LLOQ, ULOQ, LOD)

#' Add Target QC nodes
#'
#' Adds `TargetQC/Target` nodes with `QCFlag` entries and parameter children
#' (`LLOQ`, `ULOQ`, `LOD`) derived from `data$AQ$targetAQ_param` and `data$lod`.
#'
#' @param root XML root node.
#' @param data Data list from `loadNULISAseq()`; only runs when `data$AbsAssay` is TRUE.
#'
#' @return Invisibly returns NULL; mutates `root` with target QC nodes.
#' @noRd
#' @keywords internal
#'
#' @importFrom xml2 xml_add_child xml_set_attr xml_set_attrs xml_text
.add_target_qc <- function(root, data){
  if(!isTRUE(data$AbsAssay)) return(invisible(NULL))
  TargetQCNode <- xml2::xml_add_child(root, "TargetQC")
  vals <- data$qcTarget %>% dplyr::select(val)
  targetFlags <- dplyr::rename(data$qcTarget, name = "flagName", method = "normMethod", set = "status", format = "QCformat") %>%
    dplyr::select(-QCthreshold, -QCoperator)
  uniqTargets <- unique(data$qcTarget$target)

  for(tg in uniqTargets){
    targetQC <- xml2::xml_add_child(TargetQCNode, "Target")
    targetBarcode <- data$targets %>% dplyr::filter(targetName == tg) %>% dplyr::pull(targetBarcode)
    xml2::xml_set_attr(targetQC, "name", targetBarcode)

    qcflags <- targetFlags %>% dplyr::filter(target == tg) %>% dplyr::select(-val, -target)
    qcflagVals <- targetFlags %>% dplyr::filter(target == tg) %>% dplyr::pull(val)
    for(j in seq_len(nrow(qcflags))){
      qcFlagNode <- xml2::xml_add_child(targetQC, "QCFlag")
      xml2::xml_text(qcFlagNode) <- as.character(qcflagVals[j])
      attrs <- setNames(as.character(qcflags[j, ]), colnames(qcflags))
      xml2::xml_set_attrs(qcFlagNode, attrs)
    }

    .createQCnode(targetQC, "LLOQ", tg, data$AQ$targetAQ_param, c("aq_aM" = "LLOQ", "aq_pgmL" = "LLOQ_pg_ml"))
    .createQCnode(targetQC, "ULOQ", tg, data$AQ$targetAQ_param, c("aq_aM" = "ULOQ", "aq_pgmL" = "ULOQ_pg_ml"))
    .createQCnode(targetQC, "LOD",  tg, data$lod,               c("rq" = "LODNPQ", "aq_aM" = "LOD_aM", "aq_pgmL" = "LOD_pgmL"), use_named_vector = TRUE)
  }
}

# Annotate NPQ/AQ/aM/dr/aB/aR on ReadCount nodes using precomputed maps

#' Annotate ReadCount nodes
#'
#' Annotates `ReadCount` nodes under `Data/Sample` with computed attributes.
#' Attributes include `NPQ`, `aR` (raw read cutoff), `aB` (LOD flag), `AQ`, `aM`, and `dr`.
#' The `aB` flag ("above Limit of Detection") is written per target–sample
#' measurement, reflecting whether that specific measurement is above LOD.
#'
#' @param root XML root node.
#' @param data Data list from `loadNULISAseq()` providing matrices and parameters.
#' @param rawReadCutoff Integer raw read threshold for `aR`.
#'
#' @return Invisibly returns NULL; mutates `root` with attributes on `ReadCount` nodes.
#' @noRd
#' @keywords internal
#'
#' @importFrom xml2 xml_find_all xml_attr xml_set_attr xml_text
.annotate_readcounts <- function(root, data, rawReadCutoff){
  sampleNameByBarcode <- stats::setNames(data$samples$sampleName, data$samples$sampleBarcode)
  targetNameByBarcode <- stats::setNames(data$targets$targetName, data$targets$targetBarcode)
  for(sample in xml2::xml_find_all(root, './/Data//Sample')){
    sBarcode <- xml2::xml_attr(sample, "barcode")
    sName <- sampleNameByBarcode[[sBarcode]]
    if(is.null(sName) || is.na(sName)) next
    for(readcount in xml2::xml_find_all(sample, "ReadCount")){
      tBarcode <- xml2::xml_attr(readcount, "target")
      tName <- targetNameByBarcode[[tBarcode]]
      if(is.null(tName) || is.na(tName)) next

      # NPQ
      if(!is.null(rownames(data$NPQ)) && !is.null(colnames(data$NPQ)) &&
         tName %in% rownames(data$NPQ) && sName %in% colnames(data$NPQ)){
        xml2::xml_set_attr(readcount, "NPQ", as.character(data$NPQ[tName, sName]))
      }

      # aR based on raw read cutoff
      value <- as.integer(xml2::xml_text(readcount)) > rawReadCutoff
      xml2::xml_set_attr(readcount, "aR", as.character(as.integer(value)))

      # aB (LOD above flag) — per target–sample measurement
      if(!is.null(data$lod$aboveLOD) &&
         !is.null(rownames(data$lod$aboveLOD)) && !is.null(colnames(data$lod$aboveLOD)) &&
         tName %in% rownames(data$lod$aboveLOD) && sName %in% colnames(data$lod$aboveLOD)){
        xml2::xml_set_attr(readcount, "aB", as.character(as.integer(data$lod$aboveLOD[tName, sName])))
      }

      # AQ / aM / dr for absolute assay
      if(isTRUE(data$AbsAssay)){
        if(!is.null(data$AQ$Data_AQ) && tName %in% rownames(data$AQ$Data_AQ) && sName %in% colnames(data$AQ$Data_AQ)){
          xml2::xml_set_attr(readcount, "AQ", as.character(data$AQ$Data_AQ[tName, sName]))
        }
        if(!is.null(data$AQ$Data_AQ_aM) && tName %in% rownames(data$AQ$Data_AQ_aM) && sName %in% colnames(data$AQ$Data_AQ_aM)){
          xml2::xml_set_attr(readcount, "aM", as.character(data$AQ$Data_AQ_aM[tName, sName]))
        }
        if(!is.null(data$AQ$withinDR) && tName %in% rownames(data$AQ$withinDR) && sName %in% colnames(data$AQ$withinDR)){
          xml2::xml_set_attr(readcount, "dr", as.character(as.integer(data$AQ$withinDR[tName, sName])))
        }
      }
    }
  }
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
writeUpdatedXML <- function(input_XML, rawReadCutoff = 200, data = NULL){
  doc <- xml2::read_xml(input_XML)
  if(is.null(data)){
    data <- loadNULISAseq(input_XML)
  }
  doc <- .ensure_xml_version(doc)
  root <- xml2::xml_root(doc)

  .add_qc_thresholds(root, data)
  .add_plate_qc(root, data)
  .add_sample_qc(root, data)
  .add_target_qc(root, data)
  .annotate_readcounts(root, data, rawReadCutoff)

  # Preserve existing API: return XML as character string
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
