#' Write Processed NULISAseq XML
#'
#' Writes NULISAseq XML file with unnormalized / normalized data and QC flags
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
#' @return NULL
#'
#' @examples
#' writeXML('filename.xml')
#'
#' @export
#'
writeProcessedXML <- function(in_xml_file, out_xml_file){
  c(plateID, ExecutionDetails, RunSummary, targets, samples, Data) %<-%  
    readNULISAseq(in_xml_file,
                  plateID=NULL,
                  file_type='xml_no_mismatches')
  write(capture.output(targets), stdout())
}
