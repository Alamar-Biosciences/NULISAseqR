#' Read NULISAseq XML
#'
#' Reads NULISAseq XML file, where XML file is output from the Alamar
#' Biosciences Galaxy NULISAseq (Beta) tool or the NULISAseq Normalization (Alpha) tool.
#'
#' @param xml_file Path and name of the file.
#' @param file_type Type of input file, as output from Galaxy. Options include
#' xml_full_output, xml_no_mismatches (both from NULISAseq tool),
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
readNULISAseq <- function(xml_file, file_type){

}
