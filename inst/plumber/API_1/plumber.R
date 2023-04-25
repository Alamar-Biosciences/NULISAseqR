##################
# NULISAseqR API #
##################

#* @apiTitle NULISAseqR Reporting
#* @apiDescription RESTful API for NULISAseqR data analysis, reporting and related operations
#* @apiVersion 1.0.0

#* Print a log message
#*
#* @filter logger
function(req) {
    cat(as.character(Sys.time()), "-", req$REQUEST_METHOD, req$PATH_INFO, "-",
        req$HTTP_USER_AGENT, "@", req$REMOTE_ADDR, "\n")
    plumber::forward()
}

# Define Rmarkdown file path for endpoints
rmd_path <- "/workingDir/NULISAseqR/inst/rmarkdown/templates/nulisaseq/skeleton/skeleton.Rmd"

#* @param in_xml:file Character string. Path and name of the file.
#* @param IPC Name to search for Interprocess control (IPC) samples
#* @param NC Name to search for Negative control (NC) samples
#* @param IC Name to search for Internal Control (IC) targets
#* @param barcodeB:file optional BarcodeB file
#* @serializer text
#* @post /normXML
normXML <- function(in_xml, IPC=c("InterProcessControl"), NC=c("NegativeControl"), IC=c("mCherry"), barcodeB=""){
  future_promise({
    return(processXML(toString(in_xml), IPC, NC, IC, toString(barcodeB)))
  })
}

#* @param in_xml:[file] Character string vector. Path(s) and name(s) of the file(s).
#* @param IPC Name to search for Interprocess control (IPC) samples
#* @param NC Name to search for Negative control (NC) samples
#* @param IC Name to search for Internal Control (IC) targets
#* @serializer html
#* @post /xml2html
xml2html <- function(res, in_xml, IPC=c("InterProcessControl"), NC=c("NegativeControl"), IC=c("mCherry")){
  future_promise({
    UUID <- UUIDgenerate()
    tempFile <- paste0(UUID, '.Rmd')
    outFile <- paste0(UUID, ".html")
    file.copy(rmd_path, tempFile)

    # Handle multiple XML inputs
    xml_files_list <- lapply(in_xml, toString)
    xml_files_vec <- as.character(unlist(xml_files_list))

    rmarkdown::render(tempFile, output_format="html_document",
                                      output_file=outFile,
                                      params=list(
                                                  xmlFiles=xml_files_vec,
                                                  dataDir=NULL,
                                                  reportType="WebApp",
                                                  IPC=IPC,
                                                  NC=NC,
                                                  IC=IC
                                                  ))
    unlink(tempFile)
    readBin(outFile, "raw", n=file.info(outFile)$size)
  })
}

#* @param in_xml:file Character string. Path and name of the file.
#* @param IPC Name to search for Interprocess control (IPC) samples
#* @param NC Name to search for Negative control (NC) samples
#* @param IC Name to search for Internal Control (IC) targets
#* @serializer contentType list(type="application/pdf")
#* @post /xml2pdf
xml2pdf <- function(res, in_xml, IPC=c("InterProcessControl"), NC=c("NegativeControl"), IC=c("mCherry")){
  future_promise({
    UUID <- UUIDgenerate()
    tempFile <- paste0(UUID, '.Rmd')
    outFile <- paste0(UUID, ".html")
    outFilePDF <- paste0(UUID, ".pdf")
    file.copy(rmd_path, tempFile)
    rmarkdown::render(tempFile, output_format="html_document", 
                                      output_file=outFile,
                                      params=list(
                                                  xmlFiles=toString(in_xml),
                                                  dataDir=NULL,
                                                  reportType="WebApp",
                                                  IPC=IPC, 
                                                  NC=NC, 
                                                  IC=IC
                                                  ))
    unlink(tempFile)
    chrome_print(outFile, outFilePDF)
    unlink(outFile)
    readBin(outFilePDF, 'raw', n=file.info(outFilePDF)$size)
  })
}

#* To check the health of the API server, do nothing -> Returns Ok
#*
#* @get /healthCheck
#* @serializer contentType list(type="text/plain; charset=UTF-8")
function(res){
  future({
    res$status <- 200    
    return("Ok")
    })
}
