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

#* @param in_xml:file Character string. Path and name of the file.
#* @param IPC Name to search for Interprocess control (IPC) samples
#* @param NC Name to search for Negative control (NC) samples
#* @param IC Name to search for Internal Control (IC) targets
#* @param output_format Output format to return (e.g. html_document, pdf_document)
#* @serializer html
#* @post /xml2report
xml2report <- function(res, in_xml, IPC=c("InterProcessControl"), NC=c("NegativeControl"), IC=c("mCherry"), output_format="html_document"){
  future_promise({
    UUID <- UUIDgenerate()
    tempFile <- paste0(UUID, '.Rmd')
    outFile <- paste0(UUID, if(output_format == "html_document") ".html" else ".pdf")
    file.copy('skeleton.Rmd', tempFile)
    rmarkdown::render('skeleton.Rmd', output_format=output_format, 
                                      output_file=outFile,
                                      params=list(
                                                  xmlFiles=toString(in_xml),
                                                  dataDir=NULL,
                                                  type="WebApp",
                                                  IPC=IPC, 
                                                  NC=NC, 
                                                  IC=IC
                                                  ))
    unlink(tempFile)
    readBin(outFile, "raw", n=file.info(outFile)$size)
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
