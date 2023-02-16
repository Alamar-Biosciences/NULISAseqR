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
