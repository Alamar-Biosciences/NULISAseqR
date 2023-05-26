# =========================================================================
# Copyright © 2023 Alamar Biosciences, Inc.
# NULISAseqR Data Analysis and Reporting API
# =========================================================================

##################
# NULISAseqR API #
##################

# Define global vars and other options
rmd_path <- "/workingDir/NULISAseqR/inst/rmarkdown/templates/nulisaseq/skeleton/skeleton.Rmd"
options(pagedown.remote.maxattempts = 40)
options(pagedown.remote.sleeptime = 2)

#* @apiTitle NULISAseqR Analysis/Reporting
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

#* To check the health of the API server, do nothing -> Returns Ok
#*
#* @get /healthCheck
#* @serializer contentType list(type="text/plain; charset=UTF-8")
function(res) {
  future({
    res$status <- 200
    return("Ok")
    })
}

#* @param in_xml:file Character string. Path and name of the file.
#* @param IPC Name to search for Interprocess control (IPC) samples
#* @param NC Name to search for Negative control (NC) samples
#* @param IC Name to search for Internal Control (IC) targets
#* @param barcodeB:file optional BarcodeB file
#* @serializer text
#* @post /normXML
normXML <- function(in_xml,
                    IPC = c("InterProcessControl"),
                    NC = c("NegativeControl"),
                    IC = c("mCherry"),
                    barcodeB = "") {
  promises::future_promise({
    return(processXML(toString(in_xml), IPC, NC, IC, toString(barcodeB)))
  })
}

#* @param in_xml:[file] Character string vector. Path(s) and name(s) of the file(s).
#* @param IPC Name to search for Interprocess control (IPC) samples
#* @param NC Name to search for Negative control (NC) samples
#* @param IC Name to search for Internal Control (IC) targets
#* @param study_name Name of the study
#* @param assayName Name of the assay
#* @param reportType Type of the report. Options: "WebApp", "internal".
#* @param excludeSamples Sample barcodes to be excluded from analysis
#* @param outputPlots Output HTML document and required plots for slides. Returned will be a zip file.
#* @param sampleGroupCovar Covariate to retrieve sample group information
#* @post /xml2html
xml2html <- function(res,
                     in_xml,
                     IPC = c("InterProcessControl"),
                     NC = c("NegativeControl"),
                     IC = c("mCherry"),
                     study_name = "Study Name",
                     assayName = "NULISAseq 200-plex Inflammation Panel",
                     reportType = "WebApp",
                     excludeSamples = NULL,
                     outputPlots = FALSE,
                     sampleGroupCovar = NULL) {
  promises::future_promise({
    UUID <- uuid::UUIDgenerate()

    # Create required files
    out_html <- paste0(UUID, ".html")
    tempFile <- paste0(UUID, ".Rmd")
    file.copy(rmd_path, tempFile)

    # Convert to vectors: If a comma-separated string encounters
    excludeSamples <- if (!is.null(excludeSamples)) unlist(strsplit(excludeSamples, "\\s*,\\s*"))
    IPC <- unlist(strsplit(IPC, "\\s*,\\s*"))
    NC <- unlist(strsplit(NC, "\\s*,\\s*"))
    IC <- unlist(strsplit(IC, "\\s*,\\s*"))
    outputPlots <- as.logical(outputPlots)

    # Handle multiple XML inputs
    xml_files_list <- lapply(in_xml, toString)
    xml_files_vec <- as.character(unlist(xml_files_list))

    rmarkdown::render(tempFile,
                      output_format = "html_document",
                      output_file = out_html,
                      params = list(xmlFiles = xml_files_vec,
                                    dataDir = NULL,
                                    reportType = reportType,
                                    IPC = IPC,
                                    NC = NC,
                                    IC = IC,
                                    study_name = study_name,
                                    assayName = assayName,
                                    excludeSamples = excludeSamples,
                                    outputPlots = outputPlots,
                                    sampleGroupCovar = sampleGroupCovar))

    # Handle if output plots are also requested: Return a zip containing HTML+figures
    if (outputPlots) {
      out_zip <- paste0(UUID, ".zip")

      # Add plots to the zip archive
      plot_files <- list.files("./figures/", full.names = TRUE)
      for (file in plot_files) {
        zip(out_zip, file)
      }

      # Add HTML file to the zip archive
      zip(out_zip, out_html)

      # Process and prepare zip file to return
      if (file.exists(out_zip)) {
        cat(paste0("Zip file created: ", out_zip))
        bin <- readBin(out_zip, "raw", file.info(out_zip)$size)
        unlink(out_zip)
      } else {
        stop("Failed to generate the zip file!")
      }

      # Remove the figures directory: If present
      if (file.exists("./figures") && file.info("./figures")$isdir) {
        unlink("./figures", recursive = TRUE)
      }

    } else {
      bin <- readBin(out_html, "raw", n = file.info(out_html)$size)
    }

    unlink(c(tempFile, out_html))

    base64_content <- base64enc::base64encode(bin)

    return(base64_content)
  })
}

#* @param in_xml:[file] Character string. Path and name of the file.
#* @param IPC Name to search for Interprocess control (IPC) samples
#* @param NC Name to search for Negative control (NC) samples
#* @param IC Name to search for Internal Control (IC) targets
#* @param study_name Name of the study
#* @param assayName Name of the assay
#* @serializer contentType list(type="application/pdf")
#* @post /xml2pdf
xml2pdf <- function(res,
                     in_xml,
                     IPC = c("InterProcessControl"),
                     NC = c("NegativeControl"),
                     IC = c("mCherry"),
                     study_name = "Study Name",
                     assayName = "NULISAseq 200-plex Inflammation Panel") {
   promises::future_promise({
     UUID <- uuid::UUIDgenerate()
     tempFile <- paste0(UUID, ".Rmd")
     outFile <- paste0(UUID, ".html")
     outFilePDF <- paste0(UUID, ".pdf")
#     file.copy(rmd_path, tempFile)

#     # Handle multiple XML inputs
#     xml_files_list <- lapply(in_xml, toString)
#     xml_files_vec <- as.character(unlist(xml_files_list))

#     rmarkdown::render(tempFile,
#                       output_format = "html_document",
#                       output_file = outFile,
#                       params = list(xmlFiles = xml_files_vec,
#                                     dataDir = NULL,
#                                     reportType = "WebApp",
#                                     IPC = IPC,
#                                     NC = NC,
#                                     IC = IC,
#                                     study_name = study_name,
#                                     assayName = assayName))
#     unlink(tempFile)
#     chrome_path <- Sys.which("chromium-browser")
#     pagedown::chrome_print(outFile, outFilePDF, extra_args = c("--no-sandbox", "--disable-gpu"), browser = chrome_path, verbose = 2) # nolint
#     unlink(outFile)
     pdf(outFilePDF)
     plot(c(1,2,3,4,5,6,7,8))
     dev.off()

     bin <- readBin(outFilePDF, "raw", n = file.info(outFilePDF)$size)
     unlink(outFilePDF)

     base64_content <- base64enc::base64encode(bin)

     return(base64_content)
   })
 }

#* @param in_xml:[file] Character string vector. Path(s) and name(s) of the file(s).
#* @param target_info_file:file Target information file. Note that this is version specific!
#* @param PanelLotNumber Panel lot number
#* @param sample_info_file:file Sample information file. Note that this is experiment specific and has sampleName and plateID columns! # nolint
#* @param panel Name of the panel. Defaults to "200-plex Inflammation v1".
#* @param ICs Vector of string(s). Internal control names. Default is "mCherry".
#* @param SC_string Vector of character string(s) that represents SCs in the column. Default is "SC".
#* @param excludeSamples Sample barcodes to be excluded from analysis.
#* @serializer contentType list(type="text/csv")
#* @post /xml2counts
xml2counts <- function(res,
                       in_xml,
                       target_info_file,
                       PanelLotNumber,
                       sample_info_file = NULL,
                       panel = "200-plex Inflammation v1",
                       ICs = c("mCherry"),
                       SC_string = c("SC"),
                       excludeSamples = c(NULL)) {
  promises::future_promise({
    # Temp output csv file
    UUID <- uuid::UUIDgenerate()
    outFile <- paste0(UUID, ".csv")

    # Convert to vectors: If a comma-separated string encounters
    excludeSamples <- if (!is.null(excludeSamples)) unlist(strsplit(excludeSamples, "\\s*,\\s*"))
    ICs <- unlist(strsplit(ICs, "\\s*,\\s*"))
    SC_string <- unlist(strsplit(SC_string, "\\s*,\\s*"))

    # Write XML content to temporary files
    xml_files <- NULL
    for (i in seq_along(in_xml)) {
      xml_files[i] <- tempfile(tmpdir = ".", fileext = ".xml")
      writeLines(in_xml[[i]], xml_files[i])
    }

    # Write target information and sample information files to temp files
    target_info_file_local <- tempfile(tmpdir = ".", fileext = ".csv")
    writeLines(toString(target_info_file), target_info_file_local)

    if (!is.null(sample_info_file)) {
      sample_info_file_local <- tempfile(tmpdir = ".", fileext = ".csv")
      writeLines(toString(sample_info_file), sample_info_file_local)
    } else {
      sample_info_file_local <- NULL
    }

    # Generate counts output file
    NULISAseqR::writeNULISAseq(xml_files = xml_files,
                               dataDir = ".",
                               target_info_file = target_info_file_local,
                               sample_info_file = sample_info_file_local,
                               output_filename = outFile,
                               Panel = panel,
                               PanelLotNumber = PanelLotNumber,
                               ICs = ICs,
                               SC_string = SC_string,
                               excludeSamples = excludeSamples)
    bin <- readBin(outFile, "raw", n = file.info(outFile)$size)

    # Clean the workspace
    unlink(c(xml_files, target_info_file_local, sample_info_file_local, outFile))

    return(bin)
  })
}

#* @param in_xml:[file] Character string vector. Path(s) and name(s) of the file(s).
#* @param IPC Name to search for Interprocess control (IPC) samples
#* @param NC Name to search for Negative control (NC) samples
#* @param IC Name to search for Internal Control (IC) targets
#* @serializer contentType list(type="application/octet-stream")
#* @post /xml2qs
xml2qs <- function(res,
                     in_xml,
                     IPC = c("InterProcessControl"),
                     NC = c("NegativeControl"),
                     IC = c("mCherry")){
  promises::future_promise({

    UUID <- uuid::UUIDgenerate()
    tempFile <- paste0(UUID, ".qs")
    data <- loadNULISAseq(toString(in_xml), IPC, IC)
    qsave(data, tempFile)
    bin <- readBin(tempFile, "raw", n = file.info(tempFile)$size)
    
    # Clean the workspace
    unlink(c(tempFile))

    return(bin)
  })
}
