###################
# Setup Workspace #
###################

# Load other libraries required to run the API
library(logger)

# Define global vars and other options
options(pagedown.remote.maxattempts = 40)
options(pagedown.remote.sleeptime = 2)

##################
# NULISAseqR API #
##################

#* @apiTitle NULISAseqR Analysis/Reporting
#* @apiDescription RESTful API for NULISAseqR data analysis, reporting and related operations
#* @apiVersion 1.0.0
#* @apiContact list(name = "Alamar Bioscience Support", url = "https://alamarbio.com/contact-us/", email = "support@alamarbio.com") # nolint
#* @apiLicense list(name = "Apache 2.0", url = "https://www.apache.org/licenses/LICENSE-2.0.html")
#* @apiTag General General use endpoints
#* @apiTag Reporting Endpoints used to generate various encoded outputs from XML inputs

#* Print a log message
#*
#* @filter logger
function(req) {
  message(as.character(Sys.time()), "-", req$REQUEST_METHOD, req$PATH_INFO, "-",
          req$HTTP_USER_AGENT, "@", req$REMOTE_ADDR, "\n")
  plumber::forward()
}

#* To check the health of the API server, do nothing -> Returns Ok
#*
#* @get /healthCheck
#* @serializer contentType list(type="text/plain; charset=UTF-8")
#* @tag General
function(res) {
  future({
    res$status <- 200
    return("Ok")
  })
}

#######################################
# Endpoints related to data reporting #
#######################################

#* Generate XML with normalized data
#*
#* @param in_xml:file Character string. Path and name of the file.
#* @param IPC:[str] Name to search for Interprocess control (IPC) samples
#* @param NC:[str] Name to search for Negative control (NC) samples
#* @param IC:[str] Name to search for Internal Control (IC) targets
#* @param barcodeB:file optional BarcodeB file
#* @serializer text
#* @post /normXML
#* @tag Reporting
normXML <- function(in_xml,
                    IPC = c("InterProcessControl"),
                    NC = c("NegativeControl"),
                    IC = c("mCherry"),
                    barcodeB = "") {
  promises::future_promise({
    return(processXML(toString(in_xml), IPC, NC, IC, toString(barcodeB)))
  })
}

#* Generate NULISAseq HTML QC report from XML(s)
#*
#* @param in_xml:[file] Character string vector. Path(s) and name(s) of the file(s).
#* @param IPC:[str] Name to search for Interprocess control (IPC) samples
#* @param NC:[str] Name to search for Negative control (NC) samples
#* @param IC:[str] Name to search for Internal Control (IC) targets
#* @param study_name:str Name of the study
#* @param assayName:str Name of the assay
#* @param reportType:str Type of the report. Options: "WebApp", "internal".
#* @param excludeSamples:[str] Sample barcodes to be excluded from analysis
#* @param excludeTargets:[str] Targets to be excluded from analysis
#* @param outputPlots:bool Output HTML document and required plots for slides. Returned will be a zip file.
#* @param sampleGroupCovar:str Covariate to retrieve sample group information
#* @post /xml2html
#* @tag Reporting
xml2html <- function(res,
                     in_xml,
                     IPC = NULL,
                     NC = NULL,
                     IC = c("mCherry"),
                     study_name = "Study Name",
                     assayName = "NULISAseq 200-plex Inflammation Panel",
                     reportType = "WebApp",
                     excludeSamples = NULL,
                     excludeTargets = NULL,
                     outputPlots = FALSE,
                     sampleGroupCovar = NULL) {
  promises::future_promise({

    tryCatch(
      {

        # Determine the Rmd file path
        current_wd <- getwd()
        rmd_path_suffix <- "/inst/rmarkdown/templates/nulisaseq/skeleton/skeleton.Rmd"
        if (file.exists(paste0(current_wd, rmd_path_suffix))) {
          rmd_path <- paste0(current_wd, rmd_path_suffix)
        } else {
          rmd_path <- paste0(current_wd, "/NULISAseqR", rmd_path_suffix)
          rmd_path <- gsub("NULISAseqR/inst/plumber/API_1/", "", rmd_path)
        }
        log_info(paste("Rmd file found at:", rmd_path))

        # Create required temporary files
        temp_rmd <- tempfile(pattern = paste0(Sys.Date(), "-skeleton-"), fileext = ".rmd")
        out_html <- tempfile(pattern = paste0(Sys.Date(), "-qc-report-"), fileext = ".html")
        out_zip <- tempfile(pattern = paste0(Sys.Date(), "-output-files-"), fileext = ".zip")

        # Clean workspace
        temp_files_to_delete <- c(temp_rmd, out_html, out_zip)
        on.exit({
          for (temp_file in temp_files_to_delete) {
            if (file.exists(temp_file)) {
              unlink(temp_file)
            }
          }
        }, add = TRUE)

        # Make a copy of the Rmd
        file.copy(rmd_path, temp_rmd)
        log_info("Successfully generated a copy of the Rmd.")

        # Handle multiple XML inputs
        xml_files_list <- lapply(in_xml, toString)
        xml_files_vec <- as.character(unlist(xml_files_list))

        # excludeSamples and excludeTargets require repetition of the sample name/target values
        # Note that these samples/targets gets excluded from all plates
        excludeSamples <- if (!is.null(excludeSamples)) rep(list(excludeSamples), length(xml_files_vec))
        excludeTargets <- if (!is.null(excludeTargets)) rep(list(excludeTargets), length(xml_files_vec))

        # Try generating the HTML QC report + output plots/files
        rmarkdown::render(temp_rmd,
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
                                        excludeTargets = excludeTargets,
                                        outputPlots = as.logical(outputPlots),
                                        sampleGroupCovar = sampleGroupCovar))
        log_info("Successfully generated the HTML QC report.")

        # Handle if output plots/files are also requested
        # Return a zip containing HTML+figures
        # If not, just the HTML
        if (outputPlots) {

          # Add plots to the zip archive
          plot_files <- list.files("./outputFiles/", full.names = TRUE)
          for (file in plot_files) {
            zip(out_zip, file)
          }

          # Add HTML file to the zip archive
          zip(out_zip, out_html)

          # Process and prepare zip file to return
          if (file.exists(out_zip)) {
            log_info(paste0("Zip file created: ", out_zip))
            bin <- readBin(out_zip, "raw", file.info(out_zip)$size)
          } else {
            stop("Failed to generate the zip file!")
          }

          # Remove the figures/files directory: If present
          if (file.exists("./outputFiles") && file.info("./outputFiles")$isdir) {
            unlink("./outputFiles", recursive = TRUE)
          }

          log_info("Generated a zip archive with required additional plots and files.")

        } else {
          bin <- readBin(out_html, "raw", n = file.info(out_html)$size)
        }

        # Base64 encoding of the binary output
        base64_content <- base64enc::base64encode(bin)

        return(base64_content)

      },
      error = function(err) {
        log_error(paste("xml2html endpoint failed. See the error: ", err))
        res$status <- 500
        res$body <- paste("xml2html endpoint failed. See the error: ", err)
      }
    )
  })
}

#* @param in_xml:[file] Character string. Path and name of the file.
#* @param IPC:[str] Name to search for Interprocess control (IPC) samples
#* @param NC:[str] Name to search for Negative control (NC) samples
#* @param IC:[str] Name to search for Internal Control (IC) targets
#* @param study_name Name of the study
#* @param assayName Name of the assay
#* @serializer contentType list(type="application/pdf")
#* @post /xml2pdf
#* @tag Reporting
xml2pdf <- function(res,
                     in_xml,
                     IPC = c("InterProcessControl"),
                     NC = c("NegativeControl"),
                     IC = c("mCherry"),
                     study_name = "Study Name",
                     assayName = "NULISAseq 200-plex Inflammation Panel") {
   promises::future_promise({
     UUID <- uuid::UUIDgenerate()
     temp_rmd <- paste0(UUID, ".Rmd")
     outFile <- paste0(UUID, ".html")
     outFilePDF <- paste0(UUID, ".pdf")
#     file.copy(rmd_path, temp_rmd)

#     # Handle multiple XML inputs
#     xml_files_list <- lapply(in_xml, toString)
#     xml_files_vec <- as.character(unlist(xml_files_list))

#     rmarkdown::render(temp_rmd,
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
#     unlink(temp_rmd)
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

#* Generate NULISAseq NPQ counts reports from XML(s)
#*
#* @param in_xml:[file] Character string vector. Path(s) and name(s) of the file(s).
#* @param target_info_file:file Target information file. Note that this is version specific!
#* @param PanelLotNumber:str Panel lot number
#* @param sample_info_file:file Sample information file. Note that this is experiment specific and has sampleName and plateID columns!
#* @param panel:str Name of the panel. Defaults to "200-plex Inflammation v1".
#* @param ICs:[str] Vector of string(s). Internal control names. Default is "mCherry".
#* @param SC_string:[str] Vector of character string(s) that represents SCs in the column. Default is "SC".
#* @param include_IPC:bool Logical value indicating whether IPC samples be included in output. Default is `FALSE`.
#* @param include_NC:bool Logical value indicating whether NC samples be included in output. Default is `FALSE`.
#* @param include_unnorm_counts:bool Logical value indicating whether unnormalized counts be included as an additional column in output. Default is `FALSE`.
#* @param include_IC_counts:bool Logical value indicating whether IC counts be included in output. Default is `FALSE`.
#* @param interPlateNorm_transformReverse_covariateName:str Name of the covariate specifying curve quantification method. Default is `Curve_Quant`.
#* @param excludeSamples:[str] Sample barcodes to be excluded from analysis.
#* @serializer contentType list(type="text/csv")
#* @post /xml2counts
#* @tag Reporting
xml2counts <- function(res,
                       in_xml,
                       target_info_file,
                       PanelLotNumber,
                       sample_info_file = NULL,
                       panel = "200-plex Inflammation v1",
                       ICs = c("mCherry"),
                       SC_string = c("SC"),
                       include_IPC = FALSE,
                       include_NC = FALSE,
                       include_unnorm_counts = FALSE,
                       include_IC_counts = FALSE,
                       interPlateNorm_transformReverse_covariateName = "Curve_Quant",
                       excludeSamples = c(NULL)) {
  promises::future_promise({

    tryCatch(
      {
        # Create required temporary files
        out_csv <- tempfile(pattern = paste0(Sys.Date(), "-count-report-"), tmpdir = ".", fileext = ".csv")

        # Write XML content to temporary files
        xml_files <- NULL
        for (i in seq_along(in_xml)) {
          xml_files[i] <- tempfile(tmpdir = ".", fileext = ".xml")
          writeLines(in_xml[[i]], xml_files[i])
        }
        log_info("Local copies of input XMLs created.")

        # excludeSamples and excludeTargets require repetition of the sample name/target values
        # Note that these samples/targets gets excluded from all plates
        excludeSamples <- if (!is.null(excludeSamples)) rep(list(excludeSamples), length(in_xml))

        # Write target information to a temp file
        target_info_file_local <- tempfile(pattern = paste0(Sys.Date(), "-target-info-"), tmpdir = ".", fileext = ".csv")
        writeLines(toString(target_info_file), target_info_file_local)
        log_info("Local copy of target information file created.")

        # Write sample information to a temp file
        if (!is.null(sample_info_file)) {
          sample_info_file_local <- tempfile(pattern = paste0(Sys.Date(), "-sample-info-"), tmpdir = ".", fileext = ".csv")
          writeLines(toString(sample_info_file), sample_info_file_local)
          log_info("Local copy of sample information created.")
        } else {
          sample_info_file_local <- NULL
        }

        # Clean workspace
        temp_files_to_delete <- c(out_csv, xml_files, target_info_file_local, sample_info_file_local)
        on.exit({
          for (temp_file in temp_files_to_delete) {
            if (file.exists(temp_file)) {
              unlink(temp_file)
            }
          }
        }, add = TRUE)

        # Generate counts output file
        NULISAseqR::writeNULISAseq(xml_files = xml_files,
                                   dataDir = ".",
                                   target_info_file = target_info_file_local,
                                   sample_info_file = sample_info_file_local,
                                   output_filename = out_csv,
                                   Panel = panel,
                                   PanelLotNumber = PanelLotNumber,
                                   ICs = ICs,
                                   SC_string = SC_string,
                                   include_IPC = as.logical(include_IPC),
                                   include_NC = as.logical(include_NC),
                                   include_unnorm_counts = as.logical(include_unnorm_counts),
                                   include_IC_counts = as.logical(include_IC_counts),
                                   interPlateNorm_transformReverse_covariateName = interPlateNorm_transformReverse_covariateName,
                                   excludeSamples = excludeSamples)
        log_info("Count report generation is complete.")
        bin <- readBin(out_csv, "raw", n = file.info(out_csv)$size)

        return(bin)
      },
      error = function(err) {
        log_error(paste("xml2counts endpoint failed. See the error: ", err))
        res$status <- 500
        res$body <- paste("xml2counts endpoint failed. See the error: ", err)
      }
    )
  })
}

#* Convert XML to QS format
#*
#* @param in_xml:[file] Character string vector. Path(s) and name(s) of the file(s).
#* @param IPC:[str] Name to search for Interprocess control (IPC) samples
#* @param NC:[str] Name to search for Negative control (NC) samples
#* @param IC:[str] Name to search for Internal Control (IC) targets
#* @post /xml2qs
#* @tag Reporting
xml2qs <- function(res,
                   in_xml,
                   IPC = c("InterProcessControl"),
                   NC = c("NegativeControl"),
                   IC = c("mCherry")) {
  promises::future_promise({

    tryCatch(
      {
        # Create required temporary files
        out_qs <- tempfile(pattern = paste0(Sys.Date(), "-qs-file-"), fileext = ".qs")

        # Clean workspace
        on.exit({
          if (file.exists(out_qs)) {
            unlink(out_qs)
          }
        }, add = TRUE)

        data <- loadNULISAseq(toString(in_xml), IPC, IC)
        qs::qsave(data, out_qs)
        log_info("Successfully written data to qs format.")
        bin <- readBin(out_qs, "raw", n = file.info(out_qs)$size)

        return(bin)
      },
      error = function(err) {
        log_error(paste("xml2qs endpoint failed. See the error: ", err))
        res$status <- 500
        res$body <- paste("xml2qs endpoint failed. See the error: ", err)
      }
    )
  })
}
