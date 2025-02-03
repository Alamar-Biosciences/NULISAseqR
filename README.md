# NULISAseqR

NULISAseq R package v1.2

## System Requirements:
1. R (version 4.4+)

## How to install

1. Run the following R code (typical installation time < 5 minutes):
```
    install.packages('devtools')
    devtools::install_github('Alamar-Biosciences/NULISAseqR',
                              ref = 'main'
                            )
```

2. Load the package in R with `library(NULISAseqR)`.

## Demo: Loading Data (XML)

```
    # Load the NULISAseqR library
    library(NULISAseqR)

    # This is an input XML file example that is included with NULISAseqR. You could
    # alternatively change the variable inputFile to be one of your own XML files,
    # e.g. inputFile <- "filename.xml"
    inputFile <- paste0(
                         dirname(system.file("rmarkdown/templates/nulisaseq/skeleton", "skeleton.Rmd", package="NULISAseqR")),
                         "/detectability_P1_Tr03.xml"
                       )

    # Load the XML data
    dataXML <- loadNULISAseq(inputFile)

```
## Demo: Loading Data (XLSX)

```
    # Load the NULISAseqR library
    library(NULISAseqR)

    # Load the XLSX data where "XLSX inputFile" is the name of your XLSX file
    dataXLSX <- readNULISAseq("XLSX inputFile", file_type="xlsx") # example of loading an XLSX file

```

## Demo: Generating Report

Run the following to generate a QC report ( < 5 minutes)
```
  # Location of QC report template included with NULISAseqR
  template <- system.file("rmarkdown/templates/nulisaseq/skeleton", "skeleton.Rmd", package="NULISAseqR")

  # Command to create the QC report (Note that QC reports can only be created using XML files)
  rmarkdown::render(
                      template, 
                      output_file = "~/output.html",
                      params = list(
                        dataDir = dirname(template), 
                        xmlFiles = c("detectability_P1_Tr03.xml")
                      )
                   )

```

