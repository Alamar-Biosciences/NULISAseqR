# NULISAseqR

![AWS CodeBuild](https://codebuild.us-west-1.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoib045RnFOTFB4Wmo5OHBDTGJySnNJK3dtN2I3a0MwQm96UVZyMnp1anl3cGZtMWs5dVowMVl5TVlLUEw4RnNiZWlscnNTdE5KV2xQSlVyN3YrZUVvYTZRPSIsIml2UGFyYW1ldGVyU3BlYyI6InNtclNBUGloQjJEdytnMUQiLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=main)
[![Deploy SAM](https://github.com/Alamar-Biosciences/NULISAseqR/actions/workflows/deploy_sam.yml/badge.svg)](https://github.com/Alamar-Biosciences/NULISAseqR/actions/workflows/deploy_sam.yml)

NULISAseq R package

## System Requirements:
1. R (version 4.3+)

## How to install

1. Run the following R code (typical installation time < 5 minutes):
```
    install.packages('devtools')
    devtools::install_github('Alamar-Biosciences/NULISAseqR',
                              ref = 'main'
                              )
```

2. Load the package in R with `library(NULISAseqR)`.

## Demo: Loading Data

```
    library('NULISAseqR'))
    data <- loadNULISAseq('<NULISAseqR Directory>/inst/rmarkdown/templates/nulisaseq/skeleton/detectability_P1_Tr03.xml')

```

## Demo: Generating Report

```
    rmarkdown::render("<NULISAseqR Directory>/inst/rmarkdown/templates/nulisaseq/skeleton/skeleton.Rmd", params=list(dataDir="<NULISAseqR Directory>/inst/rmarkdown/templates/nulisaseq/skeleton", xmlFiles=c("detectability_P1_Tr03.xml"")))
    [Report](/inst/rmarkdown/templates/nulisaseq/skeleton/skeleton.html)
```

