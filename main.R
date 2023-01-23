#!/usr/bin/env Rscript
library(plumber)
library(optparse)

pr() %>%
  pr_mount("/", plumb("/workingDir/NULISAseqR/R/processXML.R")) %>%
  pr_run(host = "0.0.0.0", port = 8000)
