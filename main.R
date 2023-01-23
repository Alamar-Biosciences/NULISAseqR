#!/usr/bin/env Rscript
library(plumber)
library(optparse)

pr() %>%
  pr_mount("/reporting", plumb("/workingDir/processXML.R")) %>%
  pr_run(host = "0.0.0.0", port = 8000)
