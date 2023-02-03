#!/usr/bin/env Rscript
library(plumber)
library(optparse)
library(future)
library(NULISAseqR)
args = commandArgs(trailingOnly=T)
PORT <- if (length(args) == 0 ) 8000 else args[1]
plumb_api(package="NULISAseqR", name="API_1", edit=FALSE) %>%
  pr_run(host = "0.0.0.0", port = as.integer(PORT))
