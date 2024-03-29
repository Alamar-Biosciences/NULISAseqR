% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intraPlateNorm.R
\name{intraPlateNorm}
\alias{intraPlateNorm}
\title{NULISAseq Intra-Plate Normalization}
\usage{
intraPlateNorm(
  data_matrix,
  method = "IC",
  IC = NULL,
  NC_wells = NULL,
  TC_omit = NULL,
  scaleFactor = 1
)
}
\arguments{
\item{data_matrix}{The Data matrix output from readNULISAseq.R.}

\item{method}{\code{'IC'} (internal control; default) divides 
each target count by the specified IC count for a given sample well; 
if multiple ICs are specified uses geometric mean of several ICs.
\code{'TC'} (total count)  scales each sample well so data_matrix sum to 
one million (except NC wells and omitted target data_matrix).}

\item{IC}{Required for 'single' and 'geom_mean' methods. 
Vector of either the row(s) (numeric) or the 
row name(s) (character string) of the internal control(s) to be used 
for normalization. Raw data_matrix will be divided by IC data_matrix or 
by the geometric mean of IC data_matrix.}

\item{NC_wells}{Recommended for 'total_count' method when negative control 
wells are present. Vector of either 
the column indices (numeric) or the column names (character string) 
for the negative control wells. For total count method, 
NCs will be scaled such
that the (total non-NC count : total NC count) ratio is the same 
in both the un-normalized and the normalized data.}

\item{TC_omit}{Option for 'total_count' method only.
Vector of either the row indices (numeric) or the 
row names (character string) of any targets, if any, 
that should be omitted
from the total count normalization. For example, it may be desirable to
omit IC if IC counts are high. 
Omitted targets will be rescaled but will not be used in calculating
the scaling factors.}

\item{scaleFactor.}{Optional numeric value used to rescale data 
after normalizing. Default is 1. This may be desirable if there are 
many targets with read counts less than the IC counts, to avoid
normalized quantities between 0 and 1 (which will be negative
in the log scale).}

\item{normData}{A matrix of normalized count data (not log-transformed).}

\item{normFactor}{A vector of the normalization factors used for 
each sample (column) in the data matrix.}
}
\value{
A list.
}
\description{
Does intra-plate normalization either using single internal control, 
multiple internal controls (using geometric mean), or total data_matrix. 
Input is the output of readNULISAseq.R Data matrix.
}
