% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare.R
\name{prepare}
\alias{prepare}
\title{pipeline for preparing genomics data for statistical analysis}
\usage{
prepare(
  bait,
  infile,
  cols = NULL,
  impute = list(stdwidth = 0.5, shift = -1.8),
  transform = "log2",
  normalization = "median",
  filter = "HUMAN",
  raw = F,
  firstcol = "gene",
  control = "mock",
  peptide.threshold = 2,
  filter.ignore = NULL,
  verbose = F
)
}
\arguments{
\item{bait}{a vector or string containing that should be matched in the columns containing the bait.}

\item{infile}{the file path or data.frame that contains the raw data, i.e. accession numbers, intensity values, ratios etc.}

\item{cols}{optional manual entry. A vector of columns that are present in the dataset header. Follows 
the format of columns: acession, bait1, mock1, bait2, mock2, bait3, ..}

\item{impute}{how should missing data be imputed? NULL means that missing data rows are dropped. Will take a list
containing shift and stdwidth for gaussian imputation. For further details, see \code{?impute}.}

\item{transform}{charcacter. an R-command for how the data should be transformed.}

\item{normalization}{character. an R-command for how the data should be collumn-wise transformned.}

\item{filter}{character. only accession IDs of the filter specified are included.}

\item{raw}{will return the data.frame alongside the raw intensity values.}

\item{firstcol}{will change the name of the first column to the string indicated}

\item{control}{the control or references. Normally bait vs control, but it could be bait in wildtype versus 
bait in mutant for differential expression.}

\item{peptide.threshold}{how many peptides must be dectected in the ms, in order to deem it a valid observation. Default is 2.}

\item{filter.ignore}{will try to match the inputted vector or character to acession IDs. If sucessful,
it will ignore further filtering of this item. This could for instance be used, if the bait only has
one unique protein, and would therefore otherwise be filtered.}

\item{verbose}{boolean. if true, returns the table and a list with updates.}
}
\value{
a table that can be inputted to genoppi
}
\description{
Process IP intensities by normalizing and calculating logfold change.
}
