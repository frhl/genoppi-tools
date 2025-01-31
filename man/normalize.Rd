% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize.R
\name{normalize}
\alias{normalize}
\title{Normalize data}
\usage{
normalize(table = NULL, type = "median")
}
\description{
Normalize a table of label or label free data. 
The function will calculate the \code{type} e.g. median of
all numeric columns and subtract the median from the 
intensity value.
}
\author{
flassen
}
