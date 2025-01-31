% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mttest.R
\name{mttest}
\alias{mttest}
\title{Calculate a moderated t-test}
\usage{
mttest(df, keep = c("imputed"))
}
\arguments{
\item{df}{a data.frame with 'rep' in column names.}

\item{keep}{a vector of characters containing colnames that should be kept in the data.frame}
}
\description{
Will do a moderated t.test pairwise on replicate columns, and
calculate the FDR and p-value.
}
\seealso{
Other genoppi: 
\code{\link{designate}()},
\code{\link{enrichment}()},
\code{\link{interactors}()},
\code{\link{plotOverlap}()},
\code{\link{plotScatter}()},
\code{\link{plotVolcano}()}
}
\author{
April Kim / Frederik Heymann
}
\concept{genoppi}
