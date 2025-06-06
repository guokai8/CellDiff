% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{normalizeMatrix}
\alias{normalizeMatrix}
\title{Normalize Matrix by Maximum or Sum}
\usage{
normalizeMatrix(mat, method = c("max", "sum", "none"))
}
\arguments{
\item{mat}{Matrix to normalize}

\item{method}{Normalization method: "max", "sum", or "none"}
}
\value{
Normalized matrix
}
\description{
Helper function to normalize a matrix by its maximum value or sum.
}
\keyword{internal}
