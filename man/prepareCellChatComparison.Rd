% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{prepareCellChatComparison}
\alias{prepareCellChatComparison}
\title{Prepare CellChat List for Comparison}
\usage{
prepareCellChatComparison(
  object.list,
  comparison,
  slot.name = "netP",
  measure = "weight",
  thresh = 0.05
)
}
\arguments{
\item{object.list}{List of CellChat objects}

\item{comparison}{Vector of indices to compare}

\item{slot.name}{Slot name containing network information}

\item{measure}{Measure to use for comparison ("weight" or "count")}

\item{thresh}{Significance threshold}
}
\value{
List of prepared data for comparison
}
\description{
This function prepares a list of CellChat objects for comparison.
}
\keyword{internal}
