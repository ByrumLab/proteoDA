% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/limma_tables.R
\name{write_limma_excel}
\alias{write_limma_excel}
\title{Write an .xlsx file of limma results}
\usage{
write_limma_excel(
  filename,
  statlist,
  annotation,
  data,
  norm.method,
  pval_thresh,
  lfc_thresh,
  add_filter
)
}
\arguments{
\item{filename}{The filename of the Excel spreadsheet to be saved.}

\item{statlist}{A list of the per-contrast statistical results.}

\item{annotation}{A data frame of annotation date for each protein.}

\item{data}{A data frame containing the average expression data for each sample.}

\item{norm.method}{The method that was used to normalize the data for the
statistical model being output.}

\item{pval_thresh}{The p-value threshold that was used to determine significance.}

\item{lfc_thresh}{The logFC threshold that was used to determine significance.}

\item{add_filter}{Should per-column filters be added to the spreadsheet?}
}
\value{
Invisibly returns a list, where the first element is the filename
of the saved Excel spreadsheet and the second element is the openxlsx
workbook object.
}
\description{
Internal function for creating a nicely formatted Excel spreadsheet of
differential expression results.
}
\keyword{internal}
