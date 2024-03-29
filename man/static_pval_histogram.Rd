% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/limma_report.R
\name{static_pval_histogram}
\alias{static_pval_histogram}
\title{Make a p-value histogram plot}
\usage{
static_pval_histogram(data, contrast)
}
\arguments{
\item{data}{Per-contrast DE results to be plotted, as prepared by
\code{\link{prep_plot_model_data}}.}

\item{contrast}{The contrast being plotted. Used for generating the plot title.}
}
\value{
A ggplot object.
}
\description{
Internal function for plotting p-value histograms of raw- and adjusted p-values.
}
\keyword{internal}
