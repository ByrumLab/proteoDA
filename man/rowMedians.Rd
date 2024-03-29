% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{rowMedians}
\alias{rowMedians}
\title{Calculate per-row medians of a numeric array}
\usage{
rowMedians(x, ...)
}
\arguments{
\item{x}{The array for which to calculate per-row medians}

\item{...}{Additional arguments to be passed to internal functions.
Meant for na.rm.}
}
\value{
A numeric vector of appropriate length, named if input was named, with
per-row medians
}
\description{
Have done testing, gives same results as the matrixStats::rowMedians() function
it replaces, though it is much slower (the matrixStats version uses C code).
}
\examples{
# No examples yet.
}
\keyword{internal}
