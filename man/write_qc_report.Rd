% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_report.R
\name{write_qc_report}
\alias{write_qc_report}
\title{Create a quality control report}
\usage{
write_qc_report(
  DAList,
  color_column = NULL,
  label_column = NULL,
  output_dir = NULL,
  filename = NULL,
  overwrite = FALSE,
  top_proteins = 500,
  standardize = TRUE,
  pca_axes = c(1, 2),
  dist_metric = "euclidean",
  clust_method = "complete",
  show_all_proteins = F
)
}
\arguments{
\item{DAList}{A DAList.}

\item{color_column}{The name of the column in the metadata which
gives information on how to color samples in plots within the report. If not
supplied, all samples will be the same color.}

\item{label_column}{Optional. The name of column within the targets data frame
which contains labels to use for plotting figures. When not supplied,
defaults to using the column names of the data in processed_data. To ensure
good plot formatting, sample labels will be truncated to 20 characters, and
the function will give an error if the sample labels are not unique.}

\item{output_dir}{The directory in which to save the report. If the directory
does not exist, it will be created. If not provided, will default to
the current working directory.}

\item{filename}{The file name of the report to be saved. Must end in .pdf. Will
default to "QC_Report.pdf" if no filename is provided.}

\item{overwrite}{Should report file be overwritten if it already exists?
Default is FALSE.}

\item{top_proteins}{The number of most variable proteins to use
for the PCA and dendrogram clustering. Default is 500.}

\item{standardize}{Should input data be standardized to a mean of 0 and std.dev of
1 before performing PCA and dendrogram clustering? If input data are
not yet standardized, should be TRUE. Default is TRUE.}

\item{pca_axes}{A numeric vector of length 2 which lists the PC axes to plot.
Default is c(1,2), to plot the first two principal components.}

\item{dist_metric}{The metric used to define distance for dendrogram clustering.
Default is "euclidean". See \code{\link[stats:dist]{stats::dist}} for options.}

\item{clust_method}{The agglomeration method to use for dendrogram clustering.
Default is "complete", See \code{\link[stats:hclust]{stats::hclust}} for options.}

\item{show_all_proteins}{Should all proteins be shown in missing value heatmap,
of only those with missing data? Default is F (only those with missing data).}
}
\value{
If report is created successfully, invisibly returns the input DAList.
}
\description{
Saves a PDF report containing a variety of plots which provide information on the
distribution, clustering, and correlation of protein intensities
across samples. See arguments for options for customizing the report.
}
\examples{
\dontrun{
# Color samples according to group identities
# in the "treatment" column of the metadata
write_qc_report(DAList,
                color_column = "treatment")

# Change the default directory and file names
write_qc_report(DAList,
                color_column = "treatment",
                output_dir = "my/chosen/directory",
                filename = "my_report.pdf")

# Overwrite an existing report
write_qc_report(DAList,
                color_column = "treatment",
                overwrite = T)

# Customize PCA and clustering plots
write_qc_report(DAList,
                color_column = "treatment",
                top_proteins = 1000,
                pca_aces = c(2,3),
                dist_metric = "manhattan",
                clust_method = "average")

}
}
