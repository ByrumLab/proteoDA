# Create a quality control report

Saves a PDF report containing a variety of plots which provide
information on the distribution, clustering, and correlation of protein
intensities across samples. See arguments for options for customizing
the report.

## Usage

``` r
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
```

## Arguments

- DAList:

  A DAList.

- color_column:

  The name of the column in the metadata which gives information on how
  to color samples in plots within the report. If not supplied, all
  samples will be the same color.

- label_column:

  Optional. The name of column within the targets data frame which
  contains labels to use for plotting figures. When not supplied,
  defaults to using the column names of the data in processed_data. To
  ensure good plot formatting, sample labels will be truncated to 20
  characters, and the function will give an error if the sample labels
  are not unique.

- output_dir:

  Directory to save the report (created if missing). Defaults to the
  current working directory.

- filename:

  The file name of the report to be saved. Must end in .pdf. Will
  default to "QC_Report.pdf" if no filename is provided.

- overwrite:

  Overwrite if the report already exists? Default `FALSE`.

- top_proteins:

  The number of most variable proteins to use for the PCA and dendrogram
  clustering. Default is 500.

- standardize:

  Should input data be standardized to a mean of 0 and std.dev of 1
  before performing PCA and dendrogram clustering? If input data are not
  yet standardized, should be TRUE. Default is TRUE.

- pca_axes:

  A numeric vector of length 2 which lists the PC axes to plot. Default
  is c(1,2), to plot the first two principal components.

- dist_metric:

  The metric used to define distance for dendrogram clustering. Default
  is "euclidean". See [`stats::dist`](https://rdrr.io/r/stats/dist.html)
  for options.

- clust_method:

  The agglomeration method to use for dendrogram clustering. Default is
  "complete", See [`stats::hclust`](https://rdrr.io/r/stats/hclust.html)
  for options.

- show_all_proteins:

  Should all proteins be shown in missing value heatmap, of only those
  with missing data? Default is F (only those with missing data).

## Value

If report is created successfully, invisibly returns the input DAList.

## Examples

``` r
if (FALSE) { # \dontrun{
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

} # }
```
