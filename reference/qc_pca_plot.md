# PCA plot for QC report

Performs and then plots a principal component analysis of sample
intensities. Samples are colored according to the groups argument, with
colors determined by the
[`colorGroup`](https://byrumlab.github.io/proteoDA/reference/colorGroup.md)
function. By default, uses only the 500 most variable proteins for the
analysis.

## Usage

``` r
qc_pca_plot(
  data,
  groups = NULL,
  sample_labels = colnames(data),
  top_proteins = 500,
  standardize = TRUE,
  pca_axes = c(1, 2)
)
```

## Arguments

- data:

  A data frame of intensity data, likely normalized. Rows should be
  proteins and columns should be samples.

- groups:

  A character or factor vector, listing the group(s) the samples belong
  to.

- sample_labels:

  Optional, a vector of sample labels to use. If not supplied, defaults
  to using the column names in the data.

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

## Value

A ggplot object of the plot.
