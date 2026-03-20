# Dendrogram for QC report

Performs and then plots a hierarchical clustering analysis of sample
intensities across samples. Sample labels are colored according to the
"groups" argument, with colors determined by the
[`colorGroup`](https://byrumlab.github.io/proteoDA/reference/colorGroup.md)
function. By default, uses only the 500 most variable proteins for the
analysis.

## Usage

``` r
qc_dendro_plot(
  data,
  groups = NULL,
  sample_labels = NULL,
  top_proteins = 500,
  standardize = TRUE,
  dist_metric = "euclidean",
  clust_method = "complete"
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

- dist_metric:

  The metric used to define distance for dendrogram clustering. Default
  is "euclidean". See [`stats::dist`](https://rdrr.io/r/stats/dist.html)
  for options.

- clust_method:

  The agglomeration method to use for dendrogram clustering. Default is
  "complete", See [`stats::hclust`](https://rdrr.io/r/stats/hclust.html)
  for options.

## Value

A ggplot object of the plot.
