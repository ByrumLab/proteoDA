# Dendrogram for QC Report with Subgroups

Generates dendrograms for different subgroups within a dataset. Each
dendrogram highlights one subgroup while keeping others in grayscale.

## Usage

``` r
qc_dendrogram_subgroups(
  DAList,
  grouping_column,
  label_column = NULL,
  text.sizes = c(12, 4, 12),
  point.size = 3,
  group_color = NULL,
  top = 500,
  standardize = TRUE,
  dist_metric = "euclidean",
  clust_method = "complete",
  legend.position = "right",
  show.plot = TRUE
)
```

## Arguments

- DAList:

  A list containing data and metadata.

- grouping_column:

  A character string specifying the sample metadata column to use for
  subgrouping (e.g., "group","batch", etc)

- label_column:

  An optional character string specifying the sample metadata column for
  sample labels.

- text.sizes:

  A numeric vector specifying the text sizes for the title, leaf labels,
  and legend text.

- point.size:

  A numeric vector specifying the size of the node circles.

- group_color:

  A character string specifying the highlight color for the subgroup.

- top:

  The number of the most variable proteins to use in the analysis.

- standardize:

  A logical value indicating whether to standardize data (e.g., mean =
  0, SD= 1).

- dist_metric:

  A character string specifying the distance metric (e.g., "euclidean",
  "manhattan").

- clust_method:

  A character string specifying the clustering method (e.g, "complete",
  "ward.D2")

- legend.position:

  A character string specifying the legend position (e.g., "right",
  "left", "top", "bottom").

- show.plot:

  A logical value indicating whether to display the plot.

## Value

A list of dendrogram plots for each subgroup.

## Examples
