# PCA Plot for Subgroups

Generates PCA plots for different subgroups based on provided metadata.

## Usage

``` r
qc_pca_plot7_subgroups(
  DAList,
  grouping_column,
  label_column = NULL,
  group_color = NULL,
  top = 500,
  standardize = TRUE,
  pca_axes = c(1, 2),
  text.sizes = c(14, 14, 12, 10),
  point.size = 3,
  label.size = 4,
  max.labels = 100,
  legend.position = "right",
  show.plot = TRUE,
  max_pc = NULL
)
```

## Arguments

- DAList:

  A list containing data and metadata for differential analysis. Should
  include a 'data' matrix and 'metadata' data frame.

- grouping_column:

  The name of the column in the metadata that defines the grouping for
  the samples.

- label_column:

  The name of the column in the metadata for sample labels. Default:
  NULL

- group_color:

  A color to represent the groups in the plot. Default: NULL (will use
  default color).

- top:

  The number of the most variable proteins to include in the PCA
  analysis. Default: 500

- standardize:

  Logical indicating whether to standardize the input data to have mean
  0 and standard deviation 1. Default: TRUE

- pca_axes:

  A numeric vector of length 2 specifying which PCA axes to plot.
  Default: c(1, 2)

- text.sizes:

  A numeric vector specifying text sizes for different plot elements:
  title, axis titles, axis labels, and legend. Default: c(14, 14, 12,
  10)

- point.size:

  Size of points in the PCA plot. Default: 3

- label.size:

  Size of sample labels in the plot. Default: 4

- max.labels:

  Maximum number of labels to show in the plot. Default: 100

- legend.position:

  Position of the legend in the plot. Default: 'right'

- show.plot:

  Logical indicating whether to display the plot. Default: TRUE

- max_pc:

  Optional integer giving the maximum principal component index to
  display in accompanying scree plots. Passed to
  [`qc_pca_plot7()`](https://byrumlab.github.io/proteoDA/reference/qc_pca_plot7.md)
  and
  [`qc_pca_scree_plot7()`](https://byrumlab.github.io/proteoDA/reference/qc_pca_scree_plot7.md).
  Default: NULL (shows all PCs).

## Value

A list containing PCA plots for each subgroup and the PCA analysis
results.
