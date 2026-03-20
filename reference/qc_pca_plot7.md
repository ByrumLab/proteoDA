# PCA Plot

Generates PCA plots from normalized intensity data.

## Usage

``` r
qc_pca_plot7(
  data,
  groups = NULL,
  sample_labels = NULL,
  text.sizes = c(14, 14, 12, 10),
  point.size = 3,
  label.size = 4,
  max.labels = 100,
  top = 500,
  standardize = TRUE,
  pca_axes = c(1, 2),
  colors = NULL,
  title = NULL,
  max_pc = NULL,
  legend.position = "right",
  show.plot = TRUE
)
```

## Arguments

- data:

  A data frame of intensity data, likely normalized. Rows should be
  proteins and columns should be samples.

- groups:

  A character or factor vector, listing the group(s) the samples belong
  to. Default: NULL

- sample_labels:

  Optional, a vector of sample labels to use. If not supplied, defaults
  to using the column names in the data. Default: NULL

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

- top:

  The number of most variable proteins to use for the PCA analysis.
  Default: 500

- standardize:

  Logical indicating whether to standardize the input data before PCA.
  Default: TRUE

- pca_axes:

  A numeric vector of length 2 specifying which PCA axes to plot.
  Default: c(1, 2)

- colors:

  A vector of colors corresponding to groups. Default: NULL

- title:

  Title for the plot. Default: NULL

- max_pc:

  Optional integer giving the maximum principal component index to
  include in the scree plot produced by
  [`qc_pca_scree_plot7()`](https://byrumlab.github.io/proteoDA/reference/qc_pca_scree_plot7.md).
  Default: NULL (shows all PCs).

- legend.position:

  Position of the legend in the plot. Default: 'right'

- show.plot:

  Logical indicating whether to display the plot. Default: TRUE

## Value

A list containing the PCA plot, the PCA analysis results, and additional
information.

## Details

This function performs PCA on the provided data and generates a
corresponding plot.

## See also

[`arg_match`](https://rlang.r-lib.org/reference/arg_match.html),
[`cli_abort`](https://cli.r-lib.org/reference/cli_abort.html),
[`prcomp`](https://rdrr.io/r/stats/prcomp.html),
[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html),
[`aes`](https://ggplot2.tidyverse.org/reference/aes.html),
[`geom_label_repel`](https://ggrepel.slowkow.com/reference/geom_text_repel.html),
[`as.ggplot`](https://rdrr.io/pkg/ggplotify/man/as.ggplot.html)

## Examples

``` r
if (FALSE) { # \dontrun{
if(interactive()){
 # Example usage
 qc_pca_plot7(data, groups)
 }
} # }
```
