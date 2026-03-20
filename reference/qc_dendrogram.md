# Dendrogram for QC report

Performs and then plots a hierarchical clustering analysis of sample
intensities across samples. Sample labels are colored according to the
`groups` argument, with colors determined by the
[`colorGroup`](https://byrumlab.github.io/proteoDA/reference/colorGroup.md)
function. By default, uses only the 500 most variable proteins for the
analysis.

## Usage

``` r
qc_dendrogram(
  data,
  groups = NULL,
  sample_labels = NULL,
  text.sizes = c(14, 3, 12),
  point.size = 3,
  top = 500,
  standardize = TRUE,
  dist_metric = "euclidean",
  clust_method = "complete",
  colors = NULL,
  legend.position = "right",
  title = NULL,
  subtitle = NULL,
  show.plot = TRUE
)
```

## Arguments

- data:

  A numeric matrix or data frame with samples as columns and features
  (e.g. proteins) as rows.

- groups:

  A vector indicating sample groupings for color coding in the
  dendrogram.

- sample_labels:

  A vector of sample labels. If NULL, column names of `data` are used.

- text.sizes:

  A numeric vector specifying text sizes for title, leaf labels, and
  legend text.

- point.size:

  A numeric value specifying the size of the node circles.

- top:

  The number of the most variable proteins to use in the analysis.

- standardize:

  A logical value indicating whether to standardize data (mean = 0,
  SD=1).

- dist_metric:

  A character string specifying the distance metric (e.g., "euclidean",
  "manhattan").

- clust_method:

  A character string specifying the clustering method (e.g., "complete",
  "ward.D2").

- colors:

  A vector of colors corresponding to group levels.

- legend.position:

  A character string specifying legend position (e.g., "right",
  "bottom").

- title:

  A character string specifying the plot title.

- subtitle:

  A character string specifying the plot subtitle. It can be NULL.

- show.plot:

  A logical value indicating whether to display the plot.

## Value

A ggplot list containing

- `p` - a ggplot dendrogram object.

- `hc` - The hierarchical clustering object

- `data_na` - the processed data used in the clustering

- `sample_group_info` - A data frame containing sample labels and groups

- `param` - A list of parameters used in the function

## Examples

``` r
if (FALSE) { # \dontrun{
den <- qc_dendrogram(data = results$data,
                 groups        = results$metadata$group,
                 sample_labels = results$metadata$sample,
                 top           = nrow(results$data),
                 standardize   = TRUE,
                 dist_metric   = "euclidean",
                 clust_method  = "complete",
                 colors        = all_colors2$group,
                 point.size    = 3,
                 text.sizes    = c(14,3,9),
                 legend.position = "right",
                 title         = "",
                 show.plot     = FALSE)

} # }
```
