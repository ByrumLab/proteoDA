# Functions for plotting normalization metrics

A set of functions that take in normalized sample data and a grouping
factor and calculate some metric of of variability, error, etc., that we
can use to evaluate normalization methods. See
[`norm_metrics`](https://byrumlab.github.io/proteoDA/reference/norm_metrics.md)
for explanations of each metric.

## Usage

``` r
pn_plot_PCV(normList, grouping)

pn_plot_PMAD(normList, grouping)

pn_plot_PEV(normList, grouping)

pn_plot_COR(normList, grouping)

pn_plot_log2ratio(normList, grouping, zoom = F, legend = T)

pn_plot_MD(normList, grouping, use_ggrastr = F)
```

## Arguments

- normList:

  A named list of normalized data matrices.

- grouping:

  A character or factor vector, listing the group(s) the samples belong
  to.

- zoom:

  Should the plot cover the full range of log2ratios, or zoom in around
  0? Default is FALSE.

- legend:

  Should the plot include the legend? Default is TRUE.

- use_ggrastr:

  Use ggrastr in MD plots to reduce size (if installed)? Default
  `FALSE`.

## Value

A ggplot object of the plot.
