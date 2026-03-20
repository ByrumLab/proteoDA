# Make a p-value histogram plot

Internal function for plotting p-value histograms of raw- and adjusted
p-values.

## Usage

``` r
static_pval_histogram(data, contrast)
```

## Arguments

- data:

  Per-contrast DE results to be plotted, as prepared by
  [`prep_plot_model_data`](https://byrumlab.github.io/proteoDA/reference/prep_plot_model_data.md).

- contrast:

  The contrast being plotted. Used for generating the plot title.

## Value

A ggplot object.
