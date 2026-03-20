# Make a DE MD plot.

Internal function for plotting static versions of MD plots.

## Usage

``` r
static_MD_plot(data, lfc_thresh, contrast, pval_type)
```

## Arguments

- data:

  Per-contrast DE results to be plotted, as prepared by
  [`prep_plot_model_data`](https://byrumlab.github.io/proteoDA/reference/prep_plot_model_data.md).

- lfc_thresh:

  The logFC threshold used to determine significance (significant when
  \|logFC\| \> lfc.tresh). LogFC are base 2.

- contrast:

  The contrast being plotted. Used for generating the plot title.

- pval_type:

  The type of p-value to plot. Can be "raw" or "adjusted".

## Value

A ggplot object.
