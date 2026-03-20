# Generic plotting functions for normalization report

An internal set of functions that take in processed data ready for
plotting and create generic plots (e.g., means, violin plots, line
plots). These generic plots are then passed to metric-specific plotting
functions (see
[`pn_plots`](https://byrumlab.github.io/proteoDA/reference/pn_plots.md))
to create the final plot objects.

## Usage

``` r
pn_mean_plot(plotData)

pn_violin_plot(plotData)

pn_density_plot(plotData)
```

## Arguments

- plotData:

  A data frame of data to be plotted, created with the
  [`eval_pn_metric_for_plot`](https://byrumlab.github.io/proteoDA/reference/eval_pn_metric_for_plot.md)
  function.

## Value

A ggplot object of the plot.
