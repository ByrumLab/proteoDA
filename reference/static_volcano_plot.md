# Make a DE Volcano plot.

Internal function for plotting static versions of volcano plots.

## Usage

``` r
static_volcano_plot(
  data,
  lfc_thresh,
  pval_thresh,
  contrast,
  pval_type,
  control_proteins = NULL,
  anno = NULL,
  highlight_by = "uniprot_id"
)
```

## Arguments

- data:

  Per-contrast DE results to be plotted, as prepared by
  [`prep_plot_model_data`](https://byrumlab.github.io/proteoDA/reference/prep_plot_model_data.md).

- lfc_thresh:

  The logFC threshold used to determine significance (significant when
  \|logFC\| \> lfc.tresh). LogFC are base 2.

- pval_thresh:

  The p-value threshold used to determine significance (significant when
  p \< pval_thresh).

- contrast:

  The contrast being plotted. Used for generating the plot title.

- pval_type:

  The type of p-value to plot. Can be "raw" or "adjusted".

## Value

A ggplot object.
