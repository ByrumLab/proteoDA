# Make an HTMLwidget for the interactive report

This internal function is called within the .Rmd report template as part
of the user-facing
[`write_limma_plots`](https://byrumlab.github.io/proteoDA/reference/write_limma_plots.md)
function. It takes in data on statistical results, annotation, and
sample intensities (partially assembled in
[`write_limma_plots`](https://byrumlab.github.io/proteoDA/reference/write_limma_plots.md)
and then passed on the to .RMD environment), does some further
processing, packages it into the form needed for our interactive report,
and then outputs an HTMLwidget for use in the HTML file generated from
the .Rmd.

## Usage

``` r
uams_glimmaXY(
  model_data,
  counts,
  groups,
  anno,
  display.columns,
  status.cols,
  sample.cols,
  width,
  height
)
```

## Arguments

- model_data:

  The output from the
  [`prep_plot_model_data`](https://byrumlab.github.io/proteoDA/reference/prep_plot_model_data.md)
  function, containing statistical results for a single contrast.

- counts:

  A matrix or dataframe of the raw or normalized data from the data slot
  of a DAList.

- groups:

  A vector of group identities for each sample. Should have same length
  as the number of cols in counts.

- anno:

  annotation data, with added p-value columns.

- display.columns:

  A vector of columns to display in the output table.

- status.cols:

  A vector of colors to use for down regulated, nonDE, and upregulated
  proteins. Must be of length 3.

- sample.cols:

  A vector of colors for each sample. Should have the same length as
  groups.

- width:

  The width of the interactive report objects, in pixels..

- height:

  The height of the interactive report objects, in pixels.

## Value

An HTMLwidget containing our interactive plots and tables
