# Generate a QC Boxplot

This function draws a boxplot of sample intensities for quality control.
The `data` argument may be a numeric matrix or data.frame of
intensities, or a DAList object produced by the proteoDA pipeline. If a
DAList is supplied, the function will prefer the normalized data stored
in `DAList$data_per_contrast[[contrast]]` when that element exists. If
'contrast' is NULL and data_per_contrast has exactly one element, that
element will be used. Otherwise DAList\$data is used.

## Usage

``` r
qc_boxplot(
  data,
  contrast = NULL,
  groups = NULL,
  sample_labels = NULL,
  title = NULL,
  text.sizes = c(12, 10, 10, 10),
  legend.position = "right",
  boxplot_width = NULL,
  boxplot_alpha = NULL,
  plot_margin = NULL,
  colorblind_palette = NULL
)
```

## Arguments

- data:

  A matrix or data.frame of intensities, or a DAList object.

- contrast:

  Optional character string naming the contrast to use when a DAList is
  provided. If NULL and data_per_contrast has exactly one element, that
  element is used.

- groups:

  A vector indicating class membership (numeric, integer, character,
  factor, or NULL). Length must match number of samples.

- sample_labels:

  A character vector of names to display on the plot, or NULL to use
  column names from the data.

- title:

  A character string for the plot title.

- text.sizes:

  Numeric vector of length 4 specifying text sizes for title, y-axis,
  x-axis labels, and legend text, in that order.

- legend.position:

  Character string specifying legend position (default: \\right\\).

- boxplot_width:

  Numeric width for the boxplots. If NULL a dynamic value based on
  number of groups is used.

- boxplot_alpha:

  Numeric transparency for the boxplots. If NULL a dynamic value based
  on number of groups is used.

- plot_margin:

  A grid::unit object specifying the plot margins. If NULL a dynamic
  margin is used.

- colorblind_palette:

  Character vector of color values to use for group fills. If NULL a
  default colorblind friendly palette is used.

## Value

A ggplot2 object displaying the QC boxplot.
