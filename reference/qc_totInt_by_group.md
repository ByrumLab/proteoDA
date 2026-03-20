# Generate Total Intensity Barplots by Group

This function creates barplots showing total intensity values for
samples, grouped by a specified column. It includes an optional
percentile threshold to flag potential outliers and uses a color-blind
friendly palette.

## Usage

``` r
qc_totInt_by_group(
  DAList,
  label_column,
  grouping_column,
  percentile,
  colors = NULL,
  nrow = NULL,
  ncol = NULL,
  legend.position = "none",
  plot_size = NULL,
  save_path = "Intensity_barplot.png"
)
```

## Arguments

- DAList:

  A list containing `data` (matrix) and `metadata` (data frame).

- label_column:

  A character string specifying the column name in metadata for sample
  labels.

- grouping_column:

  A character string specifying the column name in metadata for
  grouping.

- percentile:

  A numeric value (0-1) for the percentile threshold to flag
  low-intensity samples.

- colors:

  A vector of colors for groups. Defaults to a color-blind friendly
  palette (max 12, \>12 uses Polychrome).

- nrow:

  Number of rows for facet wrapping.

- ncol:

  Number of columns for facet wrapping.

- legend.position:

  Position of the legend ('none', 'right', 'top', 'left', 'bottom').

- plot_size:

  A numeric value of length 2 indicating the width and height of the
  plot (in inches). Optional.

- save_path:

  A character string specifying the file path to save the plot. If NULL,
  the plot is not saved.

## Value

A list containing ggplot objects for total intensity, total number, and
total missing values barplots.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Example of barplot colored by group
  barplot1 <- qc_totInt_by_group(
    DAList          = results,
    label_column    = "sample",
    grouping_column = "group",
    percentile      = 0,
    colors          = NULL,  # or c("#E69F00", "#56B4E9", "#009E73")
    nrow            = NULL,
    ncol            = NULL,
    legend.position = "right",
    plot_size       = c(12, 6),
    save_path       = "Intensity_barplot.png"
  )

  # Save plot (example; adjust to match returned object structure)
  # ggplot2::ggsave(
  #   "total_intensity_plot.png",
  #   barplot1$plot,
  #   width  = barplot1$plot_size[1],
  #   height = barplot1$plot_size[2]
  # )
} # }
```
