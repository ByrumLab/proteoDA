# Generate a Barplot of Total Intensities by Sample Group

This function creates a barplot where samples are organized and colored
by group. It also identifies potential outliers based on a specified
percentile of total intensity values.

## Usage

``` r
qc_totInt(
  DAList,
  label_column,
  grouping_column,
  percentile,
  colors = NULL,
  legend.position = "right"
)
```

## Arguments

- DAList:

  A list containing `metadata` and `data`. The `metadata` should include
  sample labels and grouping information.

- label_column:

  A string specifying the column name in `metadata` containing sample
  labels.

- grouping_column:

  A string specifying the column name in `metadata` containing grouping
  information.

- percentile:

  A numeric value (0-1) specifying the percentile threshold for
  identifying low-intensity samples.

- colors:

  A vector of colors (one per group). If NULL, a predefined color
  palette is used.

- legend.position:

  A string specifying legend position. Options: "none", "right", "top",
  "left", "bottom".

## Value

A list containing:

- `p`: A ggplot object of the barplot.

- `bar_data`: A data frame with total intensity values and grouping
  information.

- `colors`: The color palette used.

- `percentile`: The input percentile value.

- `perc_int_thresh`: The computed intensity threshold.
