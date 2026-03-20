# Generate a Violin Plot for Intensity Data

Generate a Violin Plot for Intensity Data

## Usage

``` r
qc_violin(
  data,
  groups = NULL,
  sample_labels = NULL,
  title = "",
  text.sizes = c(12, 10, 10, 10),
  legend.position = "right",
  color_palette = NULL
)
```

## Arguments

- data:

  A matrix or data frame containing intensity values.

- groups:

  A vector (numeric, integer, character, or factor) indicating class
  membership. Defaults to NULL.

- sample_labels:

  A character vector specifying sample names to display on the plot.
  Defaults to column names of `data`.

- title:

  A character string specifying the plot title. Defaults to an empty
  string.

- text.sizes:

  A numeric vector of length 4 specifying font sizes for title, x-axis,
  y-axis, and legend text, respectively.

- legend.position:

  A character string specifying legend position (e.g., "right", "left",
  "none"). Defaults to "right".

- color_palette:

  A character vector specifying a colorblind-friendly palette. If NULL,
  a default base palette is used based on RColorBrewer.

## Value

A ggplot object representing the violin plot.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Example of violin plot colored by groups 
  violin1 <- qc_violin(
    data            = results$data,
    groups          = results$metadata$group,
    sample_labels   = results$metadata$sample,
    title           = "",
    text.sizes      = c(12, 10, 10, 10),
    legend.position = "right",
    color_palette   = c("#66c2A5", "#8DA0CB", "#FC8D62")
  )
} # }
```
