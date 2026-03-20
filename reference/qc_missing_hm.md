# Create a missing values heatmap

Makes a ComplexHeatmap object showing a heatmap of missing values in the
input data.

## Usage

``` r
qc_missing_hm(
  data,
  groups,
  sample_labels = NULL,
  column_sort = c("cluster", "group"),
  group_var_name = "",
  show_all_proteins = F
)
```

## Arguments

- data:

  A data frame of intensity data, likely normalized. Rows should be
  proteins and columns should be samples.

- groups:

  A character or factor vector, listing the group(s) the samples belong
  to.

- sample_labels:

  Optional, a vector of sample labels to use. If not supplied, defaults
  to using the column names in the data.

- column_sort:

  How should the columns of the heatmap be sorted? Options are:
  "cluster"- sort by similarity in missing values, "group"- sort samples
  by the grouping variable.

- group_var_name:

  The name of the variable being used for group sorting.

- show_all_proteins:

  Should all proteins be shown in missing value heatmap, of only those
  with missing data? Default is F (only those with missing data).

## Value

A [`grid::gTree`](https://rdrr.io/r/grid/grid.grob.html) object of the
ComplexHeatmap, which can be plotted with
[`grid::grid.draw`](https://rdrr.io/r/grid/grid.draw.html).
