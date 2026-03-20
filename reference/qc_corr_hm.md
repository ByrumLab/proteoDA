# Correlation heatmap for QC report

Plots a ComplexHeatmap object showing the pairwise correlations of
intensities across samples. Samples are grouped by similarity on the x-
and y-axis, with labels colored by group and batch.

## Usage

``` r
qc_corr_hm(data, groups, sample_labels = NULL)
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

## Value

A [`grid::gTree`](https://rdrr.io/r/grid/grid.grob.html) object of the
ComplexHeatmap, which can be plotted with
[`grid::grid.draw`](https://rdrr.io/r/grid/grid.draw.html).
