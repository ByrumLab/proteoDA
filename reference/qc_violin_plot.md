# Violin plot for QC report

Makes a violin plot of per-sample intensities. Samples are grouped by
the "groups" argument on the x-axis.

## Usage

``` r
qc_violin_plot(data, groups = NULL, sample_labels = colnames(data))
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

A ggplot object of the plot.
