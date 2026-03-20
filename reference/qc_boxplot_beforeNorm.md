# QC boxplot for data before normalization

Generate a QC boxplot from raw intensity data prior to normalization.
The input may be a numeric matrix or data.frame, or a DAList-like object
containing a top-level element named 'data'. If the data do not appear
to be on the log2 scale, a log2 transform is applied before plotting.

## Usage

``` r
qc_boxplot_beforeNorm(
  data,
  groups = NULL,
  sample_labels = NULL,
  pseudo = NULL,
  ...
)
```

## Arguments

- data:

  A numeric matrix or data.frame of intensities, or a DAList-like object
  (list) containing a 'data' matrix.

- groups:

  Optional vector indicating class membership. Length must match the
  number of samples (columns).

- sample_labels:

  Optional character vector of sample labels. Length must match the
  number of samples (columns).

- pseudo:

  Optional positive numeric value added before log2 transform to avoid
  log2(0). If NULL, a small data-driven constant is chosen.

- ...:

  Additional arguments passed to qc_boxplot.

## Value

A ggplot2 object displaying the QC boxplot of log2-scaled raw data.

## Details

Log2 scale detection uses a conservative heuristic based on data range
and value distribution. When the scale is ambiguous, the function will
apply a log2 transform to ensure comparable distributions.
