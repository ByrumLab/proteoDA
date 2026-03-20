# Write a post-normalization evaluation report for a single method

Creates a PDF report with the same panels used in
[`write_norm_report()`](https://byrumlab.github.io/proteoDA/reference/write_norm_report.md),
but **using only the already-normalized values** from your `DAList`.
This is ideal for evaluating what a chosen normalization (e.g., cyclic
loess) actually did to your data.

## Usage

``` r
write_norm_eval_report(
  DAList,
  norm_label = "cycloess",
  grouping_column,
  output_dir = NULL,
  filename = NULL,
  overwrite = FALSE,
  suppress_zoom_legend = FALSE,
  use_ggrastr = FALSE,
  input_is_log2 = TRUE,
  contrasts = NULL,
  sample_id_col = NULL,
  groups_override = NULL,
  metrics_csv = NULL
)
```

## Arguments

- DAList:

  A DAList containing already-normalized matrices (per-contrast or
  global).

- norm_label:

  Character scalar used as the list name in plots (e.g., "cycloess").

- grouping_column:

  Name of the metadata column giving sample groups. Must contain \>= 2
  groups.

- output_dir:

  Directory to save the PDF (created if missing). Default: working
  directory.

- filename:

  File name for the PDF (must end with `.pdf`). Default:
  "norm_eval_report.pdf".

- overwrite:

  Overwrite existing file? Default FALSE.

- suppress_zoom_legend:

  Remove legend from the zoomed log2 ratio plot? Default FALSE.

- use_ggrastr:

  Use ggrastr in MD plots (if installed) to reduce file size? Default
  FALSE.

- input_is_log2:

  Logical, whether the matrices are already on the log2 scale for
  metrics that expect log2 (most are visualization/relative; this flag
  is unused here but kept for symmetry).

- contrasts:

  Optional character vector: restrict to these contrasts when
  `data_per_contrast` is present.

- sample_id_col:

  Optional name of the metadata column whose values match matrix
  colnames. If `NULL`, the function will try to auto-detect a suitable
  column.

- groups_override:

  Optional named vector of group labels with names equal to **matrix
  column names**. If supplied, overrides
  `grouping_column`/`sample_id_col` matching.

- metrics_csv:

  Optional path to a CSV file where per-sample normalization metrics
  (e.g., PCV, PMAD, PEV, COR, mean intensity) will be written. If
  `NULL`, metrics are not exported to CSV.

## Value

Invisibly returns `DAList`.

## Details

If `DAList$data_per_contrast` exists, a two-page block is produced **for
each contrast**: metrics grid (PCV, PMAD, PEV, COR, Log2 ratio and zoom)
and faceted MD plots. If no per-contrast data exist, a single two-page
report is produced from `DAList$data`.
