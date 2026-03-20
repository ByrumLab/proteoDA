# Create a normalization report (optionally contrast-aware)

Saves a PDF report containing plots/metrics to compare normalization
methods.

## Usage

``` r
write_norm_report(
  DAList,
  grouping_column = NULL,
  output_dir = NULL,
  filename = NULL,
  overwrite = FALSE,
  suppress_zoom_legend = FALSE,
  use_ggrastr = FALSE,
  input_is_log2 = FALSE,
  contrasts = NULL,
  sample_id_col = NULL,
  groups_override = NULL,
  metrics_csv = NULL,
  per_contrast = FALSE,
  include_MD_plots = TRUE
)
```

## Arguments

- DAList:

  A DAList.

- grouping_column:

  Name of the metadata column that gives sample groups (used by
  within-group cyclic loess). Must contain at least two groups.

- output_dir:

  Directory to save the report (created if missing). Defaults to the
  current working directory.

- filename:

  File name of the report (must end with `.pdf`). Default
  `"normalization_report.pdf"`.

- overwrite:

  Overwrite if the report already exists? Default `FALSE`.

- suppress_zoom_legend:

  Remove the legend from the zoomed log2-ratio plot? Default `FALSE`.

- use_ggrastr:

  Use ggrastr in MD plots to reduce size (if installed)? Default
  `FALSE`.

- input_is_log2:

  Logical. If `TRUE`, indicates per-contrast (or global) input matrices
  are already on the log2 scale for methods that expect log2 input.
  Default `FALSE`.

- contrasts:

  Optional character vector: when `per_contrast = TRUE` and
  `data_per_contrast` is present, restrict the report to these
  contrasts. Default = all contrasts in `data_per_contrast`.

- sample_id_col:

  Optional name of the metadata column whose values match the matrix
  column names. If `NULL`, the function attempts to auto-detect a
  suitable column (e.g., `"sample"`).

- groups_override:

  Optional named vector of group labels with names equal to **matrix
  column names**. If supplied, this overrides the group labels derived
  from `grouping_column` and `sample_id_col`.

- metrics_csv:

  Optional path to a CSV file where per-sample normalization metrics
  (PCV, PMAD, PEV, COR, mean intensity, etc.) will be written. If
  `NULL`, metrics are not exported to CSV.

- per_contrast:

  Logical. If `TRUE` and `DAList$data_per_contrast` is present, produce
  a separate set of normalization plots per contrast. If `FALSE`
  (default), ignore `data_per_contrast` and evaluate normalization
  globally using `DAList$data`.

- include_MD_plots:

  Logical. If `TRUE` (default), include the extended MD/MA diagnostic
  pages for each contrast (or global). Set to `FALSE` to omit these
  plots, which can substantially reduce the PDF file size.

## Value

Invisibly returns the input DAList.

## Details

By default, normalization is evaluated **globally** on `DAList$data`,
which is typically the raw (unnormalized) expression matrix.

If `per_contrast = TRUE` and `DAList$data_per_contrast` exists, the
function instead evaluates normalization **separately for each
contrast**, using the matrices in that slot (leaving `DAList$data`
untouched).

The report always includes per-sample normalization metrics (PCV, PMAD,
PEV, COR, mean intensity, etc.). When `include_MD_plots = TRUE`, it
additionally shows a *set* of MD/MA-style diagnostics for each
normalization method:

- classic MD plots (points + loess trend)

- trend-only MD plots (no points, to compare bias across methods)

- 2D-binned (heatmap-style) MD plots

- Delta-trend curves: difference between each method's loess fit and the
  raw log2 trend, highlighting over-/under-correction relative to the
  input

- residual plots: distributions of loess residuals per method

- log2FC distribution plots (per method), to visualize how normalization
  compresses or preserves fold-change variation
