# Run limma analysis per contrast on filtered proteins

This function subsets the DAList for each contrast according to
`filtered_proteins_per_contrast` and fits the limma model separately. It
is useful when proteins are filtered differently per contrast. After
model fitting, the function computes rolling standard deviations (moving
SDs) of log fold changes based on intensity-sorted proteins, and
computes z-scores.

## Usage

``` r
run_filtered_limma_analysis(
  DAList,
  design_formula = ~0 + group,
  pval_thresh = 0.05,
  lfc_thresh = 1,
  adj_method = "BH",
  contrasts_file = NULL,
  binsize = "auto",
  plot_movingSD = FALSE
)
```

## Arguments

- DAList:

  A DAList object containing full data and per-contrast filtered
  proteins.

- design_formula:

  A formula for the model design (e.g., `~0 + group`, or
  `~0 + cell:treatment`).

- pval_thresh:

  P-value threshold used for defining significance. Default = 0.05.

- lfc_thresh:

  Log2 fold change threshold for significance. Default = 1.

- adj_method:

  Adjustment method for multiple testing ("BH", "BY", "none", etc).
  Default = "BH".

- contrasts_file:

  Optional CSV file with contrast definitions if not present in the
  DAList. Each row should contain a full contrast statement of the form
  `Label = expression`. The label (left-hand side) may be any string; it
  does not need to encode group names.

- binsize:

  Either an integer for the moving window size, or `"auto"` (default) to
  select an appropriate bin size automatically.

- plot_movingSD:

  Logical. If TRUE (default), plot moving SD curves for each contrast.

## Value

The input DAList, updated with per-contrast model fits and results,
moving SDs, and logFC z-scores.

## Contrast metadata (for downstream writers)

For each contrast, a sidecar list is saved at
`DAList$tags$per_contrast[[label]]$contrast_info` so downstream
functions (e.g.,
[`write_limma_tables()`](https://byrumlab.github.io/proteoDA/reference/write_limma_tables.md))
do not need to parse labels. `contrast_info` contains:

- `label`: contrast label (left-hand side)

- `contrast_expression_raw`: the RHS as written in the file

- `contrast_expression`: the RHS actually passed to limma (after
  translation)

- `design_formula`: the design formula as a character string

- `group_col`: grouping column if applicable (or `NULL` for pure
  interaction designs)

- `factors`: factor names participating in a simple two-factor
  interaction (e.g., `c("cell","treatment")`)

- `design_columns_involved`: model-matrix columns referenced by the
  translated expression

- `involved_levels`: levels from `group_col` appearing in the expression
  (if applicable)

## Examples

``` r
if (FALSE) { # \dontrun{
# Using default binsize auto-selection
filtered_DAList <- run_filtered_limma_analysis(filtered_DAList)

# With specified bin size and no plots
filtered_DAList <- run_filtered_limma_analysis(
    filtered_DAList, 
    binsize = 200, 
    plot_movingSD = FALSE)
} # }
```
