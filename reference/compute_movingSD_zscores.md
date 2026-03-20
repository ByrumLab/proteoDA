# Compute Rolling Standard Deviations and Z-scores for logFC Comparisons

This function computes a rolling (moving) standard deviation of log fold
change (logFC) values for each comparison in the `DAList$results` list.
For each contrast, proteins are ordered by mean log2 intensity (without
modifying the underlying data matrix). A moving standard deviation of
logFC is computed along this ordering and used to derive logFC Z-scores.
Results are stored both in `DAList$results[[contrast]]` and
`DAList$tags`.

## Usage

``` r
compute_movingSD_zscores(
  DAList,
  binsize = "auto",
  plot = FALSE,
  contrasts_file = NULL,
  qc_dir = NULL
)
```

## Arguments

- DAList:

  A DAList-like object containing at least:

  - `data`: A matrix/data frame with numeric expression values (log2
    intensities; proteins in rows, samples in columns).

  - `results`: A named list of per-contrast data frames, each containing
    a `logFC` column and rownames as protein IDs.

  - Optionally, `data_per_contrast`: a named list of per-contrast
    matrices/data frames. When present, `data_per_contrast[[contrast]]`
    is used for that contrast instead of the global `data`.

  - Optionally, `QC_dir`: a string path for saving QC plots.

  - Optionally, `design$contrast_vector` or
    `filtered_proteins_per_contrast` to infer contrast names.

- binsize:

  Integer or `"auto"`. The number of consecutive proteins (rows) to use
  in each rolling window, or `"auto"` to determine the best size
  automatically from the actual protein count.

- plot:

  Logical. If `TRUE`, plots the final moving SDs for each contrast and
  (if `DAList$QC_dir` is available) saves them to disk.

- contrasts_file:

  Optional CSV file with contrasts (one per row, `"name = formula"`
  syntax) used only to determine contrast names if they are not already
  available in the `DAList`.

- qc_dir:

  Optional directory path where PNG QC plots will be saved (CV vs
  binsize + per-contrast movingSD plots) when `plot = TRUE`.

## Value

The input `DAList` object with updated `tags`:

- `DAList$tags$movingSDs[[contrast]]`: rolling SD vector for each
  contrast.

- `DAList$tags$logFC_z_scores[[contrast]]`: Z-score vector for each
  contrast.

- `DAList$results[[contrast]]$movingSD`: moving SD per protein.

- `DAList$results[[contrast]]$logFC_z_scores`: logFC Z-score per
  protein.

## Details

If `binsize = "auto"`, the function derives candidate bin sizes from the
protein count, evaluates their stability using the coefficient of
variation (CV) of the moving SD curves across contrasts, and chooses a
binsize as the smallest value whose CV is within 5\\ of proteins (40)
and a maximum of 20\\

If `plot = TRUE`, two plots are generated per contrast:

1.  Moving SD vs rank index (proteins sorted by intensity).

2.  Moving SD vs mean log2 intensity (MD-style).

Plots are printed to the current device and, if `DAList$QC_dir` is set,
saved as PNG files there. When `binsize = "auto"`, a global
`movingSD_binsize_CV.png` (CV vs binsize) is also saved in `QC_dir`.
