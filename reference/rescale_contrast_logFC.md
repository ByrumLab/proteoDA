# Rescale log-fold changes for a single contrast

Rescales the `logFC` (and optional `CI.L`/`CI.R`) for one contrast
inside a DAList-like object, leaving inferential statistics (`t`,
`P.Value`, `adj.P.Val`, `B`, etc.) unchanged. Optionally renames the
contrast and updates contrast metadata/tags.

## Usage

``` r
rescale_contrast_logFC(x, label, factor = 0.5, new_label = NULL)
```

## Arguments

- x:

  A DAList-like object containing `$results`, and optionally `$tags` and
  `$design$contrast_matrix`.

- label:

  Character. Name of the contrast to rescale (must exist in
  `x$results`).

- factor:

  Numeric scalar. Multiplicative factor applied to `logFC` (and
  `CI.L`/`CI.R` if present). Default `0.5`.

- new_label:

  Optional character. If supplied and different from `label`, the
  contrast is renamed consistently across `$results`, `$tags`, and
  `$design$contrast_matrix`.

## Value

The modified object `x`.

## Details

The function is robust to several result table shapes:

- a plain `data.frame` with a `logFC` column;

- a list with a `$table` `data.frame`;

- a list with a numeric `$logFC` vector;

- a plain numeric vector (assumed to be `logFC`).

When `new_label` is provided, the function also:

- renames the entry in `x$results`;

- moves any per-contrast tags from the old label to the new one;

- renames the corresponding column in `x$design$contrast_matrix`, if
  present.

A cumulative `rescale_factor` is tracked in
`x$tags$per_contrast[[label]]$rescale_factor`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Halve logFCs for a contrast and rename it
results <- rescale_contrast_logFC(results, label = "Avg_Treat_Bio_vs_DMSO",
                                  factor = 0.5, new_label = "Avg_Treat_x0.5")
} # }
```
