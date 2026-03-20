# Build a cleaned statlist for export to Excel

Constructs a per-contrast statistics list from a `DAList` object with
full protein coverage. For each contrast, the function assembles a data
frame with one row per protein in `DAList$data` and selected columns
from the limma results plus optional moving standard deviations and
z-scores stored in `DAList$tags`.

## Usage

``` r
build_statlist(
  DAList,
  stat_cols = c("logFC", "P.Value", "adj.P.Val", "sig.PVal", "sig.FDR", "movingSDs",
    "logFC_z_scores")
)
```

## Arguments

- DAList:

  A DAList object containing:

  - `$data`: matrix/data frame of log-intensities (rows = proteins).

  - `$results`: named list of per-contrast limma result tables.

  - `$tags$movingSDs`: optional named list of moving SD vectors (one per
    contrast).

  - `$tags$logFC_z_scores`: optional named list of z-score vectors (one
    per contrast).

- stat_cols:

  Character vector of column names to include for each contrast.
  Defaults to:
  `c("logFC", "P.Value", "adj.P.Val", "sig.PVal", "sig.FDR", "movingSDs", "logFC_z_scores")`.

## Value

A named list of data frames, one per contrast, each with one row per
protein in `DAList$data` and columns given by `stat_cols`.

## Examples

``` r
if (FALSE) { # \dontrun{
statlist <- build_statlist(DAList)
names(statlist)
head(statlist[[1]])
} # }
```
