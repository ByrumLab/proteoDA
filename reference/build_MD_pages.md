# Build extended MD/MA pages for normalization report

Internal helper: given a `normList` and `groups`, constructs a list of
ggplot objects with different MD-style diagnostics:

## Usage

``` r
build_MD_pages(normList, groups, use_ggrastr = FALSE, contrast_label = NULL)
```

## Arguments

- normList:

  List of normalized matrices from
  [`apply_all_normalizations_contrast()`](https://byrumlab.github.io/proteoDA/reference/apply_all_normalizations_contrast.md).

- groups:

  Group labels used to compute log2FC.

- use_ggrastr:

  Logical, whether to use ggrastr for point-heavy plots.

- contrast_label:

  Optional label used in plot titles.

## Value

A list of ggplot objects.

## Details

1.  classic MD (points + loess trend)

2.  trend-only MD (no points, SE band)

3.  2D-binned MD (heatmap-like density + trend)

4.  Delta-trend curves vs raw log2

5.  residual distributions (density of loess residuals per method)

6.  log2FC distribution per method (violin)
