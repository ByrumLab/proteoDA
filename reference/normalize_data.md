# Normalize data in a DAList (contrast-aware)

Normalizes data using one of several methods. If
`DAList$data_per_contrast` exists and is non-empty, normalization is
applied **per-contrast** and the results are written back to
`DAList$data_per_contrast[[contrast]]`, leaving `DAList$data` unchanged.
Otherwise, it falls back to normalizing `DAList$data` (legacy global
behavior).

## Usage

``` r
normalize_data(
  DAList,
  norm_method = c("log2", "median", "mean", "vsn", "quantile", "cycloess", "rlr", "gi"),
  input_is_log2 = FALSE,
  contrasts = NULL,
  use_per_contrast_if_available = TRUE,
  ...
)
```

## Arguments

- DAList:

  A `DAList` object.

- norm_method:

  One of: "log2","median","mean","vsn","quantile","cycloess","rlr","gi".

- input_is_log2:

  Logical. If TRUE, indicates inputs are already on the log2 scale for
  methods that expect log2 input. Default FALSE.

- contrasts:

  Optional character vector to restrict which contrasts in
  `DAList$data_per_contrast` are normalized. Default = all available.

- use_per_contrast_if_available:

  Logical. Default TRUE. When TRUE and `DAList$data_per_contrast` is
  present/non-empty, normalize per-contrast and leave `DAList$data`
  unchanged.

- ...:

  Additional arguments forwarded to the underlying normalization
  function. For example, pass `groups=` to
  [`cycloessNorm()`](https://byrumlab.github.io/proteoDA/reference/norm_functions.md)
  to normalize within groups; if supplied, `groups` will be
  subset/reordered to the columns present in each per-contrast matrix.

## Value

The updated `DAList`. Tags are updated as follows:

- Per-contrast path: `DAList$tags$per_contrast_normalized = TRUE`,
  `DAList$tags$norm_method_per_contrast` (named list), and
  `DAList$normalization_per_contrast[[contrast]]$diagnostics`
  (before/after per-sample summaries) are recorded. `DAList$data` is
  unchanged.

- Global path: `DAList$tags$normalized = TRUE`,
  `DAList$tags$norm_method = norm_method`.

## Details

**Available methods**

- "log2" — Binary (base-2) log transformation of raw intensities.

- "median" — Divide each sample by its median; rescale by the mean of
  sample medians.

- "mean" — Divide each sample by its mean; rescale by the mean of sample
  means.

- "vsn" — Variance-stabilizing normalization via
  [`vsn::justvsn`](https://rdrr.io/pkg/vsn/man/justvsn.html) (expects
  raw scale).

- "quantile" — Quantile normalization via
  [`preprocessCore::normalize.quantiles`](https://rdrr.io/pkg/preprocessCore/man/normalize.quantiles.html).

- "cycloess" — Cyclic loess via
  [`limma::normalizeCyclicLoess`](https://rdrr.io/pkg/limma/man/normalizeCyclicLoess.html)
  (can be applied within groups).

- "rlr" — Robust linear regression (global) normalization.

- "gi" — Global intensity normalization (expects raw scale; returns
  log2).

**Log scale behavior** Methods "median", "mean", "quantile", "cycloess",
and "rlr" operate on log2 data. By default, this function will apply
`log2` internally before those methods. If your inputs are already
log2-transformed, set `input_is_log2 = TRUE` to avoid double logging.
