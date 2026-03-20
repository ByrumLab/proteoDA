# Write a movingSD QC report (PDF and/or PNGs)

Runs
[`compute_movingSD_zscores()`](https://byrumlab.github.io/proteoDA/reference/compute_movingSD_zscores.md)
with plotting enabled and saves either a multi-page PDF, per-contrast
PNGs, or both.

## Usage

``` r
write_movingSD_report(
  DAList,
  out_dir = "movingSD_QC",
  filename_base = "movingSD_report",
  binsize = "auto",
  contrasts_file = NULL,
  device = c("pdf", "png", "both"),
  width = 7,
  height = 5
)
```

## Arguments

- DAList:

  A DAList object with limma results.

- out_dir:

  Directory where plots will be stored. Created if needed.

- filename_base:

  Base name for the PDF report (default: "movingSD_report"). PNG files
  are named automatically as in
  [`compute_movingSD_zscores()`](https://byrumlab.github.io/proteoDA/reference/compute_movingSD_zscores.md).

- binsize:

  Integer or "auto". Passed to
  [`compute_movingSD_zscores()`](https://byrumlab.github.io/proteoDA/reference/compute_movingSD_zscores.md).

- contrasts_file:

  Optional contrasts file, passed through.

- device:

  One of "pdf", "png", or "both".

- width, height:

  PDF device size (in inches).

## Value

The updated DAList (invisibly).
