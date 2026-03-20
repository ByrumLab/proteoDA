# Write tables of limma results

Generates CSV and Excel summary tables from a limma differential
abundance analysis stored in a `DAList`.

## Usage

``` r
write_limma_tables(
  DAList,
  output_dir = NULL,
  overwrite = FALSE,
  contrasts_subdir = NULL,
  summary_csv = NULL,
  combined_file_csv = NULL,
  spreadsheet_xlsx = NULL,
  add_filter = TRUE,
  color_palette = NULL,
  add_contrast_sheets = TRUE,
  statlist = NULL
)
```

## Arguments

- DAList:

  A DAList object with statistical results in `$results`, annotation in
  `$annotation`, data in `$data`, metadata in `$metadata`, and optional
  tags in `$tags`.

- output_dir:

  Directory to output tables. Defaults to the current working directory
  if `NULL`.

- overwrite:

  Logical; overwrite existing files if they are found in `output_dir`?
  Default is `FALSE`.

- contrasts_subdir:

  Subdirectory (within `output_dir`) for per-contrast CSV files.

- summary_csv:

  Filename for the summary CSV (per-contrast counts).

- combined_file_csv:

  Filename for the combined results CSV.

- spreadsheet_xlsx:

  Filename for the Excel spreadsheet.

- add_filter:

  Logical; add Excel autofilters to columns in the main worksheet
  produced by
  [`write_limma_excel()`](https://byrumlab.github.io/proteoDA/reference/write_limma_excel.md).

- color_palette:

  Optional vector of colors used for conditional formatting in the Excel
  output.

- add_contrast_sheets:

  Logical; if `TRUE`, each per-contrast CSV is added as a separate
  worksheet to the Excel workbook.

- statlist:

  Optional list of per-contrast result tables to use instead of
  `DAList$results`. If `NULL`, `DAList$results` is used.

## Value

Invisibly returns the (validated) input `DAList`.

## Details

For each contrast, per-contrast CSV files are written using
[`write_per_contrast_csvs()`](https://byrumlab.github.io/proteoDA/reference/write_per_contrast_csvs.md),
combining annotation, data, and statistics. If
`DAList$tags$per_contrast[[label]]$contrast_info` contains a usable
`group_col` and non-empty `involved_levels`, those samples are selected
for the per-contrast CSV. Otherwise the older `"_vs_"` label parser is
used. If that also fails, all samples are included with a warning.

The function then:

- writes a summary CSV of up/down/non-significant counts per contrast;

- writes per-contrast CSV files into `contrasts_subdir`;

- writes a combined results CSV with all contrasts side by side;

- writes an Excel workbook using
  [`write_limma_excel()`](https://byrumlab.github.io/proteoDA/reference/write_limma_excel.md),
  and optionally adds each per-contrast CSV as a worksheet.

The columns included in the per-contrast and Excel tables are determined
by internal settings such as `DA_table_cols` and `stat_cols`, which are
defined in the package configuration rather than as function arguments.
