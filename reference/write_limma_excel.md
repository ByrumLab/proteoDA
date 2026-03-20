# Write an .xlsx file of limma results

Creates a nicely formatted Excel spreadsheet of differential expression
results.

## Usage

``` r
write_limma_excel(
  filename,
  statlist,
  annotation,
  data,
  norm.method,
  pval_thresh,
  lfc_thresh,
  add_filter,
  color_palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
    "#CC79A7", "#999999"),
  annot_cols = NULL,
  filt_min_reps = NULL,
  filt_min_groups = NULL
)
```

## Arguments

- filename:

  The filename of the Excel spreadsheet to be saved.

- statlist:

  A list of the per-contrast statistical results.

- annotation:

  A data frame of annotation data for each protein.

- data:

  A data frame containing the average expression data for each sample.

- norm.method:

  The normalization method used in the model.

- pval_thresh:

  The p-value threshold for significance.

- lfc_thresh:

  The log-fold change threshold for significance.

- add_filter:

  Logical. Whether to add Excel filters to columns.

- color_palette:

  A vector of colors for contrast sections (optional).

- annot_cols:

  Character vector of annotation column names to include. If `NULL`, all
  columns are used.

## Value

A list with filename and workbook object (invisibly).
