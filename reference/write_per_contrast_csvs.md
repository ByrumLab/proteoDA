# Write per-contrast results csvs

Internal function used to write per-contrast statistical results to .csv
files.

## Usage

``` r
write_per_contrast_csvs(
  annotation_df,
  data,
  results_statlist,
  output_dir,
  annotation_cols = NULL,
  metadata = NULL,
  group_col = "group",
  stat_cols = c("logFC", "P.Value", "adj.P.Val", "movingSD", "logFC_z_scores"),
  per_contrast_tags = NULL
)
```

## Arguments

- annotation_df:

  A data frame of annotation data for each gene/protein.

- data:

  A data frame of average expression data for each sample.

- results_statlist:

  A list of per-contrast DE results.

- output_dir:

  The directory in which to save the per-contrast csv files.

- annotation_cols:

  Optional character vector of annotation column names to include.

- metadata:

  Optional data.frame of sample metadata (e.g. DAList\$metadata).

- group_col:

  Character. Name of the grouping column within metadata to match
  contrast samples.

- stat_cols:

  Optional character vector of statistical columns to include from the
  results_statlist.

- per_contrast_tags:

  Optional list like DAList\$tags\$per_contrast; if present, will be
  used to determine which sample groups/levels are involved in each
  contrast.

## Value

A logical vector indicating whether each contrast file was successfully
written.
