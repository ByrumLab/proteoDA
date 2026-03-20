# Analyze Within-Group Log2FC Standard Deviations

Calculates pairwise log2 fold change (log2FC) standard deviations within
each sample group, generates diagnostic density plots, flags groups with
outlier SDs, and exports a summary table.

## Usage

``` r
analyze_within_group_SDs(
  count_data,
  sample_metadata,
  group_column = "group",
  output_dir = "QC_dir",
  outlier_method = c("IQR", "z-score"),
  z_thresh = 3
)
```

## Arguments

- count_data:

  A data frame or matrix of expression data (rows = features, columns =
  samples).

- sample_metadata:

  A data frame with at least two columns: `sample` (matching column
  names of `count_data`) and a grouping column (e.g., "group") that
  specifies sample groupings.

- group_column:

  Character string specifying the column name in `sample_metadata` used
  to group samples. Default is `"group"`.

- output_dir:

  Directory to save plots and summary table. Created if it doesn't
  exist. Default is `"QC_dir"`.

- outlier_method:

  Method for identifying SD outliers: `"IQR"` (default) or `"z-score"`.

  - `"IQR"` flags values outside 1.5 × IQR from the 1st and 3rd
    quartiles.

  - `"z-score"` flags values with Z-scores exceeding `z_thresh`.

- z_thresh:

  Z-score threshold for flagging outliers. Only used when
  `outlier_method = "z-score"`. Default is `3`.

## Value

Invisibly returns a data frame with per-group statistics:

- Group:

  Group name

- Mean_SD:

  Mean of pairwise log2FC SDs

- SD_of_SDs:

  Standard deviation of SDs across log2FC pairs

- N_Pairs:

  Number of pairwise comparisons

- Outlier:

  Logical flag indicating SD outliers based on the selected method

## Examples

``` r
if (FALSE) { # \dontrun{
  analyze_within_group_SDs(
  my_data$data, 
  my_data$metadata, 
  outlier_method = "z-score",
  z_thresh = 2.5)
} # }
```
