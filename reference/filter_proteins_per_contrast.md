# Filter proteins separately for each contrast using sample group info

This function filters proteins contrast-by-contrast based on non-missing
intensity values in a defined sample group column (e.g., "group" or
"condition").

## Usage

``` r
filter_proteins_per_contrast(
  DAList,
  contrasts_file,
  min_reps = 2,
  require_both_groups = TRUE,
  grouping_column = "group"
)
```

## Arguments

- DAList:

  A DAList object containing \$data, \$annotation, and \$metadata.

- contrasts_file:

  Path to a no-header .csv file with contrast names in the form
  "GroupA_vs_GroupB=GroupA - GroupB".

- min_reps:

  Integer. Minimum number of replicates per group to retain a protein.

- require_both_groups:

  Logical. If TRUE, protein must meet `min_reps` in both contrast
  groups. If FALSE, just one.

- grouping_column:

  Character. Name of the column in DAList\$metadata used to assign group
  membership.

## Value

A DAList object with updated \$filtered_proteins_per_contrast,
\$data_per_contrast, \$annotation_per_contrast, and a
\$tags\$retention_summary table.

## Details

Contrast names are read from a no-header CSV file where each row is of
the form: `GroupA_vs_GroupB=GroupA - GroupB`. The part before the `=` is
used as the contrast ID.

For each contrast, proteins are retained only if they have at least
`min_reps` values in either or both groups involved in the contrast (as
determined by `require_both_groups`). The function also stores
per-contrast versions of the \$data and \$annotation slots into
DAList\$data_per_contrast and DAList\$annotation_per_contrast.

## Examples

``` r
if (FALSE) { # \dontrun{
filtered_DAList <- filter_proteins_per_contrast(
  DAList = raw_DAList,
  contrasts_file = "data/contrasts.csv",
  min_reps = 2,
  require_both_groups = TRUE,
  grouping_column = "group"
)
} # }
```
