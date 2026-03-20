# Metadata Grouping Column Validator

Internal function for validating the grouping column in a metadata
`data.frame`.

## Usage

``` r
check_grouping_column(metadata, grouping_column, verbose = FALSE)
```

## Arguments

- metadata:

  A `data.frame` containing sample metadata, where rows represent
  samples.

- grouping_column:

  The name or index of the column in `metadata` that provides
  information on how to group the samples. The function checks if
  `grouping_column`:

  - Exists in `metadata`.

  - Contains non-missing and non-blank values.

  - Has unique values (if required for analysis). Additionally, if
    `verbose = TRUE`, a warning is displayed when grouping values exceed
    15 characters to encourage better formatting.

- verbose:

  Logical. If `TRUE`, checks whether values in `grouping_column` exceed
  15 characters. Default: `FALSE`.

## Value

The validated column name of `grouping_column`, if all checks pass.

## Examples

``` r
if (FALSE) { # \dontrun{
if (interactive()) {

metadata <- data.frame(
  sample = c("C1", "C2", "T1", "T2"),
  sample_id = c("C100", "C202", "T303", "T100"),
  condition = c("Control", "Control", "Treatment", "Treatment")
)

metadata2 <- data.frame(
  sample = c("C1", "C2", "T1", "T2"),
  sample = c("C100", "C202", "T303", "T100"),
  condition = c("Control", "Control", "Treatment", "Treatment"),
  check.names = FALSE
)

metadata3 <- data.frame(
  sample = c("C2", "C2", "T1", "T2"),
  sample_id = c("C100", "C202", NA, "T100"),
  sample_name = c("sample_1", "sample_2", "", "sample_4"),
  condition = c("Baseline_Saline_Control", "Baseline_Saline_Control", "Treatment", "Treatment"),
  treatment = c("", "Control", "Chemo", "Chemo"),
  group = c("Control", NA, "Treatment", ""),
  chemo = c("Control", "Saline", "Treatment", "Treatment"),
  check.names = FALSE
)

check_grouping_column(metadata = metadata, grouping_column = "condition")
check_grouping_column(metadata = metadata, grouping_column = "replicates")
check_grouping_column(metadata = metadata, grouping_column = 7)
check_grouping_column(metadata = metadata2, grouping_column = 3)
check_grouping_column(metadata = metadata2, grouping_column = "sample")
check_grouping_column(metadata = metadata3, grouping_column = "condition")
check_grouping_column(metadata = metadata3, grouping_column = 5)
check_grouping_column(metadata = metadata3, grouping_column = "group")
check_grouping_column(metadata = metadata3, grouping_column = 7)

}
} # }
```
