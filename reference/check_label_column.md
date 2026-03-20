# Metadata Label Column Validator

Internal function for validating the sample label column in a metadata
`data.frame` for a `DGEList` object.

## Usage

``` r
check_label_column(metadata, label_column, verbose = FALSE)
```

## Arguments

- metadata:

  A `data.frame` containing sample metadata, where rows represent
  samples.

- label_column:

  The name or index of the column in `metadata` that contains labels
  used to identify each sample. The function checks if `label_column`:

  - Exists in `metadata`.

  - Contains non-missing and non-blank values.

  - Has unique values.

  - Adheres to R syntax rules (a warning is issued if not).
    Additionally, if `verbose = TRUE`, a warning is displayed when label
    values exceed 15 characters to encourage better plot formatting.

- verbose:

  Logical. If `TRUE`, checks whether the label values exceed 15
  characters. Default: `FALSE`.

## Value

The validated column name of `label_column`, if all checks pass.

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
  condition = c("Control", "Control", "Treatment", "Treatment")
)

metadata4 <- data.frame(
  long_sample_column_name = c("C1", "C2", "T1", "T2"),
  sample_id = c("C100", "C202", "chemo_treatment_303", "chemo_treatment_100"),
  condition = c("Control", "Control", "Treatment", "Treatment"),
  row.names = c("C100", "C202", "T303", "T100_long_row_name")
)

check_label_column(metadata = metadata, label_column = "sample_id")
check_label_column(metadata = metadata, label_column = 2)
check_label_column(metadata = metadata, label_column = "replicates")
check_label_column(metadata = metadata2, label_column = "sample")
check_label_column(metadata = metadata3, label_column = 1)
check_label_column(metadata = metadata3, label_column = 2)
check_label_column(metadata = metadata3, label_column = "sample_name")
check_label_column(metadata = metadata4, label_column = "sample_id")

}
} # }

```
