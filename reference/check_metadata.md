# Validate Metadata Structure

This internal function checks the validity of a metadata `data.frame`,
ensuring it meets structural requirements such as unique column names,
non-empty row names, and appropriate formatting.

## Usage

``` r
check_metadata(x, verbose = FALSE)
```

## Arguments

- x:

  A `data.frame` containing sample metadata. Each row represents a
  sample.

- verbose:

  A logical value indicating whether to check if row and column names
  exceed a character length threshold. Default is `FALSE`.

## Value

Invisibly returns the validated input `data.frame`.

## Details

- Ensures `x` is a `data.frame` and not empty.

- Checks that column names are present, unique, and free from `NA` or
  blank values.

- Ensures row names are not missing, blank, or duplicated.

- If `verbose = TRUE`, warns if row or column names exceed 15
  characters.

## Examples

``` r
if (FALSE) { # \dontrun{
if(interactive()){

metadata <- data.frame(sample = c("C1", "C2", "T1", "T2"),
                       sample_id = c("C100", "C202", "T303", "T100"),
                       condition = c("Control", "Control", "Treatment",
                                     "Treatment"))

metadata2 <- data.frame(sample = c("C1", "C2", "T1", "T2"),
                        sample = c("C100", "C202", "T303", "T100"),
                        condition = c("Control", "Control", "Treatment",
                                      "Treatment"),
                        check.names = FALSE)

metadata3 <- data.frame(long_sample_column_name = c("C1", "C2", "T1", "T2"),
                        sample_id = c("C100", "C202", "T303", "T100"),
                        condition = c("Control", "Control", "Treatment",
                                      "Treatment"),
                        row.names = c("C100", "C202", "T303",
                                      "T100_long_row_name"))

metadata4 <- data.frame(sample_id = c("C1", "C2", "T1", "T2"),
                        sample = c("C100", "C202", "T303", "T100"),
                        condition = c("Control", "Control", "Treatment",
                                      "Treatment"),
                        row.names = c("C1", "C2", "T2", ""), check.rows = FALSE)


metadata  <- check_metadata(x = metadata, verbose = TRUE)
metadata2 <- check_metadata(x = metadata2, verbose = TRUE)
metadata3 <- check_metadata(x = metadata3, verbose = TRUE)
metadata4 <- check_metadata(x = metadata4, verbose = TRUE)

 }
} # }
```
