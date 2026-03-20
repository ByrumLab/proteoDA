# Check data–metadata alignment consistency

A lightweight check to verify that the columns in a data matrix and the
rows in the corresponding metadata frame refer to the same samples.
Typically used after importing data but before building a `DAList`.

## Usage

``` r
check_data_metadata_match(data, metadata)
```

## Arguments

- data:

  A numeric matrix or data frame of intensity values.

- metadata:

  A data frame of sample metadata with row names matching the sample IDs
  in `data`.

## Value

Invisibly returns TRUE if the match is valid. Stops with an informative
error if there are missing or mismatched samples.

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- data.frame(S1 = 1:3, S2 = 4:6)
meta <- data.frame(group = c("A", "B"), row.names = c("S1", "S2"))
check_data_metadata_match(dat, meta)  # passes

bad_meta <- data.frame(group = c("A", "B"), row.names = c("S1", "S3"))
check_data_metadata_match(dat, bad_meta)  # stops with descriptive error
} # }
```
