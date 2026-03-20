# Conversion of missing values and 0s

Depending on your experimental design, you may wish to treat intensity
values below the detection limit of your instrument as missing data (NA)
or as a value of 0 intensity. Though this can be done manually before
importing data into a DAList object, we also provide some simple utility
functions to convert missing values to 0 and vice-verse.
`missing_to_zero` converts missing values (NA by default, or other
user-provided values) to 0, and `zero_to_missing` converts 0s to NA.

## Usage

``` r
missing_to_zero(DAList, missing_values = c(NA))

zero_to_missing(DAList)
```

## Arguments

- DAList:

  A DAList object.

- missing_values:

  A vector of values that will be converted to 0. Default is NA, but
  users may include other numbers if they represent missing data in
  their dataset (e.g., -1).

## Value

A DAList, with missing values or 0s replaced according to the funciton
used.

## Examples

``` r
if (FALSE) { # \dontrun{

# Convert missing values to 0s
# Default is to convert NAs to 0s
no_missing <- missing_to_zero(DAList)

# But can convert other numbers as well:
no_missing <- missing_to_zero(DAList, missing = c(NA, -1))

# Convert 0s to missing values
missing_as_zero <- zero_to_missing(DAList)

} # }
```
