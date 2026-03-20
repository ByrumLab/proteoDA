# Check if a Value is a TRUE or FALSE Logical of Length 1

Internal helper function to verify if `x` is a logical value (TRUE or
FALSE) of length 1.

## Usage

``` r
check_logical(x)
```

## Arguments

- x:

  A value to be checked.

## Value

Logical. Returns `TRUE` if `x` is a logical value (`TRUE` or `FALSE`) of
length 1. Returns `FALSE` if `x` is any of the following: numeric,
character string, length greater than 1, empty string, a series of blank
spaces, `NA`, `NULL`, or non-logical values.

## Examples

``` r
if (FALSE) { # \dontrun{
if (interactive()) {
  check_logical(TRUE)    # TRUE
  check_logical(FALSE)   # TRUE
  check_logical(1)       # FALSE
  check_logical("TRUE")  # FALSE
  check_logical(c(TRUE, FALSE)) # FALSE
  check_logical(NA)      # FALSE
}
} # }
```
