# Validate Single Character String

Internal helper function to check whether `x` is a valid character
string of length 1.

## Usage

``` r
check_string(x)
```

## Arguments

- x:

  A value to be checked.

## Value

Logical. Returns `TRUE` if `x` is a character string of length 1.
Returns `FALSE` if `x` is numeric, a vector (length \> 1L), an empty
string, white spaces, `NA`, `NULL`, `Inf`, etc.

## Examples

``` r
if (FALSE) { # \dontrun{
if (interactive()) {
  check_string("ProjectID_122522") # TRUE
  check_string(123)                # FALSE
  check_string(NA)                 # FALSE
  check_string(NULL)               # FALSE
  check_string(c("oh", "snap"))    # FALSE
  check_string("     ")            # FALSE
  check_string("")                 # FALSE
  check_string(Inf)                # FALSE
}
} # }
```
