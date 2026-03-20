# Identify Positive Numeric Values

Internal helper function to verify if `x` is a valid positive numeric
value (≥ 0).

## Usage

``` r
check_num(x)
```

## Arguments

- x:

  A single value to be tested.

## Value

Logical. Returns `TRUE` if `x` is a valid positive numeric value (≥ 0).
Returns `FALSE` for non-numeric inputs, characters, multiple values,
empty strings, white spaces, `NA`, `NULL`, `Inf`, or negative numbers.

## Examples

``` r
if (FALSE) { # \dontrun{
if (interactive()) {
  check_num(1.4323)   # TRUE
  check_num(0)        # TRUE
  check_num(-45)      # FALSE
  check_num(Inf)      # FALSE
  check_num("snap")   # FALSE
  check_num("  ")     # FALSE
  check_num(NA)       # FALSE
  check_num(1e6)      # TRUE
}
} # }
```
