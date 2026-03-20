# Identify Positive Integer Values

Internal helper function to check if a value is a positive whole number
(including zero).

## Usage

``` r
check_int(x)
```

## Arguments

- x:

  A single value to be tested. Must be of length 1L.

## Value

Logical. Returns `TRUE` if `x` is a positive integer (including 0),
otherwise returns `FALSE`.

## Examples

``` r
if (FALSE) { # \dontrun{
if (interactive()) {
  check_int(5)      # TRUE
  check_int(0)      # TRUE
  check_int(-3)     # FALSE
  check_int(4.5)    # FALSE
  check_int("ABC")  # FALSE
}
} # }
```
