# Check if an object is a data.frame or matrix

Internal helper function to determine whether an object is a
`data.frame` or `matrix`.

## Usage

``` r
check_data(x)
```

## Arguments

- x:

  An object to be tested (e.g., counts, metadata, annotation).

## Value

Logical. Returns `TRUE` if `x` is a `data.frame` or `matrix`, otherwise
returns `FALSE`.

## Examples

``` r
if (FALSE) { # \dontrun{
if (interactive()) {
  check_data(data.frame(a = 1:3, b = 4:6)) # TRUE
  check_data(matrix(1:9, nrow = 3))        # TRUE
  check_data(list(a = 1, b = 2))           # FALSE
}
} # }
```
