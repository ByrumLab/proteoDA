# Check if a file exists

Internal helper function to verify if an input string corresponds to an
existing file.

## Usage

``` r
check_file(x)
```

## Arguments

- x:

  A character string representing a file path.

## Value

Logical. Returns `TRUE` if `x` is a valid file path and the file exists.
Returns `FALSE` if `x` is not a valid string or the file does not exist.

## Examples

``` r
if (FALSE) { # \dontrun{
if (interactive()) {
  check_file("path/to/existing_file.txt") # TRUE if file exists
  check_file("nonexistent_file.txt")      # FALSE
}
} # }
```
