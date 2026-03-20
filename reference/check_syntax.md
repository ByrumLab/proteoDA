# Check if a vector follows R syntax rules

Internal helper function to verify whether values in a vector conform to
R's variable naming rules.

## Usage

``` r
check_syntax(x)
```

## Arguments

- x:

  A vector of values to check.

## Value

Logical. Returns `TRUE` if all values in `x` follow R syntax rules,
otherwise returns `FALSE`.

## Details

This function checks whether values in `x` comply with R's syntactic
name rules using
[`make.names()`](https://rdrr.io/r/base/make.names.html).

- A valid name consists of letters, numbers, dots (`.`), or underscores
  (`_`).

- It must start with a letter or a dot not followed by a number.

- See: [`base::make.names()`](https://rdrr.io/r/base/make.names.html)
  and [`base::make.unique()`](https://rdrr.io/r/base/make.unique.html).

## Examples

``` r
if (FALSE) { # \dontrun{
if (interactive()) {
  check_syntax(c("C1", "C2", "T1", "T2"))       # TRUE
  check_syntax(c("C 100", "C-202", "T 303"))    # FALSE
  check_syntax(c("Var_1", "Var.2", ".Valid"))   # TRUE
  check_syntax(c(1, 2, 3, 4))                   # TRUE (coerced to strings)
}
} # }
```
