# Identify Duplicate Values

This function identifies duplicate values in a vector. It can be used to
find duplicate column names, row names, or values in any particular
column of a data frame.

## Usage

``` r
check_dup(x)
```

## Arguments

- x:

  A vector of values to be tested. This is typically a vector of column
  names, row names, or values in a particular column of a data frame.

## Value

The function returns a vector of duplicate values. If no duplicates are
found, it returns `NULL`. This allows for custom error handling if
desired.

## Examples

``` r
if (FALSE) { # \dontrun{
if(interactive()){
  x <- c("A", "A", "B", "C", "D", "D", "D", "E")
  y <- c(1, 2, 3, "A", "B", "C")
## returns vector of duplicates identified
  dup_present <- check_dup(x = x)
## returns NULL if all unique
  dup_not_present <- check_dup(x = y)

}
} # }
```
