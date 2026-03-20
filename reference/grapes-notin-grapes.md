# Not in operator

A convenience operator equivalent to `!(x %in% table)`., i.e., checks
whether elements of `x` are not in `table`.

## Usage

``` r
x %notin% table
```

## Arguments

- x:

  A vector of values to be matched.

- table:

  A vector of values to be matched against.

## Value

A logical vector indicating if elements of `x` are not in `table`.
