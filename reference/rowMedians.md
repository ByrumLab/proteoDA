# Calculate per-row medians of a numeric array

Have done testing, gives same results as the matrixStats::rowMedians()
function it replaces, though it is much slower (the matrixStats version
uses C code).

## Usage

``` r
rowMedians(x, ...)
```

## Arguments

- x:

  The array for which to calculate per-row medians

- ...:

  Additional arguments to be passed to internal functions. Meant for
  na.rm.

## Value

A numeric vector of appropriate length, named if input was named, with
per-row medians

## Examples

``` r
# No examples yet.
```
