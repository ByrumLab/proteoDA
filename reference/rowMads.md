# Calculate per-row MADs of a numeric array

Have done testing, gives same results as the matrixStats::rowMads()
function it replaces, though it is much slower (the matrixStats version
uses C code).

## Usage

``` r
rowMads(x, ...)
```

## Arguments

- x:

  The array for which to calculate per-row MADs

- ...:

  Additional arguments to be passed to internal functions. Meant for
  na.rm.

## Value

A numeric vector of appropriate length, named if input was named, with
per-row MADs

## Examples

``` r
# No examples yet.
```
