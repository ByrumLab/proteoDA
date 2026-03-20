# All pairwise differences of a vector

Taken from
<https://stackoverflow.com/questions/48445003/compute-all-pairwise-differences-of-elements-in-a-vector>,
uses a nice R matrix algebra trick, taking the outer product of the
vectors (in this case, the difference instead of the product). Then,
just return the lower triangle, so we don't double-count pw comparisons

## Usage

``` r
all_pw_diffs(vector)
```

## Arguments

- vector:

  A numeric vector of items for which you want all pairwise differences

## Value

A vector of all pairwise differences between the elements in the vector

## Examples

``` r
# No examples yet
```
