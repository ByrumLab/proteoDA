# Generate colors for Sample Groups

Assigns colors to each unique sample group in metadata. If the number of
groups is 12 or fewer, colors are assigned using a color-blind friendly
palette. For more than 12 groups, the function uses the Polychrome
package to generate a distinct color palette.

## Usage

``` r
get_colors(group)
```

## Arguments

- group:

  A vector representing the sample group classifications. If not a
  factor, it is converted into one.

## Value

A named vector of colors, where names correspond to unique group labels.

## Details

- If there are 12 or fewer unique groups, the function assigns colors
  using a predefined color-blind friendly palette.

- If there are more than 12 groups, it leverages Polychrome to generate
  distinct colors.

- If grDevices is not installed, an error is thrown when handling more
  than 12 groups.

- group must be a vector.

## Examples

``` r
group <- factor(c("A", "B", "C"))
get_colors(group)
#>         A         B         C 
#> "#E69F00" "#56B4E9" "#009E73" 
```
