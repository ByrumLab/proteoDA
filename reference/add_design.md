# Prepare limma model design matrix

A design matrix is a model matrix of explanatory variables of a set of
objects. Each row represents individual samples and the columns
represent the sample groups and factors. This function utilizes the
function
[limma::modelMatrix](https://rdrr.io/pkg/limma/man/modelMatrix.html).

## Usage

``` r
add_design(DAList, design_formula = NULL)
```

## Arguments

- DAList:

  The DAList of normalized proteins.

- design_formula:

  A string for the design matrix using intercept, no intercept,
  additive, or interaction models. See examples below.

## Value

A DAList object with a design.

## Examples

``` r
if (FALSE) { # \dontrun{
# if x is a numerical covariate, then ~x and ~0+x are different models where the
# second model assumes the expected response is zero when x is zero.
# if A is a factor, then ~A and ~0+A are the same model, just parametrized slightly different.

# Intercept Model
data <- add_design(data, "~ group")

# No Intercept Model
data <- add_design(data, "~ 0 + group")

# Additive Model
data <- add_design(data, "~ 0 + group + batch")
data <- add_design(data, "~ 0 + group + pair_sample") # fixed effects
data <- add_design(data, "~ 0 + group + (1 | xyz)")   # mixed effects (fix + random effect)

# Interaction Model
data <- add_design(data, "~ group*gender")

} # }

```
