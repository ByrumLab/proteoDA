# Identify Values Above a Specified Character Length

This function identifies values in a vector that exceed a specified
maximum number of characters. This is often used to identify long sample
names or group labels that could cause issues in visualizations or
presentations.

## Usage

``` r
check_long(x, char_thresh = 15, na.rm = FALSE)
```

## Arguments

- x:

  A vector of values. Typically column names, row names, or metadata
  values.

- char_thresh:

  A non-negative numeric value specifying the maximum allowed character
  length. Values exceeding this length are returned.

- na.rm:

  Logical; if `TRUE`, `NA` values are removed before checking. If
  `FALSE`, `NA` and `NaN` may appear in the output.

## Value

A vector of values in `x` whose character lengths exceed `char_thresh`.
Returns `NULL` if no values meet the criterion.

## Details

Behavior under different combinations of `x` contents and `na.rm`:

- **No NAs in `x`**

  - `na.rm = TRUE` or `FALSE`, and **all values are below** the
    threshold → returns `NULL`.

  - `na.rm = TRUE` or `FALSE`, and **some values exceed** the threshold
    → returns those values.

- **NAs present in `x`**

  - `na.rm = TRUE`: NAs are removed; behavior is based only on non-NA
    values.

  - `na.rm = FALSE`:

    - If all non-NA values are below threshold → returns `c(NA, NaN)` if
      present.

    - If some exceed threshold → returns the long values and any
      `NA`/`NaN` values.

## Examples

``` r
if (FALSE) { # \dontrun{
x <- c("Con_01", "Con_02", NA, "TAC_Treatment_01", "", Inf, NaN, "TAC_Treat_02")
y <- c("Con_01", "Con_02", "Treat_01", "Treat_101")
z <- c("Con_01", NA, "TAC_Treatment_01", "TAC_Treat_02", "", Inf)

# identify values with more than 15 characters
check_long(x = x, char_thresh = 15, na.rm = FALSE)

# identify values with more than 8 characters
check_long(x = y, char_thresh = 8)

# identify values with more than 6 characters
check_long(x = y, char_thresh = 6)
} # }
```
