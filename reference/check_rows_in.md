# Internal function that checks that the rownames of set of data frames have matching rownames in a reference set

Internal function that checks that the rownames of set of data frames
have matching rownames in a reference set

## Usage

``` r
check_rows_in(obj = list(), ref_rows = c())
```

## Arguments

- obj:

  A list of data frames to evaluate

- ref_rows:

  A vector of row names to which the row names of each obj data frame
  will be compared

## Value

A logical, with TRUE indicating all data frames have matching rows in
the reference set, and FALSE indicates that the one or more data frames
do not have matching rownames with the reference
