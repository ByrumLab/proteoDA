# Make hyperlinks for an Excel column

Adds hyperlinks to a column in a data frame to be exported to an Excel
file.

## Usage

``` r
make_excel_hyperlinks(data, url.col, url)
```

## Arguments

- data:

  The data frame in which to add hyperlinks.

- url.col:

  The name of the column to add hyperlinks to.

- url:

  The URL that will be prepended to the info in the column.

## Value

The original data frame, now with hyperlinks in the desired column.
