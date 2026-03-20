# Validate filename

Internal function used to check constraints on user-defined filenames.
By default checkes whether there are spaces in the filename, and whether
the extension is in the list of allowed extensions.

## Usage

``` r
validate_filename(filename, allowed_exts, check_space = T)
```

## Arguments

- filename:

  File name to check

- allowed_exts:

  Character vector of allowable extensions, with no period

- check_space:

  Check for whitespace in file name?

## Value

If filename checks out, returns invisible TRUE. Otherwise, throws error.

## Examples

``` r
# No examples yet
```
