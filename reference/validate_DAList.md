# DAList validator

Internal function for validating and normalizing a DAList-like object.
Tolerant of arbitrary result labels (custom contrast names). Returns the
(possibly re-ordered) object and ensures class "DAList".

## Usage

``` r
validate_DAList(x)
```
