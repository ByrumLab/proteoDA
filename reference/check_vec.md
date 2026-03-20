# Check if a Vector Includes a Set of Reference Values

This function checks if a reference vector (`ref`) is a subset of the
input vector (`x`). The function ensures that `x` has the same or
greater length than `ref` and that no duplicates exist in `ref`. Both
`x` and `ref` must follow R's syntax rules for valid variable names
(e.g., column names and row names).

## Usage

``` r
check_vec(x, ref)
```

## Arguments

- x:

  A vector of values to be checked. The length of `x` must be equal to
  or greater than `ref` with no duplicates. This is typically a metadata
  or annotation column.

- ref:

  A vector of allowed reference values. All values in `ref` must be
  unique, and it must contain no `NA` or empty values. This is typically
  column names or row names in a counts matrix.

## Value

Logical. Returns `TRUE` if all values in `ref` are found in `x`,
otherwise returns `FALSE`.

## Details

See [`make.names`](https://rdrr.io/r/base/make.names.html) for details
on syntactic names.

## Examples

``` r
if (FALSE) { # \dontrun{
if (interactive()) {
  # ref is a subset of x
  x <- paste0("gene_", 1:10)
  ref <- paste0("gene_", 3:7)
  check_vec(x = x, ref = ref)

  # x does not include all ref values
  x <- paste0("gene_", 1:5)
  ref <- paste0("gene_", 4:6)
  check_vec(x = x, ref = ref)

  # x does not include any ref values
  x <- paste0("gene_", 1:5)
  ref <- paste0("gene_", 20:22)
  check_vec(x = x, ref = ref)

  # x contains duplicates, NA, empty strings, or blank spaces
  ref <- c("gene_", "NA", " ")
  check_vec(x = x, ref = ref)

## x is < ref
## x contains duplicates
## x contains NA
## ref includes duplicates
## ref includes NA
## ref includes empty strings
## ref includes blank spaces

## ref contains values equal to blank spaces


samples <- data.frame(sample_id = c("C100","C202","T303","T100"),
                      sample    = c("C1","C2","T1","T2"),
                      group     = c("Con","Con","Treat","Treat"))

counts <- data.frame(C100 = c(0,10,21,3,4),
                     C202 = c(0,0,0,2,1),
                     T100 = c(15,33,21,55,42),
                     row.names = c(paste0("gene_",1:5)))

genes <- data.frame(gene_id   = c(paste0("gene_",1:5),"NA.","NA..1",""),
                    gene_name = LETTERS[1:8])


is_id <- unlist(lapply(samples,
                       function(x) {
                         check_vec(x = x, ref = colnames(counts))
                     }))

is_id <- sapply(as.list.data.frame(samples),
                function(x) {
                  check_vec(x = x, ref = colnames(counts))
                })

id_column <- colnames(samples)[is_id][1]

sort(samples[, id_column]); sort(colnames(counts))


## find id column in annotation matching counts row names
is_id <- unlist(lapply(genes,
                       function(x) {
                         check_vec(x = x, ref = row.names(counts))
                     }))


 }
} # }
```
