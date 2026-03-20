# Split a data.frame or matrix according to group membership.

Split a data.frame or matrix according group membership, creating a list
where each element consists of a data.frame/matrix for a particular
experimental group.

## Usage

``` r
split_by_groups(data, groups)
```

## Arguments

- data:

  a data.frame or matrix. For our pipeline this is typically a counts
  matrix or matrix of counts-per-million values.

- groups:

  a vector or factor giving the experimental group/condition for each
  sample (i.e. column) in data. The length of groups must equal the
  number of columns in data.

## Value

The function returns a list in which each element is a data.frame for an
individual group.

## Examples

``` r
if (FALSE) { # \dontrun{
if(interactive()){

 counts <- matrix(rpois(4*5, lambda = 5), nrow = 4,ncol = 5)
 rownames(counts) <- paste0("gene_", 1:4)
 colnames(counts) <- paste0("Sample_", 1:5)
 y <- read_rnaseq_data(counts = counts, project_id = "LastName_Date")
 y <- validate_DGEList(y = y)

 tmp <- split_by_groups(data   = y$counts,
                        groups = c("control","treat","treat","control","treat"))

 }
} # }
```
