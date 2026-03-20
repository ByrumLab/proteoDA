# Remove proteins based on annotation data

This function is used to remove proteins from a DAList, filtering using
data in the annotation data frame of the DAList. Proteins which do not
produce a value of TRUE for the supplied condition are removed from both
the data and annotation slots of the DAList. If condition evaluates to
NA, the function will return an error.

## Usage

``` r
filter_proteins_by_annotation(DAList, condition)
```

## Arguments

- DAList:

  A DAList object to be filtered.

- condition:

  An expression that returns a logical value, defined in terms of
  variables present in the annotation data frame of the supplied DAList.
  Proteins are kept if the condition is TRUE for that protein.

## Value

A DAList, with proteins that do not meet the condition removed.

## Examples

``` r
if (FALSE) { # \dontrun{
# Suppose the DAList$annotation data frame contains three columns:
# uniprot_id = An alpha-numeric ID uniquely identifying a protein
# protein_name = A character giving the protein name
# molecular_weight = A numeric value, giving the protein's MW in kDa.

# Remove a specific protein by ID
filtered <- filter_proteins_by_annotation(DAList,
                                          uniprot_id != "abc123")

# Remove any protein which contains "keratin" in the name:
filtered <- filter_proteins_by_annotation(DAList,
                                          !grepl(pattern = "keratin",
                                                 x = protein_name))

# Remove any protein with molecular weight < 30 kDa
filtered <- filter_proteins_by_annotation(DAList,
                                          molecular_weight > 30)


# Filtering functions can be chained together
filtered <- DAList |>
  filter_proteins_by_annotation(!grepl(pattern = "keratin",
                                       x = protein_name)) |>
  filter_proteins_by_annotation(molecular_weight > 30)
} # }
```
