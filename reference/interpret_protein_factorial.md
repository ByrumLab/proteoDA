# Interpret a protein across factorial contrasts

Extracts per-contrast statistics for a single protein from a DAList-like
object fitted with a factorial (e.g., BioID) design, and returns a
compact table plus a human-readable summary.

## Usage

``` r
interpret_protein_factorial(
  DA_results,
  protein,
  protein_col = NULL,
  alpha = 0.05,
  contrast_map = c(Avg_Treat_Bio_vs_DMSO = "Avg_Treat_Bio_vs_DMSO",
    Interaction_CHLA_vs_SKNF1 = "Interaction_CHLA_vs_SKNF1", TreatEffect_CHLA90 =
    "TreatEffect_CHLA90", TreatEffect_SKNF1 = "TreatEffect_SKNF1",
    CHLA90_Bio_vs_SKNF1_Bio = "CHLA90_Bio_vs_SKNF1_Bio", CHLA90_DMSO_vs_SKNF1_DMSO =
    "CHLA90_DMSO_vs_SKNF1_DMSO")
)
```

## Arguments

- DA_results:

  A DAList-like object with `$results` (per-contrast result tables) and
  optionally `$annotation`.

- protein:

  Character. Protein identifier (e.g., UniProt ID) to look up.

- protein_col:

  Optional character. Name of an explicit ID column in the per-contrast
  result tables. If `NULL`, matching falls back to row names and/or
  `DA_results$annotation`.

- alpha:

  Numeric. FDR threshold used to annotate the `effect` column. Default
  `0.05`.

- contrast_map:

  Named character vector mapping display labels (names) to keys inside
  `DA_results$results`. Defaults assume contrasts like average treatment
  effect, interaction, and per-cell effects.

## Value

A list with:

- `table`:

  A `data.frame` with columns `contrast`, `logFC`, `adj.P.Val`,
  `effect`, `t`, `B`, `P.Value`.

- `summary`:

  A single character string summarizing the key contrasts.

## Details

Row matching is flexible:

- If `protein_col` names a column in the per-contrast table, it is used
  directly.

- Otherwise, the function attempts to resolve the row using
  [`rownames()`](https://rdrr.io/r/base/colnames.html) or
  `DA_results$annotation` (preferring `annotation$uniprot_id` if
  available).

The returned table includes an `effect` label indicating direction
(up/down) and significance (`"*"` if `adj.P.Val < alpha`).

## Examples

``` r
if (FALSE) { # \dontrun{
out <- interpret_protein_factorial(DA_results, protein = "Q16666")
out$summary
out$table
} # }
```
