# Classify proteins under a factorial design

Applies a simple rule-based classifier to each protein using
per-contrast statistics extracted by
[`interpret_protein_factorial`](https://byrumlab.github.io/proteoDA/reference/interpret_protein_factorial.md).
Intended for BioID-style designs with average treatment, interaction,
and baseline (between-cell) contrasts.

## Usage

``` r
classify_factorial_batch(
  DA_results,
  id_vector = NULL,
  contrast_map = c(Avg_Treat_Bio_vs_DMSO = "Avg_Treat_Bio_vs_DMSO",
    Interaction_CHLA_vs_SKNF1 = "Interaction_CHLA_vs_SKNF1", CHLA90_Bio_vs_SKNF1_Bio =
    "CHLA90_Bio_vs_SKNF1_Bio", CHLA90_DMSO_vs_SKNF1_DMSO = "CHLA90_DMSO_vs_SKNF1_DMSO"),
  p_thresh = 0.05,
  lfc_thresh = 1
)
```

## Arguments

- DA_results:

  A DAList-like object with `$results` and optionally `$annotation`.

- id_vector:

  Optional character vector of protein IDs to classify. If `NULL`, the
  function tries `unique(DA_results$annotation$uniprot_id)`.

- contrast_map:

  Named character vector mapping display labels to keys in
  `DA_results$results`. Defaults match the typical factorial setup.

- p_thresh:

  Numeric. FDR threshold for significance tests. Default `0.05`.

- lfc_thresh:

  Numeric. Absolute log2 fold-change threshold used to decide whether
  the baseline contrast is *near zero*. Default `1`.

## Value

A `data.frame` with columns:

- `id` - protein ID,

- `class` - assigned class label,

- `base_logFC`, `base_q` - baseline contrast statistics,

- `treat_logFC`, `treat_q` - average treatment statistics,

- `int_logFC`, `int_q` - interaction statistics.

## Details

Classes returned include (non-exhaustive):

- `"Shared_Tx"` - significant average treatment effect; non-significant
  interaction.

- `"WTPreferential_Tx"` - significant average treatment; interaction
  significantly negative; baseline near zero.

- `"WTPreferential_Tx_withBaselineShift"` - as above but with
  significant baseline shift.

- `"TruncationGain"` - significant average treatment; interaction
  significantly positive; baseline near zero.

- `"TruncationGain_withBaselineShift"` - as above but with significant
  baseline shift.

- `"BaselineDifferenceOnly"` - significant baseline difference; average
  treatment not significant.

- `"Unclassified"` - none of the above rules matched.

## See also

[`interpret_protein_factorial`](https://byrumlab.github.io/proteoDA/reference/interpret_protein_factorial.md)

## Examples

``` r
if (FALSE) { # \dontrun{
cls <- classify_factorial_batch(DA_results)
table(cls$class)

# Use stricter thresholds
cls2 <- classify_factorial_batch(DA_results, p_thresh = 0.01, lfc_thresh = 1.5)
} # }
```
